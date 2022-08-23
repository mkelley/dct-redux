#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord

from calviacat import RefCat2
from calviacat.catalog import CalibrationError

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+")
parser.add_argument(
    "--fetch",
    default="default",
    choices=["default", "all", "none"],
    help="fetch catalogs as needed (default), for all images (all), or for none of them (none)",
)
parser.add_argument(
    "--filter",
    action="append",
    help="only calibrate this filter, may be used multiple times",
)
parser.add_argument(
    "--no-color-correction",
    action="store_true",
    help="do not derive a color correction",
)
parser.add_argument(
    "--plot",
    action="store_true",
    help='save a plot of the results under "zeropoints/"',
)
parser.add_argument(
    "--fit-color-thresh",
    type=int,
    default=30,
    help="number of sources required to fit the color term, or -1 to disable",
)
parser.add_argument(
    "-n",
    action="store_true",
    help="no-operation mode, FITS file is opened read only, no plots are saved",
)
args = parser.parse_args()


# setup logging
logger = logging.Logger("LMI catalog calibration")
logger.setLevel(logging.DEBUG)

# this allows logging to work when run multiple times in a python session
if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(levelname)s: %(message)s")

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler("lmi-catalog-calibrate.log")
    logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info("#" * 70)
logger.info(datetime.now().isoformat())

if args.plot and not args.n:
    logger.info("saving plots to ./zeropoints")
    if not os.path.exists("zeropoints"):
        os.mkdir("zeropoints")


LMI2REFCAT_FILTER = {
    "SDSS-G": "g",
    "SDSS-R": "r",
    "SDSS-I": "i",
    "SDSS-Z": "z",
    "VR": "r",
    "BC": "g",
    "RC": "i",
}

# do not define C if no color correction is to be used
# note, these are preliminary default color corrections, 2022-07-29
cal_kwargs = {
    "SDSS-G": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=0.162),
    "SDSS-R": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=0.011),
    "SDSS-I": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=0.020),
    "SDSS-Z": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=-0.050),
    "VR": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=-0.062),
    "BC": dict(mlim=[14, 19], gmi_lim=[0.2, 3.0], C=-0.503),
    "RC": dict(mlim=[14, 19], gmi_lim=[0.2, 1.5], C=-0.135),
}

cal_color = {
    "SDSS-G": "g-r",
    "SDSS-R": "g-r",
    "SDSS-I": "g-r",
    "SDSS-Z": "g-r",
    "VR": "g-r",
    "BC": "g-r",
    "RC": "g-r",
}

# only calibrate these filters
if args.filter is not None:
    LMI2REFCAT_FILTER = {
        filt: LMI2REFCAT_FILTER.get(filt) for filt in args.filter
    }

# track centers of the fetched catalogs
fetched = []

for f in args.files:
    basename = os.path.basename(f)[:-5]

    mode = "readonly" if args.n else "update"
    hdu = fits.open(os.path.realpath(f), mode=mode)
    h = hdu[0].header
    try:
        filt = LMI2REFCAT_FILTER[h["FILTERS"]]
    except KeyError:
        logger.info(f'{basename}: not calibrating {h["FILTERS"]}')
        hdu.close()
        continue

    if "cat" not in hdu:
        logger.info(f"{basename}: no catalog to calibrate")
        hdu.close()
        continue

    logger.info(basename)

    if "mask" in hdu:
        mask = hdu["mask"].data.astype(bool)
    else:
        mask = np.zeros_like(hdu[0].data)

    phot = Table(hdu["cat"].data)
    dct = SkyCoord(phot["ra"], phot["dec"], unit="deg")
    refcat2 = RefCat2("cat.db")
    refcat2.logger.addHandler(logger.handlers[1])

    # get a new field if more than an arcmin from any previous field
    center = SkyCoord(dct.ra.mean(), dct.dec.mean())
    distance = (
        center.separation(SkyCoord(fetched)).min()
        if len(fetched) > 0
        else 180 * u.deg
    )
    if (
        not h.get("cfetched", False)
        and (distance > 1 * u.arcmin)
        and (args.fetch != "none")
    ) or (args.fetch == "all"):
        refcat2.fetch_field(dct)
        fetched.append(center)
        h["cfetched"] = (True, "catalog fetched for this image")

    try:
        objids, distances = refcat2.xmatch(dct)
    except TypeError:
        refcat2.db.close()
        hdu.close()
        continue

    m_inst = -2.5 * np.log10(phot["flux"])
    m_err = phot["fluxerr"] / phot["flux"] * 1.0857

    # cal_args = (objids, m_inst, filt)
    kwargs = cal_kwargs[h["filters"]].copy()
    color = cal_color.get(h["filters"])
    C = None
    cal_constant = args.no_color_correction or "C" not in kwargs
    try:
        if cal_constant:
            zp, C, zp_unc, m, gmr, gmi = refcat2.cal_constant(
                objids, m_inst, filt, **kwargs
            )
        else:
            if (
                objids.compressed().size > args.fit_color_thresh
                and args.fit_color_thresh > 0
            ):
                # let C vary freely
                kwargs.pop("C")

            zp, C, zp_unc, m, gmr, gmi = refcat2.cal_color(
                objids, m_inst, filt, color, **kwargs
            )
    except CalibrationError:
        continue
    finally:
        refcat2.db.close()

    n_cal = (~m.mask).sum()

    logger.info(
        "Calibrated with %d sources: zp=%.3f C=%s unc=%.3f",
        n_cal,
        zp,
        "None" if C is None else "{:.3f}".format(C),
        zp_unc,
    )

    # save matched stars to FITS file catalog
    cat_header = hdu["cat"].header
    tab = Table(hdu["cat"].data)
    tab["rf2id"] = objids.filled(0)
    tab["m_rf2"] = m.filled(99)
    tab["g-r_rf2"] = gmr.filled(99)
    cat_header["calfilt"] = filt, "catalog filter"
    cat_header["color"] = "g-r", "catalog color"
    cat_header["nsrcs"] = n_cal, "number of sources fit"
    cat_header["magzp"] = zp, "zeropoint magnitude"
    cat_header["clrcoeff"] = C, "color coefficient"
    cat_header["mzpunc"] = zp_unc, "zeropoint uncertainty"

    del hdu["cat"]
    hdu.append(fits.BinTableHDU(tab, name="cat", header=cat_header))
    hdu.close()

    if args.plot:
        fig = plt.figure(1)
        fig.clear()
        ax = fig.gca()

        ax.scatter(gmr, m - m_inst, marker=".", color="k")

        x = np.linspace(0, 1.5)
        label = "{:.4f} + {:.4f} ($g-r$)".format(zp, C)
        ax.plot(x, C * x + zp, "r-", label=label)

        plt.setp(ax, xlabel="$g-r$ (mag)", ylabel=r"$m-m_{\rm inst}$ (mag)")
        plt.legend()
        plt.tight_layout()
        if not args.n:
            plt.savefig("zeropoints/{}.png".format(basename))
