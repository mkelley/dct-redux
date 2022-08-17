#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from astropy.io.registry import identify_format

import numpy as np
from astropy.io import ascii, fits
from astropy.io.votable import parse_single_table
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.time import Time
import astropy.units as u
from mskpy.photometry import airmass_app

parser = argparse.ArgumentParser(
    "lmi-standard-phot.py",
    description="Extract standard stars from LMI photometry catalogs.",
    epilog="Use after lmi-add-cat.py.",
)

parser.add_argument("files", nargs="+", help="files to search for photometry")
parser.add_argument(
    "--catalogs",
    default=os.path.expanduser("~/data/catalogs"),
    help="directory holding the catalog files.",
)
parser.add_argument(
    "--dmax",
    type=Angle,
    default="2 arcsec",
    help="maximum distance for catalog matches",
)
parser.add_argument(
    "--hb-dmax",
    type=Angle,
    default="20 arcsec",
    help="maximum distance for HB catalog matches",
)
parser.add_argument(
    "--keep-all",
    action="store_true",
    help="when there are multiple matches, save all measurements, otherwise only save the brightest",
)
parser.add_argument("--logfile", default="lmi-standard-phot.log")
parser.add_argument("-o", default="standard-phot.txt", help="output file name")
parser.add_argument(
    "-f", action="store_true", help="force overwrite of output file"
)

args = parser.parse_args()

if os.path.exists(args.o) and not args.f:
    raise SystemExit("Refusing to overwrite output file.")

HB_FILTERS = [
    "OH",
    "NH",
    "UC",
    "CN",
    "C3",
    "CO+",
    "BC",
    "C2",
    "GC",
    "H2O+",
    "RC",
]
elevation = 2380 * u.m  # DCT elevation


def setup_logging(logfile):
    logger = logging.Logger("lmi-standard-phot.py")
    logger.setLevel(logging.DEBUG)

    # This test allows logging to work when it is run multiple times
    # from ipython
    if len(logger.handlers) == 0:
        formatter = logging.Formatter("%(levelname)s: %(message)s")

        console = logging.StreamHandler(sys.stdout)
        console.setLevel(logging.DEBUG)
        console.setFormatter(formatter)
        logger.addHandler(console)

        logfile = logging.FileHandler(logfile)
        logfile.setLevel(logging.INFO)
        logfile.setFormatter(formatter)
        logger.addHandler(logfile)

    logger.info("#" * 70)
    logger.info(Time.now().iso + "Z")
    logger.info("Command line: " + " ".join(sys.argv[1:]))
    for handler in logger.handlers:
        if hasattr(handler, "baseFilename"):
            logger.info("Logging to " + handler.baseFilename)

    return logger


def landolt09():
    table = parse_single_table(
        args.catalogs + "/landolt09-sources.xml"
    ).to_table(use_names_over_ids=True)
    cat = dict()
    cat["names"] = table["Star"].data.data.astype(str)
    ra = Angle(
        (table["RAh"].data, table["RAm"].data, table["RAs"].data),
        unit=u.hourangle,
    )
    dec = Angle(
        (table["DEd"].data, table["DEm"].data, table["DEs"].data), unit=u.deg
    )
    dec[table["DE-"].data == b"-"] *= -1
    cat["coords"] = SkyCoord(ra, dec)
    cat["V"] = table["Vmag"].data.data
    cat["R"] = table["Vmag"].data.data - table["V-R"].data.data
    cat["R_err"] = np.sqrt(
        table["e_Vmag"].data.data ** 2 + table["e_V-R"].data.data ** 2
    )
    cat["color"] = table["V-R"].data.data
    cat["color_err"] = table["e_V-R"].data.data
    return cat


def smith02():
    table = parse_single_table(
        "/home/msk/data/catalogs/smith02-standards.xml"
    ).to_table(use_names_over_ids=True)
    cat = dict()
    cat["names"] = table["Name"].data.data.astype(str)
    ra = [s for s in table["RAJ2000"]]
    dec = [s for s in table["DEJ2000"]]
    cat["coords"] = SkyCoord(ra, dec, unit=[u.hourangle, u.deg])
    cat["SDSS-R"] = table["r'mag"].data.data
    cat["SDSS-R_err"] = table["e_r'mag"].data.data
    cat["color"] = table["g'-r'"].data.data
    cat["color_err"] = table["e_g'-r'"].data.data
    cat["SDSS-G"] = table["r'mag"].data.data + cat["color"]
    cat["SDSS-G_err"] = np.sqrt(
        table["e_r'mag"].data.data ** 2 + cat["color_err"] ** 2
    )
    return cat


def farnham00():
    table = ascii.read("/home/msk/data/catalogs/hb-standards.txt")
    cat = dict()
    cat["names"] = ["HD {}".format(n) for n in table["HD"]]
    cat["coords"] = SkyCoord(
        table["RA"], table["Dec"], unit=[u.hourangle, u.deg]
    )
    for filt in HB_FILTERS:
        cat[filt] = np.ma.filled(table[filt], np.nan).data
        cat[filt + "_err"] = np.repeat(0.01, len(table))
    cat["color"] = np.ma.filled(table["B-V"], np.nan).data
    cat["color_err"] = np.repeat(0.01, len(table))
    return cat


logger = setup_logging(args.logfile)

catalogs = {
    "Farnham et al. 2000": farnham00(),
    "Landolt 2009": landolt09(),
    "Smith et al 2002": smith02(),
}

columns = [
    "file",
    "catalog",
    "object",
    "date",
    "za",
    "airmass",
    "filter",
    "color",
    "dist",
    "m",
    "m_err",
    "m_inst",
    "m_inst_err",
]
formats = {
    "airmass": "{:.4f}",
    "color": "{:.4f}",
    "dist": "{:3f}",
    "m": "{:.4f}",
    "m_err": "{:.4f}",
    "m_inst": "{:.4f}",
    "m_inst_err": "{:.4f}",
}

rows = []
for f in sorted(args.files):
    with fits.open(f, mode="readonly") as hdu:
        h = hdu[0].header
        if h["OBSTYPE"] != "OBJECT":
            logger.info(f + ": wrong OBSTYPE")
            continue

        if "CAT" not in hdu:
            logger.info(f + ": missing photometry catalog")
            continue

        date = h["DATE-OBS"].replace("T", " ")
        za = h["ZA"]
        am = airmass_app(h["ZA"] * u.deg, elevation)
        filt = h["filters"]
        exptime = h["exptime"]

        phot = Table(hdu["CAT"].data)
        phot = phot[phot["krflux"] > 0]  # clean out bad data
        w = WCS(hdu[0])
        radec = np.array(
            w.pixel_to_world_values(list(zip(phot["x"], phot["y"])))
        )
        coords = SkyCoord(*radec.T, unit="deg")
        m = -2.5 * np.log10(phot["krflux"] / exptime)
        merr = 1.0857 * phot["krfluxerr"] / phot["krflux"]
        n = 0
        for name, cat in catalogs.items():
            if filt not in cat:
                continue
            indices, dist, d3 = coords.match_to_catalog_sky(cat["coords"])

            if "Farnham" in name:
                i = dist < args.hb_dmax
            else:
                i = dist < args.dmax

            indices, dist, d3 = indices[i], dist[i], d3[i]
            for i in set(indices):
                n += 1
                matches = np.flatnonzero((indices == i))
                brightest = matches[m[matches].argmin()]
                for j in matches:
                    if not args.keep_all and j != brightest:
                        continue
                    rows.append(
                        (
                            f,
                            name,
                            cat["names"][i],
                            date,
                            za,
                            am,
                            filt,
                            cat["color"][i],
                            dist[j].arcsec,
                            cat[filt][i],
                            cat[filt + "_err"][i],
                            m[j],
                            merr[j],
                        )
                    )

        logger.info("{}: {} standards found".format(f, n))

if len(rows) == 0:
    logger.info("No standards found in any file.")
else:
    tab = Table(rows=rows, names=columns)
    tab.meta["comments"] = [Time.now().iso]
    for c, f in formats.items():
        tab[c].format = f
    tab.sort(["object", "filter", "date"])
    tab.write(args.o, format="ascii.fixed_width_two_line", overwrite=args.f)
