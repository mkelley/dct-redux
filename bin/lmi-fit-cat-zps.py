#!/usr/bin/env python3
import os
import sys
import logging
import argparse
from glob import glob
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.modeling import models, fitting

from mskpy.photometry import airmass_app


def int_list(s):
    return [int(x) for x in s.split(",")]


def float_list(s):
    return [float(x) for x in s.split(",")]


parser = argparse.ArgumentParser()
parser.add_argument("path", help="directory with LMI data to calibrate")
parser.add_argument(
    "--nrange",
    type=int_list,
    default=[[None, None]],
    action="append",
    help="comma-separated file number ranges to calibrate, may be repeated",
)
parser.add_argument(
    "--arange",
    type=float_list,
    default=[],
    action="append",
    help="comma-separated airmass ranges to calibrate, may be repeated",
)
parser.add_argument(
    "--ignore",
    "-i",
    type=int_list,
    default=[],
    help="comma-separated file numbers to ignore",
)
args = parser.parse_args()

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
logger.info("Command-line options: %s", " ".join(sys.argv[1:]))
logger.debug("Ignore list: %s", str(args.ignore))

files = sorted(glob(f"{args.path}/lmi*fits"))
if len(files) == 0:
    logger.error(f"No lmi*.fits files at {args.path}")
    exit(1)

DCT_elevation = 2361 * u.m  # elevation of DCT

summary = []
for f in files:
    hdu = fits.open(f, mode="readonly")
    if "cat" not in hdu:
        continue

    ch = hdu["cat"].header
    if "magzp" not in ch:
        continue

    h = hdu[0].header
    filt = h["filters"]
    airmass = airmass_app(h["ZA"] * u.deg, DCT_elevation)
    zp = ch["magzp"] - 2.5 * np.log10(h["exptime"])
    summary.append(
        {
            "file": os.path.basename(f),
            "number": h["OBSERNO"],
            "date": h["DATE-OBS"],
            "filter": filt,
            "airmass": airmass,
            "exptime": h["exptime"],
            "catalog filter": ch["calfilt"],
            "catalog color": ch["color"],
            "n sources": ch.get("nsrcs", -1),
            "magzp": np.round(zp, 5),
            "clrcoeff": np.round(ch["clrcoeff"], 5),
            "mzpunc": np.round(ch["mzpunc"], 5),
            "ignore": (h["OBSERNO"] in args.ignore) or (ch["mzpunc"] == 0),
        }
    )

    if ch.get("nsrcs", -1) < 3:
        continue

if len(summary) == 0:
    raise ValueError(f"No calibrated catalogs found in {len(files)} files.")

summary = Table(summary)
summary.meta["comments"] = ["magzp includes correction for exposure time"]
summary.write(
    "catalog-cal-summary.txt",
    format="ascii.fixed_width_two_line",
    overwrite=True,
)
filters = set(summary["filter"])

analyze = summary[~summary["ignore"]]

plt.clf()
for filt in filters:
    i = analyze["filter"] == filt
    plt.errorbar(
        analyze["airmass"][i],
        analyze["magzp"][i],
        analyze["mzpunc"][i],
        ls="none",
        marker="o",
        alpha=0.5,
        label=filt,
    )


def fit_and_plot(summary):
    fit, mask = fitter(
        models.Linear1D(),
        summary["airmass"],
        summary["magzp"],
        weights=summary["mzpunc"] ** -2,
    )
    mms = sigma_clipped_stats(fit(summary["airmass"]) - summary["magzp"])

    a = np.linspace(1, max(3, max(summary["airmass"]) + 1))
    plt.plot(a, fit(a), color="k", ls=":", lw=0.5)

    # label = f"{filt}, {start}-{stop} (σ={mms[2]:.3f})"
    a = np.linspace(summary["airmass"].min(), summary["airmass"].max())
    plt.plot(a, fit(a), ls="-", lw=1.5, alpha=0.5)

    return {
        "zeropoint": round(fit.intercept.value, 5),
        "extinction": round(fit.slope.value, 5),
        "stdev": round(mms[2], 5),
    }


fitter = fitting.FittingWithOutlierRemoval(
    fitting.LinearLSQFitter(), sigma_clip, niter=3, sigma=3.0
)
extinction = []
for filt in filters:
    for start, stop in args.nrange:
        _start = min(analyze["number"]) if start is None else start
        _stop = max(analyze["number"]) if stop is None else stop
        i = (
            (analyze["number"] >= _start)
            * (analyze["number"] <= _stop)
            * (analyze["filter"] == filt)
        )
        if not any(i):
            continue

        row = {
            "file start": start,
            "file stop": stop,
            "airmass start": None,
            "airmass stop": None,
            "filter": filt,
        }
        result = fit_and_plot(analyze[i])
        row.update(result)
        extinction.append(row)

for filt in filters:
    for start, stop in args.arange:
        i = (
            (analyze["airmass"] >= start)
            * (analyze["airmass"] <= stop)
            * (analyze["filter"] == filt)
        )
        if not any(i):
            continue

        row = {
            "file start": None,
            "file stop": None,
            "airmass start": start,
            "airmass stop": stop,
            "filter": filt,
        }
        result = fit_and_plot(analyze[i])
        row.update(result)
        extinction.append(row)

tab = Table(extinction)
tab.write(
    "catalog-extinction.txt",
    format="ascii.fixed_width_two_line",
    overwrite=True,
)
plt.legend()
ax = plt.gca()
plt.setp(ax, ylabel="zeropoint (mag)", xlabel="X, airmass", xscale="log")
ax.minorticks_on()
ax.invert_yaxis()
plt.savefig("catalog-extinction.png")
plt.savefig("catalog-extinction.pdf")
