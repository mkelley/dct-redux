#!/usr/bin/env python3
import os
import sys
import logging
from glob import glob
from datetime import datetime
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.modeling import models, fitting

from mskpy.photometry import airmass_app

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


DCT_elevation = 2361 * u.m  # elevation of DCT

summary = []
airmass = defaultdict(list)
zp = defaultdict(list)
unc = defaultdict(list)
for f in sorted(glob("../ppp/lmi*fits")):
    hdu = fits.open(f, mode="readonly")
    if "cat" not in hdu:
        continue

    ch = hdu["cat"].header
    if "magzp" not in ch:
        continue

    h = hdu[0].header
    filt = h["filters"]
    summary.append(
        {
            "file": os.path.basename(f),
            "filter": filt,
            "catalog filter": ch["calfilt"],
            "catalog color": ch["color"],
            "n sources": ch.get("nsrcs", -1),
            "magzp": np.round(ch["magzp"], 5),
            "clrcoeff": np.round(ch["clrcoeff"], 5),
            "mzpunc": np.round(ch["mzpunc"], 5),
        }
    )

    if ch.get("nsrcs", -1) < 3:
        continue

    zp[filt].append(ch["magzp"] - 2.5 * np.log10(h["exptime"]))
    unc[filt].append(ch["mzpunc"])
    airmass[filt].append(airmass_app(h["ZA"] * u.deg, DCT_elevation))


Table(summary).write(
    "catalog-cal-summary.txt",
    format="ascii.fixed_width_two_line",
    overwrite=True,
)

plt.clf()
fitter = fitting.FittingWithOutlierRemoval(
    fitting.LinearLSQFitter(), sigma_clip, niter=3, sigma=3.0
)
extinction = []
for filt in zp.keys():
    fit, mask = fitter(
        models.Linear1D(),
        np.array(airmass[filt]),
        np.array(zp[filt]),
        weights=np.array(unc[filt]) ** -2,
    )

    mms = sigma_clipped_stats(fit(airmass[filt]) - np.array(zp[filt]))
    a = np.linspace(1, max(3, max(airmass[filt]) + 1))
    label = f"{filt}: {fit.intercept.value:.3f} + {fit.slope.value:.3f} X (σ={mms[2]:.3f})"
    plt.errorbar(
        airmass[filt],
        zp[filt],
        unc[filt],
        ls="none",
        marker="o",
        alpha=0.5,
        label=label,
    )
    plt.plot(a, fit(a), color="k", ls="-", lw=1)
    extinction.append([filt] + [fit.intercept.value, fit.slope.value, mms[2]])

tab = Table(
    rows=extinction, names=["filter", "zeropoint", "extinction", "stdev"]
)
tab.write(
    "catalog-extinction.txt",
    format="ascii.fixed_width_two_line",
    overwrite=True,
)
plt.legend()
ax = plt.gca()
plt.setp(ax, ylabel="zeropoint (mag)", xlabel="airmass", xscale="log")
ax.minorticks_on()
ax.invert_yaxis()
plt.savefig("lmi-zeropoints.png")
