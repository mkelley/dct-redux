#!/usr/bin/env python3
import os
import re
import argparse
import warnings
import numpy as np
from astropy.wcs import WCS, FITSFixedWarning
from astropy.io import fits
import astropy.coordinates as coords
from astropy.coordinates import SkyCoord
from astropy.table import vstack
import astropy.units as u

# from astroquery.jplhorizons import Horizons
from sbpy.data import Ephem
import mskpy

comet_pat = (
    "^("
    "([1-9]{1}[0-9]*[PD](-[A-Z]{1,2})?)"
    "|([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?)"
    ")"
)
asteroid_pat = (
    "^("
    "([1-9][0-9]*( [A-Z]{1,2}([1-9][0-9]{0,2})?)?)"
    "|(\(([1-9][0-9]*)\))"
    "|(A/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?)"
    ")"
)

######################################################################
# suppress unnecessary warnings
warnings.simplefilter("ignore", FITSFixedWarning)
######################################################################

parser = argparse.ArgumentParser(
    description="Add a WCS centered on a moving target."
)
parser.add_argument("files", nargs="+", help="FITS images to update.")
parser.add_argument(
    "--observatory", default="G37", help="Observatory location (for HORIZONS)."
)
parser.add_argument(
    "--source", choices=["jpl", "mpc"], default="jpl", help="ephemeris source"
)

args = parser.parse_args()

# first pass to get target and dates
objects = []
dates = []
for f in args.files:
    hdu = fits.open(f, mode="readonly")
    objects.append(hdu[0].header["OBJECT"])
    dates.append(hdu[0].header["DATE-OBS"])

obj = np.unique(objects)
assert len(obj) == 1, "Multiple objects found: {}".format(obj)
t = mskpy.date2time(dates)

# put files in time order
files = [args.files[i] for i in np.argsort(t)]
t = t[np.argsort(t)]

obj = obj[0]
comet_match = re.findall(comet_pat, obj)
asteroid_match = re.findall(asteroid_pat, obj)
m = comet_match + asteroid_match
if len(m) == 0:
    raise ValueError(
        "Designation is not that of a comet or asteroid: {}".format(obj)
    )
target = m[0][0]

c = []
print(obj)

if args.source == "jpl":
    if comet_match:
        opts = dict(
            closest_apparition=True, no_fragments=True, id_type="designation"
        )
    else:
        opts = {}

    if target.isdigit():
        opts["id_type"] = "smallbody"

    # only get ~30 at a time
    eph = None
    chunks = max(len(t) // 30, 1)
    for _t in np.array_split(t, chunks):
        _eph = Ephem.from_horizons(target, epochs=_t, location="G37", **opts)
        if eph is None:
            eph = _eph
        else:
            eph.table = vstack((eph.table, _eph.table))

elif args.source == "mpc":
    eph = Ephem.from_mpc(target, epochs=t, location="G37")

# `MaskedQuantity`s cause problems with SkyCoord.to_string
ra = eph["RA"].filled(np.nan)
dec = eph["Dec"].filled(np.nan)
for j in range(len(eph)):
    _c = SkyCoord(ra[j], dec[j])
    print("  ", t[j].iso, _c.to_string("hmsdms"))
    c.append(_c)

c = coords.concatenate(c)
k = t.argmin()
for j, f in enumerate(files):
    hdu = fits.open(os.path.realpath(f), mode="update")

    w0 = WCS(hdu[0].header)

    w = WCS(naxis=2)
    w.wcs.crpix = w0.wcs_world2pix(c[j].ra, c[j].dec, 1)
    w.wcs.crval = (c[k].ra.deg, c[k].dec.deg)
    try:
        w.wcs.cd = w0.wcs.cd
    except AttributeError:
        pass
    w.wcs.ctype = w0.wcs.ctype
    hdu[0].header.update(w.to_header(key="M"))

    s = "Added moving target WCS (key M)."
    for line in hdu[0].header.get("HISTORY", []):
        if s in line:
            break
    else:
        hdu[0].header.add_history(s)

    hdu.close()
