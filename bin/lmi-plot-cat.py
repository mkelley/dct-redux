#!/usr/bin/env python3
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import ascii, fits
from astropy.io.votable import parse_single_table
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import astropy.units as u

parser = argparse.ArgumentParser(description="Plot LMI image and catalog.")
parser.add_argument(
    "filename", help="FITS file to plot, must have a catalog extension"
)
parser.add_argument("-o", help="save plot to this output filename")
parser.add_argument(
    "--show",
    action="store_true",
    help="show the plot in a window (default if -o is not used)",
)
parser.add_argument(
    "--catalogs",
    default=os.path.expanduser("~/data/catalogs"),
    help="directory holding the catalog files.",
)
parser.add_argument("--dpi", default=200, help="output dots per inch")

args = parser.parse_args()

if args.filename == args.o:
    raise ValueError("Refusing to overwrite input file with a plot.")


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


catalogs = {
    "Farnham et al. 2000": farnham00(),
    "Landolt 2009": landolt09(),
    "Smith et al 2002": smith02(),
}

show = True if args.o is None else args.show

hdu = fits.open(args.filename)
im = hdu[0].data
vmin, vmax = ZScaleInterval().get_limits(im)
filt = hdu[0].header["filters"]
wcs = WCS(hdu[0].header)
cat = Table(hdu["cat"].data)
marker_range = np.sqrt(cat["krflux"].ptp())
marker_size = (
    np.sqrt(cat["krflux"] - cat["krflux"].min()) * 38 / marker_range + 2
)

fig = plt.figure(1, (8, 8), clear=True)
ax = plt.subplot(projection=wcs)
ax.imshow(im, vmin=vmin, vmax=vmax, origin="lower", cmap="gray")
ax.scatter(
    cat["x"],
    cat["y"],
    s=marker_size,
    color="tab:red",
    facecolor="none",
    lw=1,
    label="Sources",
)

for name, cat in catalogs.items():
    if filt not in cat:
        continue

    x, y = skycoord_to_pixel(cat["coords"], wcs)
    i = (x >= 0) * (y >= 0) * (x < im.shape[1]) * (y < im.shape[0])
    ax.scatter(x[i], y[i], s=10, color="tab:blue", lw=1, marker="x", label=name)

plt.setp(ax, xlabel="RA (deg)", ylabel="Dec (deg)")
plt.legend()

if show:
    plt.show()

if args.o is not None:
    plt.savefig(args.o, dpi=args.dpi)
