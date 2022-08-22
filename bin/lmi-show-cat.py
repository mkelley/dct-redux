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

parser = argparse.ArgumentParser(description="Print LMI source catalog.")
parser.add_argument(
    "filename", help="FITS file to examine, must have a catalog extension"
)
args = parser.parse_args()


hdu = fits.open(args.filename)
cat = Table(hdu["cat"].data)
cat.pprint(-1, -1)
