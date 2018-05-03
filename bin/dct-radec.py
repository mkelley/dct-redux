#!/usr/bin/env python3
import argparse
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time

parser = argparse.ArgumentParser(description='DCT RA, Dec from file accouting'
                                 ' for equinox.')
parser.add_argument('files', nargs='+', help='FITS images to inspect.')
parser.add_argument('--raunit', default='hourangle',
                    help='angular unit for RA output')
parser.add_argument('--decunit', default='deg',
                    help='angular unit for Dec output')
parser.add_argument('--sep', default=':', help='separator for sexagesimal'
                    ' representation')
parser.add_argument('--precision', help='precision for decimal'
                    ' representation')
parser.add_argument('--alwayssign', action='store_true', help='include the'
                    ' sign no matter what')

args = parser.parse_args()
opts = dict(sep=args.sep, precision=args.precision,
            alwayssign=args.alwayssign)
for f in args.files:
    h = fits.getheader(f)
    equinox = Time(h['EQUINOX'], format='decimalyear')
    c = SkyCoord(h['RA'], h['DEC'], unit=('hourangle', 'deg'),
                 frame='fk5', equinox=equinox)
    ra = c.icrs.ra.to_string(unit=args.raunit, **opts)
    dec = c.icrs.dec.to_string(unit=args.decunit, **opts)
    print(f, ra, dec)

