#!/usr/bin/env python3
import re
import argparse
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.simbad import Simbad
import astropy.coordinates as coords
from astropy.coordinates import SkyCoord
import astropy.units as u
import callhorizons
import mskpy

comet_pat = ('^(([1-9]{1}[0-9]*[PD](-[A-Z]{1,2})?)'
             '|([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?))')
aster_pat = ('^(([1-9][0-9]*( [A-Z]{1,2}([1-9][0-9]{0,2})?)?)'
             '|(\(([1-9][0-9]*)\)))')

parser = argparse.ArgumentParser(description='Add a WCS centered on a moving target.')
parser.add_argument('file', nargs='+', help='FITS images to update.')
parser.add_argument('--observatory', default='G37', help='Observatory location (for HORIZONS).')

args = parser.parse_args()
files = sorted(args.file)

# first pass to get target and dates
objects = []
dates = []
for f in files:
    hdu = fits.open(f, mode='readonly')
    objects.append(hdu[0].header['OBJECT'])
    dates.append(hdu[0].header['DATE-OBS'])

obj = np.unique(objects)
assert len(obj) == 1, "Multiple objects found: {}".format(obj)
t = mskpy.date2time(dates)

obj = obj[0]
comet_match = re.findall(comet_pat, obj)
aster_match = re.findall(aster_pat, obj)
m = comet_match + aster_match
assert len(m) != 0, 'Designation is not that of a comet or asteroid: {}'.format(obj)
target = m[0][0]

c = []
print(obj)
for i in range(len(t)):
    q = callhorizons.query(target, comet=len(comet_match) > 0)
    q.set_discreteepochs(t.jd[i])
    assert q.get_ephemerides(args.observatory) == 1, 'Error getting ephemerides.'
    _c = SkyCoord(q['RA'], q['DEC'], unit=u.deg)
    print('  ', t[i].iso, _c.to_string('hmsdms')[0])
    c.append(_c)

c = coords.concatenate(c)
i = t.argmin()
for j, f in enumerate(files):
    hdu = fits.open(f, mode='readonly')

    w0 = WCS(hdu[0].header)

    w = WCS(naxis=2)
    w.wcs.crpix = w0.wcs_world2pix(c[j].ra, c[j].dec, 1)
    w.wcs.crval = (c[i].ra.deg, c[i].dec.deg)
    w.wcs.cd = np.array([[-1, 0], [0, 1]]) * 0.24 / 3600
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hdu[0].header.update(w.to_header(key='M'))
    
    s = 'Added moving target WCS (key M).'
    for line in hdu[0].header.get('HISTORY', []):
        if s in line:
            break
    else:
        hdu[0].header.add_history(s)

    hdu.writeto(f, overwrite=True)
    hdu.close()
