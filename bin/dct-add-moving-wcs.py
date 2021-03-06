#!/usr/bin/env python3
import os
import re
import argparse
import warnings
import numpy as np
from astropy.wcs import WCS, FITSFixedWarning
from astropy.io import fits
from astroquery.simbad import Simbad
import astropy.coordinates as coords
from astropy.coordinates import SkyCoord
import astropy.units as u
#from astroquery.jplhorizons import Horizons
from sbpy.data import Ephem
import mskpy

comet_pat = ('^(([1-9]{1}[0-9]*[PD](-[A-Z]{1,2})?)'
             '|([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?))')
aster_pat = ('^(([1-9][0-9]*( [A-Z]{1,2}([1-9][0-9]{0,2})?)?)'
             '|(\(([1-9][0-9]*)\)))')

######################################################################
# suppress unnecessary warnings
warnings.simplefilter('ignore', FITSFixedWarning)
######################################################################

parser = argparse.ArgumentParser(
    description='Add a WCS centered on a moving target.')
parser.add_argument('files', nargs='+', help='FITS images to update.')
parser.add_argument('--observatory', default='G37',
                    help='Observatory location (for HORIZONS).')
parser.add_argument('--source', choices=['jpl', 'mpc'], default='jpl',
                    help='ephemeris source')

args = parser.parse_args()

# first pass to get target and dates
objects = []
dates = []
for f in args.files:
    hdu = fits.open(f, mode='readonly')
    objects.append(hdu[0].header['OBJECT'])
    dates.append(hdu[0].header['DATE-OBS'])

obj = np.unique(objects)
assert len(obj) == 1, "Multiple objects found: {}".format(obj)
t = mskpy.date2time(dates)

# put files in time order
files = [args.files[i] for i in np.argsort(t)]
t = t[np.argsort(t)]

obj = obj[0]
comet_match = re.findall(comet_pat, obj)
aster_match = re.findall(aster_pat, obj)
m = comet_match + aster_match
assert len(m) != 0, 'Designation is not that of a comet or asteroid: {}'.format(obj)
target = m[0][0]

c = []
print(obj)

if args.source == 'jpl':
    if comet_match:
        opts = dict(closest_apparition=True, no_fragments=True)
    else:
        opts = {}
    eph = Ephem.from_horizons(
        target, id_type='designation', epochs=t, location='G37',
        **opts
    )
elif args.source == 'mpc':
    eph = Ephem.from_mpc(target, epochs=t, location='G37')

for j in range(len(eph)):
    _c = SkyCoord(eph[j]['RA'], eph[j]['DEC'], unit=u.deg)
    print('  ', t[j].iso, _c.to_string('hmsdms')[0])
    c.append(_c)

c = coords.concatenate(c)
k = t.argmin()
for j, f in enumerate(files):
    hdu = fits.open(os.path.realpath(f), mode='update')

    w0 = WCS(hdu[0].header)

    w = WCS(naxis=2)
    w.wcs.crpix = w0.wcs_world2pix(c[j].ra, c[j].dec, 1)
    w.wcs.crval = (c[k].ra.deg, c[k].dec.deg)
    try:
        w.wcs.cd = w0.wcs.cd
    except AttributeError:
        pass
    w.wcs.ctype = w0.wcs.ctype
    hdu[0].header.update(w.to_header(key='M'))

    s = 'Added moving target WCS (key M).'
    for line in hdu[0].header.get('HISTORY', []):
        if s in line:
            break
    else:
        hdu[0].header.add_history(s)

    hdu.close()
