#!/usr/bin/env python3
import re
import argparse
from itertools import chain
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import callhorizons
import mskpy

comet_pat = ('^(([1-9]{1}[0-9]*[PD](-[A-Z]{1,2})?)'
             '|([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?))')
aster_pat = ('^(([1-9][0-9]*( [A-Z]{1,2}([1-9][0-9]{0,2})?)?)'
             '|(\(([1-9][0-9]*)\)))')

parser = argparse.ArgumentParser(description='Update LMI WCS using already centered targets.')
parser.add_argument('file', nargs='+', help='FITS images to update.')
parser.add_argument('--celestial', default=None, help='The aligned object is a celestial source with this "RA,Dec".  Celestial coordinates can be in any format usable by astropy.Angle (include units).')
parser.add_argument('--observatory', default='G37', help='Observatory location (for HORIZONS).')

args = parser.parse_args()
saved_coords = {}

for f in args.file:
    hdu = fits.open(f, mode='readonly')
    if 'CX' not in hdu[0].header or 'CY' not in hdu[0].header:
        print(f, ' not yet centered')
        continue

    obj = hdu[0].header['OBJECT']
    cxy = hdu[0].header['CX'], hdu[0].header['CY']

    if args.celestial is None:
        m = re.findall(comet_pat, obj) + re.findall(aster_pat, obj)
        if len(m) == 0:
            # not a comet
            i = obj.find('/')
            if i > 0:
                obj = obj[:i].strip()

            i = obj.find('(')
            if i > 0:
                obj = obj[:i].strip()

            if obj in saved_coords:
                c = saved_coords[obj]
            else:
                q = Simbad.query_object(obj)
                if q is None:
                    print(f, obj, ' not found in Simbad.')
                    continue
                c = SkyCoord(q['RA'][0], q['DEC'][0], unit=(u.hourangle, u.deg))
                saved_coords[obj] = c
        else:
            # comet or asteroid
            q = callhorizons.query(m[0][0])
            d = mskpy.date2time(hdu[0].header['DATE-OBS'])
            q.set_discreteepochs([d.jd])
            q.get_ephemerides(args.observatory)
            c = SkyCoord(q['RA'][0], q['DEC'][0], unit=u.deg)
    else:
        ra, dec = [Angle(a) for a in args.celestial.split(',')]
        c = SkyCoord(ra, dec)
        obj = 'celestial source'

    if 'EQUINOX' in hdu[0].header:
        hdu[0].header.add_history('Removed EQUINOX keyword, value: {}'.format(
            hdu[0].header['EQUINOX']))
    for k in ['EQUINOX', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1', 'RADECSYS']:
        if k in hdu[0].header:
            del hdu[0].header[k]

    w = WCS(naxis=2)
    w.wcs.crpix = cxy
    w.wcs.cd = np.array([[-1, 0], [0, 1]]) * 0.24 / 3600
    w.wcs.crval = (c.ra.deg, c.dec.deg)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hdu[0].header.update(w.to_header())

    s = 'Aligned on {} at {}'.format(obj, c.to_string(style='hmsdms'))
    print(s)
    hdu[0].header.add_history(s)

    hdu.writeto(f, overwrite=True)
    hdu.close()
    
