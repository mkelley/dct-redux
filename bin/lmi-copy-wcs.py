#!/usr/bin/env python3
import re
import argparse
import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astroquery.jplhorizons import Horizons
import mskpy

comet_pat = ('^(([1-9]{1}[0-9]*[PD](-[A-Z]{1,2})?)'
             '|([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?))')
aster_pat = ('^(([1-9][0-9]*( [A-Z]{1,2}([1-9][0-9]{0,2})?)?)'
             '|(\(([1-9][0-9]*)\)))')

parser = argparse.ArgumentParser(description='Copy LMI WCS from one filter to another.',
                                 epilog='Determines which files to pair up using file name sort order.')
parser.add_argument(
    'from_filter', help='Images with this filter have complete WCS.')
parser.add_argument(
    'to_filter', help='Images with this filter will receive a new WCS.')
parser.add_argument('file', nargs='*', help='FITS images to consider.')
parser.add_argument(
    '--offset-filter', help='Compute offset between from_filter and to_filter using distance between from_ and offset_filter and --offset_scale.')
parser.add_argument('--offset-scale', default=2.4,
                    help='Scale factor for offsets, default value works for RC,BC to OH on 20160212.')
parser.add_argument('-n', action='store_true',
                    help='No-operation mode.  No files are changed.')

args = parser.parse_args()

from_files = []
to_files = []
offset_files = []
for f in sorted(args.file):
    print(f)
    with fits.open(f, mode='readonly') as hdu:
        if hdu[0].header['FILTERS'] == args.from_filter:
            from_files.append(f)
        if hdu[0].header['FILTERS'] == args.to_filter:
            to_files.append(f)
        if hdu[0].header['FILTERS'] == args.offset_filter:
            offset_files.append(f)

print()

assert len(from_files) == len(
    to_files), "File lists have different lengths: from, to."
if args.offset_filter is None:
    offset_files = [None] * len(from_files)
else:
    assert len(from_files) == len(
        offset_files), "File lists have different lengths: from, offset."

for ff, tf, of in zip(from_files, to_files, offset_files):
    print(ff, '->', tf)
    f_hdu = fits.open(ff, mode='readonly')
    t_hdu = fits.open(tf, mode='readonly')
    if of is not None:
        o_hdu = fits.open(of, mode='readonly')

    obj = f_hdu[0].header['OBJECT']
    assert obj == t_hdu[0].header['OBJECT'], 'Objects do not match: '.format(
        obj, t_hdu[0].header['OBJECT'])
    m = re.findall(comet_pat, obj)
    comet = m[0][0]
    if len(m) != 0:
        # comet or asteroid
        comet_flag = re.findall(comet_pat, obj) is not None
        moving = True
        d = mskpy.date2time(f_hdu[0].header['DATE'])
        q = Horizons(id=comet, id_type='designation', epochs=d.jd,
                     location='G37')
        eph = q.ephemerides(closest_apparition=True, no_fragments=True)
        c_f = SkyCoord(eph['RA'][0], eph['DEC'][0], unit=u.deg)
        print(d.iso, c_f.to_string('hmsdms'))

        d = mskpy.date2time(t_hdu[0].header['DATE'])
        q = Horizons(id=comet, id_type='designation', epochs=d.jd,
                     location='G37')
        if comet:
            opts = dict(closest_apparition=True, no_fragments=True)
        else:
            opts = {}
        eph = q.ephemerides(**opts)
        c_t = SkyCoord(eph['RA'][0], eph['DEC'][0], unit=u.deg)
        print('  ->', d.iso, c_t.to_string('hmsdms'))

        dc = c_t.spherical_offsets_to(c_f)
    else:
        moving = False

    if 'EQUINOX' in t_hdu[0].header:
        t_hdu[0].header.add_history('Removed EQUINOX keyword, value: {}'.format(
            t_hdu[0].header['EQUINOX']))
    for k in ['EQUINOX', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1', 'RADECSYS']:
        if k in t_hdu[0].header:
            del t_hdu[0].header[k]

    if of is None:
        linear_offset = np.r_[0, 0]
    else:
        fxy = f_hdu[0].header['CRPIX1'], f_hdu[0].header['CRPIX2']
        oxy = o_hdu[0].header['CRPIX1'], o_hdu[0].header['CRPIX2']
        linear_offset = args.offset_scale * \
            np.r_[fxy[0] - oxy[0], fxy[1] - oxy[1]]
        print(linear_offset)

    w0 = WCS(f_hdu[0].header)
    w0.wcs.crpix = w0.wcs.crpix + linear_offset
    if moving:
        c = SkyCoord(w0.wcs.crval[0] + dc[0].deg, w0.wcs.crval[1] + dc[1].deg,
                     unit=u.deg)
        w0.wcs.crpix = skycoord_to_pixel(c, w0, 1)

        wm = WCS(f_hdu[0].header, key='M')
        wm.wcs.crpix = skycoord_to_pixel(c_t, w0, 1)
        t_hdu[0].header.update(wm.to_header(key='M'))

    import pdb
    pdb.set_trace()

    t_hdu[0].header.update(w0.to_header())
    t_hdu[0].header.add_history('Updated WCS based on {}'.format(ff))

    if not args.n:
        t_hdu.writeto(tf, overwrite=True)

    f_hdu.close()
    t_hdu.close()
