#!/usr/bin/env python3
import os
import sys
from glob import glob
import logging
import argparse
from datetime import datetime
import numpy as np
import scipy.ndimage as nd
import astropy.units as u
from astropy.io import fits, ascii
from astropy.table import Table, hstack
import mskpy
from mskpy import apphot
from mskpy.photometry import hb
import json

filter_names = {
    'SDSS-R': 'SDSSr',
    'BC': 'BC',
    'RC': 'RC',
    'OH': 'OH',
    'CN': 'CN',
    'V': 'V'
}
solar_gmr = 0.98 * 0.649 - 0.19  # Smith et al. 2002 and Colina et
                                 # al. 1996 / Bessell et al. 1998
DCT_elevation = 2361 * u.m  # elevation of DCT

parser = argparse.ArgumentParser(description='Moving target photometry.', epilog="""Configuration file format:
{
  "filters": [list, of, filters, to, consider],
  "calN": {
    "filter": [[zp, Ex, color correction (optional)], overall uncertainty],
    "OH": [[zp, toz], unc],
    "uncalibrated filter": [[0, 0], 0]
  },
  "targets": {
    "target set": {
      "observation numbers": [image, numbers, to, include, [may be, range]],
      "ignore": [file, numbers, to, skip],
      "rap": [list, of, aperture, radii, in, pixels],
      "calset": "name of calibration set to use, e.g., cal0",
      "horiz destripe": [list, of, filters, to, remove, horizontal, striping],
      "solar color": true to assume a solar color for SDSSr color correction,
      "centroid": true to centroid on the target
    }
  }
}
""", formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('files', nargs='+', help='Input data.')
parser.add_argument('--config', default='comet-phot.config', help='Configuration file.')
parser.add_argument('-o', default='comet-phot.txt', help='Output file.')
parser.add_argument('--append', action='store_true', help='Append to output file.')
parser.add_argument('--no-fixpix', action='store_false', dest='fixpix', help='Do not replace NaNs with nearby median.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite output file, if it exists.')
args = parser.parse_args()

assert not (args.append and args.overwrite), 'Only one of --append or --overwrite may be specified.'

######################################################################
# setup logging
logger = logging.Logger('LMI comet photometry')
logger.setLevel(logging.DEBUG)

# this allows logging to work when lmi-rx is run multiple times from
# ipython
if len(logger.handlers) == 0:
    formatter = logging.Formatter('%(levelname)s: %(message)s')

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler('lmi-comet-phot.log')
    logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info('#' * 70)
logger.info(datetime.now().isoformat())
logger.info('Command line: ' + ' '.join(sys.argv[1:]))

########################################################################
def correction(header, calset, solar_color=False):
    from mskpy.photometry import airmass_app, hb

    filt = header['filters']
    am = airmass_app(header['ZA'] * u.deg, DCT_elevation)
    ext = calset[filt][0]
    m0_unc = calset[filt][1]
    if filt == 'OH':
        m0 = 0
        m0_unc = 0
    else:
        m0 = ext[0] - ext[1] * am  # mag
        if solar_color and len(ext) == 3 and filt == 'SDSS-R':
            m0 -= ext[2] * solar_gmr

    return am, m0, m0_unc

class IgnoreObservation(Exception):
    pass

def get_parameters(config, obserno):
    for k, p in config['targets'].items():
        for n in p['observation numbers']:
            if n in p.get('ignore', []):
                raise IgnoreObservation
            if isinstance(n, list):
                if n[0] <= obserno and obserno <= n[1]:
                    return k, p
            else:
                if obserno == n:
                    return k, p

    raise ValueError('Observation number {} not found in config.'.format(obserno))

########################################################################
with open(args.config, 'r') as inf:
    config = json.load(inf)
logger.info('Using configuration:\n{}'.format(str(config)))

########################################################################
phot = Table(names=('file', 'target', 'filter', 'mosaic', 'date', 'time',
                    'airmass', 'ZA', 'rap', 'flux', 'f unc', 'm', 'm unc'),
             dtype=['U64', 'U64', 'U6', 'U6', 'U10', 'U8'] + [float] * 7)
for col in ['airmass', 'flux', 'f unc']:
    phot[col].format = '{:.3f}'
for col in ['m', 'm unc']:
    phot[col].format = '{:.4f}'
for col in ['rap']:
    phot[col].format = '{:.1f}'

########################################################################
# photometry in magnitudes
for f in args.files:
    hdu = fits.open(f, mode='readonly')

    logger.info(f)

    # to which target set does this belong?
    obserno = hdu[0].header['OBSERNO']
    try:
        target_set, p = get_parameters(config, obserno)
    except IgnoreObservation:
        logger.info('{} in ignore list, skipping.'.format(f))
        continue

    filt = hdu[0].header['filters']
    am, m0, m0_unc = correction(hdu[0].header, config[p['calset']],
                                p.get('solar color', False))
    hdu[0].header['AIRMASS'] = am, "Hardie 1962, with refraction corr."
    hdu[0].header['M0'] = m0, 'Magnitude zero-point, incl. extinction'
    hdu[0].header['M0UNC'] = m0_unc, 'Uncertainty on M0'
    hdu[0].header['F0'] = hb.F_0[filt].value, 'Magnitude zero-point, W/m2/um'

    cyx = hdu[0].header['CRPIX2M'], hdu[0].header['CRPIX1M']
    bgarea = hdu[0].header['bgarea']
    bgsig = hdu[0].header['bgsig']
    mosaic = hdu[0].header.get('MOSAIC', '')

    im = np.ma.MaskedArray(hdu[0].data)
    im.mask = ~np.isfinite(im)
    if args.fixpix and np.any(im.mask):
        m = nd.median_filter(im, 11)
        im[im.mask] = m
        im.mask = False
    
    area, flux = apphot(hdu[0].data, cyx, p['rap'], subsample=2)
    if area.ndim == 0:
        area = area.reshape((1,))
        flux = flux.reshape((1,))

    var = area * bgsig**2 * (1 + area / bgarea)
    i = flux > 0
    if any(i):
        var[i] += flux[i] / hdu[0].header['GAIN']

    func = np.sqrt(var)

    m = np.ones_like(flux) * np.nan
    m[i] = -2.5 * np.log10(flux[i]) + m0
    munc = np.sqrt((1.0857 * func / np.abs(flux))**2 + m0_unc**2)
    
    for i in range(len(p['rap'])):
        row = [f, hdu[0].header['OBJECT'], filter_names[filt], mosaic]
        row.extend(hdu[0].header['DATE-OBS'].split('T'))
        row.extend([am, hdu[0].header['ZA'], p['rap'][i]])
        row.extend([flux[i], func[i]])
        row.extend([m[i], munc[i]])
        phot.add_row(row)

phot.sort(['target', 'filter', 'rap', 'date', 'mosaic', 'time'])
phot.pprint(max_lines=-1, max_width=-1)

if os.path.exists(args.o) and not (args.overwrite or args.append):
    logger.warning('{} exists, not saving table'.format(args.o))
else:
    if args.append:
        saved = ascii.read(args.o)
        phot = hstack((saved, phot))
    
    mskpy.write_table(args.o, phot, {},
                      comments=['OH photometry is in instrumental magnitudes.',
                                'Units: UT, deg, pixels, DN/s, mag'],
                      overwrite=True)

