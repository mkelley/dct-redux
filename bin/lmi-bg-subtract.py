#!/usr/bin/python3
import os
import sys
from glob import glob
import logging
import argparse
from datetime import datetime
import numpy as np
import astropy.units as u
from astropy.io import fits, ascii
from astropy.table import Table, hstack
import mskpy
import json

filter_names = {
    'SDSS-R': 'SDSSr',
    'BC': 'BC',
    'RC': 'RC',
    'OH': 'OH',
    'CN': 'CN',
    'V': 'V'
}

parser = argparse.ArgumentParser(description='Background subtraction and exposure time normalization.', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('files', nargs='+', help='Data to background subtract.')
parser.add_argument('target', help='Path for background subtracted data.')
parser.add_argument('--config', default='comet-phot.config', help='Configuration file.  Used to determine which files to consider, which files to ignore, and if additional filtering is needed (e.g., destripe).')
parser.add_argument('--sections', default='bg-sections.txt', help='Name of the background section definition file.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite exsiting files.')
args = parser.parse_args()

######################################################################
# setup logging
logger = logging.Logger('LMI background subtract.')
logger.setLevel(logging.DEBUG)

# this allows logging to work when lmi-rx is run multiple times from
# ipython
if len(logger.handlers) == 0:
    formatter = logging.Formatter('%(levelname)s: %(message)s')

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler('lmi-bg-subtract.log')
    logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info('#' * 70)
logger.info(datetime.now().isoformat())
logger.info('Command line: ' + ' '.join(sys.argv[1:]))

########################################################################
class BackgroundError(Exception):
    pass

def measure_bg(im, yx, section):
    from mskpy import bgphot, meanclip

    if 'r:' in section:
        annulus = np.array(section[3:].split(','), float)
        n, bg, sig = bgphot(im, yx, annulus)
    else:
        j = section.find('[')
        s = np.zeros_like(im, bool)
        for a in section[j:].split('+'):
            s[eval('np.s_' + a.strip())] = True
        mc = meanclip(im[s], full_output=True)
        n = len(mc[2])
        bg = mc[0]
        sig = mc[1]
        if not np.isfinite(bg):
            raise BackgroundError(section)

    return float(bg), int(n), float(sig)

def horiz_destripe(hdu, p):
    import scipy.ndimage as nd
    from astropy.convolution import convolve_fft, Gaussian2DKernel
    from mskpy import imstat

    mask = hdu['mask'].data.astype(bool)
    r = mskpy.rarray(hdu[0].data.shape, cyx)
    mask += r < 2 * max(p['rap'])

    c = convolve_fft(hdu[0].data * ~mask, Gaussian2DKernel(2.5))
    mm = nd.morphological_gradient(c, 5)
    stat = imstat(mm[~mask])
    thresh = stat['scmedian'] + 5 * stat['scstdev']
    mask += nd.binary_closing(mm > thresh, iterations=3)

    masked = hdu[0].data.copy()
    masked[mask] = np.nan
    stripes = np.outer(np.nanmean(masked, 1), np.ones(hdu[0].data.shape[1]))
    return stripes

########################################################################
def add_bg_to_header(header, bg, area, sig, sect):    
    header.add_history('Background estimate based on a sigma-clipped mean.')
    if 'backgrnd' in header:
        # background already subtracted at least once, add to backgrnd keyword
        header['backgrnd'] = header['backgrnd'] + bg
    else:
        header['backgrnd'] = bg, 'Background estimate, {} per pixel.'.format(header['bunit'])
        header['bgarea'] = area, 'Number of pixels used for background estimate.'
        header['bgsig'] = sig, 'Background standard deviation.'
        header['bgsect'] = sect, 'Background area definition.'

    return header

########################################################################
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
sections = ascii.read(args.sections)
for i in range(len(sections)):
    sections[i]['file'] = os.path.basename(sections[i]['file'])

os.system('mkdir -p {}'.format(args.target))

with open(args.config, 'r') as inf:
    config = json.load(inf)
logger.info('Using configuration:\n{}'.format(str(config)))

########################################################################
# Clean the data...
for f in args.files:
    hdu = fits.open(f, mode='readonly')

    # to which target set does this belong?
    obserno = hdu[0].header['OBSERNO']
    try:
        target_set, p = get_parameters(config, obserno)
    except IgnoreObservation:
        logger.info('{} in ignore list, skipping.'.format(f))
        continue

    date, time = hdu[0].header['DATE-OBS'].split('T')
    date = date.replace('-', '')
    time = time[:8].replace(':', '')

    filt = hdu[0].header['filters']
    
    file_parts = [target_set, date, time, filter_names[filt]]
    if hdu[0].header.get('MOSAIC') is None:
        mosaic = False
    else:
        mosaic = True
        file_parts.append(hdu[0].header.get('MOSAIC'))

    outfn = '{}/{}.fits'.format(args.target, '-'.join(file_parts))

    if os.path.exists(outfn) and not args.overwrite:
        logger.info('{} exists, skipping.'.format(f))
        continue

    logger.info(f)
    # if TARGTSET is already present, then we must be working with a
    # mosaicked image and can skip the following
    if 'TARGTSET' not in hdu[0].header:
        hdu[0].header['TARGTSET'] = target_set
        hdu[0].data /= hdu[0].header['EXPTIME']
        hdu[0].header['BUNIT'] = 'DN/s'
        hdu[0].header['GAIN'] = hdu[0].header['GAIN'] * hdu[0].header['EXPTIME'], 'e-/(DN/s)'

    cyx = hdu[0].header['CRPIX2M'], hdu[0].header['CRPIX1M']
    i = np.flatnonzero(sections['file'] == os.path.basename(f))[0]
    sect = sections[i]['sections']

    bg, bgarea, bgsig = measure_bg(hdu[0].data, cyx, sect)
    hdu[0].header = add_bg_to_header(hdu[0].header, bg, bgarea, bgsig, sect)

    hdu[0].data -= bg

    if not mosaic and len(p.get('horiz destripe', [])) > 0:
        if filt in p['horiz destripe']:
            stripes = horiz_destripe(hdu, p)
            hdu[0].data -= stripes
            msg = 'Horizontal destripe, mean {:.3f} DN/s.'.format(
                np.nanmean(stripes))
            logger.info('  ' + msg)
            hdu[0].header.add_history(msg)

    hdu.writeto(outfn, overwrite=args.overwrite)
