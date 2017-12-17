#!/usr/bin/env python3
import os
import re
import sys
import shutil
import argparse
import logging
from datetime import datetime
import numpy as np
from astropy import units as u
import ccdproc
from ccdproc import CCDData, ImageFileCollection, combine

class Config:
    file_template = 'lmi_[0-9]{8}_[0-9]{4}_raw.fits'

def mode_scaler(im, s=np.s_[660:2600, 660:2600]):
    return 1 / (3 * np.ma.median(im[s]) - 2 * np.ma.mean(im[s]))

def translate_slice(s):
    """Translate a FITS slice (string) to Python slice."""
    m = re.match('\[(\d+):(\d+),(\d+):(\d+)\]', s.replace(' ', ''))
    i = [int(x) for x in m.groups()]
    return np.s_[i[2]-1:i[3], i[0]-1:i[1]]

def subarray_slice(s, ccdsum):
    """Define a Python slice from LMI subframe keywords."""
    s1 = translate_slice(s)
    xb, yb = [int(x) for x in ccdsum.split()]
    y0 = s1[0].start // yb
    y1 = s1[0].stop // yb
    x0 = s1[1].start // xb
    x1 = s1[1].stop // xb
    return np.s_[y0:y1, x0:x1]

class IntListAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        v = [int(x) for x in values.split(',')]
        setattr(namespace, self.dest, v)

parser = argparse.ArgumentParser(description='Generate partially processed LMI data.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_list', help='A list of files to process.  They must all have file names of the form: lmi_YYYYMMDD_NNNN_raw.fits.  Blank lines and lines beginning with # are ignored.')
parser.add_argument('--target', default='ppp', help='Save processed files to this target directory.')
parser.add_argument('--flat-keys', default='SKY FLAT,DOME FLAT', type=lambda s: s.split(','), help='Object type for flat fields, may be a comma-separated list.')
parser.add_argument('--no-lacosmic', dest='lacosmic', action='store_false', help='Do not run L.A.Cosmic.')
parser.add_argument('--reprocess-data', action='store_true', help='Reprocess data files.')
parser.add_argument('--reprocess-all', action='store_true', help='Recreate bias, flat, and all data files.')

args = parser.parse_args()

if args.reprocess_all:
    args.reprocess_data = True

######################################################################
# setup logging
logger = logging.Logger('LMI Reduction')
logger.setLevel(logging.DEBUG)

# this allows logging to work when lmi-rx is run multiple times from
# ipython
if len(logger.handlers) == 0:
    formatter = logging.Formatter('%(levelname)s: %(message)s')

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler('lmi-rx.log')
    logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info('#' * 70)
logger.info(datetime.now().isoformat())
logger.info('Command line: ' + ' '.join(sys.argv[1:]))

######################################################################
# copy files
if not os.path.exists(args.target):
    os.mkdir(args.target)
    logger.info('Created directory {}.'.format(args.target))

files = []
file_list = np.loadtxt(args.file_list, dtype='S')
if file_list.ndim == 0:
    file_list = file_list.reshape((1, ))

for fn in file_list:
    fn = fn.decode()
    if re.fullmatch(Config.file_template, os.path.split(fn)[1]) is None:
        logger.error('{} does not match template, lmi-YYYYMMDD-NNNN-raw.fits.'.format(fn))
        continue

    fn_new = os.path.split(fn)[1].replace('raw', 'ppp')
    files.append(fn_new)

    fn_new = os.sep.join((args.target, fn_new))
    if os.path.exists(fn_new) and not args.reprocess_data:
        continue
    
    shutil.copy(fn, fn_new)
    os.chmod(fn_new, 0o644)

    logger.debug('{} -> {}'.format(fn, fn_new))

######################################################################
# load files
location = args.target
keywords = ['obserno', 'object', 'obstype', 'date-obs', 'subarser',
            'ra', 'dec', 'airmass', 'filters', 'tempamb', 'humidity',
            'subbias', 'flatcor', 'ccdsum', 'fmdstat']
ic = ImageFileCollection(location=location, filenames=files,
                         keywords=keywords)
nbias = len(ic.files_filtered(obstype='BIAS'))
nflat = 0
for k in args.flat_keys:
    nflat += len(ic.files_filtered(obstype=k))
ndata = len(ic.files_filtered(obstype='OBJECT'))

filters = []
for h in ic.headers():
    if h['OBSTYPE'] == 'OBJECT' or h['OBSTYPE'] in args.flat_keys:
        if h['FILTERS'] not in filters:
            filters.append(h['FILTERS'])

filters.sort()

flat_breakdown = []
data_breakdown = []
for filt in filters:
    for k in args.flat_keys:
        flat_breakdown.append('{} {} {}'.format(
            len(ic.files_filtered(obstype=k, filters=filt)), filt, k))
    data_breakdown.append('{} {}'.format(
        len(ic.files_filtered(obstype='OBJECT', filters=filt)), filt))

subframes = sorted(list(set(ic.summary['subarser'])))
nsubframes = len(subframes)
subframe_breakdown = []
for i in subframes:
    subframe_breakdown.append('{} subframe {}'.format(
        len(ic.files_filtered(subarser=i)), i))

logger.info('''
  {} files:
    {} bias
    {} flats
      {}
    {} object
      {}
  {} subframes
    {}
'''.format(len(ic.files), nbias,
           nflat, '\n      '.join(flat_breakdown),
           ndata, '\n      '.join(data_breakdown),
           nsubframes, '\n    '.join(subframe_breakdown)))

######################################################################
# bias subtract
bias = dict()
logger.info('Bias frames.')
for subframe in subframes:
    fn = '{}/bias{}.fits'.format(args.target, subframe)
    if os.path.exists(fn) and not args.reprocess_all:
        logger.info('  Read subframe = {}.'.format(subframe))
        bias[subframe] = CCDData.read(fn)
    elif nbias == 0:
        logger.warning('  No bias files provided and {} not found.  Not subtracting bias for subframe == {}.'.format(fn, subframe))
        bias[subframe] = 0 * u.adu
    else:
        logger.info('  Create subframe = {}.'.format(subframe))
        files = ic.files_filtered(obstype='BIAS', subarser=subframe)
        files = [os.sep.join([ic.location, f]) for f in files]
        bias[subframe] = combine(files, method='average',
                                 clip_extrema=True, nlow=1, nhigh=1)
        bias[subframe].meta['FILENAME'] = os.path.split(fn)[1]

        n = str([int(f.split('_')[2]) for f in files])
        bias[subframe].meta.add_history(
            'Created from file numbers: {}'.format(n))

        bias[subframe].write(fn, overwrite=True)

for subframe in subframes:
    logger.info('Bias subtract and trim data.')
    i = (ic.summary['subarser'] == subframe) & ic.summary['subbias'].mask
    logger.info('  {} files to bias subtract for subframe {}'.format(
        sum(i), subframe))
    for fn in ic.summary['file'][i]:
        ccd = ccdproc.fits_ccddata_reader(os.sep.join([ic.location, fn]))
        logger.debug(fn)

        ccd = ccdproc.subtract_bias(ccd, bias[subframe])
        ccd.meta['BIASFILE'] = (bias[subframe].meta['FILENAME'],
                                'Name of the bias frame used.')
        biassec = translate_slice(ccd.meta['BIASSEC'])[1]
        ccd = ccdproc.subtract_overscan(ccd, ccd[:, biassec])

        trimsec = ccd.meta['TRIMSEC']
        ccd = ccdproc.trim_image(ccd, fits_section=trimsec)

        ccd.write(os.sep.join([ic.location, fn]), overwrite=True)

ic.refresh()

######################################################################
# flat field and correction
def style2key(style):
    """Generate a flat field key based on FITS keywords.

    filter
    Binning key: 2x2 for 2x2 binning
    NIHTS dichroic key: +D for dichroic in.

    """

    k = '{}-{}{}D'.format(
        style['filters'],
        style['CCDSUM'].strip().replace(' ', 'x'),
        '+' if style['FMDSTAT'] == 'EXTENDED' else '-'
    )

    return k

def needed_flat_styles(ic):
    """Generate a list of needed flats for this data set.

    Requires FILTERS, CCDSUM, and FMDSTAT in the image collection
    summary.

    """
    data_styles = []
    for obs in ic.summary:
        style = (obs['filters'], obs['ccdsum'], obs['fmdstat'])
        if style not in data_styles:
            data_style.append(combo)
            yield {'filters': style[0],
                   'ccdsum': style[1],
                   'fmdstat': style[2]}

logger.info('Flat fields.')
for style in needed_flat_styles(ic):
    k = style2key(style)
    flats[k] = 1
    for flat_key in args.flat_keys:
        fn = '{}/{}-{}.fits'.format(ic.location,
                                    flat_key.lower().replace(' ', ''),
                                    style2key(style))

        if os.path.exists(fn) and not args.reprocess_all:
            logger.info('  Reading {}.'.format(fn))
            flats[k] = CCDData.read(fn)
        elif len(ic.files_filtered(obstype=flat_key, **style)) == 0:
            logger.warning('No {} files provided for {} and {} not found.'.format(flat_key, filt, fn))
        else:
            logger.info('  Generating {}.'.format(fn))
            files = ic.files_filtered(obstype=flat_key, subarser=0, **style)
            files = [os.sep.join([ic.location, f]) for f in files]
            flat = combine(files, method='median', scale=mode_scaler)

            flat.mask = (flat.data > 1.2) + (flat.data < 0.8)

            n = str([int(f.split('_')[2]) for f in files])
            flat.meta.add_history('Created from file numbers: {}'.format(n))
            flat.meta['FILENAME'] = os.path.split(fn)[1]

            flat.write(fn, overwrite=True)

            if flats[k] == 1:
                logger.info('    Using {} for {}'.format(flat_key, k))
                flats[k] = flat

ic.refresh()
i = (ic.summary['obstype'] != 'BIAS') & ic.summary['flatcor'].mask & ~ic.summary['subbias'].mask
logger.info('  {} files to flat correct.'.format(sum(i)))
for fn in ic.summary['file'][i]:
    ccd = ccdproc.fits_ccddata_reader(os.sep.join([ic.location, fn]))

    flatkey = style2key(ccd.meta)
    if flats[flatkey] == 1:
        logger.info('{} skipped, no {} flatfield provided.'.format(fn, flatkey))
        continue

    logger.debug(fn)
    subframe = ccd.meta['SUBARSER']
    if subframe == 0:
        s = np.s_[:, :]
        t = np.s_[:, :]
    else:
        s = subarray_slice(ccd.meta['SUBARR{:02}'.format(subframe)],
                           ccd.meta['CCDSUM'])
        t = subarray_slice(ccd.meta['TRIMSEC'], '1 1')  # flat was not trimmed
    ccd = ccdproc.flat_correct(ccd, flats[flatkey][s][t])
    
    if args.lacosmic:
        cleaned = ccdproc.cosmicray_lacosmic(
            ccd, pssl=1100, gain=ccd.meta['gain'], readnoise=6)
        # LA cosmic converts to e- and adds the bias.  Revert back.
        ccd.data = cleaned.data / ccd.meta['gain'] - 1100
        ccd.mask += cleaned.mask
        ccd.header['LACOSMIC'] = 1, 'L.A.Cosmic processing flag.'
    else:
        ccd.header['LACOSMIC'] = 0, 'L.A.Cosmic processing flag.'

    ccd.meta['FLATFILE'] = (flats[flatkey].meta['FILENAME'],
                            'Name of the flat field correction used.')
    ccd.write(os.sep.join([ic.location, fn]), overwrite=True)
