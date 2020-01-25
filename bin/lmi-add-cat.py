#!/usr/bin/env python3
import sys
import logging
import argparse
from datetime import datetime
from glob import glob

import numpy as np
import scipy.ndimage as nd
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.table import Table
import photutils as pu
import sep


######################################################################
# setup logging
logger = logging.Logger('LMI Add Catalog')
logger.setLevel(logging.DEBUG)

# this allows logging to work when lmi-add-cat is run multiple times from
# ipython
if len(logger.handlers) == 0:
    formatter = logging.Formatter('%(levelname)s: %(message)s')

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler('lmi-add-cat.log')
    logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info('#' * 70)
logger.info(datetime.now().isoformat())
logger.info('Command line: ' + ' '.join(sys.argv[1:]))

######################################################################


def show_objects(im, objects):
    # plot background-subtracted image
    fig, ax = plt.subplots()
    m, s = np.mean(im), np.std(im)
    im = ax.imshow(im, interpolation='nearest', cmap='gray',
                   vmin=m-s, vmax=m+s, origin='lower')

    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)


parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='*', help='files to process')
parser.add_argument('--reprocess', action='store_true')
args = parser.parse_args()

for f in args.files:
    logger.info(f)
    if fits.getheader(f)['IMAGETYP'] != 'OBJECT':
        continue

    with fits.open(f, mode='update') as hdu:
        im = hdu[0].data + 0
        h = hdu[0].header
        if 'MASK' in hdu:
            mask = hdu['MASK'].data.astype(bool)
        else:
            mask = False

        if 'cat' in hdu:
            if not args.reprocess:
                continue
            else:
                del hdu['cat']

        shape = im.shape
        center = np.s_[shape[0] - 150:shape[0] + 150,
                       shape[1] - 150:shape[1] + 150]
        mms = sigma_clipped_stats(im[center])
        det = nd.binary_closing(
            (nd.uniform_filter(im, 9) > mms[1] + mms[2] / 9 * 2),
            iterations=10)

        bkg = sep.Background(im, mask=det + mask, bw=64, bh=64, fw=3, fh=3)
        if 'bg' in hdu:
            del hdu['bg']
        hdu.append(fits.ImageHDU(bkg.back(), name='bg'))
        hdu['bg'].header['bg'] = bkg.globalback
        hdu['bg'].header['rms'] = bkg.globalrms

        data = im - bkg
        data[mask] = 0
        objects, labels = sep.extract(data, 3, err=bkg.globalrms,
                                      segmentation_map=True)
        hdu[0].header['ncat'] = len(objects), 'number of objects in catalog'
        if len(objects) == 0:
            continue

        #show_objects(data, objects)

        segmap = pu.SegmentationImage(labels)
        fwhms = []
        for i in np.random.choice(len(segmap.segments), 50):
            obj = segmap.segments[i].make_cutout(data, masked_array=True)
            try:
                g = pu.fit_2dgaussian(obj)
            except:
                continue

            fwhms.append(np.mean((g.x_stddev, g.y_stddev)) * 2.35)

        fwhm = sigma_clipped_stats(fwhms)[1]

        flux, fluxerr, flag = sep.sum_circle(
            data, objects['x'], objects['y'], fwhm * 2, err=bkg.globalrms,
            gain=h['gain'])

        kronrad, krflag = sep.kron_radius(data, objects['x'], objects['y'],
                                          objects['a'], objects['b'],
                                          objects['theta'], 6.0)
        krflux, krfluxerr, _flag = sep.sum_ellipse(
            data, objects['x'], objects['y'], objects['a'], objects['b'],
            np.minimum(objects['theta'], np.pi / 2.00001),
            2.5 * kronrad, subpix=1, err=bkg.globalrms,
            gain=h['gain'])
        krflag |= _flag  # combine flags

        wcs = WCS(h)
        ra, dec = wcs.all_pix2world(objects['x'], objects['y'], 0)

        tab = Table((objects['x'], objects['y'], ra, dec, flux, fluxerr,
                     flag, objects['a'], objects['b'], objects['theta'],
                     kronrad, krflux, krfluxerr, krflag),
                    names=('x', 'y', 'ra', 'dec', 'flux', 'fluxerr', 'flag',
                           'a', 'b', 'theta', 'kronrad', 'krflux',
                           'krfluxerr', 'krflag'))
        if 'cat' in hdu:
            del hdu['cat']
        hdu.append(fits.BinTableHDU(tab, name='cat'))
        hdu['cat'].header['FWHM'] = fwhm, 'estimated median FWHM'
        hdu['cat'].header['RADIUS'] = 2 * fwhm, 'aperture photometry radius'
