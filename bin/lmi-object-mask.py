#!/usr/bin/env python3
import os
import argparse
import numpy as np
import scipy.ndimage as nd
from astropy.io import fits
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from mskpy import rarray

parser = argparse.ArgumentParser(description='Find objects and generate a mask.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files', nargs='+', help='Files to process.')
parser.add_argument('--center-keys', default='CX,CY', help='FITS header keywords to use for target x and y centers.')
parser.add_argument('--fwhm', type=float, default=6, help='FWHM for detection algorithm.')
parser.add_argument('--thresh', type=float, default=5, help='Number of background standard deviations to use as a detection threshold.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite previous files.')

args = parser.parse_args()

for f in args.files:
    print(f)
    outf = f.replace('.fits', '_mask.fits')
    if os.path.exists(outf) and not args.overwrite:
        print('  Mask exists.')
        
    with fits.open(f, mode='readonly') as hdu:
        r = rarray(hdu[0].data.shape)
        bg = nd.grey_opening(hdu[0].data, int(args.fwhm) * 3)
        im = hdu[0].data - bg
        im[r > 1500] = 0

        print('  Mean, median, stdev of background.')
        mean, median, std = sigma_clipped_stats(im[r < 1500], sigma=3, iters=5)
        im -= median
        
        print('  Find sources.')
        daofind = DAOStarFinder(fwhm=args.fwhm, threshold=args.thresh * std)
        sources = daofind(im)

        mask = np.zeros_like(im, bool)
        for source in sources:
            x = int(np.round(source['xcentroid']))
            y = int(np.round(source['ycentroid']))
            mask[y, x] = True

        # find nearest neighbors above threshold
        print('  Generate object mask.')
        det = nd.gaussian_filter(im, args.fwhm) > std

        # remove comet from detection map
        cyx = [hdu[0].header[k] for k in args.center_keys.split(',')[::-1]]
        cyx = np.round(cyx).astype(int)
        labels, n = nd.label(det)
        det[labels == labels[cyx[0], cyx[1]]] = False
        
        n = 1
        iters = 0
        while n > 0 and iters < 100:
            n = mask.sum()
            mask = nd.binary_dilation(mask, iterations=3)
            mask *= det
            n = mask.sum() - n
            iters += 1

        hdu[0].header['OBJMASK.FWHM'] = args.fwhm
        hdu[0].header['OBJMASK.THRESH'] = args.thresh
        fits.writeto(outf, mask * 1, hdu[0].header, overwrite=args.overwrite)

