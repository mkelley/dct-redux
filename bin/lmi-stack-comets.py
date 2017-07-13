#!/usr/bin/env python3
import os
import argparse
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import montage_wrapper as montage

parser = argparse.ArgumentParser(description='Stack comet images together.')
parser.add_argument('file', nargs='+', help='Files to stack by filter.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite previously saved files.')
parser.add_argument('--mask', default='_mask', help='File name suffix for object masks.')
parser.add_argument('--bgmatch', action='store_true', help='Enable background matching.')
parser.add_argument('--no-cleanup', action='store_false', dest='cleanup', help='Disable montage_wrapper cache cleanup.')

args = parser.parse_args()

# collect image and filter names
files = {}
for f in sorted(args.file):
    h = fits.getheader(f)
    if h['FILTERS'] not in files:
        files[h['FILTERS']] = []
    files[h['FILTERS']].append(f)

for filt, ff in files.items():
    print(filt)
    print(ff)

    if os.path.exists('/tmp/lmi-stack-comets/input/'):
        os.system('rm -f /tmp/lmi-stack-comets/input/*')
    else:
        os.system('mkdir -p /tmp/lmi-stack-comets/input')

    for f in ff:
        #os.system('cp {} /tmp/lmi-stack-comets/input'.format(f))
        cf = '/tmp/lmi-stack-comets/input/' + f.split('/')[-1]
        mf = f.replace('.fits', '{}.fits'.format(args.mask))
        with fits.open(f, mode='readonly') as hdu:
            if os.path.exists(mf):
                print('  Masking {}'.format(f))
                mask = fits.getdata(mf) == 1
                hdu[0].data[mask] = np.nan

            hdu.writeto(cf, overwrite=True)

    copied = sorted(glob('/tmp/lmi-stack-comets/input/*fits'))
    exptime = 0
    obsnums = []
    for f in copied:
        with fits.open(f, mode='update') as hdu:
            exptime += hdu[0].header['EXPTIME']
            obsnums.append(hdu[0].header['OBSERNO'])
            w = WCS(hdu[0].header, key='M')
            hdu[0].header.update(w.to_header(key=' '))

        if f == copied[0]:
            montage.mGetHdr(f, '/tmp/lmi-stack-comets/input/header')

    for combine in ['median', 'mean']:
        print('  -', combine, end=' ')
        outf = ff[0][:-5] + '-' + combine
        if os.path.exists(outf) and not args.ovewrite:
            print('exists, skipping.')
            continue
        else:
            print()

        if os.path.exists('/tmp/lmi-stack-comets/output/'):
            os.system('rm -f /tmp/lmi-stack-comets/output/*')
            os.system('rmdir /tmp/lmi-stack-comets/output')

        montage.mosaic('/tmp/lmi-stack-comets/input',
                       '/tmp/lmi-stack-comets/output',
                       header='/tmp/lmi-stack-comets/input/header',
                       combine=combine, bitpix=-32,
                       background_match=args.bgmatch,
                       cleanup=args.cleanup)

        os.system('cp /tmp/lmi-stack-comets/output/mosaic.fits {}.fits'.format(outf))
        os.system('cp /tmp/lmi-stack-comets/output/mosaic_area.fits {}_cov.fits'.format(outf))

        with fits.open(outf + '.fits', mode='update') as hdu:
            hdu[0].header['EXPTIME'] = exptime
            hdu[0].header['MOSAIC'] = combine
            hdu[0].header.add_history('Mosaicked: {}'.format(obsnums))
            hdu[0].header.add_history('EXPTIME updated for this mosaic.')

        print(outf, 'saved')
