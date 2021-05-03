#!/usr/bin/env python3
from glob import glob
import os
from astropy.io import fits, ascii
import argparse

allowed_keywords = ['OBSERVER', 'OBSAFFIL', 'OBJNAME', 'OBJECT',
                    'OBSTYPE', 'IMAGETYP']

parser = argparse.ArgumentParser(description='Modify FITS headers.',
                                 epilog="""File format for fixes.txt:

first frame, last frame, keyword, value
11, 26,	object|objname, 96 Her
27, 42, object|objname, HD 219188
43, 56, object|objname, C/2014 Q2 (Lovejoy)
57, 70, object|objname, 53P/Van Biesbroeck

where keyword is a "|" separated list of keywords to verify or change.

Allowed keywords:
  {}
""".format('\n  '.join(allowed_keywords)),
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('fixesfile', type=str, action='store',
                    help='Name of parameter file.')
parser.add_argument('-n', action='store_true',
                    help='No-operation mode: summarizes changes that would be made.')
args = parser.parse_args()

if not os.path.exists('original-headers'):
    if args.n:
        print('I would have created the original-headers directory to backup files.')
    else:
        os.system('mkdir original-headers')

if args.n:
    mode = 'readonly'
else:
    mode = 'update'

files_edited = 0

for row in ascii.read(args.fixesfile, format='csv', delimiter=',', fast_reader=False):
    for i in range(row['first frame'], row['last frame'] + 1):
        file_tests = ['lmi.{:04d}.fits'.format(i),
                      'lmi_????????_{:04d}_raw.fits'.format(i),
                      'lmi_????????_{:04d}_ppp.fits'.format(i),
                      '2???????.{:04d}.fits'.format(i),
                      'deveny_2???????_{:04d}_raw.fits'.format(i),
                      'deveny_2???????_{:04d}_ppp.fits'.format(i)]
        for fn in file_tests:
            try:
                fn = glob(fn)[0]
                h = fits.getheader(fn)
                break
            except IndexError:
                continue
        else:
            raise FileNotFoundError("Tried: {}".format(', '.join(file_tests)))

        test = []
        for k in row['keyword'].split('|'):
            test.append(h[k] != row['value'])

        if not any(test):
            continue
        else:
            backup = 'original-headers/{}'.format(fn).replace('.fits', '.header')
            if not os.path.exists(backup):
                if args.n:
                    print('I would have saved a backup header:'
                          ' {} -> {}'.format(fn, backup))
                else:
                    print('Backup original FITS header: {} -> {}'.format(
                        fn, backup))
                    fits.getheader(fn).tofile(backup)

            if not args.n:
                assert os.path.exists(backup)

            with fits.open(fn, mode) as hdu:
                for k in row['keyword'].split('|'):
                    assert k in hdu[0].header
                    if hdu[0].header[k] == row['value']:
                        print('{}: {} keyword verified'.format(fn, k))
                    else:
                        msg = '{} from {} to {}'.format(
                            k, hdu[0].header[k], row['value'])

                        if args.n:
                            msg = 'I would have changed {}'.format(msg)
                        else:
                            hdu[0].header[k] = row['value']
                            msg = 'Changed {}'.format(msg)
                            hdu[0].header.add_history(msg)

                        print('{}: {}'.format(fn, msg))
