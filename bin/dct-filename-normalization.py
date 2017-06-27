#!/usr/bin/env python3
import os
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description='Rename DCT data files.',
                                 epilog='Requires TELESCOP="DCT", INSTRUME="LMI", BITPIX=16, BUNIT="ADU", valid OBSERNO and DATE-OBS.  OBSERNO must be <10000.')
parser.add_argument('files', nargs='*', help='The files to rename.  Files with the appropriate FITS keywords will be modified, others will be ignored.')
parser.add_argument('-n', action='store_true', help='No-operation mode.  Just print what would have happened.')

args = parser.parse_args()

renamed = 0
skipped = 0
not_fits = 0
could_not_rename = 0
good_file_name = 0

class CouldNotRenameFile(Exception):
    pass

class FileExists(Exception):
    pass

for fn in args.files:
    try:
        with fits.open(fn, mode='readonly') as hdu:
            for k in ['TELESCOP', 'INSTRUME', 'OBSERNO', 'DATE-OBS']:
                assert k in hdu[0].header

            assert hdu[0].header['BITPIX'] == 16
            assert hdu[0].header['BUNIT'] == 'ADU'
            assert hdu[0].header['TELESCOP'] == 'DCT'
            assert hdu[0].header['INSTRUME'] == 'LMI'

            number = hdu[0].header['OBSERNO']
            if number > 9999:
                skipped += 1
                continue
            
            date = hdu[0].header['DATE-OBS'][:10].split('-')
            assert int(date[0]) > 2010
            assert int(date[1]) > 0 and int(date[1]) < 13
            assert int(date[2]) > 0 and int(date[1]) < 32
            date = ''.join(date)

        path = os.path.split(fn)[0]
        new_fn = os.path.join(path, 'lmi_{}_{:04d}_raw.fits'.format(date, number))

        if fn == new_fn:
            print('{} already a good file name.'.format(fn))
            good_file_name += 1
            continue
        
        print('{} -> {}'.format(fn, new_fn))
        if os.path.exists(new_fn):
            raise FileExists
        if not args.n:
            try:
                os.rename(fn, new_fn)
            except OSError:
                raise CouldNotRenameFile

        renamed += 1
    except AssertionError:
        print(fn, 'bad keywords or data.')
        skipped += 1
    except (IsADirectoryError, OSError):
        print(fn, 'does not seem to have FITS format.')
        not_fits += 1
    except FileExists:
        print(fn, 'new file name already exists.')
        could_not_rename += 1
    except CouldNotRenameFile:
        print(fn, 'error renaming file.')
        could_not_rename += 1

print('''
{good_file_name} files with good file names.
{renamed} files {would_have}renamed.
{skipped} files {would_have}skipped.
{not_fits} files do not seem to have FITS format.
{could_not_rename} files could not be renamed.
'''.format(good_file_name=good_file_name, renamed=renamed, skipped=skipped,
           not_fits=not_fits, could_not_rename=could_not_rename,
           would_have='would have been ' if args.n else ''))
