#!/usr/bin/python3
import re
import os
import sys
import argparse
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, Column
import pyds9
import mskpy

parser = argparse.ArgumentParser(description='Define background sections for target data.')
parser.add_argument('file', nargs='+', help='FITS files to examine.')
parser.add_argument('--overwrite', action='store_true', help='Re-define the background section.')
parser.add_argument('-o', default='bg-sections.txt', help='Save sections to this file.')

args = parser.parse_args()

number = '([0-9]+(\.[0-9]*(e[+-][0-9]+)*)*)'
boxpat = 'box\({}\)'.format(','.join([number] * 5))
annpat = 'annulus\({}\)'.format(','.join([number] * 4))

ds9 = pyds9.DS9('bgsections', wait=60)
ds9.set('cmap viridis')
ds9.set('scale zscale')
ds9.set('mode pointer')

sections = Table(names=('file', 'sections'), dtype=('U64', 'U128'))

if os.path.exists(args.o):
    tab = ascii.read(args.o)
    os.system('cp {0} {0}~'.format(args.o))
    print('Read in previous table and created backup file.')
    sections = Table(names=('file', 'sections'), dtype=('U64', 'U128'),
                     data=(tab['file'].data, tab['sections'].data))
else:
    sections = Table(names=('file', 'sections'), dtype=('U64', 'U128'))

c = None
last_object = None
for f in sorted(args.file):
    if (f in sections['file']) and (not args.overwrite):
        continue

    h = fits.getheader(f)
    if h['OBSTYPE'] != 'OBJECT':
        continue

    print(f, h['OBJECT'], h['FILTERS'])

    if f in sections['file']:
        section = sections['sections'][sections['file'] == f][0]
        ds9.set('fits {}'.format(f))
        ds9.set('zoom to 0.25')
        ds9.set('regions color red')
        ds9.set('regions system physical')
        if 's:' in section:
            j = section.find('[')
            for s in section[j:].split('+'):
                corners = [int(x) for x in re.split('[,:]', s[1:-1])]
                xc = np.array(corners[:2]).mean().astype(int) + 1
                xw = np.array(corners[:2]).ptp()
                yc = np.array(corners[2:]).mean().astype(int) + 1
                yw = np.array(corners[2:]).ptp()
                ds9.set('regions', 'box({}, {}, {}, {}, 0)'.format(xc, yc, xw, yw))
        else:
            xc, yc = h['CX'] + 1, h['CY'] + 1
            rinner, router = [int(x) for x in section[3:].split(',')]
            ds9.set('regions', 'annulus({}, {}, {}, {})'.format(rinner, router))
    else:
        if h['OBJECT'] == last_object:
            ds9.set('regions system wcs')
            regions = ds9.get('regions')
            
        ds9.set('fits {}'.format(f))
        ds9.set('zoom to 0.25')

        if h['OBJECT'] == last_object:
            ds9.set('regions system wcs')
            ds9.set('regions color green')
            regions = ds9.set('regions', regions)
            ds9.set('regions system physical')

    print("- Mark background with annulus or boxes and press enter or")
    print("  - [s]kip\n  - [w]rite and quit\n  - [q]uit")
    sys.stdin.flush()
    c = None
    while c not in ['', 's', 'w', 'q']:
        c = sys.stdin.readline().strip()
    section = None
    if c == 's':
        continue
    elif c in ['w', 'q']:
        break

    ds9.set('regions system physical')
    regions = ds9.get('regions')
    if 'box' in regions:
        boxes = []
        for box in re.findall(boxpat, regions):
            x, y, xw, yw, angle = [float(v) for v in box[::3]]
            corners = [y - yw/2, y + yw/2, x - xw/2, x + xw/2]
            boxes.append('[{:.0f}:{:.0f}, {:.0f}:{:.0f}]'.format(
                *corners))
        section = 's: {}'.format('+'.join(boxes))
    elif 'annulus' in regions:
        for annulus in re.findall(annpat, regions)[:1]:
            x, y, r0, r1  = [float(v) for v in annulus[::3]]
        section = 'r: {:.0f}, {:.0f}'.format(r0, r1)

    if section is not None:
        assert len(section) > 3, regions
        i = np.flatnonzero(sections['file'] == f)
        if len(i) == 0:
            sections.add_row((f, section))
        else:
            sections[i]['sections'] = section
        print(section)

    last_object = h['OBJECT']

if (c != 'q') and (c is not None):
    sections.write('bg-sections.txt', format='ascii', delimiter=',',
                   overwrite=True)
