#!/usr/bin/env python3
import os
import argparse
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from sbpy.data import Names, Ephem, natural_sort_key

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+')

args = parser.parse_args()


def get_ephemeris(target):
    objtype = Names.asteroid_or_comet(target)
    opts = dict(epochs=Time(h['DATE-OBS']), id_type='designation')
    if objtype == 'comet':
        parsed = Names.parse_comet(target)
        t = str(parsed.get('number', '')) + parsed['type']
        if 'desig' in parsed:
            t += '/' + parsed['desig']

        opts.update(dict(
            no_fragments=True,
            closest_apparition=True
        ))
    else:
        parsed = Names.parse_asteroid(target)
        t = parsed.get('number', parsed.get('desig'))

    if t is None:
        raise ValueError

    return Ephem.from_horizons(t, **opts)


# get all headers, sort by object and date
headers = sorted([fits.getheader(f) for f in args.files],
                 key=lambda h: (h['object'].replace(' ', ''), h['date-obs']))
headers.append(None)  # triggers summary of last target

targets = {}
last_k = None
for h in headers:
    if h is None:
        k = None
    else:
        target = h['object']
        k = target.replace(' ', '')

    if k in targets:
        filters.add(h['FILTERS'])
        airmass.append(h['airmass'])
        mjd.append(Time(h['DATE-OBS']).mjd)
        continue
    elif last_k not in [k, None] or h is None:
        print(last_k, k, 'asdf')
        targets[last_k]['filters'] = ' '.join(filters)
        date, time = Time(np.mean(mjd), format='mjd').iso.split()
        targets[last_k]['date'] = date
        targets[last_k]['time'] = time
        targets[last_k]['airmass'] = np.round(np.mean(airmass), 1)

        if h is None:
            break

    date, time = h['DATE-OBS'].split('T')
    mjd = [Time(h['DATE-OBS']).mjd]
    airmass = [h['airmass']]
    filters = set([h['FILTERS']])
    targets[k] = {
        'target': target,
        'ra': h['ra'],
        'dec': h['dec']
    }
    last_k = k

    try:
        eph = get_ephemeris(target)
    except:
        continue

    targets[k].update({
        'rh': np.round(np.sign(eph['rdot'][0]) * eph['rh'][0].value, 2),
        'delta': np.round(eph['delta'][0].value, 2),
        'phase': np.round(eph['phase'][0].value, 1)
    })

rows = sorted(list(targets.values()),
              key=lambda row: natural_sort_key(row['target']))
tab = Table(rows)[
    'target', 'date', 'time', 'ra', 'dec', 'airmass', 'filters', 'rh',
    'delta', 'phase'
]

tab.write('summary.csv')
os.system('cat summary.csv')
