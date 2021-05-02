#!/usr/bin/env python3
from collections import defaultdict
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.modeling import models, fitting
from mskpy import linefit
from mskpy.photometry import airmass_app

DCT_elevation = 2361 * u.m  # elevation of DCT

airmass = defaultdict(list)
zp = defaultdict(list)
unc = defaultdict(list)
for f in glob('../ppp/*fits'):
    hdu = fits.open(f, mode='readonly')
    if 'cat' not in hdu:
        continue

    ch = hdu['cat'].header
    if 'magzp' not in ch:
        continue

    h = hdu[0].header
    filt = h['filters']
    # exptime[filt].append(h['exptime'])
    zp[filt].append(ch['magzp'] - 2.5 * np.log10(h['exptime']))
    unc[filt].append(ch['mzpunc'])
    # airmass[filt].append(h['airmass'])
    airmass[filt].append(airmass_app(h['ZA'] * u.deg, DCT_elevation))

plt.clf()
fitter = fitting.FittingWithOutlierRemoval(
    fitting.LevMarLSQFitter(),
    sigma_clip, niter=3, sigma=3.0
)
extinction = []
for filt in zp.keys():
    fit, mask = fitter(models.Linear1D(), np.array(airmass[filt]),
                       np.array(zp[filt]),
                       weights=np.array(unc[filt])**-2)
    mms = sigma_clipped_stats(fit(airmass[filt]) - np.array(zp[filt]))
    a = np.linspace(1, max(3, max(airmass[filt]) + 1))
    label = f'{filt}: {fit.intercept.value:.3f} + {fit.slope.value:.3f} X (σ={mms[2]:.3f})'
    plt.errorbar(airmass[filt], zp[filt], unc[filt],
                 ls='none', marker='o', alpha=0.5,
                 label=label)
    plt.plot(a, fit(a), color='k', ls='-', lw=1)
    extinction.append([filt] + [fit.intercept.value, fit.slope.value, mms[2]])

tab = Table(rows=extinction, names=['filter', 'zeropoint', 'extinction',
                                    'stdev'])
tab.write('catalog-extinction.txt', format='ascii.fixed_width_two_line')
plt.legend()
ax = plt.gca()
plt.setp(ax, ylabel='zeropoint (mag)', xlabel='airmass', xscale='log')
ax.minorticks_on()
ax.invert_yaxis()
plt.savefig('lmi-zeropoints.png')
