import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import ascii
import mskpy
from mskpy import photometry, hms2dh, dh2hms, niceplot
from mskpy.photometry import hb

phot = ascii.read('standard-phot.txt', format='fixed_width_two_line')
phot = phot[np.isfinite(phot['m_inst'] * phot['m_inst_err'])]
phot.sort('date')

dm = (phot['m'] - phot['m_inst']).data
dm_unc = np.sqrt(phot['m_inst_err'].data**2 + phot['m_err'].data**2)
h = 2380 * u.m  # elevation of DCT
filters = {
    'RC': phot['filter'] == 'RC',
    'BC': phot['filter'] == 'BC',
    'UC': phot['filter'] == 'UC',
    'NH': phot['filter'] == 'NH',
    'C2': phot['filter'] == 'C2',
    'CN': phot['filter'] == 'CN',
    'OH': phot['filter'] == 'OH',
    'r\'': phot['filter'] == 'SDSS-R',
    'g\'': phot['filter'] == 'SDSS-G',
}
color_correct = ["r'", "g'"]
X = phot['airmass'].data
time = np.array(hms2dh([d.split()[1] for d in phot['date'].data]))

# Ignore G-types for OH filter (in the future, just don't observe G
# types)
HB_G_types = ['HD 11131', 'HD 25680', 'HD 28099', 'HD 29461', 'HD 30246',
              'HD 76151', 'HD 81809', 'HD 146233', 'HD 186408', 'HD 186427',
              'HD 191854', 'HD 217014']
G_types = np.prod([phot['object'] == g for g in HB_G_types], 0, bool)
# ignore_g_oh = G_types * oh

# bad_sources are excluded from the fit
bad_sources = []
good_sources = np.prod([phot['object'] != b for b in bad_sources], 0, bool)

good_times = True
# good_times = ~(bc * mskpy.between(t, [2 + 58/60, 3 + 14/60]))

good = good_sources * good_times

outf = open('cal-phot.txt', 'w')

#

A, A_unc = {}, {}
for filt, pset in filters.items():
    if sum(pset) == 0:
        continue

    if filt == 'OH':
        continue
    m_inst = phot['m_inst'][good * pset]
    m_inst_unc = phot['m_inst_err'][good * pset]
    M = phot['m'][good * pset]
    x = X[good * pset]
    color = phot['color'][good * pset]
    if filt in color_correct:
        A[filt], A_unc[filt] = photometry.cal_color_airmass(
            m_inst, m_inst_unc, M, color, x)
        m = (M - A[filt][0] + A[filt][1] * x + A[filt][2] * color).data
        res = (m - m_inst).data
    else:
        A[filt], A_unc[filt] = photometry.cal_airmass(
            m_inst, m_inst_unc, M, x)
        m = (M - A[filt][0] + A[filt][1] * x).data
        res = (m - m_inst).data

    outf.write("""
{} calibration
{}------------
""".format(filt, '-' * len(filt)))
    outf.write("""
zp = {:9.4f} +/- {:9.4f}  magnitude zeropoint [mag]
Ex = {:9.4g} +/- {:9.4g}  airmass extinction  [mag/airmass]
""".format(A[filt][0], A_unc[filt][0], A[filt][1], A_unc[filt][1]))
    if filt in color_correct:
        outf.write(
            "C  = {:9.4f} +/- {:9.4f}  color correction    [mag/color index]\n".format(A[filt][2], A_unc[filt][2]))

    outf.write("""
Residuals
  N = {}
  mean = {: 9.4g}
  stdev = {: 9.4g}
  error on the mean = {: 9.4g}
""".format(len(res),
           np.mean(res),
           np.std(res),
           np.std(res) / np.sqrt(len(res))))

if 'OH' in filters:
    assert 'BC' in A, 'OH requires BC calibration'
    i = good * filters['OH']
    m_inst = phot['m_inst'][i].data
    m_inst_unc = phot['m_inst_err'][i].data
    M = phot['m'][i].data
    z = phot['za'][i].data * u.deg

    A['OH'], A_unc['OH'] = hb.cal_oh(m_inst, m_inst_unc, M, z,
                                     'b', 'b', A['BC'][1], h)
    ext_oh = hb.ext_total_oh(A['OH'][1], z, 'b', 'b', A['BC'][1], h)
    m = M - A['OH'][0] + ext_oh
    res = m - m_inst
    outf.write("""
OH calibration
--------------

zp  = {:9.4f} +/- {:9.4f}    magnitude zeropoint [mag]
toz = {:9.4g} +/- {:9.4g}    ozone  [unk]

Residuals
  N = {}
  mean = {:9.4g}
  stdev = {:9.4g}
  error on the mean = {:9.4g}
""".format(A['OH'][0], A_unc['OH'][0], A['OH'][1], A_unc['OH'][1], len(res),
           np.mean(res), np.std(res), np.std(res) / np.sqrt(len(res))))

outf.close()
print(open('cal-phot.txt', 'r').read(-1))


#

mcolor = plt.cm.rainbow(np.linspace(0, 1, len(phot)))

fig = plt.figure(1)
fig.clear()
plt.minorticks_on()
am = np.linspace(1.0, 4)
for filt, pset in filters.items():
    if filt == 'OH' or filt not in A:
        continue
    i = good * pset
    bad = ~good * pset
    if filt in color_correct:
        y = dm + A[filt][2] * phot['color']
    else:
        y = dm

    plt.scatter(X[i], y[i], color=mcolor[i], alpha=0.5)
    if np.any(bad):
        plt.scatter(X[bad], y[bad], color=mcolor[bad], marker='x', alpha=0.5)

    plt.errorbar(X[i], y[i], dm_unc[i], color=(0, 0, 0, 0), ecolor='k')
    plt.plot(am, A[filt][0] - A[filt][1] * am, 'k-')
    plt.text(1.0, A[filt][0] - A[filt][1], filt + ' ', ha='right',
             va='center')

plt.setp(plt.gca(), xlabel=r'Airmass', ylabel=r'$\Delta m + E(X)$ (mag)')
niceplot(lw=1, mew=0.5, ms=5)
plt.draw()
plt.savefig('cal-phot-extinction.pdf')

if len(color_correct) > 0:
    fig = plt.figure(2)
    fig.clear()
    plt.minorticks_on()
    ci = np.linspace(-0.2, 1.5)
    for filt in color_correct:
        if filt not in A:
            continue

        y = dm - A[filt][0] + A[filt][1] * X
        i = good * filters[filt]
        bad = ~good * filters[filt]
        plt.scatter(phot['color'][i], y[i], color=mcolor[i], alpha=0.5)
        if np.any(bad):
            plt.scatter(phot['color'][bad], y[bad], color=mcolor[bad],
                        marker='x', alpha=0.5)
        plt.errorbar(phot['color'][i], y[i], dm_unc[i], color=(0, 0, 0, 0),
                     ecolor='0.5')
        plt.plot(ci, -A[filt][2] * ci, 'k-')
        plt.text(ci[0], -A[filt][2] * ci[0], filt + ' ', ha='right',
                 va='center')
        plt.setp(plt.gca(), xlabel='Airmass',
                 ylabel=r"$\Delta m + C (g'-r')$ (mag)")
    plt.savefig('cal-phot-colorcor.pdf')

if 'OH' in A:
    i = good * filters['OH']
    bad = ~good * filters['OH']

    fig = plt.figure(3)
    fig.clear()
    plt.minorticks_on()

    plt.scatter(X[i], dm[i], color=mcolor[i], alpha=0.5)
    if np.any(bad):
        plt.scatter(X[bad], dm[bad], color=mcolor[bad], marker='x', alpha=0.5)
    plt.errorbar(X[i], dm[i], dm_unc[i], color=(0, 0, 0, 0), ecolor='k')

    z = np.linspace(0, 81, 300) * u.deg
    x = photometry.airmass_app(z, h)
    ext_oh = hb.ext_total_oh(A['OH'][1], z, 'b', 'b', A['BC'][1], h)
    plt.plot(x, A['OH'][0] - ext_oh, 'k-')
    plt.setp(plt.gca(), xlabel='Apparent Airmass',
             ylabel=r'$\Delta OH + E(X)$ (mag)')
    plt.savefig('cal-phot-oh.pdf')
