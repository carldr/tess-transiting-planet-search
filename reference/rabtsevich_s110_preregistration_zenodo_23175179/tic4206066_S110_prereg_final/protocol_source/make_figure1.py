"""Figure 1: archival transits of A and B, sector by sector, folded on the frozen ephemeris (style of the K2 and TESS
Synergy papers: one panel per dataset), and an empty Sector 110 panel with the expected depth.
Only the transits counted by the protocol are used (calib/expected_depths.json); detrending is the frozen one."""
import json, sys, os, hashlib, numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEP, LC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
RU = len(sys.argv) > 4 and sys.argv[4] == 'ru'
cal = json.load(open(os.path.join(DEP, 'calib/expected_depths.json')))
eph = json.load(open(os.path.join(DEP, 'ephemeris_final.json')))
prov = cal['provenance']['inputs']
sys.path.insert(0, os.path.join(DEP, 'deps'))
from lcio import detrend
E = {k: dict(P=eph[k]['P'], T0=eph[k]['T0'], T14=eph[k]['T14']) for k in 'AB'}
D = {k: tuple(cal['expected'][k]['primary_S6_S32']) for k in 'AB'}   # primary expected depths (Table 3 of the full protocol)
assert round(D['A'][0], 1) == 441.9 and round(D['B'][0], 1) == 585.5
CA, CB, INK, DOT, BAND = '#0072B2', '#D55E00', '#1A1A1A', '#9DB4C8', '#FBE3D6'
HW = 0.25                                                 # half width of the panel, days (6 h)

def load(fn, col, qm):
    p = os.path.join(LC, fn)
    assert hashlib.sha256(open(p, 'rb').read()).hexdigest() == prov[fn], fn
    d = fits.getdata(p, 1)
    t = np.asarray(d['TIME'], float); f = np.asarray(d[col], float); q = np.asarray(d['QUALITY']).astype(np.int64)
    m = np.isfinite(t) & np.isfinite(f) & ((q & qm) == 0) & (f > 0)
    o = detrend(dict(t=t[m], f=f[m] / np.nanmedian(f[m])), window=0.75, eph=E)
    return o['t'], (o['fd'] - 1) * 1e6

secs = [('S6', 'PDCSAP_FLUX', 17087, 'Сектор 6 · 2018' if RU else 'Sector 6 · 2018'),
        ('S32', 'PDCSAP_FLUX', 17087, 'Сектор 32 · 2020' if RU else 'Sector 32 · 2020'),
        ('S98', 'SYS_RM_FLUX', 1610619071, 'Сектор 98 · 2025/26' if RU else 'Sector 98 · 2025/26')]
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'font.serif': ['Times New Roman'], 'mathtext.fontset': 'stix',
                     'pdf.fonttype': 42, 'axes.linewidth': 0.6, 'text.color': INK, 'axes.edgecolor': INK,
                     'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
fig, axes = plt.subplots(2, 4, figsize=(6.5, 3.3), sharex=True, sharey='row', gridspec_kw=dict(wspace=0.08, hspace=0.18))
data = {s: load(cal['per_sector'][s]['file'], col, qm) for s, col, qm, _ in secs}
for r, sig in enumerate('AB'):
    c = CA if sig == 'A' else CB
    mk = 'o' if sig == 'A' else 'D'                      # shape as well as colour (AAS graphics guide)
    P, T14 = E[sig]['P'], E[sig]['T14']
    ylo, yhi = (-1300, 1000) if sig == 'A' else (-1900, 1400)
    for k, (s, _, _, title) in enumerate(secs):
        ax = axes[r, k]
        t, f = data[s]
        tcs = [e['tc'] for e in cal['per_sector'][s][sig]['events'] if e['counted']]
        n = len(tcs)
        if n:
            dt = np.concatenate([t - tc for tc in tcs]); ff = np.concatenate([f for _ in tcs])
            w = np.abs(dt) < HW; dt, ff = dt[w], ff[w]
            ax.scatter(dt * 24, ff, s=1.0, color=DOT, lw=0, alpha=0.6, zorder=1)
            edges = np.arange(-HW, HW + 1e-9, 1 / 48)   # 30-min bins
            idx = np.digitize(dt, edges)
            sig1 = 1.4826 * np.median(np.abs(ff - np.median(ff)))
            u = [i for i in np.unique(idx)]
            bx = [dt[idx == i].mean() * 24 for i in u]
            by = [ff[idx == i].mean() for i in u]
            be = [sig1 / np.sqrt((idx == i).sum()) for i in u]
            ax.errorbar(bx, by, be, fmt=mk, ms=3, color=c, ecolor=c, elinewidth=0.8, zorder=3)
            lab = (('%d транзит' if n == 1 else '%d транзита' if n < 5 else '%d транзитов') % n) if RU else \
                  ('%d transit' % n if n == 1 else '%d transits' % n)
        else:
            # no counted transit (Sector 32, signal B): show the not-counted events that have data, as open symbols,
            # and mark the nearby transits of A that made the protocol leave them out
            evs = [e['tc'] for e in cal['per_sector'][s][sig]['events'] if not e['counted']
                   and np.sum(np.abs(t - e['tc']) < HW) > 20]
            dt = np.concatenate([t - tc for tc in evs]); ff = np.concatenate([f for _ in evs])
            w = np.abs(dt) < HW; dt, ff = dt[w], ff[w]
            ax.scatter(dt * 24, ff, s=1.0, color=DOT, lw=0, alpha=0.6, zorder=1)
            edges = np.arange(-HW, HW + 1e-9, 1 / 48); idx = np.digitize(dt, edges)
            sig1 = 1.4826 * np.median(np.abs(ff - np.median(ff)))
            u = np.unique(idx)
            ax.errorbar([dt[idx == i].mean() * 24 for i in u], [ff[idx == i].mean() for i in u],
                        [sig1 / np.sqrt((idx == i).sum()) for i in u], fmt=mk, ms=3, mfc='white', mec=c, ecolor=c,
                        elinewidth=0.8, zorder=3)
            tA = [e['tc'] for e in cal['per_sector'][s]['A']['events']]
            offs = [(min(tA, key=lambda x: abs(x - tc)) - tc) * 24 for tc in evs]
            ax.plot(offs, [yhi * 0.45] * len(offs), ls='', marker='v', color=CA, ms=4.5, zorder=4)
            ax.text(max(offs) + 0.5, yhi * 0.45, 'transit of A' if not RU else 'транзит A', ha='left', va='center',
                    fontsize=7, color=INK)
            lab = ('%d события,\nне использованы' if RU else '%d events, not used') % len(evs)
        ax.text(0.04, 0.93 if n else 0.06, lab, transform=ax.transAxes, fontsize=7.5, va='top' if n else 'bottom')
        if r == 0:
            ax.set_title(title, fontsize=8.5)
    ax = axes[r, 3]
    ax.axvspan(-T14 * 12, T14 * 12, color=BAND, zorder=0)
    ax.text(0.5, 0.06, ('ожидаемая глубина:\n%d ± %d ppm' if RU else 'expected depth:\n%d ± %d ppm') % (round(D[sig][0]), round(D[sig][1])),
            transform=ax.transAxes, ha='center', fontsize=7.5)
    ax.text(0.5, 0.93, 'данных ещё нет' if RU else 'no data yet', transform=ax.transAxes, ha='center', va='top', fontsize=7.5)
    if r == 0:
        ax.set_title('Сектор 110 · нояб. 2026' if RU else 'Sector 110 · Nov 2026', fontsize=8.5, color=INK)
    for k in range(4):
        a = axes[r, k]
        a.axhline(0, color=INK, lw=0.5, zorder=0)
        a.axhline(-D[sig][0], color=c, lw=1.3, ls=(0, (4, 2)), zorder=2)
        a.set_ylim(ylo, yhi); a.set_xlim(-HW * 24, HW * 24)
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
    axes[r, 0].set_ylabel(('%s: изменение потока, ppm' if RU else '%s: flux change, ppm') % sig)
for k in range(4):
    axes[1, k].set_xlabel('часы от транзита' if RU else 'hours from transit')
fig.savefig(OUT, bbox_inches='tight')
