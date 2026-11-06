#!/usr/bin/env python3
"""
Reference values from the finished OpenFOAM runs, the same quantities as ../check.py prints for CaNS.

  python3 extract.py slab|wavy|cylinder [mesh]     (mesh folder, default n064 / n032 / ppd016)

slab    : x-averaged T(z) against the exact steady solution, per k_s/k_f
wavy    : bulk velocity of the fluid and volume-mean T of the fluid and of the solid, per k_s/k_f
cylinder: t >= 200: U_b (flow rate per box area), C_D, C_D and C_L amplitudes, St (from the lift), mean fluid T
Only NumPy; reads the ASCII fields of the latest time (cell centres C and volumes V are written by Allrun).
"""
import re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
T0, LY, D = 300., 0.1, 1.


def info(case):
    return dict(re.findall(r'(\w+) = (\S+)', (case/'case.txt').read_text()))


def latest(case):
    return max((d for d in case.iterdir() if re.fullmatch(r'[0-9.eE+-]+', d.name) and d.name != '0'), key=lambda d: float(d.name))


def field(path):
    """internal field (uniform or nonuniform scalar/vector list) of an ASCII OpenFOAM field file"""
    txt = path.read_text()
    m = re.compile(r'internalField\s+(uniform|nonuniform\s+List<(scalar|vector)>)\s*').search(txt)
    if m.group(1) == 'uniform':
        return np.array([float(x) for x in re.compile(r'\(?([^;]*?)\)?\s*;').match(txt, m.end()).group(1).split()])
    n = re.compile(r'(\d+)\s*\(').match(txt, m.end())
    cnt, start = int(n.group(1)), n.end()
    body = txt[start:txt.index('\n)', start)].replace('(', ' ').replace(')', ' ')
    a = np.array(body.split(), dtype=float)
    return a.reshape(cnt, -1) if m.group(2) == 'vector' else a


def region(d, r, names=('T',)):
    out = dict(C=field(d/r/'C'), V=field(d/r/'V'))
    for nm in names:
        out[nm] = field(d/r/nm)
    return out


def cases(study, mesh):
    cs = sorted(p.parent for p in (HERE/study/mesh).glob('*/case.txt'))
    if not cs:
        sys.exit(f'no cases in {study}/{mesh}')
    return cs


def slab(mesh='n064'):
    for case in cases('slab', mesh):
        p = info(case); d_ = latest(case)
        H, ks, nu = float(p['H']), float(p['ks_kf']), float(p['nu'])
        LZ = 4.; d = LZ/2 - H; kf = nu
        R = {r: region(d_, r, ('T', 'U') if r == 'fluid' else ('T',)) for r in ('lowerSolid', 'fluid', 'upperSolid')}
        z = np.concatenate([R[r]['C'][:, 2] for r in R]); T = np.concatenate([R[r]['T'] for r in R]) - T0
        zu = np.unique(np.round(z, 10)); Tz = np.array([T[np.abs(z - zz) < 1e-9].mean() for zz in zu])
        y = zu - LZ/2; fl = np.abs(y) < d; Ti = d*H/(ks*kf)
        Tex = np.where(fl, Ti + (d*d - y*y)/(2*kf), Ti*np.minimum(zu, LZ - zu)/H)
        zf = R['fluid']['C'][:, 2]; yf = zf - LZ/2
        uex = float(p['bforce'])/(2*nu)*(d*d - yf*yf)
        print(f'slab {case.name}, iteration {d_.name}: T_centre = {Tz[len(Tz)//2 - 1:len(Tz)//2 + 1].mean():9.4f} '
              f'(exact {Tex.max():9.4f}), max |T - T_ex|/T_ex(centre) = {np.abs(Tz - Tex).max()/Tex.max():.2e}, '
              f'max |u - u_ex|/U_max = {np.abs(R["fluid"]["U"][:, 0] - uex).max()/uex.max():.2e}')


def wavy(mesh='n032'):
    for case in cases('wavy', mesh):
        d_ = latest(case)
        f = region(d_, 'fluid', ('T', 'U')); s = [region(d_, r) for r in ('lowerSolid', 'upperSolid')]
        mean = lambda R, a: (a*R['V']).sum()/R['V'].sum()
        Ts = sum((r['T']*r['V']).sum() for r in s)/sum(r['V'].sum() for r in s)
        print(f'wavy {case.name}, iteration {d_.name}: U_b = {mean(f, f["U"][:, 0]):.5f}, '
              f'mean T fluid {mean(f, f["T"]) - T0:9.4f}, solid {Ts - T0:9.4f}')


def cylinder(mesh='ppd016', t0=200., LX=12., LZ=4.):
    for case in cases('cylinder', mesh):
        pp = case/'postProcessing'
        fo = np.loadtxt(pp/'fluid'/'forces'/'0'/'force.dat', usecols=(0, 1, 3))
        hf = np.loadtxt((pp/'fluid'/'historyFluid'/'0'/'volFieldValue.dat').read_text().replace('(', ' ').replace(')', ' ').splitlines(),
                        usecols=(0, 1, 2))
        Vf = field(latest(case)/'fluid'/'V').sum()
        t = fo[:, 0]; m = t >= t0
        ub = np.interp(t, hf[:, 0], hf[:, 2])*Vf/(LX*LZ*LY)                 # flow rate per box area, as CaNS forcing.out
        q = 0.5*ub**2*D; cd, cl = (fo[:, 1]/LY/q)[m], (fo[:, 2]/LY/q)[m]
        tu = np.linspace(t[m][0], t[m][-1], m.sum()); s = np.interp(tu, t[m], cl - cl.mean())
        sp = np.abs(np.fft.rfft(s*np.hanning(len(s))))**2; k = np.argmax(sp[1:]) + 1
        y0, y1, y2 = np.log(sp[k - 1:k + 2]); k += 0.5*(y0 - y2)/(y0 - 2*y1 + y2)
        st = k/(tu[-1] - tu[0])*D/ub[m].mean()
        mh = hf[:, 0] >= t0
        print(f'cylinder {case.name}, t = {t0:g} ... {t[-1]:g}: U_b = {ub[m].mean():.4f}, C_D = {cd.mean():.4f}, '
              f'C_D amplitude = {0.5*(cd.max() - cd.min()):.4f}, C_L amplitude = {0.5*(cl.max() - cl.min()):.4f}, '
              f'St = {st:.4f}, mean fluid T = {hf[mh, 1].mean() - T0:.4f}')


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in ('slab', 'wavy', 'cylinder'):
        sys.exit(__doc__)
    globals()[sys.argv[1]](*sys.argv[2:3])
