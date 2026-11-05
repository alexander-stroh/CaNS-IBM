#!/usr/bin/env python3
"""
Results of the CHT examples (after the run), to compare with the values in README.md.

  python check.py slab | wavy | cylinder [case folder]

slab    : last checkpoint against the exact steady solution (Poiseuille flow; in the fluid T = T_i + q (d^2 - y^2)/(2 k_f),
          interface T_i = q d H/k_s): max |T - T_ex|/T_ex(centre) per scalar and max |u - u_ex|/U_max
wavy    : last checkpoint: bulk velocity of the fluid and mean temperature of the fluid and of the solid per scalar
cylinder: t >= 200: bulk velocity, drag and lift coefficients (C = 2 F/(U_b^2 D), the force of data/ibm_force.out minus
          the driving force on the solid volume), Strouhal number from the lift; mean T (data/forcing_s_*.out)
"""
import re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def nml(case):
    txt = (case/'input.nml').read_text()
    get = lambda k: re.search(k + r'\s*=\s*([^\n]+)', txt).group(1)
    vals = lambda k: [float(v) for v in re.split(r'[,\s]+', get(k).strip()) if v]
    return dict(ng=[int(v) for v in vals(r'ng\(1:3\)')], l=vals(r'\bl\(1:3\)'), visci=vals('visci')[0],
                bforce=vals(r'bforce\(1:3\)')[0], ks=vals(r'ks_kf\(:\)'), alphai=vals(r'alphai\(:\)')[0])


def fields(case, p, names):
    nx, ny, nz = p['ng']; N = nx*ny*nz
    out = {}
    for n in names:
        a = np.fromfile(case/'data'/f'fld_{n}.bin')
        out[n], t = a[:N].reshape(nz, ny, nx), a[N]
    return out, t


def slab(case):
    p = nml(case); nz, lz = p['ng'][2], p['l'][2]
    zc = (np.arange(nz) + 0.5)*lz/nz
    H = 1.02; d = lz/2 - H; nu = 1/p['visci']; kf = 1/p['alphai']
    F, t = fields(case, p, ['u'] + [f's_{i:03d}' for i in range(1, len(p['ks']) + 1)])
    y = zc - lz/2; fl = np.abs(y) < d
    uex = np.where(fl, p['bforce']/(2*nu)*(d*d - y*y), 0.)
    print(f'slab, t = {t:g}: max |u - u_ex|/U_max = {np.abs(F["u"].mean((1, 2)) - uex).max()/uex.max():.2e}')
    for i, ks in enumerate(p['ks']):
        Ti = d*H/(ks*kf)                                                   # q = 1
        Tex = np.where(fl, Ti + (d*d - y*y)/(2*kf), Ti*np.minimum(zc, lz - zc)/H)
        T = F[f's_{i + 1:03d}'].mean((1, 2))
        print(f'  k_s/k_f = {ks:5g}: T_centre = {T[nz//2 - 1:nz//2 + 1].mean():9.4f} (exact {Tex.max():9.4f}), '
              f'max |T - T_ex|/T_ex(centre) = {np.abs(T - Tex).max()/Tex.max():.2e}')


def wavy(case):
    p = nml(case); nx, ny, nz = p['ng']
    phi = np.fromfile(case/'data'/'sdfp.bin').reshape(nz, ny, nx)
    F, t = fields(case, p, ['u'] + [f's_{i:03d}' for i in range(1, len(p['ks']) + 1)])
    print(f'wavy, t = {t:g}: U_b = {F["u"].mean()/(phi > 0).mean():.5f} (flow rate / fluid area)')
    for i, ks in enumerate(p['ks']):
        T = F[f's_{i + 1:03d}']
        print(f'  k_s/k_f = {ks:5g}: mean T fluid {T[phi > 0].mean():9.4f}, solid {T[phi < 0].mean():9.4f}')


def cylinder(case, t0=200., D=1.):
    p = nml(case); (nx, ny, nz), (lx, ly, lz) = p['ng'], p['l']
    a, b = np.loadtxt(case/'data'/'ibm_force.out'), np.loadtxt(case/'data'/'forcing.out')
    vs = (np.fromfile(case/'data'/'sdfp.bin') < 0).mean()*lx*lz           # solid area of the x-z plane
    t = a[:, 0]; m = t >= t0
    ub = np.interp(t, b[:, 0], b[:, 4])
    fx, fz = a[:, 4]*lx*lz - p['bforce']*vs, a[:, 6]*lx*lz                  # per unit span
    q = 0.5*ub**2*D; cd, cl = (fx/q)[m], (fz/q)[m]
    tu = np.linspace(t[m][0], t[m][-1], m.sum()); s = np.interp(tu, t[m], cl - cl.mean())
    sp = np.abs(np.fft.rfft(s*np.hanning(len(s))))**2; k = np.argmax(sp[1:]) + 1
    y0, y1, y2 = np.log(sp[k - 1:k + 2]); k += 0.5*(y0 - y2)/(y0 - 2*y1 + y2)  # parabolic peak
    st = k/(tu[-1] - tu[0])*D/ub[m].mean()
    print(f'cylinder, t = {t0:g} ... {t[-1]:g}: U_b = {ub[m].mean():.4f}, Re_D = {ub[m].mean()*p["visci"]:.2f}, '
          f'C_D = {cd.mean():.4f}, C_D amplitude = {0.5*(cd.max() - cd.min()):.4f}, '
          f'C_L amplitude = {0.5*(cl.max() - cl.min()):.4f}, St = {st:.4f}')
    for i, ks in enumerate(p['ks']):
        f = np.loadtxt(case/'data'/f'forcing_s_{i + 1:03d}.out')
        print(f'  scalar {i + 1} (k_s/k_f = {ks:g}): mean T = {f[f[:, 0] >= t0, 2].mean():.4f}')


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in ('slab', 'wavy', 'cylinder'):
        sys.exit(__doc__)
    globals()[sys.argv[1]](Path(sys.argv[2]) if len(sys.argv) > 2 else HERE/sys.argv[1])
