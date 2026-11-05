#!/usr/bin/env python3
"""
Writes the level set data/sdfp.bin (and for the cylinder the initial fields data/fld_*.bin) of the CHT examples.

  python setup.py slab | wavy | cylinder      (run in this folder; the grid is read from <case>/input.nml)

phi > 0 in the fluid, float64, x fastest (array (nz, ny, nx)), at the cell centres of the uniform grid.
"""
import re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def grid(case):
    txt = (HERE/case/'input.nml').read_text()
    ng = [int(v) for v in re.search(r'ng\(1:3\)\s*=\s*([^\n]+)', txt).group(1).split(',')]
    l = [float(v) for v in re.search(r'\bl\(1:3\)\s*=\s*([^\n]+)', txt).group(1).split(',')]
    nscal = int(re.search(r'nscal\s*=\s*(\d+)', txt).group(1))
    c = [(np.arange(n) + 0.5)*L/n for n, L in zip(ng, l)]       # cell centres
    f = [(np.arange(n) + 1.0)*L/n for n, L in zip(ng, l)]       # faces i+1/2 (u at xf, w at zf)
    return ng, l, nscal, c, f


def write(case, name, a):
    d = HERE/case/'data'; d.mkdir(exist_ok=True)
    np.ascontiguousarray(a, dtype='<f8').tofile(d/name)


# flat solid slabs of thickness H at both walls (interfaces inside cells: H/dz = 16.32 for nz = 64)
def slab():
    H = 1.02
    (nx, ny, nz), (lx, ly, lz), _, (xc, yc, zc), _ = grid('slab')
    phi = np.minimum(zc - H, lz - H - zc)
    write('slab', 'sdfp.bin', np.broadcast_to(phi[:, None, None], (nz, ny, nx)))


# sinusoidal solid walls z_w = H0 + A sin(2 pi x/lx) and lz - z_w: exact distance by Newton from the nearest sample
def wavy():
    H0, A = 1.0, 0.5
    (nx, ny, nz), (lx, ly, lz), _, (xc, yc, zc), _ = grid('wavy')
    k = 2*np.pi/lx
    f, fp, fpp = (lambda s: H0 + A*np.sin(k*s)), (lambda s: A*k*np.cos(k*s)), (lambda s: -A*k*k*np.sin(k*s))
    def dist(x, z):                                               # signed distance to z = f(x), > 0 above
        s0 = np.linspace(0, lx, 4096, endpoint=False)
        ds = (s0[None, :] - x[:, None] + lx/2) % lx - lx/2
        s = x + ds[np.arange(len(x)), np.argmin(ds**2 + (f(s0)[None, :] - z[:, None])**2, axis=1)]
        for _ in range(20):
            g, dg = (s - x) + (f(s) - z)*fp(s), 1 + fp(s)**2 + (f(s) - z)*fpp(s)
            s -= g/np.where(np.abs(dg) > 1e-12, dg, 1e-12)
        return np.sign(z - f(x))*np.hypot(s - x, f(s) - z)
    X, Z = (a.ravel() for a in np.meshgrid(xc, zc))
    phi = np.minimum(dist(X, Z), dist(X, lz - Z)).reshape(nz, nx)
    write('wavy', 'sdfp.bin', np.broadcast_to(phi[:, None, :], (nz, ny, nx)))


# cylinder D = 1 at (3, 2) in a 12 x 4 channel; initial state: Poiseuille flow with U_b = 1, zero in the cylinder,
# a small asymmetric w bump behind it to start the shedding, T = 0 (restart file at t = 0)
def cylinder():
    D, XC, ZC = 1.0, 3.0, 2.0
    (nx, ny, nz), (lx, ly, lz), nscal, (xc, yc, zc), (xf, yf, zf) = grid('cylinder')
    dist = lambda x, z: np.hypot(x[None, :] - XC, z[:, None] - ZC) - D/2       # (nz, nx)
    b = lambda a: np.broadcast_to(a[:, None, :], (nz, ny, nx))
    write('cylinder', 'sdfp.bin', b(dist(xc, zc)))
    u = 1.5*(1 - (2*zc/lz - 1)**2)[:, None]*(dist(xf, zc) > 0)
    w = 0.05*np.exp(-((xc[None, :] - XC - 1.5)**2 + (zf[:, None] - ZC - 0.3)**2)/0.25)*(dist(xc, zf) > 0)
    zero = np.zeros((nz, ny, nx))
    for n, a in [('u', b(u)), ('v', zero), ('w', b(w)), ('p', zero)] + [(f's_{i:03d}', zero) for i in range(1, nscal + 1)]:
        write('cylinder', f'fld_{n}.bin', np.concatenate([np.ravel(a), [0., 0.]]))   # + time, istep


if __name__ == '__main__':
    case = sys.argv[1] if len(sys.argv) > 1 else ''
    if case not in ('slab', 'wavy', 'cylinder'):
        sys.exit(__doc__)
    globals()[case]()
    print(f'wrote {case}/data')
