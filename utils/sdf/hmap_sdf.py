#!/usr/bin/env python3
"""
Signed-distance level set of rough walls given by a height map, at the cell centres of a CaNS grid: the
level-set input of CaNS-IBM (`&ibm geometry = 'sdf', sdf_file = ...`).

  python hmap_sdf.py <case> <height map> [options]

<case>/input.nml gives ng and l; the z grid is read from <case>/grid_zf.in (gtype = 9; n or n+1 face
coordinates), else from <case>/data/grid.bin (written by CaNS at start-up), else it is uniform.
The height map is a periodic 2D array h[j, i] (rows = y, columns = x, no duplicated end row/column) as .npy
or plain text; it spans the whole domain lx x ly unless --lmap gives its own extent (then it is repeated).

Lower wall:  z_w(x, y) = zm + h(x, y) - <h>    (zm = mean wall position, --zm)
Upper wall:  --upper mirror: lz - z_w(x, y)   (reflection z -> lz - z)
             --upper shift : lz - z_w(x + lx/2, y + ly/2)   (the convention of &hmap in CaNS-IBM)
             --upper none  : no upper wall
phi = min(phi_lower, phi_upper) > 0 in the fluid, float64, x fastest (C order (nz, ny, nx)).

Method: periodic bicubic spline of the map; within |z - z_w| < B (B = --band cells x max(dx, dy)) the exact
Euclidean distance to the surface by Newton projection from nine start points around the vertical foot (the
smallest converged distance wins, sign from z - z_w); beyond the band the slope-corrected vertical distance
(z - z_w)/sqrt(1 + |grad z_w|^2), which is smooth and keeps the sign. The solver needs the correct sign
everywhere and an accurate, smooth zero level; |grad phi| = 1 matters only for the per-point wall outputs.
"""
import argparse, re
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline


def read_nml(case):
    txt = re.sub(r'!.*', '', (Path(case)/'input.nml').read_text())
    def get(name, cast):
        m = re.search(name + r'\s*\(1:3\)\s*=\s*([^\n/]+)', txt, flags=re.I) or \
            re.search(name + r'\s*=\s*([^\n/]+)', txt, flags=re.I)
        return [cast(v) for v in re.split(r'[,\s]+', m.group(1).strip()) if v][:3]
    return get('ng', int), get('l', lambda v: float(v.replace('d', 'e')))


def z_faces(case, nz, lz):
    case = Path(case)
    if (case/'grid_zf.in').exists():
        zf = np.loadtxt(case/'grid_zf.in')
        zf = zf - zf[0] if len(zf) == nz + 1 else np.concatenate([[0.], zf])
    elif (case/'data'/'grid.bin').exists():
        g = np.fromfile(case/'data'/'grid.bin').reshape(4, nz)  # dzc, dzf, zc, zf(1:n) (CaNS save_grid)
        zf = np.concatenate([[0.], g[3]])
    else:
        zf = np.linspace(0., lz, nz + 1)
    assert len(zf) == nz + 1 and abs(zf[-1] - lz) < 1e-8*lz, 'z grid does not match ng(3), l(3)'
    return zf


def height_spline(h, lxm, lym, zm):
    ny, nx = h.shape
    dx, dy = lxm/nx, lym/ny
    p = 6                                                   # periodic padding for the spline
    hp = np.pad(h, p, mode='wrap')
    xs, ys = np.arange(-p, nx + p)*dx, np.arange(-p, ny + p)*dy
    spl = RectBivariateSpline(ys, xs, hp, kx=3, ky=3, s=0)
    off = zm - h.mean()
    def zw(xx, yy, dxo=0, dyo=0):
        return spl.ev(np.mod(yy, lym), np.mod(xx, lxm), dx=dyo, dy=dxo) + (off if (dxo == 0 and dyo == 0) else 0.)
    return zw


def project(zw, px, py, pz, x0, y0, it=30):
    """distance from p to the closest point of z = zw(x, y), Newton on |(s, t, zw(s, t)) - p|^2 from (x0, y0)"""
    s, t = x0.copy(), y0.copy()
    for _ in range(it):
        h = zw(s, t); hx = zw(s, t, 1, 0); hy = zw(s, t, 0, 1)
        hxx = zw(s, t, 2, 0); hxy = zw(s, t, 1, 1); hyy = zw(s, t, 0, 2)
        r = h - pz
        gs = (s - px) + r*hx; gt = (t - py) + r*hy
        a11 = 1 + hx*hx + r*hxx; a12 = hx*hy + r*hxy; a22 = 1 + hy*hy + r*hyy
        det = a11*a22 - a12*a12
        ok = det > 1e-12                                      # Newton only where the Hessian is positive definite
        ds = np.where(ok, (a22*gs - a12*gt)/np.where(ok, det, 1.), 0.5*gs)
        dt = np.where(ok, (a11*gt - a12*gs)/np.where(ok, det, 1.), 0.5*gt)
        stp = np.maximum(1., np.sqrt(ds*ds + dt*dt)/0.02)     # limit the step
        s -= ds/stp; t -= dt/stp
    return np.sqrt((s - px)**2 + (t - py)**2 + (zw(s, t) - pz)**2)


def sdf_wall(zw, X, Y, Z, B):
    """signed distance to the wall z = zw(x, y) (> 0 above), exact within |z - z_w| < B"""
    h = zw(X, Y); hx = zw(X, Y, 1, 0); hy = zw(X, Y, 0, 1)
    vert = Z - h
    phi = vert/np.sqrt(1 + hx*hx + hy*hy)
    band = np.abs(vert) < B
    px, py, pz = X[band], Y[band], Z[band]
    best = np.abs(vert[band])                                 # the vertical distance bounds the distance
    for ox in (-0.5, 0., 0.5):
        for oy in (-0.5, 0., 0.5):
            best = np.minimum(best, project(zw, px, py, pz, px + ox*best, py + oy*best))
    phi[band] = np.sign(vert[band])*best
    return phi


_G = {}
def _init(h, lxm, lym, zm, xc, yc, lx, ly, lz, B, upper):
    Y2, X2 = np.meshgrid(yc, xc, indexing='ij')
    _G.update(zw=height_spline(h, lxm, lym, zm), X2=X2, Y2=Y2, lx=lx, ly=ly, lz=lz, B=B, upper=upper)


def _level(z):
    g = _G; X2, Y2 = g['X2'], g['Y2']
    phi = sdf_wall(g['zw'], X2, Y2, np.full_like(X2, z), g['B'])
    if g['upper'] == 'mirror':
        phi = np.minimum(phi, sdf_wall(g['zw'], X2, Y2, np.full_like(X2, g['lz'] - z), g['B']))
    elif g['upper'] == 'shift':
        phi = np.minimum(phi, sdf_wall(g['zw'], X2 + 0.5*g['lx'], Y2 + 0.5*g['ly'], np.full_like(X2, g['lz'] - z), g['B']))
    return phi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case'); ap.add_argument('hmap')
    ap.add_argument('--zm', type=float, required=True, help='mean position of the lower wall')
    ap.add_argument('--lmap', type=float, nargs=2, help='extent lx_map ly_map of the map (default: lx ly)')
    ap.add_argument('--upper', choices=['mirror', 'shift', 'none'], default='mirror')
    ap.add_argument('--band', type=float, default=5., help='exact-distance band in cells of max(dx, dy)')
    ap.add_argument('--np', type=int, default=1, help='processes')
    ap.add_argument('-o', '--out', help='output file (default <case>/data/sdfp.bin)')
    a = ap.parse_args()
    (nx, ny, nz), (lx, ly, lz) = read_nml(a.case)
    h = np.load(a.hmap) if a.hmap.endswith('.npy') else np.loadtxt(a.hmap)
    lxm, lym = a.lmap if a.lmap else (lx, ly)
    zf = z_faces(a.case, nz, lz); zc = 0.5*(zf[1:] + zf[:-1])
    dx, dy = lx/nx, ly/ny; B = a.band*max(dx, dy)
    xc, yc = (np.arange(nx) + 0.5)*dx, (np.arange(ny) + 0.5)*dy
    with Pool(a.np, initializer=_init, initargs=(h, lxm, lym, a.zm, xc, yc, lx, ly, lz, B, a.upper)) as pool:
        phi = np.stack(pool.map(_level, list(zc), chunksize=1))
    out = Path(a.out) if a.out else Path(a.case)/'data'/'sdfp.bin'
    out.parent.mkdir(parents=True, exist_ok=True)
    np.ascontiguousarray(phi).tofile(out)
    gz = np.gradient(phi, zc, axis=0)
    gy = (np.roll(phi, -1, 1) - np.roll(phi, 1, 1))/(2*dy)
    gx = (np.roll(phi, -1, 2) - np.roll(phi, 1, 2))/(2*dx)
    g = np.sqrt(gx**2 + gy**2 + gz**2)[np.abs(phi) < 0.8*B]
    vol = np.sum((phi <= 0)*np.diff(zf)[:, None, None])*dx*dy/(lx*ly)
    print(f'{nx} x {ny} x {nz}, band B = {B:.4g}, map {h.shape[1]} x {h.shape[0]} on {lxm:g} x {lym:g}, '
          f'h: mean {h.mean():.4g}, range {h.min():.4g} ... {h.max():.4g}')
    print(f'|grad phi| for |phi| < 0.8 B: 1/50/99 % = {np.quantile(g, 0.01):.3f} / {np.quantile(g, 0.5):.3f} / '
          f'{np.quantile(g, 0.99):.3f};  solid height per unit area (cell signs) {vol:.4f}')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
