#!/usr/bin/env python3
"""
Gate 1d of the port: conjugate-linear exactness of the tangential face flux (`is_cht_tang`).

A planar interface inclined at 30 deg in the y-z plane, z = Z0 + tan(30 deg) (y - 0.5), solid below, in a box
0.25 x 1 x 0.75 (cubic cells, h = 1/N). The level set (the signed distance, linear) makes the code's link
fractions and face fractions exact, so for the conjugate-linear field

    T = a + s1 sigma1 + s2 sigma2 + q n / k_m       (s1 = x, s2, n: along t2 = (0, cos, sin), n = (0, -sin, cos))

(continuous, with a continuous normal flux) every face flux of the scheme,

    F = k_face (T_N - T_P)/h + sum_m c_m T_m ,

must equal the exact one, k_area (sigma . t-part along the face normal) + q n_d, to round-off on the y and z
faces. On the x faces the face fraction comes from midpoint columns across a kinked integrand, so they are
exact to the quadrature only. The series flux alone misses each face by (k_area - k_face) (sigma t)_d.
Faces whose frame window reaches the y seam of the periodic box (where the plane jumps) are left out.

    python3 plane_exactness.py --cans <cans> [--kappas 2,12] [--N 24]
"""
import argparse, json, math, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cht_order as co

ANG = math.radians(30.)
Z0 = 0.35
NV = np.array([0., -math.sin(ANG), math.cos(ANG)])      # into the fluid (above)
T2 = np.array([0., math.cos(ANG), math.sin(ANG)])
P0 = np.array([0., 0.5, Z0])                            # a point of the plane


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cans', required=True)
    ap.add_argument('--work', default=str(HERE.parent.parent/'run'/'plane_exactness'))
    ap.add_argument('--kappas', default='2,12')
    ap.add_argument('--N', type=int, default=24)
    a = ap.parse_args()
    N = a.N; h = 1./N; nx, ny, nz = N//4, N, int(round(0.75*N))
    out = []
    for kap in [float(v) for v in a.kappas.split(',')]:
        d = Path(a.work)/f'k{kap:g}'; (d/'data').mkdir(parents=True, exist_ok=True)
        nml = co.NML.format(nx=nx, ny=ny, nz=nz, lz=0.75, dt=1.e-12, nstep=1, restart='F', kappa=kap, cs=1.0,
                            ssource=0., ssource_s=0., bc0=0., bc1=0., geometry='sdf', analytic='',
                            scheme='is_cht_rkl2(:) = T\nis_cht_tang(:) = T\n')
        nml = nml.replace('l(1:3) = 1., 1., 0.75', 'l(1:3) = %r, 1., 0.75' % (nx*h))
        (d/'input.nml').write_text(nml)
        xc = (np.arange(nx) + 0.5)*h; yc = (np.arange(ny) + 0.5)*h; zc = (np.arange(nz) + 0.5)*h
        X, Y, Z = np.meshgrid(xc, yc, zc, indexing='ij')
        ps = (Y - P0[1])*NV[1] + (Z - P0[2])*NV[2]
        np.ascontiguousarray(ps.transpose(2, 1, 0), dtype='<f8').tofile(d/'data'/'sdfp.bin')
        r = subprocess.run(['mpirun', '-np', '1', a.cans], cwd=d, stdout=open(d/'run.log', 'w'), stderr=subprocess.STDOUT)
        assert r.returncode == 0, 'CaNS failed in %s' % d
        n = nx*ny*nz
        g = np.fromfile(d/'data'/'ibm_geometry.bin', '<f8')
        th = [g[(4 + m)*n:(5 + m)*n].reshape(nz, ny, nx).transpose(2, 1, 0) for m in range(4)]
        raw = (d/'data'/'cht_tang_001.bin').read_bytes()
        nf = int(np.frombuffer(raw[:4], '<i4')[0])
        rec = np.dtype([('ijkd', '<i4', 4), ('c', '<f8', 18)])
        faces = np.frombuffer(raw[16:], rec, count=nf)
        # the conjugate-linear field at the cell centres (periodic images in x, y by index)
        sig1, sig2, qq, aa = 0.7, 1.3, -0.9, 0.4
        def T(i, j, k):
            x = np.array([(i - 0.5)*h, (j - 0.5)*h, (k - 0.5)*h]) - P0
            nn = x @ NV
            km = 1.0 if nn > 0 else kap
            return aa + sig1*x[0] + sig2*(x @ T2) + qq*nn/km
        tvec = sig1*np.array([1., 0., 0.]) + sig2*T2                 # the tangential gradient
        kf = lambda t: 1./(t + (1 - t)/kap)
        errs = {1: [], 2: [], 3: []}; miss = {1: [], 2: [], 3: []}
        for f in faces:
            i, j, k, ax = (int(v) for v in f['ijkd'])
            if not (4 <= j <= ny - 4):                                # away from the y seam
                continue
            lo = (i, j, k); up = (i + (ax == 1), j + (ax == 2), k + (ax == 3))
            # exact face fraction of the plane over the face rectangle: the fluid part of every line along z is
            # exact (the plane is invariant in x), midpoint rule over 20000 lines across it on the x faces
            cen = np.array([(lo[0] - 0.5)*h, (lo[1] - 0.5)*h, (lo[2] - 0.5)*h]); cen[ax - 1] += 0.5*h
            zw = lambda y: P0[2] + math.tan(ANG)*(y - P0[1])            # fluid above the plane
            zfl = lambda y, z0, z1: min(max(z1 - max(zw(y), z0), 0.), z1 - z0)/(z1 - z0)
            if ax == 1:                                                # y-z face: lines along z at many y
                ys = cen[1] + ((np.arange(20000) + 0.5)/20000 - 0.5)*h
                ff = float(np.mean([zfl(y, cen[2] - 0.5*h, cen[2] + 0.5*h) for y in ys]))
            elif ax == 2:                                              # x-z face at fixed y: one line along z
                ff = zfl(cen[1], cen[2] - 0.5*h, cen[2] + 0.5*h)
            else:                                                      # x-y face at fixed z: fluid where zw(y) < z
                y0, y1 = cen[1] - 0.5*h, cen[1] + 0.5*h
                yb = P0[1] + (cen[2] - P0[2])/math.tan(ANG)             # zw(yb) = z
                ff = min(max(yb - y0, 0.), y1 - y0)/(y1 - y0)
            karea = ff*1.0 + (1 - ff)*kap
            Fex = karea*tvec[ax - 1] + qq*NV[ax - 1]
            t = [th[0], th[1], th[2]][ax - 1][(i - 1) % nx, (j - 1) % ny, k - 1]
            Fs1 = kf(t)*(T(*up) - T(*lo))/h
            st = sum(f['c'][m - 1]*T(*(np.array(lo) + np.array(co.block_offset(m, ax)))) for m in range(1, 19))
            errs[ax].append(abs(Fs1 + st - Fex)/max(abs(Fex), 1.))
            miss[ax].append(abs(Fs1 - Fex)/max(abs(Fex), 1.))
        row = dict(kappa=kap, N=N, faces={ax: len(errs[ax]) for ax in errs},
                   max_err_tang={ax: (max(errs[ax]) if errs[ax] else None) for ax in errs},
                   max_err_step1={ax: (max(miss[ax]) if miss[ax] else None) for ax in miss})
        out.append(row)
        print('kappa %g: faces (x, y, z) %s; max relative flux error, tangential: x %.1e  y %.1e  z %.1e;'
              ' series alone: x %.1e  y %.1e  z %.1e' % (kap, [len(errs[q]) for q in (1, 2, 3)],
              *(row['max_err_tang'][q] or 0 for q in (1, 2, 3)), *(row['max_err_step1'][q] or 0 for q in (1, 2, 3))))
    json.dump(out, open(Path(a.work)/'plane_exactness.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
