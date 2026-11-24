#!/usr/bin/env python3
"""
Gate 1e of the port: the tangential face flux of CaNS (`is_cht_tang`) against AMPHIBIOUS on the node of
"The example in numbers" (cht_method_results.tex; reference values in
`worked_example.json` of the AMPHIBIOUS project, passed with --ref).

AMPHIBIOUS: node (j,k) = (5,7) at (y,z) = (0.2, 0.3), wall y_w(z) = 0.225 + 0.1 cos(2 pi z) (solid below),
h = 0.05, kappa = C_s = 2, x-invariant. Its nodes sit at (j-1)h, (k-1)h; CaNS's cell centres at (k-1/2)h. In
CaNS axes (z wall-normal, y spanwise) the same node is cell (i, 7, 5) on the wall z_w(y) = 0.25 + 0.1 cos 2 pi
(y - 0.025), box 0.2 x 1 x 1, 4 x 20 x 20 cells.

Compared per face (S: the z face below the cell, W/E: the y faces): the series conductance k_face and the
coefficients c = (k_area - k_face) w of the tangential gradient on the block, summed over the three x-planes
(the 2-D code fits on the 2x3 in-plane block; with the frame's first tangent along x the 3-D fit is that fit
on each x-plane, a third of the weight each). The geometry differs by design: AMPHIBIOUS bisects the exact
wall, CaNS works from the level set at the cell centres (quadratic link fractions, Taylor-model face
fractions), so agreement to the accuracy of that path is expected, not to round-off.

    python3 worked_example_check.py --cans <cans> --ref <worked_example.json> [--work DIR]
"""
import argparse, json, math, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cht_order as co

H = 0.05


class Wall:
    """the wall z = g(y) in the form cht_order.sdf expects (functions of x, y)"""
    def __init__(self):
        w = 2*math.pi
        self.g = lambda x, y: 0.25 + 0.1*np.cos(w*(y - 0.025)) + 0*x
        self.gx = lambda x, y: 0*x + 0*y
        self.gy = lambda x, y: -0.1*w*np.sin(w*(y - 0.025)) + 0*x
        self.gxx = lambda x, y: 0*x + 0*y
        self.gxy = lambda x, y: 0*x + 0*y
        self.gyy = lambda x, y: -0.1*w*w*np.cos(w*(y - 0.025)) + 0*x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cans', required=True)
    ap.add_argument('--ref', required=True, help='worked_example.json, the AMPHIBIOUS reference values')
    ap.add_argument('--work', default=str(HERE.parent.parent/'run'/'worked_example'))
    a = ap.parse_args()
    d = Path(a.work); (d/'data').mkdir(parents=True, exist_ok=True)
    nml = co.NML.format(nx=4, ny=20, nz=20, lz=1.0, dt=1.e-12, nstep=1, restart='F', kappa=2.0, cs=2.0,
                        ssource=0., ssource_s=0., bc0=0., bc1=0., geometry='sdf', analytic='',
                        scheme='is_cht_rkl2(:) = T\nis_cht_tang(:) = T\n')
    nml = nml.replace('l(1:3) = 1., 1., 1.0', 'l(1:3) = 0.2, 1., 1.')
    assert 'l(1:3) = 0.2, 1., 1.' in nml
    (d/'input.nml').write_text(nml)
    xc = (np.arange(4) + 0.5)*H; yc = (np.arange(20) + 0.5)*H; zc = (np.arange(20) + 0.5)*H
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing='ij')
    ps = co.sdf(Wall(), X, Y, Z)
    np.ascontiguousarray(ps.transpose(2, 1, 0), dtype='<f8').tofile(d/'data'/'sdfp.bin')
    log = open(d/'run.log', 'w')
    r = subprocess.run(['mpirun', '-np', '1', a.cans], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    assert r.returncode == 0, 'CaNS failed, see run.log'
    # CaNS: link fractions and the face list
    n = 4*20*20
    g = np.fromfile(d/'data'/'ibm_geometry.bin', '<f8')
    th = [g[(4 + m)*n:(5 + m)*n].reshape(20, 20, 4).transpose(2, 1, 0) for m in range(4)]
    raw = (d/'data'/'cht_tang_001.bin').read_bytes()
    nf = int(np.frombuffer(raw[:4], '<i4')[0])
    rec = np.dtype([('ijkd', '<i4', 4), ('c', '<f8', 18)])
    faces = {tuple(int(v) for v in f['ijkd']): f['c'] for f in np.frombuffer(raw[16:], rec, count=nf)}
    kface = lambda t: 1./(t + (1 - t)/2.0)
    ref = json.load(open(a.ref))
    out = []
    # (name, AMPHIBIOUS face id, CaNS lower cell (j,k) and axis)
    for name, fid, (j, k, ax) in (('S', '1', (7, 4, 3)), ('W', '3', (6, 5, 2)), ('E', '4', (7, 5, 2))):
        rf = ref['faces'][fid]
        P = ref['P']
        # AMPHIBIOUS coefficients by offset (jj, kk) from P: contrib = scale*wS, scale = fS*fac*d1
        scale = float(eval(rf['scale'])) if isinstance(rf.get('scale'), str) else rf['scale']
        contrib = eval(rf['contrib']) if isinstance(rf['contrib'], str) else rf['contrib']
        fS = rf['fS']; d1 = ref['case']['d1y']
        amph = {(c['jj'], c['kk']): c['value']/(fS*d1) for c in contrib}      # = fac*wS
        # CaNS: the face at every i (x-invariant: all equal), coefficients summed over the x offsets and
        # mapped to (jj, kk): jj = z offset of the block cell from the CaNS cell (i,7,5), kk = its y offset
        cs = []
        for i in range(1, 5):
            c = faces.get((i, j, k, ax))
            assert c is not None, ('no term on face', name, i)
            acc = {}
            for m in range(1, 19):
                ib, jb, kb = co.block_offset(m, ax)
                key = (k + kb - 5, j + jb - 7)
                acc[key] = acc.get(key, 0.) + c[m - 1]
            cs.append(acc)
        spread = max(abs(cs[0][kk] - c_[kk]) for c_ in cs[1:] for kk in cs[0])
        cans = cs[0]
        tval = (th[2][0, j - 1, k - 1] if ax == 3 else th[1][0, j - 1, k - 1])
        kf_c, kf_a = kface(tval), rf['kF']
        keys = sorted(amph)
        assert set(keys) == set(cans), (keys, sorted(cans))
        num = max(abs(cans[q] - amph[q]) for q in keys); den = max(abs(amph[q]) for q in keys)
        out.append(dict(face=name, kF_cans=kf_c, kF_amph=kf_a, kF_rel=abs(kf_c - kf_a)/kf_a,
                        fac_amph=rf['fac'], coeff_amph={str(q): amph[q] for q in keys},
                        coeff_cans={str(q): cans[q] for q in keys}, coeff_maxdiff_rel=num/den,
                        x_spread=spread))
        print('face %s: k_face CaNS %.6f  AMPHIBIOUS %.6f  (rel. %.1e); (k_area-k_face) w by (jj,kk): max |diff| / max |c| = %.2e; '
              'spread over x %.1e' % (name, kf_c, kf_a, abs(kf_c - kf_a)/kf_a, num/den, spread))
        for q in keys:
            print('      %-8s CaNS %+.6f  AMPHIBIOUS %+.6f' % (q, cans[q], amph[q]))
    json.dump(out, open(d/'worked_example_check.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
