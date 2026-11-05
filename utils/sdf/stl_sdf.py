#!/usr/bin/env python3
"""
Signed distance at the cell centres of a CaNS grid from an STL surface: the level-set input of CaNS-IBM
(`&ibm geometry = 'sdf', sdf_file = ...`).

  python stl_sdf.py <case> wall.stl [more.stl ...] [options]

<case>/input.nml gives ng and l; the z grid is taken from <case>/data/grid.bin (written by CaNS at start-up,
before the level set is read), <case>/grid_zf.in (gtype = 9, face positions) or is uniform (--uniform).
Output: <case>/data/sdfp.bin (or -o), float64, x fastest (C order (nz, ny, nx)), phi > 0 in the fluid.

The STL describes the boundary of the SOLID with outward normals (pointing into the fluid). It can be
 - a closed solid: faces on the domain boundary planes x = 0, lx, y = 0, ly, z = 0, lz are dropped
   (they are not fluid-solid interfaces),
 - an open surface sheet that is periodic in x and y (e.g. a rough wall as a height-field triangulation),
 - any number of objects in one or several files; the geometry is repeated periodically in x and y.
--flip if the normals point into the solid, --mirror adds the upper wall by the reflection z -> lz - z,
--scale/--shift transform the STL coordinates to the domain [0, lx] x [0, ly] x [0, lz].

Method:
 1. weld the triangle soup (exact duplicates, then within --tol), drop degenerate and domain-boundary faces,
    add periodic copies of the faces near the x/y boundaries and weld again (stitches the seams)
 2. band test: lower bound of the distance from a k-d tree of surface samples (spacing <= the cell size);
    only cells that may lie within |phi| < B go to the exact search (B = --band cells x max(dx, dy, dz))
 3. exact distance in the band: closest point on the mesh (libigl AABB tree) in parallel over z slabs;
    sign from the angle-weighted pseudonormal of the closest face/edge/vertex (Baerentzen & Aanaes 2005),
    or from the generalised winding number of the closed input mesh (--sign winding)
 4. outside the band phi = +-B (only the sign matters there); the sign is set by a flood fill of the
    connected regions outside the band (periodic in x, y), voted by the adjacent band cells
"""
import argparse, os, re, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import igl

# ------------------------------------------------------------------------------------------------ input

def read_stl(fname):
    """triangle soup (nf, 3, 3) float64 from a binary or ASCII STL"""
    size = os.path.getsize(fname)
    with open(fname, 'rb') as f:
        head = f.read(84)
    nf = int.from_bytes(head[80:84], 'little') if size >= 84 else -1
    if size == 84 + 50*nf:  # binary
        dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
        m = np.memmap(fname, dtype=dt, mode='r', offset=84, shape=(nf,))
        return np.asarray(m['v'], dtype=np.float64)
    V, F, _ = igl.readSTL(str(fname))  # ASCII
    return V[F]

def weld(tri, tol):
    """vertices and faces of a triangle soup: exact duplicates merged by sorting, then vertices closer than tol"""
    P = tri.reshape(-1, 3)
    key = np.ascontiguousarray(P).view(np.uint64)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    new = np.ones(len(P), bool)
    new[1:] = np.any(key[order[1:]] != key[order[:-1]], axis=1)
    gid = np.cumsum(new) - 1
    inv = np.empty(len(P), np.int64); inv[order] = gid
    V = P[order[new]]
    F = inv.reshape(-1, 3)
    return weld_near(V, F, tol)

def weld_near(V, F, tol):
    pairs = cKDTree(V).query_pairs(tol, output_type='ndarray')
    if len(pairs):
        g = sparse.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(V), len(V)))
        _, lab = connected_components(g, directed=False)
        first = np.full(lab.max() + 1, -1); first[lab[::-1]] = np.arange(len(V))[::-1]
        rep = first[lab]
        keep, inv = np.unique(rep, return_inverse=True)
        V, F = V[keep], inv[F]
    F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 2] != F[:, 0])]  # collapsed faces
    return V, F

def read_grid(case, nz, lz, uniform):
    case = Path(case)
    if uniform:
        zf = np.linspace(0., lz, nz + 1); return 0.5*(zf[1:] + zf[:-1]), np.diff(zf), 'uniform'
    if (case/'data'/'grid.bin').exists():
        g = np.fromfile(case/'data'/'grid.bin').reshape(4, nz)  # dzc, dzf, zc, zf (CaNS save_grid)
        return g[2], g[1], 'data/grid.bin'
    if (case/'grid_zf.in').exists():
        zf = np.concatenate([[0.], np.loadtxt(case/'grid_zf.in')])
        return 0.5*(zf[1:] + zf[:-1]), np.diff(zf), 'grid_zf.in'
    sys.exit('no z grid: run CaNS once to write data/grid.bin, provide grid_zf.in, or use --uniform')

def read_nml(case):
    txt = (Path(case)/'input.nml').read_text()
    num = lambda s: [float(v.replace('d', 'e').replace('D', 'e')) for v in re.split(r'[,\s]+', s.strip()) if v]
    ng = num(re.search(r'^\s*ng\s*(?:\(1:3\))?\s*=\s*([^!\n/]*)', txt, re.M | re.I).group(1))[:3]
    l = num(re.search(r'^\s*l\s*(?:\(1:3\))?\s*=\s*([^!\n/]*)', txt, re.M | re.I).group(1))[:3]
    return [int(v) for v in ng], l

# ------------------------------------------------------------------------------------------------ geometry

def boundary_faces(V, F, L, tol):
    """faces lying in one of the domain boundary planes"""
    drop = np.zeros(len(F), bool)
    for ax in range(3):
        c = V[F, ax]
        for p in (0., L[ax]):
            drop |= np.all(np.abs(c - p) < tol, axis=1)
    return drop

def periodic_copies(V, F, L, B):
    """faces shifted by (+-lx, +-ly) whose bounding box reaches into [-B, lx + B] x [-B, ly + B]"""
    T = V[F]; lo, hi = T.min(1), T.max(1)
    Vs, Fs = [V], [F]
    nv = len(V)
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            if sx == 0 and sy == 0: continue
            s = np.array([sx*L[0], sy*L[1], 0.])
            sel = np.all((lo[:, :2] + s[:2] < np.array(L[:2]) + B) & (hi[:, :2] + s[:2] > -B), axis=1)
            if not sel.any(): continue
            used, inv = np.unique(F[sel], return_inverse=True)
            Vs.append(V[used] + s); Fs.append(inv.reshape(-1, 3) + nv); nv += len(used)
    return np.concatenate(Vs), np.concatenate(Fs)

def pseudonormals(V, F):
    """face normals, angle-weighted vertex normals, edge normals (sum of the adjacent faces), face-edge map"""
    e1, e2 = V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]
    FN = np.cross(e1, e2); FN /= np.linalg.norm(FN, axis=1)[:, None]
    VN = np.zeros_like(V)
    for c in range(3):
        a = V[F[:, (c + 1) % 3]] - V[F[:, c]]; b = V[F[:, (c + 2) % 3]] - V[F[:, c]]
        ang = np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.einsum('ij,ij->i', a, b))
        for d in range(3):
            VN[:, d] += np.bincount(F[:, c], ang*FN[:, d], len(V))
    # edge c of a face joins its vertices c and c+1
    a = np.stack([F[:, c] for c in range(3)], 1); b = np.stack([F[:, (c + 1) % 3] for c in range(3)], 1)
    key = np.minimum(a, b).astype(np.int64)*len(V) + np.maximum(a, b)
    _, EMAP = np.unique(key.ravel(), return_inverse=True)
    EMAP = EMAP.reshape(-1, 3)
    EN = np.zeros((EMAP.max() + 1, 3))
    for d in range(3):
        EN[:, d] = np.bincount(EMAP.ravel(), np.repeat(FN[:, d], 3), len(EN))
    return FN, VN, EN, EMAP

def pn_sign(P, C, f, V, F, FN, VN, EN, EMAP, eps=1e-6):
    """sign of P - C along the pseudonormal of the feature of face f that contains the closest point C"""
    A, B_, C_ = V[F[f, 0]], V[F[f, 1]], V[F[f, 2]]
    v0, v1, v2 = B_ - A, C_ - A, C - A
    d00 = np.einsum('ij,ij->i', v0, v0); d01 = np.einsum('ij,ij->i', v0, v1); d11 = np.einsum('ij,ij->i', v1, v1)
    d20 = np.einsum('ij,ij->i', v2, v0); d21 = np.einsum('ij,ij->i', v2, v1)
    den = d00*d11 - d01*d01
    b1 = (d11*d20 - d01*d21)/den; b2 = (d00*d21 - d01*d20)/den
    bc = np.stack([1. - b1 - b2, b1, b2], 1)
    n = FN[f].copy()
    small = bc < eps
    nsmall = small.sum(1)
    e = nsmall == 1                                   # on the edge opposite vertex k: edge (k+1)
    if e.any():
        k = np.argmax(small[e], 1); n[e] = EN[EMAP[f[e], (k + 1) % 3]]
    v = nsmall >= 2                                   # at vertex k
    if v.any():
        k = np.argmax(bc[v], 1); n[v] = VN[F[f[v], k]]
    return np.sign(np.einsum('ij,ij->i', P - C, n))

def surface_samples(V, F, h):
    """points on every face with spacing <= h, and the covering radius r (every surface point is within r)"""
    T = V[F]
    emax = np.max(np.linalg.norm(T - np.roll(T, 1, axis=1), axis=2), axis=1)
    m = np.maximum(1, np.ceil(emax/h)).astype(int)
    pts, r = [T.mean(1)], np.max(np.linalg.norm(T - T.mean(1)[:, None], axis=2)[m == 1], initial=0.)
    for mm in np.unique(m[m > 1]):
        i, j = np.meshgrid(np.arange(mm + 1), np.arange(mm + 1), indexing='ij')
        ok = i + j <= mm
        bc = np.stack([i[ok], j[ok], mm - i[ok] - j[ok]], 1)/mm      # barycentric lattice
        Tm = T[m == mm]
        pts.append(np.einsum('kc,fcd->fkd', bc, Tm).reshape(-1, 3))
        r = max(r, emax[m == mm].max()/mm)
    return np.concatenate(pts), r

# ------------------------------------------------------------------------------------------------ distance

G = {}

def _slab(ks):
    g = G
    X, Y = np.meshgrid(g['xc'], g['yc'], indexing='xy')                  # (ny, nx)
    P = np.stack([np.broadcast_to(X, (len(ks),) + X.shape), np.broadcast_to(Y, (len(ks),) + X.shape),
                  np.broadcast_to(g['zc'][ks][:, None, None], (len(ks),) + X.shape)], -1).reshape(-1, 3)
    out = np.full(len(P), np.nan)                                        # nan = outside the band
    if g['B'] < np.inf:
        dl, _ = g['tree'].query(P, distance_upper_bound=g['B'] + g['r'])
        cand = np.nonzero(dl < g['B'] + g['r'])[0]
    else:
        cand = np.arange(len(P))
    if len(cand) == 0:
        return ks, out
    Pc = P[cand]
    lo, hi = Pc.min(0) - g['B'], Pc.max(0) + g['B']
    crop = np.nonzero(np.all((g['flo'] <= hi) & (g['fhi'] >= lo), axis=1))[0]
    if len(crop) == 0:
        return ks, out
    Fc = g['F'][crop]
    tree = igl.AABB(); tree.init(g['V'], Fc)
    d2, I, C = tree.squared_distance(g['V'], Fc, Pc)
    d = np.sqrt(d2)
    band = d < g['B']
    f = crop[I[band]]
    if g['sign'] == 'pseudonormal':
        s = pn_sign(Pc[band], C[band], f, g['V'], g['F'], g['FN'], g['VN'], g['EN'], g['EMAP'])*g['flip']
    else:                       # independent of the orientation; periodic images of the input mesh included
        Q, w = Pc[band], np.zeros(band.sum())
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                s = np.array([sx*g['L'][0], sy*g['L'][1], 0.])
                if np.all((Q.min(0) + s <= g['bb0'][1]) & (Q.max(0) + s >= g['bb0'][0])):
                    w = np.maximum(w, np.abs(igl.fast_winding_number(g['V0'], g['F0'], Q + s)))
        s = np.where(w > 0.5, -1., 1.)
    out[cand[band]] = s*d[band]
    return ks, out

def flood_sign(phi, B):
    """phi = +-B outside the band: sign of the connected regions (periodic in x, y) voted by the band cells"""
    out = np.isnan(phi)
    if not out.any():
        return phi, 0, 0.
    lab, nl = ndimage.label(out)                                        # 6-connected
    # periodic merge in x and y
    i, j = [], []
    for a, b in ((lab[:, :, 0], lab[:, :, -1]), (lab[:, 0, :], lab[:, -1, :])):
        m = (a > 0) & (b > 0); i.append(a[m]); j.append(b[m])
    i, j = np.concatenate(i), np.concatenate(j)
    g = sparse.coo_matrix((np.ones(len(i)), (i, j)), shape=(nl + 1, nl + 1))
    nc, comp = connected_components(g, directed=False)
    s = np.where(out, 0, np.sign(phi)).astype(np.int8)
    s[(~out) & (phi == 0)] = -1                                         # on the surface: solid
    vp, vm = np.zeros(nc), np.zeros(nc)
    for ax, per in ((2, True), (1, True), (0, False)):
        for sh in (1, -1):
            nb = np.roll(s, sh, ax)
            if not per:
                idx = [slice(None)]*3; idx[ax] = 0 if sh == 1 else -1; nb[tuple(idx)] = 0
            m = out & (nb != 0)
            c = comp[lab[m]]
            vp += np.bincount(c, nb[m] > 0, nc); vm += np.bincount(c, nb[m] < 0, nc)
    if np.any((vp + vm)[np.unique(comp[lab[out]])] == 0):
        sys.exit('a region outside the band has no band neighbour (no surface in the domain?)')
    sgn = np.where(vp >= vm, 1., -1.)
    conflict = np.minimum(vp, vm).sum()/max((vp + vm).sum(), 1)
    phi[out] = sgn[comp[lab[out]]]*B
    return phi, nc, conflict

# ------------------------------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case'); ap.add_argument('stl', nargs='+')
    ap.add_argument('-o', '--out', help='output file (default <case>/data/sdfp.bin)')
    ap.add_argument('--uniform', action='store_true', help='uniform z grid')
    ap.add_argument('--scale', type=float, default=1.); ap.add_argument('--shift', type=float, nargs=3, default=[0., 0., 0.])
    ap.add_argument('--mirror', action='store_true', help='add the reflection z -> lz - z (upper wall)')
    ap.add_argument('--flip', action='store_true', help='normals point into the solid')
    ap.add_argument('--keep-boundary-faces', action='store_true', help='do not drop faces on the domain boundary')
    ap.add_argument('--band', type=float, default=6., help='exact distance for |phi| < band x max(dx, dy, dz); 0 = everywhere')
    ap.add_argument('--sign', choices=['pseudonormal', 'winding'], default='pseudonormal',
                    help='winding: generalised winding number of the (closed) input mesh, robust to holes')
    ap.add_argument('--tol', type=float, default=1e-6, help='weld/boundary tolerance relative to max(l)')
    ap.add_argument('--np', type=int, default=os.cpu_count())
    ap.add_argument('--slab', type=int, default=2_000_000, help='grid points per parallel task')
    a = ap.parse_args()
    t0 = time.time()
    (nx, ny, nz), L = read_nml(a.case)
    zc, dzf, gsrc = read_grid(a.case, nz, L[2], a.uniform)
    dx, dy = L[0]/nx, L[1]/ny
    xc, yc = (np.arange(nx) + 0.5)*dx, (np.arange(ny) + 0.5)*dy
    tol = a.tol*max(L)
    tri = np.concatenate([read_stl(f) for f in a.stl])*a.scale + np.array(a.shift)
    if a.mirror:
        m = tri[:, ::-1].copy(); m[..., 2] = L[2] - m[..., 2]; tri = np.concatenate([tri, m])
    zmin, zmax = tri[..., 2].min(), tri[..., 2].max()
    print(f'{len(tri)} triangles, bounding box {tri.reshape(-1, 3).min(0)} - {tri.reshape(-1, 3).max(0)} '
          f'({time.time() - t0:.1f} s)')
    # band: cells up to `band` x the largest spacing near the surface (dz where the surface lies, +- band x dx)
    hxy = max(dx, dy)
    near = (zc > zmin - a.band*hxy) & (zc < zmax + a.band*hxy)
    h = max(hxy, dzf[near].max() if near.any() else dzf.max())
    B = a.band*h if a.band > 0 else np.inf
    print(f'grid {nx} x {ny} x {nz} ({gsrc}), l = {L}, band B = {B:.4g} ({a.band:g} x {h:.4g})')
    V0, F0 = weld(tri, tol); del tri
    F0 = F0[np.linalg.norm(np.cross(V0[F0[:, 1]] - V0[F0[:, 0]], V0[F0[:, 2]] - V0[F0[:, 0]]), axis=1) > (tol*tol)]
    drop = np.zeros(len(F0), bool) if a.keep_boundary_faces else boundary_faces(V0, F0, L, tol)
    V, F = periodic_copies(V0, F0[~drop], L, B if B < np.inf else max(L))
    V, F = weld_near(V, F, tol)
    print(f'{len(V0)} vertices, {len(F0)} faces, {drop.sum()} on the domain boundary dropped; with periodic copies '
          f'{len(F)} faces ({time.time() - t0:.1f} s)')
    FN, VN, EN, EMAP = pseudonormals(V, F)
    T = V[F]
    G.update(V=V, F=F, flo=T.min(1), fhi=T.max(1), FN=FN, VN=VN, EN=EN, EMAP=EMAP, xc=xc, yc=yc, zc=zc, B=B,
             sign=a.sign, flip=-1. if a.flip else 1., V0=V0, F0=F0, L=L, bb0=(V0.min(0), V0.max(0)))
    if B < np.inf:
        S, r = surface_samples(V, F, B/6)
        G.update(tree=cKDTree(S), r=r)
        print(f'{len(S)} surface samples for the band test (covering radius {r:.3g}) ({time.time() - t0:.1f} s)')
    del T
    kper = max(1, a.slab//(nx*ny))
    tasks = [np.arange(k, min(k + kper, nz)) for k in range(0, nz, kper)]
    phi = np.empty((nz, ny, nx))
    with Pool(a.np) as pool:
        for ks, o in pool.imap_unordered(_slab, tasks):
            phi[ks] = o.reshape(len(ks), ny, nx)
    print(f'distance in the band: {np.sum(~np.isnan(phi))} cells ({time.time() - t0:.1f} s)')
    phi, nc, conflict = flood_sign(phi, B)
    out = Path(a.out) if a.out else Path(a.case)/'data'/'sdfp.bin'
    out.parent.mkdir(parents=True, exist_ok=True)
    np.ascontiguousarray(phi).tofile(out)
    # checks
    gz = np.gradient(phi, zc, axis=0)
    gy = (np.roll(phi, -1, 1) - np.roll(phi, 1, 1))/(2*dy); gx = (np.roll(phi, -1, 2) - np.roll(phi, 1, 2))/(2*dx)
    gn = np.sqrt(gx**2 + gy**2 + gz**2)
    near = np.abs(phi) < 0.5*min(B, 1e300)
    near[[0, -1]] = False
    q = np.quantile(gn[near], [0.01, 0.5, 0.99]) if near.any() else [np.nan]*3
    solid = (phi <= 0).mean(axis=(1, 2))
    print(f'regions outside the band: {nc}, sign votes against the majority: {100*conflict:.3g} %')
    print(f'|grad phi| for |phi| < B/2: 1/50/99 % quantiles {q[0]:.3f} / {q[1]:.4f} / {q[2]:.3f}')
    print(f'solid fraction: first cell {solid[0]:.3f}, last cell {solid[-1]:.3f}, total {(phi <= 0).mean():.4f}')
    print(f'wrote {out} ({time.time() - t0:.1f} s)')

if __name__ == '__main__':
    main()
