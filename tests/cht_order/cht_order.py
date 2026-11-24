#!/usr/bin/env python3
"""
Order of accuracy of the steady conjugate conduction operator of the CHT scalars (`scal_cht`).

Each case is a small box, periodic in x and y over [0,1]^2, z in [0,LZ], with ONE interface
z = g(x,y): solid (k_s = ks_kf) below, fluid (k_f = 1) above. CaNS itself writes the sub-cell
geometry it uses (`&ibm dump_ibm_geometry = T`: the link fractions thx, thy, thz and the
volume fractions phi, after one time step); this script assembles from them the operator that
`scal_cht` applies,

    L T = sum over the six faces of +/- k_face (T_N - T_P)/h_c / h_f,
    k_face = 1/(th/k_f + (1-th)/k_s),                     (scal_cht, cht.f90 `kface`)

and solves the steady problem  L T + phi q_f + (1-phi) q_s = 0  (the capacity cancels) against a
manufactured conjugate solution (as theory/test-code/tangential/cht3d.py of AMPHIBIOUS):

    T_s = F,   T_f = F + c(x,y) P + d(x,y) P^2,   P = z - g(x,y),
    c = (kappa - 1) grad F . grad P / |grad P|^2 at z = g,   d = 2 + sin(2 pi x),

continuous with a continuous normal flux at the interface, varying along it, with sources in both
materials (q_s = -kappa lap F, q_f = -lap T_f) and a fluid-only jump of d2T/dn2. The box faces
z = 0 and z = LZ carry the exact solution in the ghost cells (Dirichlet data; harness choice, so
that the error is that of the interior operator). The error is measured at every cell centre,
fluid and solid, with the material of the code's mask (level set > 0: fluid).

Families (wall g):  flat-on   z_w = 0.25, on a cell face of every grid
                    flat-off  z_w = 0.25 + dz/3, 1/3 of a cell above a face (cut fixed under refinement)
                    corr      0.3 + 0.05 sin(2 pi x)           (corrugated, y-invariant)
                    curved    0.3 + 0.06 cos(2 pi x) cos(2 pi y)
                    curved2   curved + 0.03 sin(2 pi (x + y))  (generic; 'sdf' only, the wall of cht3d.py)
Geometry paths:     sdf       exact signed distance in data/sdfp.bin (Newton projection onto the wall)
                    analytic  the same wall through &ibm ibm_shape ('sinusoid' / 'macdonald'), bisection
Grids:              aspect dx/dz = 1: N = 16, 24, 32, 48, 64;  aspect 4: N = 8, 12, 16, 24, 32  (N = nx = ny,
                    nz = LZ*aspect*N)

Usage (from anywhere; cans defaults to <repo>/build/cans, work dir to <repo>/run/cht_order):

    python3 cht_order.py run      [--paths sdf,analytic] [--families ...] [--kappas 0.5,2,10]
                                  [--aspects 1,4] [--jobs 12] [--quick] [--cans PATH] [--work DIR]
    python3 cht_order.py fidelity [--cans PATH] [--work DIR]
    python3 cht_order.py report   [--work DIR]

`run` writes <work>/results.json and prints L2/Linf errors and orders; `fidelity` checks that the
operator assembled here IS the one the code applies: one explicit RK3 step of CaNS from a given
field (zero velocity, sources in both materials, capacity ratio 2) against the same step taken here.

--scheme tang (default step1): the scalar runs with `is_cht_rkl2 = T` and `is_cht_tang = T`; CaNS then also
writes its tangential face list (data/cht_tang_001.bin: per face the lower cell, the axis and the 18
coefficients of (k_area-k_face)*s_t on the 2x3x3 block), which is added here to the series operator exactly
as `cht_tang_apply` adds it. With `fidelity` the CaNS step is RKL2 (two half steps around the convection,
which vanishes at zero velocity), taken here with the same stage counts.
"""
import argparse, json, math, os, re, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import sympy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
LZ = 0.75                       # box height; x, y in [0, 1]
ZREF = 0.4                      # centre of the manufactured F
GRIDS = {1: [16, 24, 32, 48, 64], 4: [8, 12, 16, 24, 32]}
QUICK = {1: [16, 24, 32], 4: [8, 12, 16]}
FAMILIES = ['flat-on', 'flat-off', 'corr', 'curved', 'curved2']
SCHEME = 'step1'                # or 'tang' (--scheme)
TP = 2*math.pi


# ----------------------------------------------------------------------------- the wall
def wall(family, dz):
    """sympy expression of z = g(x,y) (dz: the cell height, which places the flat walls)"""
    x, y = sp.symbols('x y', real=True)
    if family == 'flat-on':                 # 0.25 lies on a cell face of every grid below
        assert abs(0.25/dz - round(0.25/dz)) < 1e-12
        return sp.Float(0.25, 17) + 0*x
    if family == 'flat-off':                # the cut held at 1/3 of a cell; the wall moves by dz/3
        return sp.Float(0.25 + dz/3, 17) + 0*x
    if family == 'corr':
        return sp.Rational(3, 10) + sp.Rational(5, 100)*sp.sin(2*sp.pi*x)
    if family == 'curved':
        return sp.Rational(3, 10) + sp.Rational(6, 100)*sp.cos(2*sp.pi*x)*sp.cos(2*sp.pi*y)
    if family == 'curved2':
        return (sp.Rational(3, 10) + sp.Rational(6, 100)*sp.cos(2*sp.pi*x)*sp.cos(2*sp.pi*y)
                + sp.Rational(3, 100)*sp.sin(2*sp.pi*(x + y)))
    raise ValueError(family)


def analytic_ibm(family, dz):
    """the &ibm entries that make the analytic geometry of CaNS represent the same wall"""
    if family in ('flat-on', 'flat-off'):
        z0 = float(wall(family, dz).subs({}))
        return dict(ibm_shape='sinusoid', l0=z0, amp=(0., 0.), nw=(0, 0), phase=(0., 0.))
    if family == 'corr':      # l0 + A/2 (1 + sin(2 pi x)) = 0.3 + 0.05 sin(2 pi x)
        return dict(ibm_shape='sinusoid', l0=0.25, amp=(0.1, 0.), nw=(1, 0), phase=(0., 0.))
    if family == 'curved':    # l0 + A (1 + cos cos)
        return dict(ibm_shape='macdonald', l0=0.24, amp=(0.06, 0.), nw=(1, 1), phase=(0., 0.))
    return None


# ----------------------------------------------------------------------------- exact solution
class Exact:
    def __init__(self, g, kappa):
        x, y, z = sp.symbols('x y z', real=True)
        eta = z - ZREF
        F = ((1 + 2*eta**2)*sp.sin(2*sp.pi*y + sp.Rational(3, 10))*sp.cos(2*sp.pi*x) + 3*eta**2
             + sp.Rational(1, 2)*sp.cos(4*sp.pi*x)*eta**2)
        P = z - g
        grad = lambda f: sp.Matrix([sp.diff(f, v) for v in (x, y, z)])
        c = ((kappa - 1)*grad(F).dot(grad(P))/grad(P).dot(grad(P))).subs(z, g)
        d = 2 + sp.sin(2*sp.pi*x)
        Tf = F + c*P + d*P**2
        lap = lambda f: sum(sp.diff(f, v, 2) for v in (x, y, z))
        L = lambda e: sp.lambdify((x, y, z), e, 'numpy')
        self.Ts, self.Tf = L(F), L(Tf)
        self.qs, self.qf = L(-kappa*lap(F)), L(-lap(Tf))
        self.g = sp.lambdify((x, y), g, 'numpy')
        self.gx = sp.lambdify((x, y), sp.diff(g, x), 'numpy')
        self.gy = sp.lambdify((x, y), sp.diff(g, y), 'numpy')
        self.gxx = sp.lambdify((x, y), sp.diff(g, x, 2), 'numpy')
        self.gxy = sp.lambdify((x, y), sp.diff(g, x, y), 'numpy')
        self.gyy = sp.lambdify((x, y), sp.diff(g, y, 2), 'numpy')
        # the two interface conditions, on a sample of interface points
        xs, ys = np.random.default_rng(0).random((2, 64))
        zs = self.g(xs, ys) + 0*xs
        jT = L(Tf - F)(xs, ys, zs)
        jq = L(grad(Tf).dot(grad(P)) - kappa*grad(F).dot(grad(P)))(xs, ys, zs)
        assert np.abs(jT).max() < 1e-12 and np.abs(jq).max() < 1e-9, (np.abs(jT).max(), np.abs(jq).max())

    def T(self, x, y, z, fluid):
        return np.where(fluid, self.Tf(x, y, z), self.Ts(x, y, z)) + 0*x

    def q(self, x, y, z, phi):
        return phi*self.qf(x, y, z) + (1 - phi)*self.qs(x, y, z) + 0*x


def sdf(ex, X, Y, Z):
    """exact signed distance to z = g(x,y), > 0 above (fluid): Newton projection onto the wall from
    the best of a local sample (the vertical distance bounds the true one, so the foot lies within
    that radius); g periodic with period 1"""
    gv = lambda u, v: ex.g(u, v) + 0*u
    dv = Z - gv(X, Y)
    smp = np.linspace(0, 1, 7)
    flat = all(np.all(f(smp[:, None], smp[None, :]) + 0*np.ones((7, 7)) == 0) for f in (ex.gx, ex.gy))
    if flat:
        return dv
    x, y, z = X.ravel(), Y.ravel(), Z.ravel()
    r = np.abs(dv.ravel())
    m = 9
    s = np.linspace(-1, 1, m)
    best = np.full(x.size, np.inf); u0 = x.copy(); v0 = y.copy()
    for a in s:
        for b in s:
            u, v = x + a*r, y + b*r
            dd = (u - x)**2 + (v - y)**2 + (gv(u, v) - z)**2
            k = dd < best
            best[k] = dd[k]; u0[k] = u[k]; v0[k] = v[k]
    u, v = u0, v0
    D = lambda u, v: 0.5*((u - x)**2 + (v - y)**2 + (gv(u, v) - z)**2)
    for it in range(200):           # Newton on the squared distance, backtracking where it does not descend
        G = gv(u, v) - z
        gx, gy = ex.gx(u, v) + 0*u, ex.gy(u, v) + 0*u
        gxx, gxy, gyy = ex.gxx(u, v) + 0*u, ex.gxy(u, v) + 0*u, ex.gyy(u, v) + 0*u
        r1 = (u - x) + G*gx
        r2 = (v - y) + G*gy
        a11 = 1 + gx*gx + G*gxx; a12 = gx*gy + G*gxy; a22 = 1 + gy*gy + G*gyy
        det = a11*a22 - a12*a12
        pd = (det > 0) & (a11 > 0)
        du = np.where(pd, (a22*r1 - a12*r2)/np.where(pd, det, 1.), r1)
        dw = np.where(pd, (a11*r2 - a12*r1)/np.where(pd, det, 1.), r2)
        D0, lam = D(u, v), np.ones_like(u)
        for _ in range(40):
            worse = D(u - lam*du, v - lam*dw) > D0*(1 + 1e-12) + 1e-30   # round-off of D: accept
            if not worse.any():
                break
            lam = np.where(worse, 0.5*lam, lam)
        u, v = u - lam*du, v - lam*dw
        if max(np.abs(lam*du).max(), np.abs(lam*dw).max()) < 1e-16:
            break
    G = gv(u, v) - z
    res = np.maximum(np.abs((u - x) + G*ex.gx(u, v)), np.abs((v - y) + G*ex.gy(u, v)))
    dist = np.sqrt((u - x)**2 + (v - y)**2 + G**2)
    # the code reads the level set beyond its sign only within ~3 half-diagonals of the wall (< 0.3 here)
    assert np.all((res < 1e-11) | (dist > 0.3)), 'Newton projection did not converge near the wall'
    assert np.all(dist <= r + 1e-14)
    return (np.sign(dv.ravel())*dist).reshape(X.shape)


# ----------------------------------------------------------------------------- the CaNS case
NML = """&dns
ng(1:3) = {nx}, {ny}, {nz}
l(1:3) = 1., 1., {lz!r}
gtype = 1, gr = 0.
cfl = 0.95, dtmax = 1.e5, dt_f = {dt!r}
visci = 1.
inivel = 'zer'
is_wallturb = F
nstep = {nstep}, time_max = 1.e9, tw_max = 1.e9
stop_type(1:3) = T, F, F
restart = {restart}, is_overwrite_save = T, nsaves_max = 0
icheck = 1000, iout0d = 1000, iout1d = 100000000, iout2d = 0, iout3d = 0, isave = 100000000
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3)   =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = F, F, F
velf(1:3) = 0., 0., 0.
iaverage = 0, average_start_time = 1.e9
gacc(1:3) = 0., 0., 0.
nscal = 1
dims(1:2) = 1, 1, ipencil_axis = 3
/

&scalar
beta = 0.
iniscal(:)  = 'zer'
alphai(:)   = 1.
cbcscal(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
bcscal(0:1,1:3,1)  =  0.,0.,   0.,0.,   {bc0!r},{bc1!r}
ssource(:)    = {ssource!r}
ssource_s(:)  = {ssource_s!r}
is_sforced(:) = F
scalf(:)      = 0.
is_cht(:)     = T
ks_kf(:)      = {kappa!r}
cs_cf(:)      = {cs!r}
cbcscal_ibm(:) = 'D'
bcscal_ibm(:)  = 0.
const_scal = T
is_boussinesq_buoyancy = F
{scheme}/

&cudecomp
cudecomp_t_comm_backend = 1, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = F
cudecomp_h_comm_backend = 2, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = F
/

&numerics
is_impdiff = F, is_impdiff_1d = F
is_poisson_pcr_tdma = F
is_fast_mom_kernels = F
/

&other_options
is_debug = F, is_timing = F
/

&io
io_backend = 'mpiio'
/

&ibm
is_ibm = T
geometry = '{geometry}'
sdf_file = 'data/sdfp.bin'
ibm_iout_wall = 0
dump_ibm_geometry = T
{analytic}/
"""


class Case:
    def __init__(self, path, family, kappa, aspect, N, work):
        self.path, self.family, self.kappa, self.aspect, self.N = path, family, kappa, aspect, N
        self.nx = self.ny = N
        self.nz = int(round(LZ*aspect*N))
        self.dx, self.dy, self.dz = 1/N, 1/N, LZ/self.nz
        assert abs(self.dz*aspect - self.dx) < 1e-14
        self.g = wall(family, self.dz)
        self.scheme = SCHEME
        pre = '' if SCHEME == 'step1' else SCHEME + '-'
        self.dir = Path(work)/f'{pre}{path}-{family}-k{kappa:g}-a{aspect}-n{N}'
        self.xc = (np.arange(self.nx) + 0.5)*self.dx
        self.yc = (np.arange(self.ny) + 0.5)*self.dy
        self.zc = (np.arange(self.nz) + 0.5)*self.dz        # checked against data/grid.out

    def mesh(self):
        return np.meshgrid(self.xc, self.yc, self.zc, indexing='ij')   # (nx, ny, nz)

    def write(self, ex, nstep=1, restart='F', dt=1.e-12, ssource=0., ssource_s=0., cs=1., bc=(0., 0.)):
        d = self.dir
        (d/'data').mkdir(parents=True, exist_ok=True)
        an = ''
        if self.path == 'analytic':
            a = analytic_ibm(self.family, self.dz)
            an = ("ibm_shape = '{s}'\nibm_direction(0,3) = T\nl_0(0,3) = {l0!r}\n"
                  "amp_l(0,1) = {a0!r}, amp_l(0,2) = {a1!r}\nn_wave(0,1) = {n0}, n_wave(0,2) = {n1}\n"
                  "phase_l(0,1) = {p0!r}, phase_l(0,2) = {p1!r}\n").format(
                      s=a['ibm_shape'], l0=a['l0'], a0=a['amp'][0], a1=a['amp'][1], n0=a['nw'][0],
                      n1=a['nw'][1], p0=a['phase'][0], p1=a['phase'][1])
        (d/'input.nml').write_text(NML.format(nx=self.nx, ny=self.ny, nz=self.nz, lz=LZ, dt=dt, nstep=nstep,
                                              restart=restart, kappa=float(self.kappa), cs=float(cs),
                                              ssource=float(ssource), ssource_s=float(ssource_s),
                                              bc0=float(bc[0]), bc1=float(bc[1]),
                                              geometry=self.path, analytic=an,
                                              scheme='' if self.scheme == 'step1' else
                                              'is_cht_rkl2(:) = T\nis_cht_tang(:) = T\n'))
        X, Y, Z = self.mesh()
        self.ps = sdf(ex, X, Y, Z)                       # also the material of the analytic path
        if self.path == 'sdf':
            np.ascontiguousarray(self.ps.transpose(2, 1, 0), dtype='<f8').tofile(d/'data'/'sdfp.bin')

    def run(self, cans):
        log = open(self.dir/'run.log', 'w')
        r = subprocess.run(['mpirun', '-np', '1', str(cans)], cwd=self.dir, stdout=log, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            raise RuntimeError(f'CaNS failed in {self.dir} (see run.log)')

    def geometry(self):
        """thx, thy, thz, phi as CaNS wrote them, arrays (nx, ny, nz)"""
        n = self.nx*self.ny*self.nz
        a = np.fromfile(self.dir/'data'/'ibm_geometry.bin', dtype='<f8')
        assert a.size == 8*n, (a.size, n)
        f = [a[(4 + m)*n:(5 + m)*n].reshape(self.nz, self.ny, self.nx).transpose(2, 1, 0) for m in range(4)]
        g = np.loadtxt(self.dir/'data'/'grid.out')
        assert np.allclose(g[1:self.nz + 1, 2], self.zc, rtol=0, atol=1e-14), 'unexpected z grid'
        return f


    def tang_faces(self):
        """the tangential face list CaNS wrote (scheme 'tang'): [(i, j, k, axis, c[18])], 1-based cells"""
        p = self.dir/'data'/'cht_tang_001.bin'
        if self.scheme == 'step1' or not p.exists():
            return []
        raw = p.read_bytes()
        nf = int(np.frombuffer(raw[:4], '<i4')[0])
        rec = np.dtype([('ijkd', '<i4', 4), ('c', '<f8', 18)])
        a = np.frombuffer(raw[16:], rec, count=nf)
        return [(*map(int, r['ijkd']), r['c'].copy()) for r in a]


def block_offset(m, d):
    """the m-th cell (m = 1..18) of the block of a face along axis d (1..3), as `block_offset` in cht_tang.f90"""
    a = (m - 1)//9; r = (m - 1) % 9; b = r//3 - 1; cc = r % 3 - 1
    return {1: (a, b, cc), 2: (b, a, cc), 3: (b, cc, a)}[d]


def tang_matrix(c, faces):
    """the tangential term as `cht_tang_apply` adds it (without the capacity): a face's
    st = sum_m c_m T(block_m) enters its lower cell with +st/W and its upper cell with -st/W, each only
    where that cell is local (a periodic face appears once for each of its two cells)"""
    nx, ny, nz = c.nx, c.ny, c.nz
    idx = np.arange(nx*ny*nz).reshape(nx, ny, nz)
    W = {1: c.dx, 2: c.dy, 3: c.dz}
    rows, cols, vals = [], [], []
    for (i, j, k, d, cf) in faces:
        cl = [idx[(i + o[0] - 1) % nx, (j + o[1] - 1) % ny, k + o[2] - 1] for o in (block_offset(m, d) for m in range(1, 19))]
        for side, sgn in ((0, 1.), (1, -1.)):
            ci, cj, ck = i + (side if d == 1 else 0), j + (side if d == 2 else 0), k + (side if d == 3 else 0)
            if not (1 <= ci <= nx and 1 <= cj <= ny and 1 <= ck <= nz):
                continue
            r = idx[ci - 1, cj - 1, ck - 1]
            rows.extend([r]*18); cols.extend(cl); vals.extend(sgn*cf/W[d])
    return sps.csr_matrix((vals, (rows, cols)), shape=(idx.size, idx.size))


# ----------------------------------------------------------------------------- the operator
def kface(th, kappa):
    """the code's face conductivity, alphaf = 1 (scal_cht: 1./(th*rf + (1.-th)*rs))"""
    return 1./(th*1. + (1. - th)*(1./kappa))


def assemble(c, th, kappa):
    """L (sparse, nx*ny*nz square, interior couplings only) and, for the ghost cells below k = 0
    and above k = nz-1, the coefficients that multiply the ghost values; the face conductivities of
    the two box faces are those of the adjacent cells' material (asserted: no interface there)."""
    thx, thy, thz, phi = th
    nx, ny, nz = c.nx, c.ny, c.nz
    idx = np.arange(nx*ny*nz).reshape(nx, ny, nz)
    rows, cols, vals = [], [], []
    def link(kf, h_c, h_f, a, b):          # flux kf*(T_b - T_a)/h_c, into a with +, out of b with -
        w = kf/(h_c*h_f)
        rows.extend([a.ravel(), a.ravel(), b.ravel(), b.ravel()])
        cols.extend([b.ravel(), a.ravel(), a.ravel(), b.ravel()])
        vals.extend([w.ravel(), -w.ravel(), w.ravel(), -w.ravel()])
    link(kface(thx, kappa), c.dx, c.dx, idx, np.roll(idx, -1, 0))
    link(kface(thy, kappa), c.dy, c.dy, idx, np.roll(idx, -1, 1))
    link(kface(thz[:, :, :-1], kappa), c.dz, c.dz, idx[:, :, :-1], idx[:, :, 1:])
    L = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(idx.size, idx.size))
    fl = c.ps > 0
    assert not fl[:, :, 0].any() and fl[:, :, -1].all(), 'the interface must not reach the box faces'
    assert np.all(thz[:, :, -1] == 1.0), 'top boundary face not in the fluid'
    kb, kt = kappa, 1.0                  # box faces: in the solid (bottom), in the fluid (top)
    wb, wt = kb/(c.dz*c.dz), kt/(c.dz*c.dz)
    dg = np.zeros(idx.size)
    dg[idx[:, :, 0].ravel()] -= wb
    dg[idx[:, :, -1].ravel()] -= wt
    L = (L + sps.diags(dg)).tocsr()
    faces = c.tang_faces()
    if faces:
        L = (L + tang_matrix(c, faces)).tocsr()
    return L, wb, wt


def check_links(c, th):
    """consistency of the dumped link fractions with the material of the cells (sign of the level set):
    a link is cut (0 < th < 1) iff its two cells differ; returns the number of mismatches"""
    fl = c.ps > 0
    bad = 0
    for ax, t in enumerate(th[:3]):
        if ax == 2:
            a, b, t = fl[:, :, :-1], fl[:, :, 1:], t[:, :, :-1]
        else:
            a, b = fl, np.roll(fl, -1, ax)
        cut = (t > 0) & (t < 1)
        bad += int(np.count_nonzero(cut != (a != b)))
        bad += int(np.count_nonzero((a == b) & (t != np.where(a, 1., 0.))))
    return bad


def linsolve(A, b):
    """direct (minimum degree on A^T+A) below 30000 unknowns; above, BiCGSTAB with an incomplete LU,
    to a relative residual of 1e-13, falling back to the direct solve"""
    A = A.tocsc()
    if A.shape[0] > 30000:
        ilu = spla.spilu(A, drop_tol=1e-4, fill_factor=10)
        x, info = spla.bicgstab(A, b, M=spla.LinearOperator(A.shape, ilu.solve), tol=1e-13, atol=0, maxiter=2000)
        if info == 0 and np.abs(A @ x - b).max() <= 1e-11*np.abs(b).max():
            return x
    return spla.spsolve(A, b, permc_spec='MMD_AT_PLUS_A')


def solve_case(c, ex):
    th = c.geometry()
    nbad = check_links(c, th)
    L, wb, wt = assemble(c, th, c.kappa)
    X, Y, Z = c.mesh()
    fl = c.ps > 0
    Te = ex.T(X, Y, Z, fl)
    rhs = -ex.q(X, Y, Z, th[3]).ravel()
    zg0, zg1 = -0.5*c.dz, LZ + 0.5*c.dz
    Xb, Yb = X[:, :, 0], Y[:, :, 0]
    idx = np.arange(c.nx*c.ny*c.nz).reshape(c.nx, c.ny, c.nz)
    rhs[idx[:, :, 0].ravel()] -= wb*ex.Ts(Xb, Yb, zg0 + 0*Xb).ravel()
    rhs[idx[:, :, -1].ravel()] -= wt*ex.Tf(Xb, Yb, zg1 + 0*Xb).ravel()
    T = linsolve(L, rhs).reshape(Te.shape)
    e = T - Te
    cut = (th[3] > 1e-9) & (th[3] < 1 - 1e-9)        # cut cells (phi of a face-aligned wall is 1 - O(1e-16))
    # conservation: the volume-weighted operator has zero column sums (except at the box faces)
    cs = np.abs(np.asarray(L.sum(axis=0)).ravel()).reshape(Te.shape)[:, :, 1:-1].max()/np.abs(L.diagonal()).max()
    return dict(path=c.path, family=c.family, kappa=c.kappa, aspect=c.aspect, N=c.N, nz=c.nz, scheme=c.scheme,
                tang_faces=len(c.tang_faces()),
                h=max(c.dx, c.dz), L2=float(np.sqrt(np.mean(e**2))), Linf=float(np.abs(e).max()),
                Linf_fluid=float(np.abs(e[fl]).max()), Linf_solid=float(np.abs(e[~fl]).max()),
                cut_cells=int(cut.sum()), link_mismatch=nbad, colsum=float(cs))


def one(args):
    path, family, kappa, aspect, N, work, cans, scheme = args
    global SCHEME
    SCHEME = scheme
    c = Case(path, family, kappa, aspect, N, work)
    ex = Exact(c.g, kappa)
    t0 = time.time()
    if not (c.dir/'data'/'ibm_geometry.bin').exists():
        c.write(ex)
        c.run(cans)
    else:
        X, Y, Z = c.mesh(); c.ps = sdf(ex, X, Y, Z)
    r = solve_case(c, ex)
    r['secs'] = time.time() - t0
    return r


def orders(res):
    out = {}
    for r in res:
        out.setdefault((r['path'], r['family'], r['kappa'], r['aspect']), []).append(r)
    for k, v in out.items():
        v.sort(key=lambda r: r['N'])
        for a, b in zip(v[:-1], v[1:]):
            q = math.log(a['h']/b['h'])
            for n in ('L2', 'Linf'):
                b['order_' + n] = math.log(a[n]/b[n])/q if a[n] > 0 and b[n] > 0 else float('nan')
    return out


def report(res, f=sys.stdout):
    groups = orders(res)
    print('%-9s %-9s %5s %2s | %-40s | %-23s | %s' % ('path', 'family', 'kappa', 'a', 'L2 by grid', 'orders L2/Linf',
                                                        'finest L2/Linf'), file=f)
    for k in sorted(groups, key=lambda k: (k[0], FAMILIES.index(k[1]), k[3], k[2])):
        v = groups[k]
        print('%-9s %-9s %5g %2d | %-40s | %-23s | %.2f/%.2f' % (
            k[0], k[1], k[2], k[3], ' '.join('%.2e' % r['L2'] for r in v),
            ' '.join('%.2f/%.2f' % (r['order_L2'], r['order_Linf']) for r in v[1:]),
            v[-1].get('order_L2', float('nan')), v[-1].get('order_Linf', float('nan'))), file=f)
    bad = sum(r['link_mismatch'] for r in res)
    print('link/material mismatches: %d;  worst relative column sum: %.1e' % (bad, max(r['colsum'] for r in res)), file=f)


# ----------------------------------------------------------------------------- fidelity check
def fidelity(cans, work, path='sdf', family='curved', kappa=2.0, N=12, aspect=1):
    """one explicit RK3 step of CaNS (zero velocity) against the same step on the operator assembled here"""
    c = Case(path, family, kappa, aspect, N, Path(work)/'fidelity')
    if c.dir.exists():
        shutil.rmtree(c.dir)
    ex = Exact(c.g, kappa)
    dt, ss, sss, cs, bc = 2.e-5, 0.7, -0.3, 2.0, (0.1, -0.2)
    c.write(ex, nstep=1, restart='T', dt=dt, ssource=ss, ssource_s=sss, cs=cs, bc=bc)
    X, Y, Z = c.mesh()
    T0 = ex.T(X, Y, Z, c.ps > 0)
    z = np.zeros((c.nz, c.ny, c.nx))
    for name, a in [('u', z), ('v', z), ('w', z), ('p', z), ('s_001', T0.transpose(2, 1, 0))]:
        np.concatenate([np.ravel(a), [0., 0.]]).astype('<f8').tofile(c.dir/'data'/f'fld_{name}.bin')
    c.run(cans)
    T1c = np.fromfile(c.dir/'data'/'fld_s_001.bin')[:T0.size].reshape(c.nz, c.ny, c.nx).transpose(2, 1, 0)
    u1 = np.fromfile(c.dir/'data'/'fld_u.bin')[:T0.size]
    th = c.geometry()
    L, wb, wt = assemble(c, th, kappa)
    phi = th[3].ravel()
    capi = 1./(1. + (cs - 1.)*(1. - phi))
    src = (ss*phi + sss*(1 - phi))*capi
    idx = np.arange(T0.size).reshape(T0.shape)
    b = np.zeros(T0.size)
    def slope(T):                  # capi*(L T) with the code's Dirichlet ghosts T_g = 2 bc - T_1
        g = np.zeros(T0.size)
        g[idx[:, :, 0].ravel()] = wb*(2*bc[0] - T[idx[:, :, 0].ravel()])
        g[idx[:, :, -1].ravel()] = wt*(2*bc[1] - T[idx[:, :, -1].ravel()])
        return capi*(L @ T + g)
    T = T0.ravel().copy()
    if c.scheme == 'step1':
        rk = [(32/60, 0.), (25/60, -17/60), (45/60, -25/60)]
        old = np.zeros(T0.size)
        for a, bb in rk:
            new = slope(T)
            T = T + a*dt*new + bb*dt*old + (a + bb)*dt*src
            old = new
    else:
        # RKL2, two half steps (the convection between them vanishes at zero velocity); the stage count
        # from the Gershgorin bound the code printed
        G = float(re.search(r'Gershgorin bound of the conduction operator\s+(\S+)', (c.dir/'run.log').read_text()).group(1))
        s = 2
        while 0.5*(s*s + s - 2) < 0.5*dt*G:
            s += 1
        M = lambda Y: slope(Y) + src
        b = lambda j: 1/3 if j < 2 else (j*j + j - 2)/(2*j*(j + 1))
        w1 = 4/(s*s + s - 2); h = 0.5*dt
        for half in range(2):
            Y0 = T.copy(); M0 = M(Y0)
            Ym2, Ym1 = Y0, Y0 + w1/3*h*M0
            for j in range(2, s + 1):
                mu = (2*j - 1)/j*b(j)/b(j - 1); nu = -(j - 1)/j*b(j)/b(j - 2)
                mut = mu*w1; gat = -(1 - b(j - 1))*mut
                Yj = mu*Ym1 + nu*Ym2 + (1 - mu - nu)*Y0 + mut*h*M(Ym1) + gat*h*M0
                Ym2, Ym1 = Ym1, Yj
            T = Ym1
    T1h = T.reshape(T0.shape)
    d = np.abs(T1h - T1c).max()
    step = np.abs(T1h - T0).max()
    r = dict(case=c.dir.name, max_abs_diff=float(d), max_change=float(step), rel_to_change=float(d/step),
             rel_to_T=float(d/np.abs(T0).max()), max_u_after=float(np.abs(u1).max()))
    return r


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['run', 'fidelity', 'report'])
    ap.add_argument('--cans', default=str(REPO/'build'/'cans'))
    ap.add_argument('--work', default=str(REPO/'run'/'cht_order'))
    ap.add_argument('--paths', default='sdf,analytic')
    ap.add_argument('--families', default=','.join(FAMILIES))
    ap.add_argument('--kappas', default='0.5,2,10')
    ap.add_argument('--aspects', default='1,4')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--quick', action='store_true', help='three coarse grids per case')
    ap.add_argument('--out', default='results.json')
    ap.add_argument('--scheme', default='step1', choices=['step1', 'tang'])
    a = ap.parse_args()
    global SCHEME
    SCHEME = a.scheme
    work = Path(a.work); work.mkdir(parents=True, exist_ok=True)
    if a.mode == 'fidelity':
        res = [fidelity(a.cans, work, path=p) for p in a.paths.split(',')]
        for r in res:
            print('fidelity %-30s max|T_harness - T_CaNS| = %.2e  (%.1e of the step change, %.1e of max|T|); '
                  'max|u| after the step %.1e' % (r['case'], r['max_abs_diff'], r['rel_to_change'], r['rel_to_T'],
                                                 r['max_u_after']))
        json.dump(res, open(work/'fidelity.json', 'w'), indent=1)
        return
    if a.mode == 'report':
        report(json.load(open(work/a.out)))
        return
    grids = QUICK if a.quick else GRIDS
    jobs = []
    for p in a.paths.split(','):
        for fam in a.families.split(','):
            if p == 'analytic' and fam == 'curved2':        # not representable by ibm_shape
                continue
            for k in [float(v) for v in a.kappas.split(',')]:
                for asp in [int(v) for v in a.aspects.split(',')]:
                    for N in grids[asp]:
                        jobs.append((p, fam, k, asp, N, work, a.cans, a.scheme))
    t0 = time.time()
    res = []
    with ProcessPoolExecutor(a.jobs) as pool:
        for r in pool.map(one, jobs):
            res.append(r)
            print('%-9s %-9s k %-4g a%d N %3d  L2 %.3e  Linf %.3e  cut cells %6d  mismatches %d  (%.0fs)' % (
                r['path'], r['family'], r['kappa'], r['aspect'], r['N'], r['L2'], r['Linf'], r['cut_cells'],
                r['link_mismatch'], r['secs']), flush=True)
    orders(res)
    json.dump(res, open(work/a.out, 'w'), indent=1)
    report(res)
    print('total %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
