#!/usr/bin/env python3
"""
Body-fitted OpenFOAM (v2412, ESI) reference cases of the three CHT examples: one region per material
(chtMultiRegionSimpleFoam for the steady slab and wavy cases, chtMultiRegionFoam for the cylinder), one run per
conductivity/capacity pair, in <study>/<mesh>/K<k_s/k_f>_C<C_s/C_f>/ with an Allrun script.

  python3 make_cases.py slab     [n]     n = cells in z (default 64, as the CaNS example)
  python3 make_cases.py wavy     [n]     n = cells in x and in z (default 32)
  python3 make_cases.py cylinder [ppd]   ppd = D/dx far from the cylinder (default 16); ppd > 16 runs in parallel

slab, wavy: steady, so only k_s/k_f matters (four runs, C = 1). cylinder: the six (k_s/k_f, C_s/C_f) pairs of the
CaNS example, started from the same initial flow (set_u0.py), t = 0 ... 400, averages from t = 200.
Same non-dimensional problem as CaNS: rho = cp = 1, mu = nu, Pr = 1, constant pressure gradient and heat source
q = 1 in the fluid, T = 0 (OpenFOAM: T0 = 300) at z = 0 and 4. The mesh is 2D (one cell in y, empty).
"""
import sys
from pathlib import Path
import numpy as np
import ofgen as G

HERE = Path(__file__).resolve().parent
LY = 0.1
PAIRS_STEADY = [(1., 1.), (10., 1.), (100., 1.), (0.1, 1.)]
PAIRS_CYL = [(1., 1.), (10., 1.), (100., 1.), (0.1, 1.), (1., 0.1), (1., 10.)]
tag = lambda K, C: f'K{K:g}_C{C:g}'


def layered_mesh(nx, lx, nsol, nflu, zlow, zup, lz, edges_fn=None):
    """3 blocks (lowerSolid, fluid, upperSolid) between z = 0, zlow(x), zup(x), lz"""
    v = [(0., 0.), (lx, 0.), (lx, zlow(lx)), (0., zlow(0.)), (lx, zup(lx)), (0., zup(0.)), (lx, lz), (0., lz)]
    blocks = [('lowerSolid', (0, 1, 2, 3), (nx, nsol), (1, 1)),
              ('fluid', (3, 2, 4, 5), (nx, nflu), (1, 1)),
              ('upperSolid', (5, 4, 6, 7), (nx, nsol), (1, 1))]
    pat = {'left': ('cyclic', [(0, 3), (3, 5), (5, 7)], '        neighbourPatch right;\n'),
           'right': ('cyclic', [(1, 2), (2, 4), (4, 6)], '        neighbourPatch left;\n'),
           'bottom': ('wall', [(0, 1)], ''), 'top': ('wall', [(7, 6)], '')}
    return G.blockmesh(v, blocks, edges_fn(v) if edges_fn else [], pat, LY)


def gen_layered(study, mname, bmd, nu, bforce, nnonorth, info):
    for K, C in PAIRS_STEADY:
        case = HERE/study/mname/tag(K, C)
        G.write_common(case, nu, bforce, {'lowerSolid': (K*nu, C), 'upperSolid': (K*nu, C)}, True, nnonorth,
                       G.control('chtMultiRegionSimpleFoam', 100000, 1, 5000))
        G.put(case/'system'/'blockMeshDict', 'dictionary', 'blockMeshDict', bmd)
        G.fields_fluid(case, None)
        for r in ('lowerSolid', 'upperSolid'):
            G.fields_solid(case, r, '"(bottom|top)"')
        G.write_allrun(case, 'chtMultiRegionSimpleFoam', ['fluid', 'lowerSolid', 'upperSolid'])
        (case/'case.txt').write_text(f'study = {study}\nmesh = {mname}\nks_kf = {K}\ncs_cf = {C}\nbforce = {bforce}\nnu = {nu}\n{info}')
    return HERE/study/mname


def slab(n=64):
    """flat slabs of thickness H = 1.02 (as the CaNS example), U_max = 1, Re = 100 based on the half gap d = 0.98"""
    H, LZ, nu = 1.02, 4., 0.01
    d = LZ/2 - H; bforce = 2*nu*1.0/d**2
    nsol = int(round(n*H/LZ)); nflu = n - 2*nsol
    bmd = layered_mesh(2, 0.5, nsol, nflu, lambda x: H, lambda x: LZ - H, LZ)
    return gen_layered('slab', f'n{n:03d}', bmd, nu, bforce, 0, f'H = {H}\ncells: solid {nsol}, fluid {nflu}\n')


def wavy(n=32):
    """z_w = 1 + 0.5 sin(2 pi x/4), upper wall mirrored; n cells in x, n/4 per solid and n/2 in the fluid in z
    (vertical mesh lines, spline interfaces)"""
    LX, LZ, H0, A, nu = 4., 4., 1., 0.5, 0.01
    bforce = 2*nu*1.0
    zl = lambda x: H0 + A*np.sin(2*np.pi*x/LX)
    zu = lambda x: LZ - zl(x)
    def edges(v):
        xs = np.linspace(0, LX, 81)[1:-1]
        return [(3, 2, 'spline', [(x, zl(x)) for x in xs]), (5, 4, 'spline', [(x, zu(x)) for x in xs])]
    bmd = layered_mesh(n, LX, n//4, n//2, zl, zu, LZ, edges)
    return gen_layered('wavy', f'n{n:03d}', bmd, nu, bforce, 2, '')


def cylinder_mesh(ppd):
    """O-grid: box 2D x 2D around the cylinder (D = 1 at (3, 2)) with 2 ppd cells per side, ppd radial cells graded
    towards the wall; outside the box dx = D/ppd as in CaNS; the solid has a square core and a ring"""
    XC, ZC, R, b, s = 3., 2., 0.5, 1.0, 0.2
    nq, nr, nrs = 2*ppd, ppd, max(ppd//2, 4)
    c45 = np.sqrt(0.5)
    v = []
    def add(x, z): v.append((x, z)); return len(v) - 1
    sq = [add(XC + sx*s, ZC + sz*s) for sx, sz in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ci = [add(XC + sx*R*c45, ZC + sz*R*c45) for sx, sz in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    bx = [add(XC + sx*b, ZC + sz*b) for sx, sz in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    corners = [(2., 1.), (4., 1.), (4., 3.), (2., 3.)]
    g = {}
    for x in (0., 2., 4., 12.):
        for z in (0., 1., 3., 4.):
            g[x, z] = bx[corners.index((x, z))] if (x, z) in corners else add(x, z)
    blocks = [('cylinder', (sq[0], sq[1], sq[2], sq[3]), (nq, nq), (1, 1)),
              ('cylinder', (ci[0], ci[1], sq[1], sq[0]), (nq, nrs), (1, 2.)),
              ('cylinder', (sq[1], ci[1], ci[2], sq[2]), (nrs, nq), (0.5, 1)),
              ('cylinder', (sq[3], sq[2], ci[2], ci[3]), (nq, nrs), (1, 0.5)),
              ('cylinder', (ci[0], sq[0], sq[3], ci[3]), (nrs, nq), (2., 1))]
    gf = 6.
    blocks += [('fluid', (bx[0], bx[1], ci[1], ci[0]), (nq, nr), (1, 1/gf)),
               ('fluid', (ci[1], bx[1], bx[2], ci[2]), (nr, nq), (gf, 1)),
               ('fluid', (ci[3], ci[2], bx[2], bx[3]), (nq, nr), (1, gf)),
               ('fluid', (bx[0], ci[0], ci[3], bx[3]), (nr, nq), (1/gf, 1))]
    ang = lambda a: (XC + R*np.cos(np.radians(a)), ZC + R*np.sin(np.radians(a)))
    edges = [(ci[0], ci[1], 'arc', [ang(270)]), (ci[1], ci[2], 'arc', [ang(0)]),
             (ci[2], ci[3], 'arc', [ang(90)]), (ci[3], ci[0], 'arc', [ang(180)])]
    ncx = {(0., 2.): 2*ppd, (2., 4.): nq, (4., 12.): 8*ppd}
    ncz = {(0., 1.): ppd, (1., 3.): nq, (3., 4.): ppd}
    for (x0, x1), nx_ in ncx.items():
        for (z0, z1), nz_ in ncz.items():
            if (x0, z0) != (2., 1.):
                blocks.append(('fluid', (g[x0, z0], g[x1, z0], g[x1, z1], g[x0, z1]), (nx_, nz_), (1, 1)))
    col = lambda x: [(g[x, 0.], g[x, 1.]), (g[x, 1.], g[x, 3.]), (g[x, 3.], g[x, 4.])]
    row = lambda z: [(g[0., z], g[2., z]), (g[2., z], g[4., z]), (g[4., z], g[12., z])]
    pat = {'left': ('cyclic', col(0.), '        neighbourPatch right;\n'),
           'right': ('cyclic', col(12.), '        neighbourPatch left;\n'),
           'bottom': ('wall', row(0.), ''), 'top': ('wall', row(4.), '')}
    return G.blockmesh(v, blocks, edges, pat, LY)


FUNCS_CYL = """    forces
    {
        type            forces;
        libs            (forces);
        region          fluid;
        patches         (fluid_to_cylinder);
        rho             rhoInf;
        rhoInf          1;
        CofR            (3 0 2);
        writeControl    timeStep;
        writeInterval   1;
    }
    average
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        region          fluid;
        timeStart       200;
        writeControl    writeTime;
        fields          ( U { mean on; prime2Mean on; base time; } T { mean on; prime2Mean on; base time; } );
    }
    averageSolid
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        region          cylinder;
        timeStart       200;
        writeControl    writeTime;
        fields          ( T { mean on; prime2Mean on; base time; } );
    }
    historyFluid
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        region          fluid;
        operation       volAverage;
        fields          (T U);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   10;
        log             false;
    }
    historyCylinder
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        region          cylinder;
        operation       volAverage;
        fields          (T);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   10;
        log             false;
    }
"""


def cylinder(ppd=16):
    """Re_D ~ 80: nu = 1/80, bforce = 0.054 (as the CaNS example), initial flow of the CaNS example"""
    nu, bforce = 1./80, 0.054
    bmd = cylinder_mesh(ppd)
    np_ = {16: 1, 32: 4}.get(ppd, 8)
    for K, C in PAIRS_CYL:
        case = HERE/'cylinder'/f'ppd{ppd:03d}'/tag(K, C)
        ctrl = G.control('chtMultiRegionFoam', 400, 1e-3, 100, adjust=True, maxco=0.5, functions=FUNCS_CYL)
        G.write_common(case, nu, bforce, {'cylinder': (K*nu, C)}, False, 1, ctrl)
        G.put(case/'system'/'blockMeshDict', 'dictionary', 'blockMeshDict', bmd)
        # no pressure work in the energy equation (as in the CaNS scalar equation; with the non-solenoidal initial
        # field the first pressure solve would otherwise put a large dp/dt into the temperature)
        tp = case/'constant'/'fluid'/'thermophysicalProperties'
        tp.write_text(tp.read_text().replace('thermoType\n', 'dpdt            no;\n\nthermoType\n', 1))
        G.fields_fluid(case, '"(bottom|top)"')
        G.fields_solid(case, 'cylinder', None)
        G.write_allrun(case, 'chtMultiRegionFoam', ['fluid', 'cylinder'], np_)
        ar = case/'Allrun'; t = ar.read_text()
        key = 'splitMeshRegions -cellZones -overwrite > log.splitMeshRegions 2>&1 || exit 1\n'
        assert key in t
        ar.write_text(t.replace(key, key + 'postProcess -func writeCellCentres -region fluid -time 0 > log.cc0.fluid 2>&1 || exit 1\n'
                                'python3 ../../../set_u0.py . > log.set_u0 2>&1 || exit 1\n'))
        (case/'case.txt').write_text(f'study = cylinder\nppd = {ppd}\nks_kf = {K}\ncs_cf = {C}\nbforce = {bforce}\nnu = {nu}\n')
    return HERE/'cylinder'/f'ppd{ppd:03d}'


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in ('slab', 'wavy', 'cylinder'):
        sys.exit(__doc__)
    d = globals()[sys.argv[1]](*(int(a) for a in sys.argv[2:3]))
    print(f'wrote {d.relative_to(HERE)}/: ' + ' '.join(sorted(p.name for p in d.iterdir())))
