# Order of accuracy of the CHT conduction operator

`cht_order.py` measures the steady order of accuracy of the conduction operator that `scal_cht` applies to CHT
scalars (`&scalar is_cht = T`), against a manufactured conjugate solution, for both geometry paths
(`&ibm geometry = 'sdf'` and `'analytic'`).

## What it does

For each case, CaNS writes the sub-cell geometry it uses (`&ibm dump_ibm_geometry = T`: link fractions `thx`,
`thy`, `thz` and volume fractions `phi`, one time step, one task). The script assembles from them the operator of
`scal_cht`, a sum of face fluxes `k_face (T_N - T_P)/h` with `k_face = 1/(th/k_f + (1-th)/k_s)`. It then solves
the steady problem `L T + phi q_f + (1-phi) q_s = 0` with a sparse solver.

* **Box.** Periodic in x and y over [0,1]², z in [0, 0.75], one interface z = g(x,y): solid (k_s = `ks_kf`)
  below, fluid (k_f = 1) above. The exact solution sits in the ghost cells at z = 0 and z = 0.75.
* **Exact solution.** T_s = F and T_f = F + c P + d P², with P = z − g and c(x,y) chosen for flux continuity
  (as `cht3d.py` of the AMPHIBIOUS project). It is continuous, has a continuous normal flux, varies along the
  interface, has sources in both materials and a fluid-only jump of d²T/dn².
* **Error.** Measured at every cell centre, fluid and solid, with the code's material (level set > 0: fluid).

| family | wall z = g(x,y) |
|---|---|
| `flat-on` | 0.25: on a cell face of every grid |
| `flat-off` | 0.25 + dz/3: the cut held at 1/3 of a cell above a face |
| `corr` | 0.3 + 0.05 sin(2πx) |
| `curved` | 0.3 + 0.06 cos(2πx) cos(2πy) |
| `curved2` | `curved` + 0.03 sin(2π(x+y)) (`sdf` only; not an `ibm_shape`) |

The `sdf` path reads the exact signed distance (`data/sdfp.bin`, Newton projection onto the wall). The
`analytic` path describes the same wall with `ibm_shape = 'sinusoid'` (flat, `corr`) or `'macdonald'` (`curved`).

The parameters are κ = `ks_kf` ∈ {0.5, 2, 10} and the cell aspect ratio dx/dz = dy/dz ∈ {1, 4}.

| aspect | N = nx = ny | nz |
|---|---|---|
| 1 | 16, 24, 32, 48, 64 | 0.75 N |
| 4 | 8, 12, 16, 24, 32 | 3 N |

Orders are computed between successive grids, from h = dx.

## Usage

```
python3 cht_order.py run      [--paths sdf,analytic] [--families ...] [--kappas 0.5,2,10] [--aspects 1,4]
                              [--jobs 8] [--quick] [--cans ../../build/cans] [--work ../../run/cht_order]
python3 cht_order.py fidelity      # the assembled operator vs one RK3 step of CaNS (zero velocity)
python3 cht_order.py report        # reprint <work>/results.json
./testit.sh                        # quick CI check (fidelity, reduced order set, pytest test.py)
```

The full matrix (2 paths, 5 families, 3 κ, 2 aspects, 5 grids) takes about 5 minutes with 8 jobs. It needs NumPy,
SciPy and SymPy. Cases already run are reused: the script skips the CaNS run when `data/ibm_geometry.bin`
exists.

`fidelity` takes one explicit RK3 step of CaNS from a manufactured field. The step has zero velocity, sources in
both materials, `cs_cf = 2` and non-zero Dirichlet values. It compares the result with the same step on the
assembled operator. The two agree to round-off (1e-16), on both paths.

## Regression of the examples

`regress_examples.sh <cans> <refdir> <rundir> [case ...] [-- key=value ...]` runs `examples/cht-levelset/slab` and
`wavy`. If `<refdir>/<case>.md5` is missing, it saves the checksums of the outputs there. Otherwise it checks them,
bit for bit. Entries after `--` are set in `input.nml`.
