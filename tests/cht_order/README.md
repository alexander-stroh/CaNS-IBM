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

## The tangential face flux (`--scheme tang`)

`--scheme tang` runs every case with `&scalar is_cht_rkl2 = T, is_cht_tang = T`. CaNS then also writes its tangential
face list (`data/cht_tang_001.bin`: per face the lower cell, the axis and the 18 coefficients of `(k_area-k_face)*s_t`
on the 2x3x3 block), which the script adds to the series operator exactly as `cht_tang_apply` does (`sdf` path only;
the tangential term needs the level set). `fidelity --scheme tang` takes one RKL2 step of CaNS (two half steps, zero
velocity) against the same step on the assembled operator.

Results (2026-11-24, `sdf`, 150 cases, `run/cht_order_dev/tang_sdf.json`):

* **Fidelity:** 4e-16 (RKL2 step of CaNS against the harness).
* **Order:** L2 2.00-2.06 on the finest pair in all 30 groups; Linf 1.97-2.03, except three groups at κ = 10
  (1.84-1.93: `flat-off` aspect 1 and 4, `curved2` aspect 4).
* **Against the series flux** (`run/cht_order/phase0_step1.json`): at aspect 1 the L2 error is 1.4-12 times and the
  Linf error 3-19 times smaller on the cut families; at aspect 4 the errors are similar on these grids (the series
  flux's first-order defect is masked there, see the port log), but the series orders keep falling while the
  tangential ones stay at 2. `flat-on` is the series operator itself (no face gets a term).
* **Conservation:** column sums of the volume-weighted operator below 3.4e-16 relative; 0 link/material mismatches.

Two further checks:

* `plane_exactness.py --cans <cans>`: an interface inclined at 30 deg and the conjugate-linear field
  `a + sigma.s + q n/k`. The face flux of the scheme is exact to round-off on the y and z faces (4e-15 to 3.5e-14
  relative at κ = 2 and 12), and to the 20-column face quadrature on the x faces (1.2e-4, 1.2e-3); the series flux
  alone misses the faces by 13-89 %.
* `worked_example_check.py --cans <cans> --ref <worked_example.json>`: the node of "The example in numbers" of the
  AMPHIBIOUS note (reference values in `worked_example.json` of the AMPHIBIOUS project, not in this repository). The series conductances agree to 1e-4 (S) and 1e-3 (W),
  the coefficients `(k_area-k_face)*w` of the three faces with a term to 1.5-2.3 % on S and E; on W the weights agree to
  0.5 % and `k_area-k_face`, small there, to 15 %. The differences come from the geometry: CaNS takes the interface
  from the level set at the cell centres, AMPHIBIOUS bisects the exact wall.
