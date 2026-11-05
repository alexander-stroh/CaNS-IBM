# Conjugate heat transfer with the level-set IBM: laminar test cases

Three small 2D laminar cases (one cell layer of 4 points in y) for the conjugate heat transfer (CHT) with the
level-set geometry (`&ibm geometry = 'sdf'`). They run in 1 to 35 minutes on one CPU core and come from a validation of
the method against exact solutions and body-fitted OpenFOAM (`chtMultiRegionFoam`) simulations.

| case | geometry | grid | reference |
|---|---|---|---|
| `slab` | flat solid slabs of thickness 1.02 at both walls, fluid gap 1.96 | 4 × 4 × 64 | exact steady solution |
| `wavy` | sinusoidal solid walls z<sub>w</sub> = 1 + 0.5 sin(2πx/4), upper wall mirrored | 32 × 4 × 32 | OpenFOAM |
| `cylinder` | cylinder D = 1 at (3, 2) in a 12 × 4 channel, Re<sub>D</sub> ≈ 80, vortex shedding | 192 × 4 × 64 (D/Δx = 16) | OpenFOAM |

Common set-up: periodic in x and y, no-slip walls with T = 0 at z = 0 and z = 4, flow driven by a constant pressure
gradient `bforce`, Pr = 1, uniform heat source `ssource = 1` in the fluid. Every case carries six CHT scalars with the
conductivity and heat-capacity ratios

| scalar | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| `ks_kf` | 1 | 10 | 100 | 0.1 | 1 | 1 |
| `cs_cf` | 1 | 1 | 1 | 1 | 0.1 | 10 |

(`wavy` uses `cs_cf = 0.1` for all six: the steady state does not depend on it, and the small solid capacity shortens
the slow solid time scale.) The interfaces lie inside the cells, not on cell faces. Diffusion is implicit in z
(`is_impdiff_1d = T`, which needs `ipencil_axis = 3`) and explicit in x and y.

## Running

```
python setup.py slab               # writes slab/data/sdfp.bin (the level set, > 0 in the fluid)
cd slab && mpirun -np 1 ../../../build/cans && cd ..
python check.py slab
```

and the same for `wavy` and `cylinder`. For the cylinder `setup.py` also writes the initial state (`restart = T`):
Poiseuille flow with U<sub>b</sub> = 1, zero velocity in the cylinder, a small asymmetric w perturbation behind it to start
the shedding, and T = 0. `setup.py` needs only NumPy; it reads the grid from `<case>/input.nml`.

| case | end time | CPU, 1 core |
|---|---|---|
| `slab` | t = 20000 (14 diffusion times of the slowest scalar, k<sub>s</sub>/k<sub>f</sub> = 0.1) | 1 min |
| `wavy` | t = 3000 | 1 min |
| `cylinder` | t = 400 (statistics over t = 200 ... 400) | ~35 min (4 min on one A6000 GPU) |

## Expected results

`check.py` prints the numbers below. They are from the GPU validation runs with this code; CPU runs of `slab` and `wavy` reproduce them to all digits shown, and a CPU run of the cylinder
agreed with them to 10<sup>-13</sup> over the first steps.

**slab** (exact solution: Poiseuille flow, T = T<sub>i</sub> + q(d² − y²)/(2k<sub>f</sub>) in the fluid, linear in the slabs,
interface T<sub>i</sub> = q d H/k<sub>s</sub>):

```
slab, t = 20000.4: max |u - u_ex|/U_max = 6.01e-04
  k_s/k_f =     1: T_centre =  147.9375 (exact  147.9312), max |T - T_ex|/T_ex(centre) = 4.28e-05
  k_s/k_f =    10: T_centre =   57.9735 (exact   57.9672), max |T - T_ex|/T_ex(centre) = 1.09e-04
  k_s/k_f =   100: T_centre =   48.9771 (exact   48.9708), max |T - T_ex|/T_ex(centre) = 1.29e-04
  k_s/k_f =   0.1: T_centre = 1047.5769 (exact 1047.5712), max |T - T_ex|/T_ex(centre) = 5.51e-06
```

(scalars 5 and 6 give the same steady state as scalar 1).

**wavy** (no exact solution; at n = 128 the bulk velocity and the interface and maximum temperatures agree with
OpenFOAM within 0.14 % after Richardson extrapolation; the case here is the coarsest grid, n = 32):

```
wavy, t = 3000.09: U_b = 0.19831 (flow rate / fluid area)
  k_s/k_f =     1: mean T fluid  123.7524, solid   52.5718
  k_s/k_f =    10: mean T fluid   34.5714, solid    6.2432
  k_s/k_f =   100: mean T fluid   23.3078, solid    0.6522
  k_s/k_f =   0.1: mean T fluid  832.6560, solid  435.0163
```

At t = 3000 the k<sub>s</sub>/k<sub>f</sub> = 0.1 scalar is still slowly rising; it is steady after t ≈ 9000.

**cylinder** (D/Δx = 16; C = 2F/(U<sub>b</sub>² D) with the IB force of `ibm_force.out` minus the driving force on the
cylinder volume):

```
cylinder, t = 200 ... 400.021: U_b = 1.0007, Re_D = 80.06, C_D = 3.4504, C_D amplitude = 0.2711, C_L amplitude = 3.1511, St = 0.2773
  scalar 1 (k_s/k_f = 1): mean T = 55.3024
  scalar 2 (k_s/k_f = 10): mean T = 55.2980
  scalar 3 (k_s/k_f = 100): mean T = 55.2960
  scalar 4 (k_s/k_f = 0.1): mean T = 55.2823
  scalar 5 (k_s/k_f = 1): mean T = 55.3605
  scalar 6 (k_s/k_f = 1): mean T = 54.1560
```

At D/Δx = 64 the same set-up agrees with OpenFOAM at D/Δx = 64 within 0.5 % in C<sub>D</sub>, the C<sub>D</sub> and C<sub>L</sub>
amplitudes, St and U<sub>b</sub>; the mean temperatures are 0.7–0.9 % higher than in OpenFOAM at every resolution
(the mean T is still rising slowly at t = 400). To run it finer, change `ng` in `cylinder/input.nml`
(D/Δx = 32: `384, 4, 128`, `l(2) = 0.125`) and rerun `setup.py cylinder`.

## Outputs used

* `data/fld_*.bin` (checkpoint at the end of the run), `data/forcing.out` (bulk velocity, column 5),
  `data/forcing_s_<n>.out` (mean of scalar n, column 3), `data/ibm_force.out` (IB force per unit box volume, x in
  column 5, z in column 7): see [`docs/INFO_INPUT.md`](../../docs/INFO_INPUT.md).
* The per-cell wall outputs (`FV-fric`, `FV-heat`, `wall_heat.out`) are off (`ibm_iout_wall = 0`); set
  `ibm_iout_wall` > 0 to write them.
