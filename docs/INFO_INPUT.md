# About the input file `input.nml`

Consider the following input file as an example (corresponding to a turbulent plane channel flow). `&dns` defines a Fortran namelist containing all the necessary physical and computational parameters to set a case.


```fortran
&dns
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
visci = 5640.
inivel = 'poi'
is_wallturb = T
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
velf(1:3) = 1., 0., 0.
dims(1:2) = 2, 2, ipencil_axis = 1
/
```
<details>

<summary>Tip for vim/nvim users</summary>
Consider adding the following lines in your `.vimrc` file for syntax highlighting of the namelist file:

```vim
if has("autocmd")
  au BufNewFile,BufRead *.nml set filetype=fortran
  au BufNewFile,BufRead *.namelist set filetype=fortran
endif
```

</details>

---
---

```fortran
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
```

These lines set the computational grid.

`ng(1:3)` and `l(1:3)` are the **number of points**  and **domain length** in each direction.

`gtype` and `gr` are the **grid stretching type** and **grid stretching parameter** that tweak the non-uniform grid in the third direction; zero `gr` implies no stretching. See `initgrid.f90` for more details. The following options are available for `gtype`:

* `1`: grid clustered towards both ends (default)
* `2`: grid clustered towards the lower end
* `3`: grid clustered towards the upper end
* `4`: grid clustered towards the middle
* `9`: face coordinates read from the file `grid_zf.in` in the run directory: `ng(3)` values `zf(1:n)` with `zf(n) = l(3)`, or `ng(3)+1` values `zf(0:n)` with `zf(n) - zf(0) = l(3)`, in which case all z coordinates are shifted so that the domain starts at `zf(0)`; `gr` is ignored

---

```fortran
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
```

This line controls the simulation time step size.

The time step size is set to be equal to `min(cfl*dt_cfl,dtmax)` if `dt_f < 0`, and to `dt_f` otherwise. In the former case, the code prescribes the minimum value between `dtmax` and `cfl` times the maximum allowable time step `dt_cfl` (computed every `icheck` time steps; see below). `dtmax` is therefore used when a constant time step, smaller than `cfl*dt_cfl`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dt_cfl`. Finally, a constant time step size may be forced, irrespective of the temporal stability evaluation through `dt_f`.

---

```fortran
visci = 5640.
```

This line defines the inverse of the fluid viscosity, `visci`, meaning that the viscosity is `visc = visci**(-1)`. Note that, for a setup defined with unit reference length and velocity scales, `visci` has the same value as the flow Reynolds number.

---

```fortran
inivel = 'poi'
is_wallturb = T
```

These lines set the initial velocity field.

`initvel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`                                     ; streamwise direction in `x`
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `tbl`: temporal boundary layer profile with wall velocity `uref`                  ; streamwise direction in `x`
* `pdc`: plane Poiseuille flow profile with constant pressure gradient              ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `hdc`: half plane Poiseuille flow profile with constant pressure gradient         ; streamwise direction in `x`
* `tgv`: three-dimensional Taylor-Green vortex
* `tgw`: two-dimensional   Taylor-Green vortex
* `ant`: three-dimensional Antuono vortex

`is_wallturb`, if true, **superimposes a high amplitude disturbance on the initial velocity field** that effectively triggers transition to turbulence in a wall-bounded shear flow.

See `initflow.f90` for more details.

---

```fortran
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
```

These lines set the simulation termination criteria and whether the simulation should be restarted from a checkpoint file.

`nstep` is the **total number of time steps**.

`time_max` is the **maximum physical time**.

`tw_max` is the **maximum total simulation wall-clock time**.

`stop_type` sets which criteria for terminating the simulation are to be used (more than one can be selected, and at least one of them must be `T`)

* `stop_type(1)`, if true (`T`), the simulation will terminate after `nstep` time steps have been simulated
* `stop_type(2)`, if true (`T`), the simulation will terminate after `time_max` physical time units have been reached
* `stop_type(3)`, if true (`T`), the simulation will terminate after `tw_max` simulation wall-clock time (in hours) has been reached

`restart`, if true, **restarts the simulation** from a previously saved checkpoint file, named `fld.bin` (or `fld.h5` / `fld.bp` for the other supported backends).

`is_overwrite_save`, if true, overwrites the checkpoint file `fld.bin` at every save; if false, a symbolic link is created which makes `fld.bin` point to the last checkpoint file with name `fld_???????.bin` (with `???????` denoting the corresponding time step number). In the latter case, to restart a run from a different checkpoint one just has to point the file `fld.bin` to the right file, e.g.: ` ln -sf fld_0000100.bin fld.bin`.

`nsaves_max` limits the number of saved checkpoint files, if `is_over_write_save` is false; a value of `0` or any negative integer corresponds to no limit, and the code uses the file format described above; otherwise, only `nsaves_max` checkpoint files are saved, with the oldest save being overwritten when the number of saved checkpoints exceeds this threshold; in this case, files with a format `fld_????.bin` are saved (with `????` denoting the saved file number), with `fld.bin` pointing to the last checkpoint file as described above; moreover, a file `log_checkpoints.out` records information about the time step number and physical time corresponding to each saved file number.

Finally, note that a checkpoint file will be saved before the simulation ends.

---

```fortran
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
```

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the stability criterion and `cfl` (above)
* every `iout0d` time steps **history files with global scalar variables** are appended; currently the forcing pressure gradient and time step history are reported
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file
* every `iout3d` time steps **3d scalar fields** are written to a file
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation

1d, 2d, and 3d outputs can be tweaked by modifying files `out?d.h90` and recompiling the source. See also `output.f90` for more details. _Set any of these variables to `0` to skip the corresponding operation._

---

```fortran
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
```

These lines set the boundary conditions (BC).

The **type** (BC) for each field variable is set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1`, where:

* `X0` `X1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `z`

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The **last four rows** follow the same logic, but now for the BC **values** (dummy for a periodic direction).

---

```fortran
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
velf(1:3) = 1., 0., 0.
```

These lines set the flow forcing.

`bforce` is a constant **body force density term** in the direction in question (e.g., the negative of a constant pressure gradient) that can be added to the right-hand side of the momentum equation. The three values correspond to the three domain directions.

`is_forced`, if true in the direction in question, **forces the flow** with a pressure gradient that balances the total wall shear (e.g., for a pressure-driven channel with zero net acceleration/constant bulk velocity). The three boolean values correspond to the three domain directions.

`velf` is the **target bulk velocity** in the direction in question (where `is_forced` is true). The three values correspond to the three domain directions.

With the immersed-boundary method (`&ibm is_ibm = T`), the bulk velocity is the mean over the fluid points only, and the forcing impulse is added as the exact IB scheme integrates it (divided by the IB coefficient `A` near the wall, and not at all at the solid points), so the fluid bulk velocity matches `velf` exactly before the pressure projection. `&ibm bulk_vel_def` selects the normalization: `'superficial'` (default; divided by the box volume) or `'intrinsic'` (divided by the fluid volume). The same applies to `scalf` for the scalars. The output `forcing.out` reports this fluid-only mean.

---

```fortran
dims(1:2) = 2, 2, ipencil_axis = 1
```

This line sets the domain decomposition and orientation of the computational subdomains.

`dims` is the **processor grid**, the number of domain partitions along the first and second decomposed directions (which depend on the selected default pencil orientation). `dims(1)*dims(2)` corresponds therefore to the total number of computational subdomains. Setting `dims(:) = [0,0]` will trigger a runtime autotuning step to find the processor grid that minimizes transpose times. Note, however, that other components of the algorithm (e.g., collective I/O) may also be affected by the choice of processor grid.

`ipencil_axis` sets the **orientation of the computational subdomains** (or pencils), being one of `[1,2,3]` for `[X,Y,Z]`-aligned pencils. X-aligned is the default if this option is not set, and should be optimal for all cases except for Z-implicit diffusion, where Z-pencils are recommended if `dims(2) > 1` in the input file; see the description of the `&numerics` namelist below.

---

```fortran
iaverage = 0, average_start_time = 0.
```
**3D time averages** (`&dns`): every `iaverage` steps (0: off) for `time >= average_start_time`, the running means of u, v, w (at the cell centres), p, uu, vv, ww, uv, uw, vw and, per scalar, s, s², u s and w s are updated. They are written at every checkpoint as `data/<name>_mean.bin` (for visualization; `vex`, `vey`, `vez`, `p`, `u2`, `v2`, `w2`, `uv`, `uw`, `vw`, `s_<nnn>`, `s2_<nnn>`, `us_<nnn>`, `ws_<nnn>`), `data/<name>_mean_res.bin` (for restarting, with the time appended) and `data/count_avg.out` (number of samples). On restart the averages continue from these files if present.

# About the `&cudecomp` namelist under `input.nml`

In addition to the `&dns` namelist in the input file, there is an **optional** namelist to set some runtime configurations for the *cuDecomp* library. Consider the following `&cudecomp` namelist, which corresponds to the default options in case the file is not provided:

```fortran
&cudecomp
cudecomp_t_comm_backend = 0, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = T
cudecomp_h_comm_backend = 0, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = T
```

The first line sets the configuration for the transpose communication backend autotuning. Here `cudecomp_t_comm_backend` can be one of:

* `1` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P`
* `2` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL`
* `3` -> `CUDECOMP_TRANSPOSE_COMM_MPI_A2A`
* `4` -> `CUDECOMP_TRANSPOSE_COMM_NCCL`
* `5` -> `CUDECOMP_TRANSPOSE_COMM_NCCL_PL`
* `6` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM`
* `7` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL`
* `8` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM_SM`
* any other value -> enable runtime transpose backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_t_enable_nccl`) and NVSHMEM (`cudecomp_is_t_enable_nvshmem`) options for *transpose* communication backend autotuning.

The second line is analogous to the first one, but for halo communication backend autotuning. Here `cudecomp_h_comm_backend` can be one of:

* `1` -> `CUDECOMP_HALO_COMM_MPI`
* `2` -> `CUDECOMP_HALO_COMM_MPI_BLOCKING`
* `3` -> `CUDECOMP_HALO_COMM_NCCL`
* `4` -> `CUDECOMP_HALO_COMM_NVSHMEM`
* `5` -> `CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING`
* any other value -> enable runtime halo backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_h_enable_nccl`) and NVSHMEM (`cudecomp_is_h_enable_nvshmem`) options for *halo* communication backend autotuning.

Finally, it is worth recalling that passing `dims(1:2) = [0,0]` under `&dns` will trigger the *processor grid* autotuning, so there is no need to provide that option in the `&cudecomp` namelist.

# About the `&scalar` namelist under `input.nml`

This namelist is **optional** and defines the parameters needed to solve the transport equations associated with an arbitrary number of scalars.

The **number of scalars** `nscal` has to be set in the `&dns` namelist, e.g., for a single scalar:
```fortran
&dns
! (...)
nscal = 1
\
```
The default value of `nscal` is zero.

---

The following example namelist defines the parameters needed for the differentially heated cavity case that may be found under `examples/`.
```fortran
&scalar
alphai(1)              = 842.61498
beta                   = 1
iniscal(1)             = 'dhc'
cbcscal(0:1,1:3,1)     = 'D'  ,'D' ,  'P','P',  'N','N'
bcscal(0:1,1:3,1)      =  -0.5,0.5 ,   0.,0. ,   0.,0.
ssource(1)             = 0.
is_sforced(1)          = F
scalf(1)               = 0.
is_boussinesq_buoyancy = T
/
```

---

```fortran
alphai(1)        = 842.61498
beta             = 1
```
These lines define the **inverse of the diffusivity** `alphai` of the scalar(s) and, in this case, the **thermal expansion coefficient** `beta`.

`alphai` is an array with size `nscal`. To define several scalars, the inverse diffusivities can be defined as:
```fortran
alphai(:)        = alphai_1, alphai_2, alphai_3, ..., alphai_nscal
```

---

```fortran
iniscal(1)       = 'dhc'
```
This line sets the **initial scalar field(s)**. `iniscal` is an array with size `nscal`. The following options are available:
* `zer`: uniform scalar field equal to zero
* `uni`: uniform scalar field equal to `sref`
* `cou`: linearly varying scalar field from bottom (`z = 0`) to top   (`z = l(3)`)
* `dhc`: linearly varying scalar field from left   (`x = 0`) to right (`x = l(1)`)
* `tbl`: temporal boundary layer profile with scalar bc `sref` at the bottom wall

---

```fortran
cbcscal(0:1,1:3,1) = 'D'  ,'D' ,  'P','P',  'N','N'
bcscal(0:1,1:3,1)  =  -0.5,0.5 ,   0.,0. ,   0.,0.
```
These lines set the **boundary conditions** for each scalar, just like the velocity components. The last dimension corresponds to the scalar index.

---

```fortran
ssource(1)       = 0.
is_sforced(1)    = F
scalf(1)         = 0.
```
These lines set the **scalar forcing**.

`ssource` is an array of size `nscal` defining a **uniform volumetric source term** added to the righ-hand-side of each scalar equation (analogous to `source` in momentum).

`is_sforced` is a logical array of size `nscal` that **triggers the bulk forcing** of the corresponding scalar field (analogous to `is_forced` in momentum).

`scalf` is an array of size `nscal` defining **the target bulk mean** for each scalar (where `is_sforced` is true; analogous to `velf` in momentum).

```fortran
const_scal        = T
scal_mean_grad(1) = -1.
```
If `const_scal` is false (default true), each scalar `s` is the fluctuation about an imposed mean streamwise gradient, `S = scal_mean_grad*x + s`, and the advection of that mean gradient, `-u*scal_mean_grad`, is added to the right-hand side of the scalar equation (together with `ssource`). It is part of the explicit right-hand side, so it is integrated with the same Runge-Kutta weights as the advection term. `scal_mean_grad` is an array of size `nscal`, default `-1.`. With conjugate heat transfer (`is_cht = T`) the term is added in the same way (per unit heat capacity of the cell); it vanishes in the solid, where u = 0, so that, e.g., the mixed-type heating of Kasagi et al. (1992) with a conducting wall (Tiselj et al. 2001, Flageul et al. 2015) is `scal_mean_grad = -1/(Re_b Pr)` with the wall heat flux imposed at the outer boundary of the solid (`cbcscal = 'N'`).

---

**Conjugate heat transfer with the immersed boundary** (`is_ibm = T`)

```fortran
is_cht(1)     = T
ks_kf(1)      = 10.
cs_cf(1)      = 1.
ssource(1)    = 1.
ssource_s(1)  = 0.
qw_ibm(1)     = 0.
qw_ibm_eps    = 1.5
bcscal_ibm_top(1) = 0.
```
With `is_cht = T` (array of size `nscal`, default `F`) the scalar is also solved inside the IBM solid, in one domain:

`C dT/dt + u.grad(T) = div(k grad(T)) + ssource*phi + ssource_s*(1-phi) + qw_ibm*a`,

with all properties relative to the fluid: `k = alpha` in the fluid and `alpha*ks_kf` in the solid, `C = 1 + (cs_cf-1)*(1-phi)`, `phi` the fluid volume fraction of the cell, and `a` the interface area per unit volume. `ks_kf` and `cs_cf` are the solid-to-fluid **conductivity** and **volumetric heat-capacity** ratios (default 1). The conductivity at a cell face comes from the series resistance of the segment between the two cell centres, `k_face = 1/(th/k_f + (1-th)/k_s)`, where `th` is the fluid fraction of the segment taken at the exact interface position, so the heat flux crossing the interface is continuous there; the interface conditions (continuity of T and of the normal heat flux) are not imposed separately. `cbcscal_ibm`/`bcscal_ibm` are not used for CHT scalars.

`cht_kface` (in `&scalar`, default `'series'`) selects the face conductivity at the faces near the interface. The series form above is exact for heat flowing across the interface, but along it a face takes the conductivity of the material of its two cell centres, however much of the face lies in the other material (a cut cell whose centre is in the fluid conducts along the wall like pure fluid). With `cht_kface = 'tensor'` (requires `geometry = 'sdf'`) a face near the interface gets the diagonal element of the effective conductivity tensor of a locally planar interface, `k_face = w*k_series(th) + (1-w)*(a*k_f + (1-a)*k_s)`, with `a` the fluid area fraction of the face (cut by the zero level of the quadratic expansions of the level set about the two adjacent cell centres) and `w = n_d**2` the squared component of the interface normal along the face direction. It is exact for a planar interface both across it and along it when the interface is aligned with the grid, and reduces to the series form away from the interface and for `ks_kf = 1`. It is applied per CHT scalar (through an equivalent segment fraction), in the explicit and the implicit diffusion and in the wall heat flux output (`FV-heat`). Steady conduction along a solid strip whose interfaces lie inside the cells (`ks_kf` = 0.01 and 100): the error of the uncut cells drops from 14 % (series) to 0.19 % (aligned grid: 0.10 %), and the solid-mean temperature error at `ks_kf = 100` from -8 % to -0.002 %. Memory: three extra arrays per CHT scalar (plus their pencil copies with implicit diffusion). Like the series form, it is limited by the stability of the factorized implicit diffusion (below): the wavy-wall example with `ks_kf = 100`, `cs_cf = 0.1` (diffusivity ratio 1000) at its laminar CFL time step (dt = 0.2, implicit factors of order 100 in all three directions) is unstable with `'tensor'` (a y odd-even mode in a cut cell) and stable at dt = 0.1.

`'tensor'` uses only the diagonal of that tensor, which keeps the 7-point stencil of the line solvers but is exact only for interfaces aligned with the grid. For an inclined interface the off-diagonal elements `K_de = (k_series(th) - k_parallel(a))*n_d*n_e` carry the tangential flux (along the interface) through faces normal to another direction; without them the face conductivity of a highly conducting solid leaks into the fluid. With `cht_kface = 'tensor_full'` they are added as an explicit, conservative term (`scal_cht_cross`): the flux through, e.g., an x face gets `K_xy*dT/dy + K_xz*dT/dz`, with the tangential derivatives at the face centre from the four surrounding cell centres; it is nonzero only at the faces near the interface (and zero through the domain boundaries along z). By default (`cht_kcross_scheme = 'explicit'`) it is part of the Runge-Kutta slope, also with implicit diffusion, where the diagonal stays implicit. Steady conduction along a solid strip inclined at 30 deg to the grid (16x32 and 32x64 cells, errors of the uncut cells against a converged fine-grid reference): `ks_kf = 0.01`, rms error 8.1/3.6 % (series), 3.4/1.8 % (`'tensor'`), 0.18/0.06 % (`'tensor_full'`); `ks_kf = 100`, rms 1.8/0.8 %, 7.3/3.8 %, 1.6/0.6 % and fluid-mean error +2.8/+1.2 %, -14.2/-7.8 %, -2.9/-1.1 %; for a grid-aligned strip `'tensor_full'` is identical to `'tensor'`. The wall heat flux output (`FV-heat`, `wall_heat.out`) includes the tangential part of the flux through the faces between fluid and solid cell centres (with the same stencil), so that it is the heat that actually enters the solid cells: for a conducting cylinder (`ks_kf = 100`, 16x32 cells) heated inside, the net heat flow out of the solid cells equals their source to round-off with all three forms, whereas without the tangential part the `'tensor_full'` value would be 7.8 times too large. Memory: six more arrays per CHT scalar. Stability: the explicit off-diagonal term limits the time step independently of the implicit diagonal. At start-up the code prints, per CHT scalar, its largest explicit rate B (the sum of the moduli of its stencil coefficients in a cell, 1/time; B ~ |k_series - k_parallel|/(C*dx*dz), i.e. it grows with the grid refinement like 1/h**2), and the time step is limited to `dt <= cht_kcross_cfl/B` (`&scalar`, default 2; with a fixed `dt_f` the limit is not applied, only a warning is printed). On the 30 deg strip (`ks_kf = 100`, 16x32 cells, implicit diffusion, 20000 steps) the explicit form was stable up to dt*B = 5 (`cs_cf = 1`, B = 18406; unstable at 6) and dt*B = 12 (`cs_cf = 100`, B = 1613; unstable at 16), with the same steady state as at small dt to 2e-13. `cht_kcross_scheme = 'douglas'` (the earlier form) puts the term into the explicit predictor of the factorized (Douglas) scheme, forward Euler per substep, without a time-step limit; it became unstable for dt*B of about 10-25 (on the strip: `cs_cf = 100` stable at dt*B = 6, unstable at 25 after 780 steps; `cs_cf = 1` stable at 2.3, unstable at 9.4 after 6400 steps and at 17 after 320). `'craig-sneyd'` adds the Craig-Sneyd correction (a second line solve per substep, theta = 1/2) and was less stable than `'douglas'` (unstable after 420 and 220 steps in the last two cases). `cht_kcross_stab` (default 0, `'douglas'` and `'craig-sneyd'` only) adds a bound of the off-diagonal term to the face conductivities of the line solvers (left-hand side only, so the steady state is unchanged); with values 0.5-2 it delayed the instability (e.g. `cs_cf = 1`, dt*B = 17: from 320 to 470-2760 steps) but did not remove it, except at `cs_cf = 100`, dt*B = 25 (stable for 6000 steps).

`'series'`, `'tensor'` and `'tensor_full'` are all first-order accurate at an inclined interface: the conductivities of `'tensor_full'` are those of the homogenized medium, which relate the mean flux to the mean gradient of a volume, whereas the scheme multiplies them with difference quotients between cell centres, and these do not see the kink of the temperature at the interface in the right proportion. `cht_kface = 'tensor_o2'` (requires `geometry = 'sdf'`) chooses the face conductivities so that the face flux is exact for a piecewise linear temperature at a planar interface of any orientation, `T = T0 + g.(x-x0) + (r-1)*(n.g)*min(psi,0)` (`g`: gradient on the fluid side, `n`: unit normal into the fluid, `psi`: signed distance, `r = k_f/k_s`; continuous temperature, tangential gradient and normal flux). The difference quotients of the flux stencil of this field are `G_d = g_d + (r-1)*(n.g)*n_d*s_d`, with `s_d` the solid fraction of the stencil segments along `d` (the two adjacent cell centres along the face normal, the four surrounding segments of `scal_cht_cross` along the tangential directions; from the segment fractions `th` of the series form, so the exact interface position enters). Inverting this for `n.g` and inserting it into the exact flux through a face of fluid area fraction `a` gives, for an x face (y and z alike), `F_x = K_xx*G_x + K_xy*G_y + K_xz*G_z` with `K_xx = k_a + c*n_x`, `K_xy = c*n_y`, `K_xz = c*n_z`, `k_a = a*k_f + (1-a)*k_s`, `c = n_x*((1-a)*(k_f-k_s) - k_a*(r-1)*s_x)/sigma`, `sigma = 1 + (r-1)*sum_d(n_d**2*s_d)`. It reduces to the series form for a grid-aligned interface and to `k_f` or `k_s` where the stencil is not cut; the diagonal element stayed between `k_f` and `k_s` in all tests (it is clipped to this range, and the number of clipped faces is printed at start-up). The implementation is that of `'tensor_full'` with other coefficients: diagonal through the equivalent segment fraction (explicit and implicit diffusion, `FV-heat`), off-diagonal elements as the explicit term of `scal_cht_cross` with the same inputs (`cht_kcross_scheme`, `cht_kcross_cfl`) and the same memory; the faces on the domain boundaries along z keep the diagonal element only. Steady conduction, solid strip inclined at 30 deg in the x-z plane (16x32 to 128x256 cells, rms error of the uncut cells in % of the maximum temperature, against an interface-conforming finite-element reference): `ks_kf = 100`: 1.77/0.77/0.35/0.18 (series, order 1.2-0.9), 1.58/0.58/0.22/0.089 (`'tensor_full'`, order 1.4-1.3), 1.09/0.29/0.071/0.017 (`'tensor_o2'`, order 1.9-2.1); `ks_kf = 0.01`: 8.1/3.6/1.7/0.89, 0.18/0.056/0.029/0.0091 and 0.33/0.076/0.028/0.012 (order 2.1-1.3; here the nearly insulating strip meets the isothermal wall at an obtuse fluid angle, and the corner singularity of the solution, not the interface treatment, limits the order; on the finest grid the error is also of the size of the uncertainty of the reference, 0.016 % in the maximum norm). The same strip in the y-z plane gives the x-z result to round-off. A strip inclined in the x-y plane (interface parallel to z, tan = 0.5, 16x32x32 to 64x128x128 cells) gives 1.05/0.27/0.068 % (`ks_kf = 100`) and 0.16/0.041/0.010 % (0.01), order 2.0 (maximum error of the uncut cells 1.9/0.66/0.13 and 0.24/0.072/0.015 %), against 1.74/0.43 and 0.21/0.052 % with `'tensor_full'` on the two coarser grids; in this configuration the heat flows along the strip in a grid direction and only crosses the inclined interface, which the series form already treats well (0.34/0.10 and 1.12/0.37 %). In a Python prototype of the same discretization, a conducting cylinder in a box with the exact solution as boundary values converges with second order in the rms and the maximum norm for `ks_kf` = 0.01, 0.1, 10 and 100 from 16 to 512 cells per box length (`ks_kf = 100`: rms 0.24 % at 16 cells, 0.0030 % at 128; series and `'tensor_full'`: order 1 to 1.5). Limits: the gain needs a resolved, smooth interface; at the coarsest resolutions (16 cells per strip period or cylinder box) `'tensor_o2'` is not more accurate than `'tensor_full'` or the series form, and corners of the geometry or of the solution reduce the order locally. The time-step limit of the explicit term is that of `'tensor_full'` with a rate B about 13 % larger: on the 30 deg strip (16x32 cells, 20000 steps) it was stable up to dt*B = 4 for `cs_cf = 1` (`ks_kf` = 100 and 0.01; unstable at 5 and 8) and up to dt*B = 16 for `cs_cf = ks_kf` (`ks_kf = 100`: unstable at 24), so the default `cht_kcross_cfl = 2` holds. Tested for steady conduction only, serial, with 2x2 MPI ranks (identical steady state) and on one GPU (OpenACC build without further changes; same rate B, and the steady states of the x-z strip with 16x32 and 128x256 cells and of the x-y strip with 16x32x32 cells equal those of the CPU build to round-off, 1e-16 to 3e-16; the 64x128x128 x-y case was run on the GPU only); not yet tested with flow, with a genuinely three-dimensional interface against a reference, or on several GPUs.

* `ssource` acts in the fluid and `ssource_s` (default 0) in the solid, both weighted with the volume fractions of the cut cells.
* `qw_ibm` (default 0) is a **heat flux released at the interface** per unit area (a jump of the normal heat flux, `-k_s dT/dn|_s + qw_ibm = -k_f dT/dn|_f`, n into the fluid), deposited with a smoothed delta function of half-width `qw_ibm_eps` grid spacings about the zero level of the level set; it requires `geometry = 'sdf'`.
* Diffusion is explicit by default, with the time step limited by the largest diffusivity `max(alpha, alpha*ks_kf/cs_cf)`. With `&numerics is_impdiff = T` it is integrated implicitly in all three directions (Crank-Nicolson in delta form, approximately factorized, Douglas ADI); this is stable for solid-to-fluid diffusivity ratios `ks_kf/cs_cf` up to about 100, not at 1000.
* A non-finite scalar value stops the run with exit code 1.

`bcscal_ibm_top` (non-CHT scalars, default `bcscal_ibm`) is the IB value of the scalar above the domain mid-plane (z > l(3)/2), e.g. a different temperature on the upper wall.

---

**Buoyancy effects**

The input value `beta` above is the thermal expansion coefficient that relates density to temperature in the so-called _Boussinesq approximation_, thereby triggering buoyancy-induced momentum transport. Its default value is zero. When buoyancy is activated, **the first scalar is always associated with the temperature difference**. Moreover, in that case, the value for the **gravitational acceleration vector** `gacc` should be provided in the `&dns` namelist:
```fortran
&dns
! (...)
gacc(1:3) = 0., 0., -1.
nscal = 1
/
```
which, in this example, defines a negative gravitational acceleration along `z` with unit magnitude. Its default value is `[0., 0., 0.]`.

Finally, **activate the buoyancy term** with `is_boussinesq_buoyancy = T` in the `&scalar` namelist.

# About the `&numerics` namelist under `input.nml`

This namelist defines parameters related to the numerical discretization and computational method. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&numerics
is_impdiff = F, is_impdiff_1d = F
is_poisson_pcr_tdma = F
adv_form = 'vv'
/
```

In these lines, `is_impdiff` and `is_impdiff_1d` enable the (semi-) **implicit temporal integration of diffusion terms**:

* `is_impdiff`, if `.true.`, the diffusion term of the Navier-Stokes and scalar equations is integrated in time implicitly, which may improve the stability of the numerical algorithm for viscous-dominated flows.
* `is_impdiff_1d` is similar to `is_impdiff`, but with implicit diffusion *only* along Z, which may be advantageous when the grid along Z is much finer than along the other directions; *for optimal parallel performance, the domain should not be decomposed along Z* (`ipencil_axis=3`, or `ipencil_axis = 1/2` with `dims(2) = 1`)
* With the IBM (`&ibm is_ibm = T`), implicit diffusion is supported only with `is_impdiff_1d = T`. Each Z column is then solved directly, with the IB inside the operator: the solid points are held at the wall value by a penalty, and the IB term keeps its exact (exponential) time integration, while the Z diffusion is Crank-Nicolson. This requires a domain that is not decomposed along Z (`ipencil_axis = 3`, or `dims(2) = 1`) and non-periodic Z boundaries. The steady state is the same as with explicit diffusion. `ibm_force.out` stays exact, but the IBM bulk forcing no longer sets the bulk velocity exactly, because the forcing is applied before the implicit solve (as without the IBM). With `is_impdiff = T, is_impdiff_1d = F` the 3D FFT-based solve cannot include the IB, so this combination is rejected.

Finally, `is_poisson_pcr_tdma`, if `.true.`, allows for solving the Poisson/Helmholtz equations along Z with a parallel cyclic reduction--tridiagonal matrix algorithm (PCR-TDMA) method. This approach may result in major gains in scalability for pencil-distributed simulations at scale, on many GPUs.

`adv_form` selects the **discretization of the advection term** of the momentum equation, both in the separate and in the fused (`is_fast_mom_kernels = T`) kernels:

* `'vv'` (default): skew-symmetric form of Verstappen & Veldman; it conserves kinetic energy also on non-uniform grids, provided the velocity field is discretely divergence-free.
* `'div'`: divergence (flux) form; it conserves momentum in the strict flux sense also where the discrete divergence is not zero (e.g., in near-wall cells with the immersed-boundary method), but it is not energy-conserving on a grid stretched along Z.

The two forms differ only by the term `u*div(u)/2`, so they coincide for a discretely divergence-free field. Before this option existed, `is_ibm = T` forced `'div'` when `is_fast_mom_kernels = F`.

```fortran
is_ibm_damp = F, ibm_damp_eps = 0.02, ibm_damp_d = 3., 6.
```
`is_ibm_damp` adds an optional fourth-difference hyperviscosity in x near the immersed boundary, `-ibm_damp_eps*sigma*|u|*dx^3*d4/dx4`, to u, v, w and the scalars, against odd-even (2 dx) oscillations that sharp IB edges can trigger. `sigma = 1` where the level set is below `ibm_damp_d(1)*dx`, with a cosine ramp to 0 at `ibm_damp_d(2)*dx`. It requires `geometry = 'sdf'` and an x direction that is not decomposed. The x stencil also reaches into the solid next to the wall and adds a small drag (about 1 % in the skin friction at `ibm_damp_eps = 0.04` on a truncated-cone surface); it is off by default.

# About the `&ibm` namelist under `input.nml` (selected entries)

```fortran
&ibm
ibm_shape    = 'macdonald'
ibm_hch      = -1.
bulk_vel_def = 'superficial'
ibm_proj     = 'none'
ibm_iout_wall = -1
geometry     = 'analytic'
sdf_file     = 'data/sdfp.bin'
dump_ibm_geometry = F
mproj_tol    = 1.e-3, mproj_maxit = 50
ibm_updiss_mode = 'none'
ibm_updiss   = 1.
ibm_updiss_ncell = 2
ibm_updiss_ntaper = 0
/
```
`ibm_shape` selects the analytic wall shape (with `&hmap use_hmap = F`). For a wall normal to direction `t` on side `s` (`ibm_direction(s,t) = T`), with in-plane directions `i = mod(t,3)+1` and `ii = mod(t+1,3)+1`, and the bottom wall at `l_0(0,t) + h`, the top wall at `l_0(1,t) - h`:
* `'macdonald'` (default): `h = amp_l(s,i)*(1 +/- cos(2 pi n_wave(s,i) x_i/l_i)*cos(2 pi n_wave(s,ii) x_ii/l_ii))` (`+` for the bottom wall); `amp_l(s,ii)` and `phase_l` are ignored, and a warning is printed if they are set.
* `'sinusoid'`: `h = sum over d = i, ii of amp_l(s,d)/2*(1 + sin(2 pi n_wave(s,d) x_d/l_d + phase_l(s,d)))`.
* `'cones'`: staggered truncated cones as in Theobald et al. (WCCM-ECCOMAS 2020): height `amp_l(s,i)`, base radius `1.25 amp_l(s,i)`, flank slope 55 degrees; `n_wave(s,i) x n_wave(s,ii)` units of size `l_i/n_wave(s,i) x l_ii/n_wave(s,ii)`, each with two cones centred at (1/4, 1/4) and (3/4, 3/4) of the unit.

The roughness height range printed at start-up (`hmin`, `hmax`) is sampled from the actual wall shape on a grid 4 times finer than the mesh in x and y.

`ibm_hch` is the reference half height (`Re_tau = ibm_hch/visc`) of the rough-wall grid `gtype = 5`, whose uniform fine region extends from `z = 0` to `hmax`. It is required for `gtype = 5` (with a height map, the value from the map's info file is used if `ibm_hch` is not set).

`geometry` selects how the immersed boundary is described:
* `'analytic'` (default): the wall shapes of `ibm_shape` or the height map of `&hmap`.
* `'sdf'`: a **level set** (signed distance, > 0 in the fluid) at the cell centres, read from `sdf_file`: raw float64, `ng(1)*ng(2)*ng(3)` values, x fastest, no header. It is interpolated linearly to the staggered velocity positions; the wall distance on each grid line comes from a quadratic crossing through three points, and the fluid volume fraction of a cut cell from the second-order Taylor expansion of the level set about the cell centre. Any geometry can be used (several bodies, overhangs, closed bodies); the solver needs the correct sign everywhere and a smooth, accurate zero level within about two cells of the wall, while `|grad phi| = 1` matters only for the per-point wall outputs and `is_ibm_damp`. Generators: `utils/sdf/hmap_sdf.py` (height maps) and `utils/sdf/stl_sdf.py` (STL surfaces). At least two cells per task along z are needed.

`dump_ibm_geometry = T` writes `data/ibm_geometry.bin` (serial runs only): the IB coefficients `lap_u, lap_v, lap_w, lap_s` and, with CHT scalars, the face fractions `thx, thy, thz` and the volume fractions `phi`, each `n(1)*n(2)*n(3)` doubles.

`bulk_vel_def` is described with `velf` above.

`ibm_proj` selects the velocity correction after the pressure Poisson equation, which is solved with the standard (not IB-aware) FFT solver:
* `'none'` (default): standard correction everywhere. The correction leaks into the solid points (O(dt^2) velocities there), which are zeroed again at the start of the next RK stage.
* `'damped'`: the correction is divided by the IB coefficient `A` at the non-solid points, as the exact IB scheme integrates a pressure-gradient impulse, and the solid points are not corrected, so their velocity stays exactly zero. The velocity is then divergence free only up to a residual in the cut cells, of the same size as the divergence that `'none'` creates when it zeroes the solid points, so the divergence check does not abort on it. The fluid bulk velocity is no longer conserved by the correction; with bulk forcing it is restored at the next stage. The cost is that of `'none'`.

In tests (turbulent channel over the MacDonald wall, 64^3) the two options give the same errors and temporal order (2 at the band points) to three digits; `'damped'` removes only the solid leakage. A lagged pre-compensation of the predictor with the previous pressure correction (Kim-Moin type) was tested and is unstable: the zeroing of the solid points at each stage start already acts as such a pre-compensation.

* `'masked'`: an IB-consistent projection (Seo & Mittal, J. Comput. Phys. 230 (2011) 7347, with binary face apertures M, 0 at the solid velocity points): `D(M G p) = D(u*)/dt`, `u = M (u* - dt G p)`, solved by conjugate gradients preconditioned with the FFT Poisson solver, to the relative residual `mproj_tol` or at most `mproj_maxit` iterations per RK stage. The velocity is then zero at the solid points and divergence free in every cell, including the cut cells. It costs about 15 iterations per stage (several times the cost of a time step); on a rough-wall test it left the statistics unchanged.

`ibm_updiss_mode` adds, near the immersed body, the dissipative part of third-order upwind advection, `-d/dx(c phi |u| dx^3/12 d^3q/dx^3) - d/dy(c phi |v| dy^3/12 d^3q/dy^3) - d/dz(c phi |w| dz^3/12 d^3q/dz^3)` for each velocity component `q`, in flux form (momentum conservative; in z with undivided differences on the stretched grid). The central advection is non-dissipative, and at cell Peclet numbers well above 2 the sharp IB makes it develop grid-scale (2 dx) wiggles upstream of the body, visible in the mean wall-normal velocity and in `u'w'`; the term damps them. `c = ibm_updiss`; for `c > 1.087` the advective time step is reduced by `1.087/c` (RK3 stability). `ibm_updiss_mode` is:

* `'none'`: no term (default); nothing is allocated or computed.
* `'body'`: the term acts with the weight `phi` in [0,1] of a layer around the immersed wall, grown from the IB mask of each velocity component: level 1 are the fluid points with a solid neighbour (the points carrying IB coefficients), level `m` the fluid points next to level `m-1`, neighbours counted in all 26 directions so that the layer is isotropic in grid cells. `phi = 1` up to level `ibm_updiss_ncell`, then a cosine taper over a further `ibm_updiss_ntaper` levels (solid points: 1). The layer depends neither on the geometry nor on the wall orientation, and shrinks with the grid. The term is computed only on the layer points.

With scalars (`nscal > 0`) the same term, with the layer grown from the scalar mask and the face velocities as advection velocities, is added to the scalar equations. Damping only the momentum would leave the scalar with the central advection next to the sharp IB and lower the heat transfer: in the MacDonald et al. (2019) case (`examples/macdonald-eggcarton-roughness`) `Ch` drops by 0.5 % with the momentum term alone and is unchanged (-0.1 %, within the statistical uncertainty) with both; drag and `Delta U+` are unchanged in both cases. The default `ibm_updiss_ncell = 2` damps the wiggles upstream of the truncated cones almost as well as 3 (high-pass rms of the mean wall-normal velocity 0.069 vs 0.066, 0.122 without the term).

Scalars with conjugate heat transfer (`is_cht = T`, solved in the fluid and the solid) get the same term on their own layer, grown from the cell-centre mask like that of the Dirichlet scalars, which also includes the solid cells next to the fluid (weight 1 there). The dissipative flux through a face with nonzero velocity between a fluid and a solid cell (with the level set, the velocity point of such a face can lie in the fluid) is thus added to both cells, so the term is conservative across the interface like the central advection of `scal_cht`, and it is divided by the heat capacity of the cell, as that advection. In the solid the face velocities are zero, so the term vanishes away from the interface. The second differences are evaluated on the layer and the solid cells next to it only (not the whole solid). `ibm_updiss_mode` and `is_ibm_damp` are independent and can be combined.

`ibm_iout_wall` is the interval, in time steps, of the per-point wall outputs `fric-*`, `pressure*`, `FV-fric_*` and `ibm_wall_stress_*` below. A negative value (default) uses `iout1d`; `0` switches them off. The plane statistics `turb_stats_*` stay at `iout1d`. On fine grids the per-cell files are large (one line per cut cell), so a larger interval than `iout1d` keeps the statistics frequent without them.

IB force output (with `is_ibm = T`):
* `data/ibm_force.out`, every `iout0d` steps: the force the IB exerts on the fluid in that time step, per unit box volume (positive = drag), computed exactly from the momentum the IB removes from the plain explicit update. Columns: time; driving force `bforce - dpdl` (x, y, z); total IB force (x, y, z); its parts in the lower and upper halves of the domain (x, y, z each); the parts from the fluid (band) points, i.e. the `lambda` term, from the solid points (predictor zeroed), and from the projection leakage into the solid (x, y, z each). With `ibm_proj = 'damped'` there are three more columns, the part of the damped correction at the fluid (band) points; its solid-point part is included in the leakage columns. At steady state, total IB force + stresses on the domain walls = driving force. The drag on a z-wall is the total minus the driving force on the solid volume (`pi*h` per unit plan area). The z component includes the pressure inside the solid and on the domain boundary, so it is not a wall force.
* `data/ibm_wall_stress_<step>.out`, every `iout1d` steps (lower z-wall only): plan-view map of the wall stress, `tau(i,j) = F(i,j)/(dx dy) - pi h(x,y)`, where `F(i,j)` is the IB force integrated over the lower half of the column and `h` is the wall height. Columns: i, j, tau_x at (xf(i),yc(j)), tau_y at (xc(i),yf(j)). Its mean is the drag per unit plan area.
* `data/turb_stats_ibm_<step>.out`, every `iout1d` steps (z-walls), next to `turb_stats_<step>.out` (unchanged, superficial over all points including the solid): plane statistics over the fluid points only. Columns: z; fluid fractions `phi` of the u-, v- and cell-centre points; intrinsic means of u, v, w; intrinsic covariances u'u', v'v', w'w', u'w', u'v', v'w' (w and the products with w at the cell centres, over the fluid cells, with the solid faces counted as 0). Superficial means are `phi` times the intrinsic ones. In a fully solid plane the values are 0.
* `data/FV-fric_<rank>_s_<step>.out`, every `iout1d` steps, one file per task: the wall force per cell in that time step, for every cell that contains a piece of the wall. It is the exact IB force of `ibm_force.out` distributed over the cells, so its sum is the drag of `ibm_wall_stress_*.out` (up to the solid volume of the cut cells, sampled with 20^3 points per cell). Each velocity node next to the wall passes its IB impulse to the cell in which its wall point lies. For a band node this is the `lambda` part, split over the wall directions in proportion to `lambda`. For a solid node it is the local momentum flux from its fluid neighbours. The pressure part of the solid nodes is telescoped onto the cut cells: a cell gets `p - p_mean` times its face area for each solid face, with the sign of the face, so the pressure inside the solid cancels. The driving force on the solid part of the cell is subtracted. Columns: i, j, k (global); wall area; unit normal out of the solid (x, y, z); force of the fluid on the wall F (x, y, z); its normal part `F.n` and vector; its tangential part (magnitude, vector); the same divided by the wall area (normal stress `F.n/A`, i.e. minus the wall pressure relative to the fluid-mean pressure, and shear stress). The values are local, but a cell holds the force of the nodes next to its piece of wall, not a pointwise traction. With `geometry = 'sdf'` the wall points, the solid volume, the wall area and the normal of each cell come from the level set (the area and normal from the signed projected areas of the sub-cell solid/fluid faces).
* `data/FV-heat_<rank>_s_<step>.out`, with scalars, whenever the wall outputs are on (`ibm_iout_wall /= 0`), same cells and order as `FV-fric`: the wall heat flux per cell, from the diffusive face fluxes of the segments between a fluid and a solid cell centre (CHT: the series resistance of the face, the flux the solver uses; non-CHT: a linear profile from the fluid node to the wall value at the distance `th*h`), assigned to the cell that contains the interface point. Columns: i, j, k (global); wall area; unit normal (x, y, z); per scalar the heat flow Q from the fluid into the wall, Q divided by the wall area, and the interface temperature (mean over the crossed faces).
* `data/wall_heat.out`, at the same steps: time; the total heat flow into the wall of each scalar (`nscal` columns); the part of it in the output cells of each scalar (`nscal` columns).

# About the `&other_options` namelist under `input.nml`

This namelist defines other parameters related to the monitoring, debugging, or benchmarking of a computation. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&other_options
is_debug = T, is_timing = T
/
```

In these lines, `is_debug` performs some **sanity checks for debugging purposes** that do not introduce computational overhead, and `is_timing` reports the wall-clock time per step.

Note: other parameters for `&numerics` and `&other_options` are not exposed here, as they are meant for very specific, advanced use or developers. They can be found in `param.f90`.

This optional namelist selects the **I/O backend** used for checkpointing and field output.

# About the `&io` namelist under `input.nml`

```fortran
&io
io_backend = 'mpiio'
/
```

The following options are currently available for `io_backend`:

* `mpiio`: raw binary output based on MPI-I/O (default)
* `hdf5`: HDF5 output for checkpoint and visualization files (requires building with `USE_HDF5=1`)
* `adios2`: ADIOS2 output for checkpoint and visualization files (requires building with `USE_ADIOS2=1`)

If `&io` is not provided, *CaNS* defaults to `io_backend = 'mpiio'`. Note also that selecting `hdf5` or `adios2` requires compiling the code with support for the corresponding library; otherwise, *CaNS* falls back to `mpiio` with a warning at runtime.
