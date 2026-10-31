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
If `const_scal` is false (default true), each scalar `s` is the fluctuation about an imposed mean streamwise gradient, `S = scal_mean_grad*x + s`, and the advection of that mean gradient, `-u*scal_mean_grad`, is added to the right-hand side of the scalar equation (together with `ssource`). It is part of the explicit right-hand side, so it is integrated with the same Runge-Kutta weights as the advection term. `scal_mean_grad` is an array of size `nscal`, default `-1.`.

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

# About the `&ibm` namelist under `input.nml` (selected entries)

```fortran
&ibm
ibm_shape    = 'macdonald'
ibm_hch      = -1.
bulk_vel_def = 'superficial'
ibm_proj     = 'none'
ibm_iout_wall = -1
ibm_updiss_mode = 'none'
ibm_updiss   = 1.
ibm_updiss_band = 2.
ibm_updiss_ncell = 2
ibm_updiss_ntaper = 0
/
```
`ibm_shape` selects the analytic wall shape (with `&hmap use_hmap = F`). For a wall normal to direction `t` on side `s` (`ibm_direction(s,t) = T`), with in-plane directions `i = mod(t,3)+1` and `ii = mod(t+1,3)+1`, and the bottom wall at `l_0(0,t) + h`, the top wall at `l_0(1,t) - h`:
* `'macdonald'` (default): `h = amp_l(s,i)*(1 +/- cos(2 pi n_wave(s,i) x_i/l_i)*cos(2 pi n_wave(s,ii) x_ii/l_ii))` (`+` for the bottom wall); `amp_l(s,ii)` and `phase_l` are ignored, and a warning is printed if they are set.
* `'sinusoid'`: `h = sum over d = i, ii of amp_l(s,d)/2*(1 + sin(2 pi n_wave(s,d) x_d/l_d + phase_l(s,d)))`.

The roughness height range printed at start-up (`hmin`, `hmax`) is sampled from the actual wall shape on a grid 4 times finer than the mesh in x and y.

`ibm_hch` is the reference half height (`Re_tau = ibm_hch/visc`) of the rough-wall grid `gtype = 5`, whose uniform fine region extends from `z = 0` to `hmax`. It is required for `gtype = 5` (with a height map, the value from the map's info file is used if `ibm_hch` is not set).

`bulk_vel_def` is described with `velf` above.

`ibm_proj` selects the velocity correction after the pressure Poisson equation, which is solved with the standard (not IB-aware) FFT solver:
* `'none'` (default): standard correction everywhere. The correction leaks into the solid points (O(dt^2) velocities there), which are zeroed again at the start of the next RK stage.
* `'damped'`: the correction is divided by the IB coefficient `A` at the non-solid points, as the exact IB scheme integrates a pressure-gradient impulse, and the solid points are not corrected, so their velocity stays exactly zero. The velocity is then divergence free only up to a residual in the cut cells, of the same size as the divergence that `'none'` creates when it zeroes the solid points, so the divergence check does not abort on it. The fluid bulk velocity is no longer conserved by the correction; with bulk forcing it is restored at the next stage. The cost is that of `'none'`.

In tests (turbulent channel over the MacDonald wall, 64^3) the two options give the same errors and temporal order (2 at the band points) to three digits; `'damped'` removes only the solid leakage. A lagged pre-compensation of the predictor with the previous pressure correction (Kim-Moin type) was tested and is unstable: the zeroing of the solid points at each stage start already acts as such a pre-compensation.

`ibm_updiss_mode` adds, near the immersed body, the dissipative part of third-order upwind advection, `-d/dx(c phi |u| dx^3/12 d^3q/dx^3) - d/dy(c phi |v| dy^3/12 d^3q/dy^3) - d/dz(c phi |w| dz^3/12 d^3q/dz^3)` for each velocity component `q`, in flux form (momentum conservative; in z with undivided differences on the stretched grid). The central advection is non-dissipative, and at cell Peclet numbers well above 2 the sharp IB makes it develop grid-scale (2 dx) wiggles upstream of the body, visible in the mean wall-normal velocity and in `u'w'`; the term damps them. `c = ibm_updiss`; for `c > 1.087` the advective time step is reduced by `1.087/c` (RK3 stability). The weight `phi` in [0,1] is:

* `'none'`: no term (default); nothing is allocated or computed.
* `'slab'`: x and y terms only, `phi` a function of z only: 1 up to `ibm_updiss_band*hmax` from a rough z-wall, then a cosine taper to 0 over a further `hmax`. Requires rough z-walls.
* `'body'`: a layer around the immersed wall, grown from the IB mask of each velocity component: level 1 are the fluid points with a solid neighbour (the points carrying IB coefficients), level `m` the fluid points next to level `m-1`, neighbours counted in all 26 directions so that the layer is isotropic in grid cells. `phi = 1` up to level `ibm_updiss_ncell`, then a cosine taper over a further `ibm_updiss_ntaper` levels (solid points: 1). The layer depends neither on the geometry nor on the wall orientation, and shrinks with the grid. The term is computed only on the layer points.

`ibm_iout_wall` is the interval, in time steps, of the per-point wall outputs `fric-*`, `pressure*`, `FV-fric_*` and `ibm_wall_stress_*` below. A negative value (default) uses `iout1d`; `0` switches them off. The plane statistics `turb_stats_*` stay at `iout1d`. On fine grids the per-cell files are large (one line per cut cell), so a larger interval than `iout1d` keeps the statistics frequent without them.

IB force output (with `is_ibm = T`):
* `data/ibm_force.out`, every `iout0d` steps: the force the IB exerts on the fluid in that time step, per unit box volume (positive = drag), computed exactly from the momentum the IB removes from the plain explicit update. Columns: time; driving force `bforce - dpdl` (x, y, z); total IB force (x, y, z); its parts in the lower and upper halves of the domain (x, y, z each); the parts from the fluid (band) points, i.e. the `lambda` term, from the solid points (predictor zeroed), and from the projection leakage into the solid (x, y, z each). With `ibm_proj = 'damped'` there are three more columns, the part of the damped correction at the fluid (band) points; its solid-point part is included in the leakage columns. At steady state, total IB force + stresses on the domain walls = driving force. The drag on a z-wall is the total minus the driving force on the solid volume (`pi*h` per unit plan area). The z component includes the pressure inside the solid and on the domain boundary, so it is not a wall force.
* `data/ibm_wall_stress_<step>.out`, every `iout1d` steps (lower z-wall only): plan-view map of the wall stress, `tau(i,j) = F(i,j)/(dx dy) - pi h(x,y)`, where `F(i,j)` is the IB force integrated over the lower half of the column and `h` is the wall height. Columns: i, j, tau_x at (xf(i),yc(j)), tau_y at (xc(i),yf(j)). Its mean is the drag per unit plan area.
* `data/turb_stats_ibm_<step>.out`, every `iout1d` steps (z-walls), next to `turb_stats_<step>.out` (unchanged, superficial over all points including the solid): plane statistics over the fluid points only. Columns: z; fluid fractions `phi` of the u-, v- and cell-centre points; intrinsic means of u, v, w; intrinsic covariances u'u', v'v', w'w', u'w', u'v', v'w' (w and the products with w at the cell centres, over the fluid cells, with the solid faces counted as 0). Superficial means are `phi` times the intrinsic ones. In a fully solid plane the values are 0.
* `data/FV-fric_<rank>_s_<step>.out`, every `iout1d` steps, one file per task: the wall force per cell in that time step, for every cell that contains a piece of the wall. It is the exact IB force of `ibm_force.out` distributed over the cells, so its sum is the drag of `ibm_wall_stress_*.out` (up to the solid volume of the cut cells, sampled with 20^3 points per cell). Each velocity node next to the wall passes its IB impulse to the cell in which its wall point lies. For a band node this is the `lambda` part, split over the wall directions in proportion to `lambda`. For a solid node it is the local momentum flux from its fluid neighbours. The pressure part of the solid nodes is telescoped onto the cut cells: a cell gets `p - p_mean` times its face area for each solid face, with the sign of the face, so the pressure inside the solid cancels. The driving force on the solid part of the cell is subtracted. Columns: i, j, k (global); wall area; unit normal out of the solid (x, y, z); force of the fluid on the wall F (x, y, z); its normal part `F.n` and vector; its tangential part (magnitude, vector); the same divided by the wall area (normal stress `F.n/A`, i.e. minus the wall pressure relative to the fluid-mean pressure, and shear stress). The values are local, but a cell holds the force of the nodes next to its piece of wall, not a pointwise traction.

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
