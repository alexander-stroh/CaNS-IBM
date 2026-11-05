! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
#if defined(_OPENACC) && !defined(_USE_DIEZDECOMP)
use cudecomp
#endif
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
#if !defined(_EPS_EXACT_ZERO) /* recommended */
real(rp), parameter :: eps = epsilon(1._rp)
#else
real(rp), parameter :: eps = 0._rp
#endif
real(rp), parameter :: small = epsilon(1._rp)*10**(precision(1._rp)/2)
real(rp), parameter :: ibm_pen = 1._rp/epsilon(1._rp) ! solid-point penalty of the z-implicit IBM solve
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
! variables to be determined from the input file
!
integer , protected, dimension(3) :: ng
real(rp), protected, dimension(3) :: l
integer , protected :: gtype
real(rp), protected :: gr
real(rp), protected :: cfl,dtmax,dt_f
real(rp), protected :: visci
!rough input
real(rp), protected :: dzp_fine
real(rp), protected :: dzp_coarse
!
character(len=100), protected :: inivel
logical, protected :: is_wallturb
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save
integer , protected :: nsaves_max
integer , protected :: icheck,iout0d,iout1d,iout2d,iout3d,isave
integer , protected :: iaverage           ! 3D time averages every iaverage steps for time >= average_start_time (0: off)
real(rp), protected :: average_start_time
!
integer , dimension(2) :: dims ! not protected -> overwritten when autotuning
integer , protected :: ipencil_axis = 1
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
!
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(3) :: is_forced
real(rp), protected, dimension(3) :: velf
!
real(rp), protected, dimension(3) :: dl,dli
real(rp), protected :: visc
!
! scalar input parameters
!
integer , protected :: nscal ! number of transported scalars
real(rp), protected, allocatable, dimension(:) :: alphai
real(rp), protected :: beta
real(rp), protected, dimension(3) :: gacc
character(len=100), protected, allocatable, dimension(:)     :: iniscal
character(len=1)  , protected, allocatable, dimension(:,:,:) :: cbcscal ! size (0:1,3,nscal)
real(rp)          , protected, allocatable, dimension(:,:,:) ::  bcscal ! size (0:1,3,nscal)
real(rp), protected, allocatable, dimension(:) :: ssource
logical,protected :: const_scal ! if false, the mean-gradient source -u*scal_mean_grad is added
real(rp), protected, allocatable, dimension(:) :: scal_mean_grad
logical , protected, allocatable, dimension(:) :: is_sforced
real(rp), protected, allocatable, dimension(:) :: scalf
logical , protected :: is_boussinesq_buoyancy = .false.
real(rp), protected :: alpha_max
! scalar-ibm
character(len=1)  , protected, allocatable, dimension(:) :: cbcscal_ibm ! size (nscal)
real(rp)          , protected, allocatable, dimension(:) ::  bcscal_ibm ! size (nscal)
real(rp)          , protected, allocatable, dimension(:) ::  bcscal_ibm_top ! size (nscal), IB value above the domain mid-plane
!
! conjugate heat transfer (CHT): scalar also solved inside the IBM solid, with
! ks_kf = solid-to-fluid conductivity ratio and cs_cf = solid-to-fluid volumetric heat capacity ratio
!
logical , protected, allocatable, dimension(:) :: is_cht ! size (nscal)
real(rp), protected, allocatable, dimension(:) :: ks_kf  ! size (nscal)
real(rp), protected, allocatable, dimension(:) :: cs_cf  ! size (nscal)
!
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
logical, protected :: is_use_diezdecomp = .false.,is_diezdecomp_x2z_z2x_transposes = .false.
#endif
!
! other options: debugging/benchmarking
!
logical, protected :: is_debug = .true., is_debug_poisson = .false., &
                      is_timing = .true., &
                      is_mask_divergence_check = .false.
!
! other options: numerics
!
logical, protected :: is_impdiff = .false., is_impdiff_1d = .false., &
                      is_poisson_pcr_tdma = .false., &
                      is_fast_mom_kernels = .true., &
                      is_gridpoint_natural_channel = .false.
character(len=3), protected :: adv_form = 'vv' ! advection form: 'vv' (skew-symmetric, Verstappen & Veldman) or 'div' (divergence)
!
! optional damping of the x odd-even (2 dx) mode near the immersed boundary (level-set geometry only, off by default):
! rhs += -ibm_damp_eps*sigma*|u|*dx^3*d4/dx4 for u, v, w and the scalars, sigma = 1 for phi <= ibm_damp_d(1)*dx,
! cosine ramp to 0 at phi = ibm_damp_d(2)*dx (phi: level set, > 0 in the fluid)
!
logical , protected :: is_ibm_damp = .false.
real(rp), protected :: ibm_damp_eps = 0.02_rp, ibm_damp_d(2) = [3._rp,6._rp]
!
! i/o backend parameters
!
character(len=16), protected :: io_backend = 'mpiio'
character(len=6) , protected :: io_ext = '.bin'
logical          , protected :: is_use_compression = .false.
! IBM
logical,protected     :: is_ibm=.false.
logical,protected     :: do_richardson=.false.
logical, protected    :: ibm_direction(0:1,3)=.false.
real(rp), protected   :: l_0(0:1,3)=0._rp
integer, protected    :: n_wave(0:1,3)=0._rp
real(rp), protected   :: amp_l(0:1,3)=0._rp
real(rp), protected   :: phase_l(0:1,3)=0._rp
character(len=16), protected :: ibm_shape = 'macdonald' ! analytic wall shape: 'macdonald' or 'sinusoid'
real(rp), protected   :: ibm_hch = -1._rp ! reference half height for the rough-wall grid (gtype = 5)
character(len=16), protected :: bulk_vel_def = 'superficial' ! IBM bulk velocity: 'superficial' or 'intrinsic'
character(len=16), protected :: ibm_proj = 'none' ! IB-aware velocity correction: 'none' or 'damped'
integer , protected   :: ibm_iout_wall = -1 ! interval of the per-point wall outputs (< 0: iout1d, 0: off)
! geometry of the IBM: 'analytic' (`isInbody`, wavy walls/height map) or 'sdf' (signed-distance level set
! at the cell centres read from `sdf_file`, see `mod_ibm_sdf`)
character(len=16) , protected :: geometry = 'analytic'
character(len=100), protected :: sdf_file = 'data/sdfp.bin'
logical           , protected :: dump_ibm_geometry = .false. ! write the IBM/CHT geometry fields (serial runs)
!hmap
logical,protected     :: use_hmap=.false.
character(len=25), protected :: hmap_mode="fit"
logical,protected     :: override_grid=.false.
character(len=25), protected :: hmap_name="thakkar"
!
contains
  subroutine read_input(myid)
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use mpi
    implicit none
    character(len=255) :: input_file,filename
    integer, intent(in) :: myid
    integer :: iunit,ierr,numarg
    character(len=1024) :: c_iomsg
    logical :: is_io_fallback
    namelist /dns/ &
                  ng, &
                  l, &
                  gtype,gr, &
                  dzp_fine,dzp_coarse,&
                  cfl,dtmax,dt_f, &
                  visci, &
                  inivel, &
                  is_wallturb, &
                  nstep,time_max,tw_max, &
                  stop_type, &
                  restart,is_overwrite_save,nsaves_max, &
                  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                  iaverage,average_start_time, &
                  cbcvel,cbcpre,bcvel,bcpre, &
                  bforce, &
                  is_forced, &
                  velf, &
                  gacc, &
                  nscal, &
                  dims,ipencil_axis
    namelist /io/  &
                  io_backend, &
                  is_use_compression
#if defined(_OPENACC) && !defined(_USE_DIEZDECOMP)
    namelist /cudecomp/ &
                       cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem, &
                       cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
#endif
    namelist /scalar/ &
                     alphai,beta, &
                     iniscal, &
                     cbcscal,bcscal, &
                     cbcscal_ibm,bcscal_ibm,bcscal_ibm_top,&
                     const_scal,scal_mean_grad,&
                     is_cht,ks_kf,cs_cf,&
                     ssource, &
                     is_sforced, &
                     scalf, &
                     is_boussinesq_buoyancy
    namelist /numerics/ &
                       is_impdiff,is_impdiff_1d, &
                       is_poisson_pcr_tdma, &
                       is_fast_mom_kernels, &
                       is_ibm_damp,ibm_damp_eps,ibm_damp_d, &
                       is_gridpoint_natural_channel, &
                       adv_form
    namelist /other_options/ &
                            is_debug,is_debug_poisson, &
                            is_timing, &
                            is_mask_divergence_check
    namelist /ibm/ &
                            is_ibm,           &
                            do_richardson,    &
                            ibm_direction,    &
                            l_0,              &
                            n_wave,           &
                            amp_l,            &
                            phase_l,          &
                            ibm_shape,        &
                            ibm_hch,          &
                            bulk_vel_def,     &
                            ibm_proj,         &
                            ibm_iout_wall,    &
                            geometry,         &
                            sdf_file,         &
                            dump_ibm_geometry
    namelist /hmap/ &
                            use_hmap,         &
                            hmap_mode,        &        
                            override_grid,    &
                            hmap_name        
         

    !
    ! defaults
    !
    dt_f = -1.
    gacc(:) = 0.
    nscal = 0
    iaverage = 0
    average_start_time = 0.
    !
    ! use a custom input file if one was provided on the command line
    !
    input_file = 'input.nml'
    numarg = command_argument_count()
    if(numarg == 1) then
      call get_command_argument(1,filename)
      input_file = trim(filename)
    end if
    open(newunit=iunit,file=input_file,status='old',action='read',iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading the input file: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read `dns` namelist
      !
      rewind(iunit)
      read(iunit,nml=dns,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading `dns` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      dl(:) = l(:)/(1.*ng(:))
      dli(:) = dl(:)**(-1)
      visc = visci**(-1)
      if(all([1,2,3] /= ipencil_axis)) then
        ipencil_axis = 1 ! default to one
        if(myid == 0) print*, 'Warning: prescribed value of `ipencil_axis` is different from 1/2/3.'
        if(myid == 0) print*, 'Defaulting to `ipencil_axis = 1` (x-aligned pencils)...'
      end if
#if defined(_OPENACC)
#if !defined(_USE_DIEZDECOMP)
      !
      ! reading cuDecomp parameters, if these are set
      !
      ! defaults
      !
      cudecomp_is_t_comm_autotune  = .true.
      cudecomp_is_h_comm_autotune  = .true.
      cudecomp_is_t_enable_nccl    = .true.
      cudecomp_is_h_enable_nccl    = .true.
      cudecomp_is_t_enable_nvshmem = .true.
      cudecomp_is_h_enable_nvshmem = .true.
      rewind(iunit)
      read(iunit,nml=cudecomp,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading `cudecomp` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 8) then
        cudecomp_is_t_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_t_comm_backend)
        case(1)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        case(2)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
        case(3)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
        case(4)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
        case(5)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
        case(6)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
        case(7)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
        case(8)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_SM
        case default
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        end select
      end if
      if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 5) then
        cudecomp_is_h_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_h_comm_backend)
        case(1)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        case(2)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI_BLOCKING
        case(3)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NCCL
        case(4)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM
        case(5)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING
        case default
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        end select
      end if
      !
      ! manually set cuDecomp out-of-place transposes by default
      !
      cudecomp_is_t_in_place = .false.
#else
      is_use_diezdecomp = .true.
      is_diezdecomp_x2z_z2x_transposes = .false.
#endif
#endif
      !
      ! reading scalar transport parameters, if these are set
      !
      if(nscal > 0) then
        !
        ! allocate memory
        !
        allocate(alphai(nscal),iniscal(nscal), &
                 cbcscal(0:1,3,nscal),bcscal(0:1,3,nscal), &
                 cbcscal_ibm(nscal),bcscal_ibm(nscal),bcscal_ibm_top(nscal),&
                 is_cht(nscal),ks_kf(nscal),cs_cf(nscal),&
                 ssource(nscal),is_sforced(nscal),scalf(nscal),scal_mean_grad(nscal))
        !
        ! set default values
        !
        beta          = 0.
        ssource(:)    = 0.
        is_sforced(:) = .false.
        scalf(:)      = 0.
        cbcscal_ibm(:) = "D"
        bcscal_ibm(:)  = 0._rp
        bcscal_ibm_top(:) = huge(1._rp) ! unset: same as bcscal_ibm
        const_scal     = .true.
        scal_mean_grad(:) = -1._rp
        is_cht(:)      = .false.
        ks_kf(:)       = 1._rp
        cs_cf(:)       = 1._rp
        !
        ! read `scalar` namelist
        !
        rewind(iunit)
        read(iunit,nml=scalar,iostat=ierr,iomsg=c_iomsg)
        if(ierr /= 0) then
          if(myid == 0) print*, 'ERROR: reading `scalar` namelist: ', trim(c_iomsg)
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          close(iunit)
          error stop
        end if
        where(bcscal_ibm_top(:) == huge(1._rp)) bcscal_ibm_top(:) = bcscal_ibm(:)
      else
        nscal = 0 ! negative values equivalent to nscal = 0
      end if
      alpha_max = huge(1._rp)
      alpha_max = minval(alphai(1:nscal))
      alpha_max = alpha_max**(-1)
      if(is_boussinesq_buoyancy) then
        if(nscal == 0) then
          if(myid == 0) print*, 'ERROR: reading the input file: `is_boussinesq_buoyancy = T` requires `nscal > 0`.'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      end if
      !
      ! read `numerics` namelist
      !
      rewind(iunit)
      read(iunit,nml=numerics,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `numerics` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      if(is_impdiff_1d) is_impdiff = .true.
      !
      ! read `other_options` namelist
      !
      rewind(iunit)
      read(iunit,nml=other_options,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `other_options` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read `io` namelist
      !
      rewind(iunit)
      read(iunit,nml=io,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `io` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      is_io_fallback = .false.
      select case(trim(io_backend))
      case('mpiio')
        io_ext = '.bin'
      case('hdf5')
#if defined(_USE_HDF5)
        io_ext = '.h5'
#else
        if(myid == 0) print*, 'Warning: `io_backend = hdf5` requires an HDF5 build.'
        is_io_fallback = .true.
#endif
      case('adios2')
#if defined(_USE_ADIOS2)
        io_ext = '.bp'
#else
        if(myid == 0) print*, 'Warning: `io_backend = adios2` requires an ADIOS2 build.'
        is_io_fallback = .true.
#endif
      case default
        if(myid == 0) print*, 'Warning: unknown `io_backend = ', trim(io_backend), '`.'
        is_io_fallback = .true.
      end select
      if(is_io_fallback) then
        if(myid == 0) print*, 'Defaulting to `io_backend = mpiio`...'
        io_backend = 'mpiio'
        io_ext = '.bin'
      end if
      if(is_use_compression .and. trim(io_backend) == 'mpiio') then
        if(myid == 0) print*, 'Warning: compression is ignored for `io_backend = mpiio`.'
      end if
    rewind(iunit)
    !
    ! read "IBM" namelist
    !
    read(iunit,nml=ibm,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `ibm` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
    rewind(iunit)
    !
    ! read hmap namelist
    !
    read(iunit,nml=hmap,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `hmap` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
    close(iunit)
    !
    ! checks of the CHT input (the `ibm` namelist is read after `scalar`)
    !
    if(nscal > 0) then
      if(any(is_cht(:))) then
        if(.not.is_ibm) then
          if(myid == 0) print*, 'ERROR: reading the input file: `is_cht = T` requires `is_ibm = T`.'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
        if(any(ks_kf(:) <= 0._rp).or.any(cs_cf(:) <= 0._rp)) then
          if(myid == 0) print*, 'ERROR: reading the input file: `ks_kf` and `cs_cf` must be positive.'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      end if
    end if
  end subroutine read_input
end module mod_param
