! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
!-------------------------------------------------------------------------------------
program cans
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,bounduvw_norm,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr,dinfo_ptdma
  use mod_correc         , only: correc,correc_ibm
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_initflow       , only: initflow,initscal
  use mod_initgrid       , only: initgrid,save_grid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load_one
  use mod_mom            , only: bulk_forcing
  use mod_rk             , only: rk,rk_scal,calc_a_b
  use mod_timer          , only: timer_tic,timer_toc
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out1d_chan_ibm,out2d,out3d, &
                                 write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: ng,l,dl,dli, &
                                 gtype,gr, &
                                 dzp_fine,dzp_coarse,&
                                 cfl,dtmax,dt_f, &
                                 visc,alpha_max, &
                                 inivel,is_wallturb, &
                                 nstep,time_max,tw_max,stop_type, &
                                 restart,is_overwrite_save,nsaves_max, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 iaverage,average_start_time, &
                                 cbcvel,bcvel,cbcpre,bcpre, &
                                 is_forced,bforce,velf, &
                                 gacc,nscal,beta, &
                                 alphai,iniscal,cbcscal,bcscal,&
                                 cbcscal_ibm,bcscal_ibm,bcscal_ibm_top,&
                                 ssource,is_sforced,scalf,scal_mean_grad,&
                                 is_boussinesq_buoyancy,&
                                 dims, &
                                 nb,is_bound, &
                                 rkcoeff,small, &
                                 datadir,io_ext, &
                                 read_input, &
                                 is_debug,is_debug_poisson, &
                                 is_timing, &
                                 is_impdiff,is_impdiff_1d, &
                                 is_poisson_pcr_tdma,      &
                                 is_mask_divergence_check, &
                                 is_ibm,           &
                                 do_richardson,    &
                                 ibm_direction,    &
                                 l_0,              &
                                 n_wave,           &
                                 amp_l,            &
                                 phase_l,          &
                                 bulk_vel_def,     &
                                 ibm_proj,ibm_iout_wall, &
                                 ibm_hch,          &
                                 geometry,sdf_file,dump_ibm_geometry, &
                                 use_hmap,         &
                                 override_grid,     &
                                 hmap_name          
  use mod_sanity         , only: test_sanity_input,test_sanity_solver
  use mod_scal           , only: scalar,initialize_scalars,bulk_forcing_s
  use mod_stats3d        , only: stats3d_init,stats3d_update,stats3d_save
  use mod_solve_helmholtz, only: solve_helmholtz,solve_helmholtz_ibm,rhs_bound
  use mod_cht            , only: cht_init
  use mod_ibm_sdf        , only: sdf_read,sdf_staggered,sdf_mask,set_ibm_2nd_sdf,calc_grad_dist_sdf, &
                                 set_cht_geometry_sdf,sdf_interface_density
  use mod_ibm_damp       , only: ibm_damp_init
  use mod_param          , only: is_ibm_damp,ibm_damp_d
  use mod_param          , only: qw_ibm,qw_ibm_eps
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
  use mod_common_cudecomp, only: istream_acc_queue_1,ap_z_ptdma
#endif
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean,calc_mean_flow,&
                                 calc_mean_flow_easy,chk_ibm, &
                                 bulk_mean_ibm,bulk_forcing_ibm,ibm_impulse,ibm_impulse_proj
#if defined(_OPENACC)
  use mod_utils          , only: device_memory_footprint
#endif
  use mod_types
  use omp_lib
!*****IBM******!
  use mod_ibm
  use mod_thakkar           
!**************!
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), target, allocatable, dimension(:,:,:) :: u,v,w,p,pp
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(0:1,3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  integer , dimension(3) :: n_z_d
  real(rp), allocatable, dimension(:,:,:) :: ap_d,cp_d
  logical :: is_ptdma_update_p
  real(rp) :: normfftp
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
#else
  integer    , dimension(2,2) :: arrplanu,arrplanv,arrplanw
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
  !
  real(rp) :: dt,dti,dt_cfl,time,dtrk,dtrki,divtot,divmax
  real(rp) :: velmax_solid(3),divtot_m,divmax_m,ubulk_fl
  logical  :: is_ibm_force,is_ibm_map,is_first_stage,is_zero_vel
  integer  :: iout_wall = 0 ! interval of the per-point IB wall outputs
  character(len=6) :: impart
  real(rp), allocatable, dimension(:,:,:) :: imp_map ! column-integrated IB force (x,y) of the lower wall
  real(rp) :: imp_ibm(2,2,3),imp_leak(2,2,3),pi_drive(3) ! IB force: (lower/upper half, fluid/solid, component)
  real(rp) :: dti_step ! 1/dt of the last time step (dti may change before the output)
  integer :: irk,istep,iunit
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         xc_g,xf_g,yc_g,yf_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !
  type(scalar), target, allocatable, dimension(:) :: scalars
  type(scalar), pointer :: s
  real(rp) :: meanscal
  real(rp) :: gain
  logical  :: is_intrinsic
  real(rp), allocatable, dimension(:) :: fs
  integer :: iscal,is,nscalibm
  real(rp), allocatable, dimension(:,:) :: swz ! IB wall value of each scalar per local level k: bcscal_ibm_top above the mid-plane
  real(rp), allocatable, dimension(:)   :: swz0 ! zero wall value (velocity)
  !
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
  !
  real(rp) :: dt12,dt12av,dt12min,dt12max
  real(rp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=3  ) :: scalnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename
  type :: arr_ptr
    real(rp), pointer, contiguous , dimension(:,:,:) :: arr
  end type arr_ptr
  type(arr_ptr)    , allocatable, dimension(:) ::   io_vars
  character(len=10), allocatable, dimension(:) :: c_io_vars
  integer :: k,kk
  logical :: is_done,kill
!********IBM*******!
  integer               :: ibm_i,np_mpi
  integer               :: size_u,size_v,size_w,size_s
  integer,allocatable   :: fri_u_ga(:),fri_v_ga(:),fri_w_ga(:),sca_ga(:)
  integer,allocatable   :: displ_u(:),displ_v(:),displ_w(:),displ_s(:)
  logical,allocatable   :: mask_u(:,:,:)
  logical,allocatable   :: mask_v(:,:,:)
  logical,allocatable   :: mask_w(:,:,:)
  logical,allocatable   :: mask_s(:,:,:)!(i,j,k)
  real(rp),allocatable  :: lap_u(:,:,:)
  real(rp),allocatable  :: lap_v(:,:,:)
  real(rp),allocatable  :: lap_w(:,:,:)
  real(rp),allocatable  :: lap_s(:,:,:)!(i,j,k)
  integer ,allocatable, dimension(:,:) :: idx_u,idx_v,idx_w,idx_s ! band points, see `set_ibm_band`
  real(rp),allocatable, dimension(:)   :: lapb_u,lapb_v,lapb_w,lapb_s

  real(rp),allocatable  :: A_u(:,:,:)
  real(rp),allocatable  :: A_v(:,:,:)
  real(rp),allocatable  :: A_w(:,:,:)
  real(rp),allocatable  :: A_s(:,:,:)!(i,j,k)
  real(rp),allocatable  :: B_u(:,:,:)
  real(rp),allocatable  :: B_v(:,:,:)
  real(rp),allocatable  :: B_w(:,:,:)
  real(rp),allocatable  :: B_s(:,:,:)!(i,j,k)

  real(rp)              :: mean_u,mean_v,mean_w
  ! fric-pre arrays
  logical,allocatable   :: band_u(:,:,:)
  logical,allocatable   :: band_v(:,:,:)
  logical,allocatable   :: band_w(:,:,:)
  logical,allocatable   :: band_s(:,:,:)

  real(rp),allocatable  :: fri_u(:,:)
  real(rp),allocatable  :: fri_v(:,:)
  real(rp),allocatable  :: fri_w(:,:)
  real(rp),allocatable  :: p_grad(:,:)
  real(rp),allocatable,target  :: fri_u_g(:,:),fri_v_g(:,:),fri_w_g(:,:),p_grad_g(:,:)
  real(rp),allocatable,target  :: dum_g(:,:)
  real(rp),pointer      :: globu(:,:),globv(:,:),globw(:,:),globp(:,:)
  real(rp),allocatable :: frifv(:,:,:,:),frifvout(:,:)
  character (len=7)   :: mpiranknum 
  real(rp) :: dtrk_last = 0._rp,p_wf,pp_wf,fac_wf ! exact wall force per cell: previous stage's dt_rk, pressure levels
  logical  :: is_wall_heat = .false.,is_per(3)
  real(rp), allocatable :: frhout(:,:),qwall(:),ze_sdf(:) ! wall heat flux per cell; level-set z coordinates
  real(rp) :: ibm_idx_off(12)
!*******HMAP********!,
  !for init grid
  real(rp)  :: hmax_g,hch_g,hmin_g
!******************!
!*******CHT*********!
  ! sub-cell geometry of the solid for the scalars with conjugate heat transfer (see `set_cht_geometry`)
  real(rp), allocatable, dimension(:,:,:) :: thx,thy,thz,phi_cht
  logical  :: is_cht_any
  real(rp) :: alpha_max_eff
  ! interfacial heat flux of CHT scalars (`qw_ibm`): interface area per unit volume at the cell centres
  real(rp), allocatable, dimension(:,:,:) :: adens
  real(rp) :: area_ib
  ! level-set geometry (`geometry = 'sdf'`): level set at the cell centres (two halo layers) and at the
  ! cell centres/velocity positions (one halo layer)
  logical  :: is_sdf
  real(rp), allocatable, dimension(:,:,:) :: ps_sdf,phi_s_sdf,phi_u_sdf,phi_v_sdf,phi_w_sdf
!******************!
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np_mpi,ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  iout_wall = merge(iout1d,ibm_iout_wall,ibm_iout_wall < 0)
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  if(is_poisson_pcr_tdma) then
#if defined(_OPENACC)
    n_z_d(:) = ap_z_ptdma%shape(:)
#else
    n_z_d(:) = dinfo_ptdma%zsz(:)
#endif
    allocate(ap_d(n_z_d(1),n_z_d(2),n_z_d(3)), &
             cp_d(n_z_d(1),n_z_d(2),n_z_d(3)))
  end if
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(xc_g(0:ng(1)+1),xf_g(0:ng(1)+1), &
           yc_g(0:ng(2)+1),yf_g(0:ng(2)+1))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  if(is_impdiff) then
    allocate(lambdaxyu(n_z(1),n_z(2)), &
             lambdaxyv(n_z(1),n_z(2)), &
             lambdaxyw(n_z(1),n_z(2)))
    allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
             av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
             aw(n_z(3)),bw(n_z(3)),cw(n_z(3)))
    allocate(rhsbu%x(n(2),n(3),0:1), &
             rhsbu%y(n(1),n(3),0:1), &
             rhsbu%z(n(1),n(2),0:1), &
             rhsbv%x(n(2),n(3),0:1), &
             rhsbv%y(n(1),n(3),0:1), &
             rhsbv%z(n(1),n(2),0:1), &
             rhsbw%x(n(2),n(3),0:1), &
             rhsbw%y(n(1),n(3),0:1), &
             rhsbw%z(n(1),n(2),0:1))
  end if
  !
  allocate(scalars(nscal))
  call initialize_scalars(scalars,nscal,n,n_z)
  allocate(fs(nscal))
  !$acc enter data copyin(scalars(:))
  is_cht_any = .false.
  if(nscal > 0) is_cht_any = any(scalars(:)%is_cht)
  is_sdf = trim(geometry) == 'sdf'
  if(is_ibm_damp .and. .not.(is_ibm .and. is_sdf)) then
    if(myid == 0) print*, 'ERROR: is_ibm_damp needs the level-set geometry (is_ibm = T, geometry = ''sdf'').'
    call MPI_FINALIZE(ierr)
    error stop
  end if
  !
  ! largest scalar diffusivity for the time step limit of explicit diffusion
  ! (inside the solid of CHT scalars: conductivity over heat capacity)
  !
  alpha_max_eff = alpha_max
  if(nscal > 0 .and. .not.is_impdiff) alpha_max_eff = max(alpha_max,maxval(scalars(:)%alpha_max))
  !
  if(use_hmap)then
    call read_thakkar_nfo(myid)
    call read_thakkar_bin(myid)
    print*, "Height Map is on use..." 
  endif
  if(is_debug) then
    if(myid == 0) print*, '*** CaNS build information ***'
    if(myid == 0) print*, 'Compiler version: ', compiler_version()
    if(myid == 0) print*, 'Compiler options: ', compiler_options()
    block
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
      integer :: ilen
      call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
      if(myid == 0) print*, 'MPI version: ', trim(mpi_version)
    end block
    if(myid == 0) print*, ''
  end if
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  hmax_g=0._rp; hch_g=0._rp; hmin_g=0._rp
  if(is_ibm)then
    ! roughness height range of the z-walls, from the actual wall shape
    call wall_height_range(ng,l,ibm_direction,amp_l,n_wave,l_0,phase_l,hmin_g,hmax_g, &
                           hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha)
    if(use_hmap .and. ibm_direction(0,3)) hmax_g = max(hmax_g,l_0(0,3)+hmap_max_tha) ! exact max of the bilinear map
    if(myid==0)print*,"IBM roughness height range (z-walls): hmin, hmax = ",hmin_g,hmax_g
  endif
  if(gtype/=5)then !5 is the rough wall IBM aware cluster two end function 
    call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g,cbcpre(0,3)//cbcpre(1,3) == 'PP')
  else
    if(ibm_hch > 0._rp)then
      hch_g=ibm_hch
    elseif(use_hmap)then
      hch_g=hch ! from the height-map info file
    else
      if(myid==0)print*,"ERROR: the rough-wall grid (gtype = 5) requires &ibm ibm_hch > 0."
      call MPI_FINALIZE(ierr)
      error stop
    endif
    if(myid==0)print*,"rough-wall grid: hmax, hch = ",hmax_g,hch_g
    if(myid==0 .and. hmax_g >= l(3)/2._rp) &
      print*,"Warning: the roughness reaches the channel centre (hmax >= l(3)/2)."
    if(myid==0 .and. hch_g > l(3)/2._rp) &
      print*,"Warning: ibm_hch is larger than half the domain height."
    call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g,cbcpre(0,3)//cbcpre(1,3) == 'PP',&
                        visc,dzp_fine,dzp_coarse,hch_g,hmax_g,myid)
  endif
  do kk=0,ng(1)+1
    xc_g(kk) = (kk-0.5_rp)*dl(1)
    xf_g(kk) = (kk-1.0_rp)*dl(1)
  end do
  do kk=0,ng(2)+1
    yc_g(kk) = (kk-0.5_rp)*dl(2)
    yf_g(kk) = (kk-1.0_rp)*dl(2)
  end do
  if(myid == 0) then
    open(newunit=iunit,file=trim(datadir)//'geometry.out',status='replace')
    write(iunit,*) ng(1),ng(2),ng(3)
    write(iunit,*) l(1),l(2),l(3)
    close(iunit)
    call save_grid(datadir,'grid',ng,zc_g,zf_g,dzc_g,dzf_g)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,dl,dli,l) async
  !$acc enter data copyin(gacc) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc data copy(ng) async
  !$acc parallel loop default(present) async
  do k=0,ng(3)+1
    dzci_g(k) = dzc_g(k)**(-1)
    dzfi_g(k) = dzf_g(k)**(-1)
  end do
  !$acc end data
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc parallel loop default(present) async
  do k=0,n(3)+1
    grid_vol_ratio_c(k) = dl(1)*dl(2)*dzc(k)/(l(1)*l(2)*l(3))
    grid_vol_ratio_f(k) = dl(1)*dl(2)*dzf(k)/(l(1)*l(2)*l(3))
  end do
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc wait
  allocate(swz(0:n(3)+1,max(nscal,1)),swz0(0:n(3)+1))
  swz(:,:) = 0._rp; swz0(:) = 0._rp
  do is=1,nscal
    do k=0,n(3)+1
      swz(k,is) = merge(bcscal_ibm_top(is),bcscal_ibm(is),zc(k) > zf_g(0)+l(3)/2) ! above the domain mid-plane
    end do
  end do
  !$acc enter data copyin(swz,swz0) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  if(is_poisson_pcr_tdma) then
    !$acc enter data create(ap_d,cp_d) async
  end if
  !$acc wait
  if(is_impdiff) then
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
                    lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu,rhsbu%x,rhsbu%y,rhsbu%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
                    lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv,rhsbv%x,rhsbv%y,rhsbv%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
                    lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw,rhsbw%x,rhsbw%y,rhsbw%z)
    do iscal=1,nscal
      s => scalars(iscal)
      call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,s%cbc,s%bc, &
                      s%lambdaxy,['c','c','c'],s%a,s%b,s%c,s%arrplan,s%normfft, &
                      s%rhsb%x,s%rhsb%y,s%rhsb%z)
    end do
    if(is_impdiff_1d) then
      deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
      call fftend(arrplanu)
      call fftend(arrplanv)
      call fftend(arrplanw)
      deallocate(rhsbu%x,rhsbu%y,rhsbv%x,rhsbv%y,rhsbw%x,rhsbw%y)
      do iscal=1,nscal
        s => scalars(iscal)
        deallocate(s%rhsb%x,s%rhsb%y)
        deallocate(s%lambdaxy)
        call fftend(s%arrplan)
      end do
    end if
    !$acc enter data copyin(lambdaxyu,au,bu,cu,lambdaxyv,av,bv,cv,lambdaxyw,aw,bw,cw) async
    !$acc enter data copyin(rhsbu,rhsbu%x,rhsbu%y,rhsbu%z) async
    !$acc enter data copyin(rhsbv,rhsbv%x,rhsbv%y,rhsbv%z) async
    !$acc enter data copyin(rhsbw,rhsbw%x,rhsbw%y,rhsbw%z) async
    do iscal=1,nscal
      s => scalars(iscal)
      !$acc enter data copyin(s%lambdaxy,s%a,s%b,s%c) async
      !$acc enter data copyin(s%rhsb,s%rhsb%x,s%rhsb%y,s%rhsb%z) async
    end do
    !$acc wait
  end if
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call set_cufft_wspace(pack(arrplanu,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanv,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanw,.true.),istream_acc_queue_1)
  end if
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z,nscal)/(1._sp*1024**3), ' ***'
#endif
  if(is_debug_poisson) then
    call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,l,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                            nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
  end if
  !
  allocate(io_vars(4+nscal),c_io_vars(4+nscal)) ! u,v,w,p,scalars
  io_vars(1)%arr => u; c_io_vars(1) = 'u'
  io_vars(2)%arr => v; c_io_vars(2) = 'v'
  io_vars(3)%arr => w; c_io_vars(3) = 'w'
  io_vars(4)%arr => p; c_io_vars(4) = 'p'
  do iscal=1,nscal
    write(scalnum,'(i3.3)') iscal
    io_vars(4+iscal)%arr => scalars(iscal)%val; c_io_vars(4+iscal) = 's_'//scalnum
  end do
  !
  is_ptdma_update_p = .true.
  !
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,visc,is_forced,velf,bforce,is_wallturb,u,v,w,p)
    do iscal=1,nscal
      s => scalars(iscal)
      call initscal(s%ini,s%bc,ng,lo,l,dl,zc,dzf,s%alpha,s%is_forced,s%scalf,s%val)
    end do
    if(myid == 0) print*, '*** Initial condition successfully set ***'
  else
    do is=1,4+nscal
      call load_one('r',trim(datadir)//'fld_'//trim(c_io_vars(is))//trim(io_ext), &
                    MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep)
    end do
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !$acc enter data copyin(u,v,w,p,dudtrko,dvdtrko,dwdtrko) create(pp)
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  do iscal=1,nscal
    s => scalars(iscal)
    !$acc enter data copyin(s%val,s%dsdtrko) async(1)
    call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
  end do
  !$acc wait
  !
  ! 3D time averages
  !
  if(iaverage > 0) call stats3d_init(n,nscal,restart,trim(datadir),ng,lo,hi,myid)
  !
  call chkdt(n,dl,dzci,dzfi,visc,alpha_max_eff,u,v,w,dt_cfl)
  dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
  dti = 1./dt
  kill = .false.
!**********IBM-PART***********!
  allocate(fri_u_ga(np_mpi),fri_v_ga(np_mpi),fri_w_ga(np_mpi),sca_ga(np_mpi))
  allocate(displ_u(np_mpi),displ_v(np_mpi),displ_w(np_mpi),displ_s(np_mpi))
  fri_u_ga(:)=0;fri_v_ga(:)=0;fri_w_ga(:)=0;sca_ga(:)=0;
  displ_u(:)=0;displ_v(:)=0;displ_w(:)=0;displ_s(:)=0;
  if(is_ibm)then
    ! we initalize the ibm coef. here
    allocate(mask_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(mask_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(mask_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(mask_s(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    ! isInbody ---> false 
    ! so they start w/ all fluid
    print*, "***IBM coefficients are initializing***"
    mask_u = .false.
    mask_v = .false.
    mask_w = .false.
    mask_s = .false.
  endif
  !
  if(is_ibm)then
    allocate(lap_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(lap_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(lap_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(lap_s(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    lap_s=0._rp
    lap_u=0._rp
    lap_v=0._rp
    lap_w=0._rp
    allocate(band_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(band_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(band_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(band_s(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    band_u=.false.
    band_v=.false.
    band_w=.false.
    band_s=.false.
  endif
  !
  allocate(A_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(A_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(A_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(A_s(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(B_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(B_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(B_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(B_s(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  A_u=1._rp
  A_v=1._rp
  A_w=1._rp
  A_s=1._rp
  B_u=1._rp
  B_v=1._rp
  B_w=1._rp
  B_s=1._rp
  !
#if defined(_OPENACC)
  !$acc enter data copyin(A_u,A_v,A_w,A_s,B_u,B_v,B_w,B_s)
#endif  
  
  if(is_ibm .and. is_sdf)then
    !
    ! level-set geometry: masks from the sign of the level set at the respective positions
    !
    if(myid == 0) print*, '*** IBM geometry from the level set `'//trim(sdf_file)//'` ***'
    allocate(ps_sdf(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
    allocate(phi_s_sdf(0:n(1)+1,0:n(2)+1,0:n(3)+1),phi_u_sdf(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             phi_v_sdf(0:n(1)+1,0:n(2)+1,0:n(3)+1),phi_w_sdf(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    call sdf_read(sdf_file,ng,lo,n,zc,ps_sdf)
    call sdf_staggered(n,dzc,dzf,ps_sdf,phi_s_sdf,phi_u_sdf,phi_v_sdf,phi_w_sdf)
    call sdf_mask(phi_u_sdf,mask_u)
    call sdf_mask(phi_v_sdf,mask_v)
    call sdf_mask(phi_w_sdf,mask_w)
    call sdf_mask(phi_s_sdf,mask_s)
    if(is_ibm_damp) call ibm_damp_init(n,ng,dl(1),ibm_damp_d,phi_s_sdf,phi_u_sdf,phi_v_sdf,phi_w_sdf,myid)
#if defined (_OPENACC)
    !$acc enter data copyin(mask_u,mask_v,mask_w,mask_s)
#endif
  else if(is_ibm)then
    ! we fill the ibm masks here
    print*, "***IBM masks are deploying***"
    call set_ibm_staircase(lo,mask_u,1,0,0,n,l,dl,&
    ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf)
    call set_ibm_staircase(lo,mask_v,0,1,0,n,l,dl,&
    ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf)
    call set_ibm_staircase(lo,mask_w,0,0,1,n,l,dl,&
    ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf)
    call set_ibm_staircase(lo,mask_s,0,0,0,n,l,dl,&
    ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf)
#if defined (_OPENACC)
    !$acc enter data copyin(mask_u,mask_v,mask_w,mask_s)
#endif
  endif
  if(is_ibm)then
    ! column-integrated IB force, also needed by the IB force output (`iout0d`) with `geometry = 'sdf'`
    allocate(imp_map(n(1),n(2),2))
    !$acc enter data create(imp_map)
  endif
  if(is_ibm)then
    print*, "***IBM coefficients are deploying***"
    if(is_sdf)then
    call set_ibm_2nd_sdf(n,dl,dzc,dzf,0,phi_u_sdf,mask_u,lap_u,band_u)
    call set_ibm_2nd_sdf(n,dl,dzc,dzf,0,phi_v_sdf,mask_v,lap_v,band_v)
    call set_ibm_2nd_sdf(n,dl,dzc,dzf,1,phi_w_sdf,mask_w,lap_w,band_w)
    call set_ibm_2nd_sdf(n,dl,dzc,dzf,0,phi_s_sdf,mask_s,lap_s,band_s)
    else
    call set_ibm_2nd(lo,mask_u,lap_u,1,0,0&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf,dzc,dzf,use_hmap,band_u)
    call set_ibm_2nd(lo,mask_v,lap_v,0,1,0&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf,dzc,dzf,use_hmap,band_v)
    call set_ibm_2nd(lo,mask_w,lap_w,0,0,1&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf,dzc,dzf,use_hmap,band_w)
    call set_ibm_2nd(lo,mask_s,lap_s,0,0,0,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,&
                            hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,zc,zf,dzc,dzf,use_hmap,band_s)
    endif
    call set_ibm_band(lap_u,idx_u,lapb_u)
    call set_ibm_band(lap_v,idx_v,lapb_v)
    call set_ibm_band(lap_w,idx_w,lapb_w)
    call set_ibm_band(lap_s,idx_s,lapb_s)
    !$acc enter data copyin(idx_u,idx_v,idx_w,idx_s,lapb_u,lapb_v,lapb_w,lapb_s)
  endif
  if(is_ibm)then
    allocate(fri_u(count(band_u(1:n(1),1:n(2),1:n(3))),12))! we do this way so we exlude ghost cells and only focus on the inner cells
    allocate(fri_v(count(band_v(1:n(1),1:n(2),1:n(3))),12))! init the local arrays
    allocate(fri_w(count(band_w(1:n(1),1:n(2),1:n(3))),12))
    allocate(p_grad(count(band_s(1:n(1),1:n(2),1:n(3))),12))
    allocate(dum_g(1,12))! dummy array 
    ! columns 7-9 of the force output hold rank-local indices; written as global ones
    ibm_idx_off(:) = 0._rp
    ibm_idx_off(7:9) = lo(1:3)-1
    ! now lets init the global arrays
    size_u=size(fri_u,1);size_v=size(fri_v,1);size_w=size(fri_w,1);size_s=size(p_grad,1)
    call MPI_GATHER(size_u,1,MPI_INT,fri_u_ga,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(size_v,1,MPI_INT,fri_v_ga,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(size_w,1,MPI_INT,fri_w_ga,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(size_s,1,MPI_INT,sca_ga,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    if(myid==0)then
      allocate(fri_u_g(sum(fri_u_ga),12),fri_v_g(sum(fri_v_ga),12),fri_w_g(sum(fri_w_ga),12),p_grad_g(sum(sca_ga),12))
      do ibm_i=2,np_mpi
        displ_u(ibm_i)=displ_u(ibm_i-1)+fri_u_ga(ibm_i-1)
        displ_v(ibm_i)=displ_v(ibm_i-1)+fri_v_ga(ibm_i-1)
        displ_w(ibm_i)=displ_w(ibm_i-1)+fri_w_ga(ibm_i-1)
        displ_s(ibm_i)=displ_s(ibm_i-1)+sca_ga(ibm_i-1)
      end do
    endif
    !
    ! set-up of the per-point wall outputs (every `iout_wall` steps; none with ibm_iout_wall = 0)
    !
    if(iout_wall > 0) then
    if(is_sdf)then
    call calc_grad_dist_sdf(fri_u ,lo,n,dl,zc,zf,band_u,1,0,0,phi_u_sdf)
    call calc_grad_dist_sdf(fri_v ,lo,n,dl,zc,zf,band_v,0,1,0,phi_v_sdf)
    call calc_grad_dist_sdf(fri_w ,lo,n,dl,zc,zf,band_w,0,0,1,phi_w_sdf)
    call calc_grad_dist_sdf(p_grad,lo,n,dl,zc,zf,band_s,0,0,0,phi_s_sdf)
    else
    call calc_grad_dist(fri_u,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_u,1,0,0,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,dzf,dzc)
    call calc_grad_dist(fri_v,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_v,0,1,0,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,dzf,dzc)
    call calc_grad_dist(fri_w,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_w,0,0,1,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,dzf,dzc)
    call calc_grad_dist(p_grad,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_s,0,0,0,hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,dzf,dzc)                    
    endif
    endif
    !
    ! exact wall force and wall heat flux per cell (`FV-fric_*`, `FV-heat_*`): with the level-set
    ! geometry, the wall points and the sub-cell geometry come from the level set (ps_sdf, ze_sdf)
    !
    allocate(frifv(0:n(1),0:n(2),0:n(3),19))
    frifv(:,:,:,:) = 0._rp
    if(iout_wall > 0) then
    if(.not.is_sdf)then
    call init_wall_force(n,lo,l,dl,zc,zf,dzc,dzf,mask_u,mask_v,mask_w,mask_s,ibm_direction,amp_l,n_wave,l_0,phase_l, &
                         hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha)
    call init_fric_cubes(wf_cell,mask_s,lo,dl,zc,zf,dzc,&
                                        dzf,xf_g,yc_g,0,0,0,ibm_direction,amp_l,&
                                        n_wave,l_0,phase_l,n,l,hmap_tha,lx_tha,ly_tha,&
                                        nx_hmap_tha,ny_hmap_tha,frifv)
    else
      allocate(ze_sdf(-1:n(3)+2))
      ze_sdf(0:n(3)+1) = zc(0:n(3)+1)
      ze_sdf(-1) = 2._rp*ze_sdf(0)-ze_sdf(1)
      ze_sdf(n(3)+2) = 2._rp*ze_sdf(n(3)+1)-ze_sdf(n(3))
      call init_wall_force(n,lo,l,dl,zc,zf,dzc,dzf,mask_u,mask_v,mask_w,mask_s,ibm_direction,amp_l,n_wave,l_0,phase_l, &
                           ps=ps_sdf,ze=ze_sdf)
      call init_fric_cubes(wf_cell,mask_s,lo,dl,zc,zf,dzc,&
                                          dzf,xf_g,yc_g,0,0,0,ibm_direction,amp_l,&
                                          n_wave,l_0,phase_l,n,l,output=frifv,ps=ps_sdf,ze=ze_sdf)
    endif
    endif
    allocate(frifvout(count(frifv(1:n(1),1:n(2),1:n(3),10)>tiny(1._rp)),26))
    allocate(frhout(size(frifvout,1),7+3*nscal),qwall(2*nscal))
    is_wall_heat = nscal > 0 .and. iout_wall > 0
    is_per(:) = cbcpre(0,:)//cbcpre(1,:) == 'PP'
    !
#if defined (_OPENACC)
    !$acc enter data copyin(fri_u,fri_v,fri_w,p_grad)
    !$acc enter data copyin(band_s,band_u)
#endif
    write(mpiranknum,'(i7.7)') myid !for naming the rank based output
  endif
  if(is_ibm .and. is_impdiff_1d)then
    ! the IB-aware z-implicit solve (`solve_helmholtz_ibm`) needs whole z columns on each task
    block
      logical :: is_ok
      is_ok = n(3) == ng(3)
      call MPI_ALLREDUCE(MPI_IN_PLACE,is_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
      if(.not.is_ok) then
        if(myid == 0) print*, "ERROR: is_ibm with is_impdiff_1d requires a domain not decomposed along z &
                              &(ipencil_axis = 3, or dims(2) = 1)."
        call MPI_FINALIZE(ierr)
        error stop
      end if
      is_ok = .not.any(cbcvel(0,3,:)//cbcvel(1,3,:) == 'PP')
      do iscal=1,nscal
        is_ok = is_ok.and.scalars(iscal)%cbc(0,3)//scalars(iscal)%cbc(1,3) /= 'PP'
      end do
      if(.not.is_ok) then
        if(myid == 0) print*, "ERROR: is_ibm with is_impdiff_1d requires non-periodic z boundaries."
        call MPI_FINALIZE(ierr)
        error stop
      end if
    end block
  endif
  !
  ! conjugate heat transfer: sub-cell geometry of the solid (interface position along the cell
  ! centre-to-centre segments and fluid volume fractions, exact with the 2nd-order IBM, on the cell
  ! faces with the staircase IBM) and the implicit diffusion solver
  !
  if(is_cht_any .or. (is_ibm .and. is_wall_heat))then ! also for the wall heat flux of scalars without CHT
    allocate(thx(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             thy(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             thz(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             phi_cht(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    if(is_sdf)then
      call set_cht_geometry_sdf(n,dl,zc,dzc,dzf,ps_sdf,.true.,thx,thy,thz,phi_cht)
    else
      call set_cht_geometry(lo,n,l,dl,zc,dzc,dzf,mask_s,.true.,ibm_direction,amp_l,n_wave,l_0,phase_l,&
                            hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha,thx,thy,thz,phi_cht)
    endif
    !$acc enter data copyin(thx,thy,thz,phi_cht)
    if(is_cht_any) then
    call cht_init(n,ng,l,dzci_g,dzfi_g,thx,thy,thz,phi_cht,scalars)
    if(myid == 0) print*, '*** CHT geometry set up ***'
    end if
  endif
  !
  ! interfacial heat flux of CHT scalars (`qw_ibm`, heat per unit interface area): interface area density
  ! from the level set, a smoothed delta of half-width `qw_ibm_eps` grid spacings (see `sdf_interface_density`)
  !
  if(nscal > 0) then
    if(any(qw_ibm(:) /= 0._rp)) then
      if(.not.is_sdf) then
        if(myid == 0) print*, 'ERROR: `qw_ibm /= 0` requires `geometry = ''sdf''`.'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
      allocate(adens(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call sdf_interface_density(n,dl,zc,dzf,ps_sdf,qw_ibm_eps,adens,area_ib)
      if(myid == 0) print*, '*** interfacial heat flux: interface area from the level set = ', area_ib, &
                            ', delta half-width ', qw_ibm_eps, ' grid spacings ***'
    end if
  end if
  if(.not.allocated(adens)) then
    allocate(adens(0:0,0:0,0:0)); adens(:,:,:) = 0._rp ! unused
  end if
  !$acc enter data copyin(adens)
  !
  ! optional output of the IBM/CHT geometry fields (interior cells, x fastest, serial runs only)
  !
  if(dump_ibm_geometry .and. is_ibm)then
    if(np_mpi == 1)then
      open(newunit=iunit,file=trim(datadir)//'ibm_geometry.bin',access='stream',form='unformatted',status='replace')
      write(iunit) lap_u(1:n(1),1:n(2),1:n(3)),lap_v(1:n(1),1:n(2),1:n(3)), &
                   lap_w(1:n(1),1:n(2),1:n(3)),lap_s(1:n(1),1:n(2),1:n(3))
      if(is_cht_any) write(iunit) thx(1:n(1),1:n(2),1:n(3)),thy(1:n(1),1:n(2),1:n(3)), &
                                  thz(1:n(1),1:n(2),1:n(3)),phi_cht(1:n(1),1:n(2),1:n(3))
      close(iunit)
      if(myid == 0) print*, '*** IBM geometry written to data/ibm_geometry.bin (lap_u,v,w,s; thx,y,z,phi) ***'
    else
      if(myid == 0) print*, 'Warning: `dump_ibm_geometry` is only supported for serial runs.'
    endif
  endif
  if(is_ibm) deallocate(lap_u,lap_v,lap_w,lap_s) ! only the band list (`set_ibm_band`) is kept
!********IBM-END*************!
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait ! not needed but to prevent possible future issues
  !$acc update self(u,v,w,p)
  do iscal=1,nscal
    !$acc update self(scalars(iscal)%val)
  end do
  ! write the inital condition
  if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
    include 'out1d.h90'
    call out1d_chan(trim(datadir)//'turb_stats_'//fldnum//'.out',ng,lo,hi,3,l,dl,zc_g,u,v,w)
    if(is_ibm) &
      call out1d_chan_ibm(trim(datadir)//'turb_stats_ibm_'//fldnum//'.out',ng,lo,hi,l,dl,zc_g,u,v,w, &
                          mask_u,mask_v,mask_w,mask_s)
  end if
  if(is_ibm.and.iout_wall > 0.and.mod(istep,max(iout_wall,1)) == 0) then
      call calc_shear_st(fri_u,u,visc)
      call calc_shear_st(fri_v,v,visc)
      call calc_shear_st(fri_w,w,visc)
      call get_wall_pres(p,p_grad,lo,dl,zc,zf,mask_s,band_s)
      if(myid==0)then
        globu=>fri_u_g
        globv=>fri_v_g
        globw=>fri_w_g
        globp=>p_grad_g
      else
        globu=>dum_g
        globv=>dum_g
        globw=>dum_g
        globp=>dum_g
      endif
      !$acc update self(fri_u,fri_v,fri_w,p_grad)
      do ibm_i=1,12
        call MPI_GATHERV(fri_u(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_u,MPI_REAL_RP,globu(:,ibm_i),fri_u_ga,displ_u,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(fri_v(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_v,MPI_REAL_RP,globv(:,ibm_i),fri_v_ga,displ_v,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(fri_w(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_w,MPI_REAL_RP,globw(:,ibm_i),fri_w_ga,displ_w,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(p_grad(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_s,MPI_REAL_RP,globp(:,ibm_i),sca_ga,displ_s,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      end do
      if(myid==0)then
        call write_data(trim(datadir)//"fric-u_"//fldnum//".out",myid,fri_u_g)
        call write_data(trim(datadir)//"fric-v_"//fldnum//".out",myid,fri_v_g)
        call write_data(trim(datadir)//"fric-w_"//fldnum//".out",myid,fri_w_g)
        call write_data(trim(datadir)//"pressure"//fldnum//".out",myid,p_grad_g)
      endif
  end if
  if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
    include 'out2d.h90'
  end if
  if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
    include 'out3d.h90'
  end if
!*****************************
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  is_first_stage = .true.
  do while(.not.is_done)
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()
    end if
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    tauxo(:,:) = 0.; tauyo(:,:) = 0.; tauzo(:,:) = 0.
    dpdl(:)     = 0.
    fs(1:nscal) = 0.
    is_ibm_map   = is_ibm.and.iout_wall > 0.and.mod(istep,max(iout_wall,1)) == 0
    is_ibm_force = is_ibm.and.iout0d > 0.and.mod(istep,max(iout0d,1)) == 0
    is_ibm_force = is_ibm_force.or.is_ibm_map
    call timer_tic(trim(merge('step_ibm_out','step        ',is_ibm_force)),nvtx_color='g',nvtx_gpu_stream=-1)
    if(is_ibm_force) then
      imp_ibm(:,:,:)  = 0.
      imp_leak(:,:,:) = 0.
      imp_map(:,:,:)  = 0.
      !$acc update device(imp_map)
    end if
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      !****//IBM\\****!
      is_zero_vel = is_ibm .and. (is_first_stage .or. trim(ibm_proj) == 'none')
      if(is_ibm)then
        call timer_tic('ibm_stage_start',nvtx_color='r',nvtx_gpu_stream=-1)
        if(is_ibm_force) then ! velocity the projection left in the solid
          call ibm_impulse('leak',n,lo,ng,grid_vol_ratio_f,mask_u,A_u,B_u,0._rp,u,imp_leak(:,:,1),imp_map(:,:,1))
          call ibm_impulse('leak',n,lo,ng,grid_vol_ratio_f,mask_v,A_v,B_v,0._rp,v,imp_leak(:,:,2),imp_map(:,:,2))
          call ibm_impulse('leak',n,lo,ng,grid_vol_ratio_c,mask_w,A_w,B_w,0._rp,w,imp_leak(:,:,3))
        end if
        if(is_ibm_map) then ! leakage = -dt_rk*grad(pp) of the previous stage ('none' only)
          fac_wf = merge(dtrk_last,0._rp,trim(ibm_proj) == 'none')
          if(fac_wf > 0._rp) then
            call calc_mean_flow(grid_vol_ratio_f,mask_s,pp,pp_wf,n)
            call wf_gather(1,n,dli,dzfi,grid_vol_ratio_f,mask_u,fac_wf,pp,pp_wf,0._rp)
            call wf_gather(2,n,dli,dzfi,grid_vol_ratio_f,mask_v,fac_wf,pp,pp_wf,0._rp)
            call wf_gather(3,n,dli,dzfi,grid_vol_ratio_f,mask_w,fac_wf,pp,pp_wf,0._rp)
          end if
          call wf_scatter('leak',1,dli,dzci,grid_vol_ratio_f,mask_u,A_u,B_u,fac_wf,0._rp,pp,u)
          call wf_scatter('leak',2,dli,dzci,grid_vol_ratio_f,mask_v,A_v,B_v,fac_wf,0._rp,pp,v)
          call wf_scatter('leak',3,dli,dzci,grid_vol_ratio_c,mask_w,A_w,B_w,fac_wf,0._rp,pp,w)
        end if
        !
        ! zero the solid points; after the first stage the scalar still holds the wall value set
        ! at the end of the previous stage, and the velocity is still zero with 'damped', since
        ! `correc_ibm` does not correct the solid points
        !
        if(is_zero_vel) then
          call apply_ibm_staircase(u,mask_u,dtrk)
          call apply_ibm_staircase(v,mask_v,dtrk)
          call apply_ibm_staircase(w,mask_w,dtrk)
        end if
        if(is_first_stage) then
          do nscalibm=1,nscal
            if(scalars(nscalibm)%is_cht) cycle ! CHT: the scalar is solved inside the solid
            call apply_ibm_staircase_scalar(scalars(nscalibm)%val,&
                mask_s,swz(:,nscalibm))
          end do
        end if
        call timer_toc('ibm_stage_start',nvtx_gpu_stream=-1)
      endif
      !
      ! boundary conditions: the halos are still up to date from the end of the previous stage,
      ! unless the solid points were zeroed
      !
      if(is_zero_vel) then
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      else
        call bounduvw_norm(cbcvel,n,bcvel,is_bound,dl,dzc,dzf,u,v,w)
      end if
      is_first_stage = .false.
      !*******************!
      !****//IBM\\****!
      if(is_ibm)then
        call timer_tic('ibm_calc_a_b',nvtx_color='r',nvtx_gpu_stream=-1)
        call calc_a_b(idx_u,lapb_u,A_u,B_u,dtrk,visc)
        call calc_a_b(idx_v,lapb_v,A_v,B_v,dtrk,visc)
        call calc_a_b(idx_w,lapb_w,A_w,B_w,dtrk,visc)
        call timer_toc('ibm_calc_a_b',nvtx_gpu_stream=-1)
      endif
      if(is_ibm_force) call timer_tic('ibm_force',nvtx_color='m',nvtx_gpu_stream=-1)
      if(is_ibm_force) then
        call ibm_impulse('pre',n,lo,ng,grid_vol_ratio_f,mask_u,A_u,B_u,0._rp,u,imp_ibm(:,:,1),imp_map(:,:,1))
        call ibm_impulse('pre',n,lo,ng,grid_vol_ratio_f,mask_v,A_v,B_v,0._rp,v,imp_ibm(:,:,2),imp_map(:,:,2))
        call ibm_impulse('pre',n,lo,ng,grid_vol_ratio_c,mask_w,A_w,B_w,0._rp,w,imp_ibm(:,:,3))
      end if
      if(is_ibm_map) then
        call wf_scatter('pre',1,dli,dzci,grid_vol_ratio_f,mask_u,A_u,B_u,0._rp,0._rp,p,u)
        call wf_scatter('pre',2,dli,dzci,grid_vol_ratio_f,mask_v,A_v,B_v,0._rp,0._rp,p,v)
        call wf_scatter('pre',3,dli,dzci,grid_vol_ratio_c,mask_w,A_w,B_w,0._rp,0._rp,p,w)
      end if
      if(is_ibm_force) call timer_toc('ibm_force',nvtx_gpu_stream=-1)
      !*******************!
      do iscal=1,nscal
        s => scalars(iscal)
        if(s%is_cht) then
          !
          ! conjugate heat transfer: no wall condition at the interface (A_s/B_s and the wall value unused)
          !
          call rk_scal(rkcoeff(:,irk),n,dli,l,dzci,dzfi,grid_vol_ratio_f,s%alpha,dt,is_bound,u,v,w, &
                       s%is_forced,s%scalf,s%source,s%fluxo,s%dsdtrko,s%val,s%f,A_s,B_s, &
                       swz(:,iscal),scal_mean_grad(iscal),mask_s, &
                       is_cht=s%is_cht,alphas=s%alphas,cs_cf=s%cs_cf,ssource_s=s%source_s,cbcz=s%cbc(:,3), &
                       thx=thx,thy=thy,thz=thz,phi=phi_cht,qw=qw_ibm(iscal),adens=adens)
        else
        if(is_ibm) then
          call timer_tic('ibm_calc_a_b',nvtx_color='r',nvtx_gpu_stream=-1)
          call calc_a_b(idx_s,lapb_s,A_s,B_s,dtrk,s%alpha) ! A_s, B_s depend on the diffusivity
          call timer_toc('ibm_calc_a_b',nvtx_gpu_stream=-1)
        end if
        call rk_scal(rkcoeff(:,irk),n,dli,l,dzci,dzfi,grid_vol_ratio_f,s%alpha,dt,is_bound,u,v,w, &
                     s%is_forced,s%scalf,s%source,s%fluxo,s%dsdtrko,s%val,s%f,A_s,B_s, &
                     swz(:,iscal),scal_mean_grad(iscal),mask_s)
        end if
        if(is_ibm) then
          if(s%is_forced) call bulk_forcing_ibm(n,mask_s,A_s,s%f,s%val)
        else
          call bulk_forcing_s(n,s%is_forced,s%f,s%val)
        end if
        fs(iscal) = fs(iscal) + s%f
        if(is_impdiff .and. .not.s%is_cht) then ! CHT: diffusion already treated implicitly in `rk_scal`
          if(is_ibm .and. is_impdiff_1d) then
            call solve_helmholtz_ibm(n,-0.5*s%alpha*dtrk,s%a,s%b,s%c,s%rhsb%x,s%rhsb%y,s%rhsb%z,is_bound, &
                                     s%cbc,['c','c','c'],mask_s,A_s,swz(:,iscal),.false.,s%val)
          else
            call solve_helmholtz(n,ng,hi,s%arrplan,s%normfft,-0.5*s%alpha*dtrk, &
                                 s%lambdaxy,s%a,s%b,s%c,s%rhsb%x,s%rhsb%y,s%rhsb%z,is_bound,s%cbc,['c','c','c'],s%val)
          end if
        end if
        call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
      end do
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,dt,visc,p, &
              is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f,&
              A_u,A_v,A_w,B_u,B_v,B_w,mask_u,mask_v,mask_w)
      if(is_ibm) then
        call timer_tic('ibm_bulk_forcing',nvtx_color='r',nvtx_gpu_stream=-1)
        if(is_forced(1)) call bulk_forcing_ibm(n,mask_u,A_u,f(1),u)
        if(is_forced(2)) call bulk_forcing_ibm(n,mask_v,A_v,f(2),v)
        if(is_forced(3)) call bulk_forcing_ibm(n,mask_w,A_w,f(3),w)
        call timer_toc('ibm_bulk_forcing',nvtx_gpu_stream=-1)
      else
        call bulk_forcing(n,is_forced,f,u,v,w)
      end if
      if(is_ibm_force) call timer_tic('ibm_force',nvtx_color='m',nvtx_gpu_stream=-1)
      if(is_ibm_force) then
        if(.not.is_impdiff_1d) then
          impart = 'post'
        else
          impart = 'post_s' ! the fluid part follows the implicit solve ('impl')
        end if
        call ibm_impulse(impart,n,lo,ng,grid_vol_ratio_f,mask_u,A_u,B_u,merge(f(1),0._rp,is_forced(1)),u,imp_ibm(:,:,1), &
                         imp_map(:,:,1))
        call ibm_impulse(impart,n,lo,ng,grid_vol_ratio_f,mask_v,A_v,B_v,merge(f(2),0._rp,is_forced(2)),v,imp_ibm(:,:,2), &
                         imp_map(:,:,2))
        call ibm_impulse(impart,n,lo,ng,grid_vol_ratio_c,mask_w,A_w,B_w,merge(f(3),0._rp,is_forced(3)),w,imp_ibm(:,:,3))
      end if
      if(is_ibm_map) then
        call calc_mean_flow(grid_vol_ratio_f,mask_s,p,p_wf,n)
        call wf_gather(1,n,dli,dzfi,grid_vol_ratio_f,mask_u,dtrk,p,p_wf,dtrk*bforce(1)+merge(f(1),0._rp,is_forced(1)))
        call wf_gather(2,n,dli,dzfi,grid_vol_ratio_f,mask_v,dtrk,p,p_wf,dtrk*bforce(2)+merge(f(2),0._rp,is_forced(2)))
        call wf_gather(3,n,dli,dzfi,grid_vol_ratio_f,mask_w,dtrk,p,p_wf,dtrk*bforce(3)+merge(f(3),0._rp,is_forced(3)))
        call wf_scatter(impart,1,dli,dzci,grid_vol_ratio_f,mask_u,A_u,B_u,dtrk,bforce(1),p,u)
        call wf_scatter(impart,2,dli,dzci,grid_vol_ratio_f,mask_v,A_v,B_v,dtrk,bforce(2),p,v)
        call wf_scatter(impart,3,dli,dzci,grid_vol_ratio_c,mask_w,A_w,B_w,dtrk,bforce(3),p,w)
      end if
      if(is_ibm_force) call timer_toc('ibm_force',nvtx_gpu_stream=-1)
      dpdl(:) = dpdl(:) + f(:)
      if(is_impdiff) then
        alpha = -.5*visc*dtrk
        if(is_ibm .and. is_impdiff_1d) then
          call solve_helmholtz_ibm(n,alpha,au,bu,cu,rhsbu%x,rhsbu%y,rhsbu%z,is_bound,cbcvel(:,:,1),['f','c','c'], &
                                   mask_u,A_u,swz0,.not.is_ibm_force,u)
          call solve_helmholtz_ibm(n,alpha,av,bv,cv,rhsbv%x,rhsbv%y,rhsbv%z,is_bound,cbcvel(:,:,2),['c','f','c'], &
                                   mask_v,A_v,swz0,.not.is_ibm_force,v)
          call solve_helmholtz_ibm(n,alpha,aw,bw,cw,rhsbw%x,rhsbw%y,rhsbw%z,is_bound,cbcvel(:,:,3),['c','c','f'], &
                                   mask_w,A_w,swz0,.not.is_ibm_force,w)
        else
          call solve_helmholtz(n,ng,hi,arrplanu,normfftu,alpha, &
                               lambdaxyu,au,bu,cu,rhsbu%x,rhsbu%y,rhsbu%z,is_bound,cbcvel(:,:,1),['f','c','c'],u)
          call solve_helmholtz(n,ng,hi,arrplanv,normfftv,alpha, &
                               lambdaxyv,av,bv,cv,rhsbv%x,rhsbv%y,rhsbv%z,is_bound,cbcvel(:,:,2),['c','f','c'],v)
          call solve_helmholtz(n,ng,hi,arrplanw,normfftw,alpha, &
                               lambdaxyw,aw,bw,cw,rhsbw%x,rhsbw%y,rhsbw%z,is_bound,cbcvel(:,:,3),['c','c','f'],w)
        end if
        if(is_ibm_force) call timer_tic('ibm_force',nvtx_color='m',nvtx_gpu_stream=-1)
        if(is_ibm_force .and. is_impdiff_1d) then
          call ibm_impulse('impl',n,lo,ng,grid_vol_ratio_f,mask_u,A_u,B_u,0._rp,u,imp_ibm(:,:,1),imp_map(:,:,1))
          call ibm_impulse('impl',n,lo,ng,grid_vol_ratio_f,mask_v,A_v,B_v,0._rp,v,imp_ibm(:,:,2),imp_map(:,:,2))
          call ibm_impulse('impl',n,lo,ng,grid_vol_ratio_c,mask_w,A_w,B_w,0._rp,w,imp_ibm(:,:,3))
        end if
        if(is_ibm_map .and. is_impdiff_1d) then
          call wf_scatter('impl',1,dli,dzci,grid_vol_ratio_f,mask_u,A_u,B_u,0._rp,0._rp,p,u)
          call wf_scatter('impl',2,dli,dzci,grid_vol_ratio_f,mask_v,A_v,B_v,0._rp,0._rp,p,v)
          call wf_scatter('impl',3,dli,dzci,grid_vol_ratio_c,mask_w,A_w,B_w,0._rp,0._rp,p,w)
        end if
        if(is_ibm_force) call timer_toc('ibm_force',nvtx_gpu_stream=-1)
      end if
      !****//IBM\\****!
      if(is_ibm)then
        call timer_tic('ibm_staircase',nvtx_color='r',nvtx_gpu_stream=-1)
        if(is_ibm_force .or. .not.is_impdiff) then ! else done in `solve_helmholtz_ibm`
          call apply_ibm_staircase(u,mask_u,dtrk)
          call apply_ibm_staircase(v,mask_v,dtrk)
          call apply_ibm_staircase(w,mask_w,dtrk)
        end if
        do nscalibm=1,nscal
          if(scalars(nscalibm)%is_cht) cycle ! CHT: the scalar is solved inside the solid
          call apply_ibm_staircase_scalar(scalars(nscalibm)%val,&
              mask_s,swz(:,nscalibm))
        end do 
        call timer_toc('ibm_staircase',nvtx_gpu_stream=-1)
      endif
      !*******************!
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      !*******************!
      call fillps(n,dli,dzfi,dtrki,u,v,w,pp)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],pp,is_ptdma_update_p,ap_d,cp_d)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,pp)
      if(trim(ibm_proj) == 'damped') then
        call correc_ibm(n,dli,dzci,dtrk,pp,mask_u,mask_v,mask_w,A_u,A_v,A_w,u,v,w)
        if(is_ibm_force) call timer_tic('ibm_force',nvtx_color='m',nvtx_gpu_stream=-1)
        if(is_ibm_force) then
          call ibm_impulse_proj(1,n,lo,ng,dli,dzci,grid_vol_ratio_f,dtrk,mask_u,A_u,pp,imp_leak(:,:,1),imp_map(:,:,1))
          call ibm_impulse_proj(2,n,lo,ng,dli,dzci,grid_vol_ratio_f,dtrk,mask_v,A_v,pp,imp_leak(:,:,2),imp_map(:,:,2))
          call ibm_impulse_proj(3,n,lo,ng,dli,dzci,grid_vol_ratio_c,dtrk,mask_w,A_w,pp,imp_leak(:,:,3))
        end if
        if(is_ibm_map) then
          call calc_mean_flow(grid_vol_ratio_f,mask_s,pp,pp_wf,n)
          call wf_gather(1,n,dli,dzfi,grid_vol_ratio_f,mask_u,dtrk,pp,pp_wf,0._rp)
          call wf_gather(2,n,dli,dzfi,grid_vol_ratio_f,mask_v,dtrk,pp,pp_wf,0._rp)
          call wf_gather(3,n,dli,dzfi,grid_vol_ratio_f,mask_w,dtrk,pp,pp_wf,0._rp)
          call wf_scatter('proj',1,dli,dzci,grid_vol_ratio_f,mask_u,A_u,B_u,dtrk,0._rp,pp,u)
          call wf_scatter('proj',2,dli,dzci,grid_vol_ratio_f,mask_v,A_v,B_v,dtrk,0._rp,pp,v)
          call wf_scatter('proj',3,dli,dzci,grid_vol_ratio_c,mask_w,A_w,B_w,dtrk,0._rp,pp,w)
        end if
        if(is_ibm_force) call timer_toc('ibm_force',nvtx_gpu_stream=-1)
      else
        call correc(n,dli,dzci,dtrk,pp,u,v,w)
      end if
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
      dtrk_last = dtrk
      call updatep(n,dli,dzci,dzfi,alpha,pp,p)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
    end do
    dti_step    = dti
    dpdl(:)     = -dpdl(:)*dti
    fs(1:nscal) = fs(1:nscal)*dti
    if(is_ibm_force) then
      call timer_tic('ibm_out',nvtx_color='m',nvtx_gpu_stream=-1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,imp_ibm ,size(imp_ibm ),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,imp_leak,size(imp_leak),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      imp_ibm(:,:,:)  = imp_ibm(:,:,:)*dti
      imp_leak(:,:,:) = imp_leak(:,:,:)*dti
      !$acc update self(imp_map)
      imp_map(:,:,:)  = imp_map(:,:,:)*dti
      call timer_toc('ibm_out',nvtx_gpu_stream=-1)
    end if
    pi_drive(:) = bforce(:) - dpdl(:) ! total driving force per unit volume in this step

    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    if(icheck > 0.and.mod(istep,max(icheck,1)) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,alpha_max_eff,u,v,w,dt_cfl)
      dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
      if(dt_cfl < small) then
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,l,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Velocity divergence norms: Mean = ', divtot, '| Maximum = ', divmax
      if(is_ibm) then
        call chk_ibm(n,l,dli,dzfi,grid_vol_ratio_f,mask_u,mask_v,mask_w,u,v,w, &
                     velmax_solid,divtot_m,divmax_m,ubulk_fl)
        if(myid == 0) print*, 'IBM: max|u,v,w| in solid = ', velmax_solid(:)
        if(myid == 0) print*, 'IBM: masked divergence norms: Mean = ', divtot_m, '| Maximum = ', divmax_m
        if(myid == 0) print*, 'IBM: fluid-only bulk velocity = ', ubulk_fl
        var(1) = 1.*istep
        var(2) = time
        var(3:5) = velmax_solid(:)
        var(6) = divtot_m
        var(7) = divmax_m
        var(8) = ubulk_fl
        call out0d(trim(datadir)//'ibm_chk.out',8,var)
      end if
      if(.not.is_mask_divergence_check) then
        if((divmax > small.and.trim(ibm_proj) /= 'damped').or.is_nan(divtot)) then
          if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
          if(myid == 0) print*, 'Aborting...'
          is_done = .true.
          kill = .true.
        end if
      end if
      !
      ! scalars: abort on a non-finite value (a NaN/Inf anywhere makes the volume mean non-finite)
      !
      do iscal=1,nscal
        call bulk_mean(n,grid_vol_ratio_f,scalars(iscal)%val,meanscal)
        if(is_nan(meanscal).or.abs(meanscal) > huge(meanscal)) then
          if(myid == 0) print*, 'ERROR: scalar ', iscal, ' is not finite.'
          if(myid == 0) print*, 'Aborting...'
          is_done = .true.
          kill = .true.
        end if
      end do
    end if
    !
    ! 3D time averages
    !
    if(iaverage > 0.and.mod(istep,max(iaverage,1)) == 0.and.time >= average_start_time) &
      call stats3d_update(n,nscal,u,v,w,p,scalars)
    !
    ! output routines below
    !
    if(iout0d > 0.and.mod(istep,max(iout0d,1)) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)) > 0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_ibm) then ! fluid-only mean, as controlled by the IBM bulk forcing
          is_intrinsic = trim(bulk_vel_def) == 'intrinsic'
          if(is_forced(1).or.abs(bforce(1)) > 0.) then
            call bulk_mean_ibm(n,grid_vol_ratio_f,mask_u,A_u,u,is_intrinsic,meanvelu,gain)
          end if
          if(is_forced(2).or.abs(bforce(2)) > 0.) then
            call bulk_mean_ibm(n,grid_vol_ratio_f,mask_v,A_v,v,is_intrinsic,meanvelv,gain)
          end if
          if(is_forced(3).or.abs(bforce(3)) > 0.) then
            call bulk_mean_ibm(n,grid_vol_ratio_c,mask_w,A_w,w,is_intrinsic,meanvelw,gain)
          end if
        else
          if(is_forced(1).or.abs(bforce(1)) > 0.) then
            call bulk_mean(n,grid_vol_ratio_f,u,meanvelu)
          end if
          if(is_forced(2).or.abs(bforce(2)) > 0.) then
            call bulk_mean(n,grid_vol_ratio_f,v,meanvelv)
          end if
          if(is_forced(3).or.abs(bforce(3)) > 0.) then
            call bulk_mean(n,grid_vol_ratio_c,w,meanvelw)
          end if
        end if
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = [meanvelu,meanvelv,meanvelw]
        call out0d(trim(datadir)//'forcing.out',7,var)
      end if
      !
      do iscal=1,nscal
        s => scalars(iscal)
        write(scalnum,'(i3.3)') iscal
        if(s%is_forced.or.abs(s%source) > 0.) then
          meanscal = 0.
          if(is_ibm) then
            call bulk_mean_ibm(n,grid_vol_ratio_f,mask_s,A_s,s%val,trim(bulk_vel_def) == 'intrinsic',meanscal,gain)
          else
            call bulk_mean(n,grid_vol_ratio_f,s%val,meanscal)
          end if
          if(.not.s%is_forced) fs(:) = s%source
          var(1:3) = [time,fs(iscal),meanscal]
          call out0d(trim(datadir)//'forcing_s_'//scalnum//'.out',3,var)
        end if
      end do
      if(is_ibm_force) then
        !
        ! force the IB exerts on the fluid in this time step, per unit box volume (sign: drag > 0),
        ! from the momentum the IB removes from the plain explicit update (see `ibm_impulse`);
        ! at steady state, total + wall stresses at the domain boundaries = pi_drive
        !
        var(1)     = time
        var(2:4)   = pi_drive(:)
        var(5:7)   = sum(sum(imp_ibm(:,:,:)+imp_leak(:,:,:),1),1) ! total
        var(8:10)  = sum(imp_ibm(1,:,:)+imp_leak(1,:,:),1)       ! lower half of the domain
        var(11:13) = sum(imp_ibm(2,:,:)+imp_leak(2,:,:),1)       ! upper half of the domain
        var(14:16) = sum(imp_ibm(:,1,:),1)                       ! fluid (band) points, lambda term
        var(17:19) = sum(imp_ibm(:,2,:),1)                       ! solid points, predictor
        var(20:22) = sum(imp_leak(:,2,:),1)                      ! solid points, projection leakage
        if(trim(ibm_proj) == 'none') then
          call out0d(trim(datadir)//'ibm_force.out',22,var)
        else
          var(23:25) = sum(imp_leak(:,1,:),1)                    ! fluid (band) points, damped correction
          call out0d(trim(datadir)//'ibm_force.out',25,var)
        end if
      end if
    end if
    write(fldnum,'(i7.7)') istep
    ! turb_statistics+ friction statistics
    if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      call out1d_chan(trim(datadir)//'turb_stats_'//fldnum//'.out',ng,lo,hi,3,l,dl,zc_g,u,v,w)
      if(is_ibm) &
        call out1d_chan_ibm(trim(datadir)//'turb_stats_ibm_'//fldnum//'.out',ng,lo,hi,l,dl,zc_g,u,v,w, &
                            mask_u,mask_v,mask_w,mask_s)
    endif
    !
    ! per-point wall outputs, every `iout_wall` steps
    !
    if(is_ibm_map) then
      call timer_tic('ibm_out',nvtx_color='m',nvtx_gpu_stream=-1)
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      call calc_shear_st(fri_u,u,visc)
      call calc_shear_st(fri_v,v,visc)
      call calc_shear_st(fri_w,w,visc)
      call get_wall_pres(p,p_grad,lo,dl,zc,zf,mask_s,band_s)
      if(myid==0)then
        globu=>fri_u_g
        globv=>fri_v_g
        globw=>fri_w_g
        globp=>p_grad_g
      else
        globu=>dum_g
        globv=>dum_g
        globw=>dum_g
        globp=>dum_g
      endif
      !$acc update self(fri_u,fri_v,fri_w,p_grad)
      do ibm_i=1,12
        call MPI_GATHERV(fri_u(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_u,MPI_REAL_RP,globu(:,ibm_i),fri_u_ga,displ_u,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(fri_v(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_v,MPI_REAL_RP,globv(:,ibm_i),fri_v_ga,displ_v,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(fri_w(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_w,MPI_REAL_RP,globw(:,ibm_i),fri_w_ga,displ_w,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(p_grad(:,ibm_i)+ibm_idx_off(ibm_i),&
                         size_s,MPI_REAL_RP,globp(:,ibm_i),sca_ga,displ_s,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      end do
      if(myid==0)then
        call write_data(trim(datadir)//"fric-u_"//fldnum//".out",myid,fri_u_g)
        call write_data(trim(datadir)//"fric-v_"//fldnum//".out",myid,fri_v_g)
        call write_data(trim(datadir)//"fric-w_"//fldnum//".out",myid,fri_w_g)
        call write_data(trim(datadir)//"pressure"//fldnum//".out",myid,p_grad_g)
      endif
      !
      ! exact wall force per cell of this time step (FV method; analytic or level-set geometry)
      !
      call calc_wall_force(n,lo,l,dti_step,pi_drive,mask_s,frifv,frifvout)
      call writeFVfric(trim(datadir)//'FV-fric_'//trim(mpiranknum)//'_s_'//fldnum//'.out',frifvout)
      ! wall stress map of the lower wall from the exact IB force of this time step (analytic wall shape)
      if(.not.is_sdf .and. ibm_direction(0,3)) &
        call out_ibm_wall_stress(trim(datadir)//'ibm_wall_stress_'//fldnum//'.out',n,lo,ng,l,dl,pi_drive,imp_map)
      !
      ! wall heat flux per cell, same cells and order as FV-fric (instantaneous, end of the time step);
      ! wall_heat.out: time, total heat flow into the wall per scalar, and that in the output cells
      !
      if(is_wall_heat) then
        do iscal=1,nscal
          s => scalars(iscal)
          call calc_wall_heat(iscal,nscal,n,lo,ng,is_per,dl,dzc,dzf,mask_s,thx,thy,thz,s%val, &
                              s%alpha,s%alphas,s%is_cht,swz(:,iscal))
        end do
        call wall_heat_out(n,lo,nscal,mask_s,frifv,frhout,qwall)
        call write_fv_rows(trim(datadir)//'FV-heat_'//trim(mpiranknum)//'_s_'//fldnum//'.out',frhout)
        var(1) = time
        var(2:1+2*nscal) = qwall(:)
        call out0d(trim(datadir)//'wall_heat.out',1+2*nscal,var)
      end if
      call timer_toc('ibm_out',nvtx_gpu_stream=-1)
    endif
    !
    if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out1d.h90'
    end if
    if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out2d.h90'
    end if
    if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out3d.h90'
    end if
    if(isave > 0.and.((mod(istep,max(isave,1)) == 0).or.(is_done.and..not.kill))) then
      if(is_overwrite_save) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      do is=1,4+nscal
        select case(trim(c_io_vars(is)))
        case('u')
          call load_one('w',trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                        MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xf_g,yc_g,zc_g)
        case('v')
          call load_one('w',trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                        MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yf_g,zc_g)
        case('w')
          call load_one('w',trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                        MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yc_g,zf_g)
        case default
          call load_one('w',trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                        MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yc_g,zc_g)
        end select
      end do
      if(.not.is_overwrite_save) then
        !
        ! fld_*.<ext> -> last checkpoint file (symbolic link)
        !
        do is=1,4+nscal
          call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                         'fld_'//trim(c_io_vars(is))//trim(io_ext))
        end do
      end if
      if(iaverage > 0) call stats3d_save(trim(datadir),ng,lo,hi,nscal,time,istep,xc_g,yc_g,zc_g,myid)
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
    call timer_toc(trim(merge('step_ibm_out','step        ',is_ibm_force)),nvtx_gpu_stream=-1)
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Average, minimum & maximum elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
    end if
  end do
  !***//calc-mean-flow\\***!
  if(do_richardson .and. is_ibm)then ! fluid-only means; the masks exist only with is_ibm
    !call calc_mean_flow_easy(mask_u,u,mean_u,n,dl,1)
    !call calc_mean_flow_easy(mask_v,v,mean_v,n,dl,2)
    !call calc_mean_flow_easy(mask_w,w,mean_w,n,dl,3)
    call calc_mean_flow(grid_vol_ratio_f,mask_u,u,mean_u,n)
    call calc_mean_flow(grid_vol_ratio_f,mask_v,v,mean_v,n)
    call calc_mean_flow(grid_vol_ratio_c,mask_w,w,mean_w,n)
    if(myid==0)then
      print*,"mean flow for each dir: x, y, z",mean_u,mean_v,mean_w
    endif
    do iscal=1,nscal
      call calc_mean_flow(grid_vol_ratio_f,mask_s,scalars(iscal)%val,meanscal,n)
      if(myid==0) print*,"mean scalar ",iscal,": ",meanscal
    end do
  endif
  !************************!
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call fftend(arrplanu)
    call fftend(arrplanv)
    call fftend(arrplanw)
  end if
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize
#endif
  call MPI_FINALIZE(ierr)
  if(kill) error stop 1 ! aborted run (time step too small, divergence, non-finite scalar): non-zero exit status
end program cans
