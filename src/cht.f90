! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_cht
  !
  ! implicit treatment of the diffusion of scalars with conjugate heat transfer (CHT), i.e., with a
  ! thermal conductivity and heat capacity that differ between the fluid and the IBM solid
  !
  ! the diffusion operator L = L_x + L_y + L_z, e.g. along z
  !
  !   L_z T(k) = az(k)*(T(k-1)-T(k)) + cz(k)*(T(k+1)-T(k)),
  !   az(k) = kzm(k)/(dzc(k-1)*dzf(k)*C(k)), cz(k) = kzp(k)/(dzc(k)*dzf(k)*C(k)),
  !
  ! with kzp/kzm the conductivity of the upper/lower cell face, kface = 1/(th/kf + (1-th)/ks), th the
  ! fluid fraction of the centre-to-centre segment (see `scal_cht`), and C = 1 + (cs_cf-1)*(1-phi) the
  ! (relative) heat capacity, phi the fluid volume fraction; it is integrated in time with the delta
  ! form of the Crank-Nicolson scheme with approximate factorization (Douglas ADI), per Runge-Kutta
  ! substep:
  !
  !   (1 - g*L_x)(1 - g*L_y)(1 - g*L_z) dT = rhs,   T^{n+1} = T^n + dT,   g = dt_rk/2,
  !
  ! where rhs contains the explicit terms and 2*g*L*T^n; the scheme is second-order accurate in time,
  ! and its steady state is that of the unfactored discretisation; it is unconditionally stable only
  ! when L_x, L_y and L_z commute: with the strongly varying k/C of an interface the three-factor
  ! product is not, and large time steps can amplify short waves (observed on a laminar wavy-wall
  ! case: stable with solid/fluid diffusivity ratios (k_s/k_f)/(C_s/C_f) <= 100, unstable at 1000
  ! for dt*(k/C)_s*4*(1/dx**2+1/dy**2+1/dz**2) >~ 1000); a non-finite scalar aborts the run (main);
  ! the tridiagonal systems along x and y (periodic) and z (Dirichlet/Neumann) are solved per line,
  ! after transposing dT to the pencil aligned with the line direction; the coefficients are
  ! evaluated on the fly from the face fractions and fluid volume fractions, which are stored once
  ! per pencil orientation
  !
  ! on GPUs, the x and y pencils have the line direction along the first (contiguous) array axis; a
  ! solve with one thread per line along that axis has uncoalesced memory accesses, so these two
  ! pencils are permuted locally to (second axis, first axis, third axis) around the line solve
  ! (`solve_swapped`), with the line direction along the second axis and the threads along the
  ! contiguous one; the face and volume fractions of these pencils are stored in the permuted layout
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_param     , only: ipencil_axis,is_impdiff
  use mod_types
#if !defined(_OPENACC)
  use decomp_2d
#else
  use cudecomp
  use mod_param          , only: cudecomp_is_t_in_place
  use mod_common_cudecomp, only: dtype_rp => cudecomp_real_rp, &
                                 solver_buf_0,solver_buf_1,work, &
                                 ap_x => ap_x_poi,ap_y => ap_y_poi,ap_z => ap_z_poi, &
                                 ch => handle,gd => gd_poi, &
                                 istream => istream_acc_queue_1
#endif
  implicit none
  private
  public cht_init,cht_solve
  !
  ! face fractions along the line direction and fluid volume fractions in the x-, y- and z-aligned
  ! pencils (memory layout of the respective transposes)
  !
  real(rp), allocatable, dimension(:,:,:) :: thx_x,thy_y,thz_z,phi_x,phi_y,phi_z
  real(rp), allocatable, dimension(:,:,:) :: cwork,uwork ! work arrays of the line solvers
  real(rp), allocatable, dimension(:,:,:) :: pwork       ! permuted x/y pencil of the GPU line solves
  real(rp), allocatable, dimension(:) :: dzci_c,dzfi_c   ! global grid spacings along z
  integer , dimension(3) :: ngl                          ! global grid size
  integer , dimension(3) :: nx_p,ny_p,nz_p               ! shapes of the x-, y- and z-pencil arrays
  integer :: iy_axis                                     ! axis along y in the y-pencil array (1 on GPUs, 2 on CPUs)
  real(rp) :: dxi2,dyi2
  contains
  subroutine cht_init(n,ng,l,dzci_g,dzfi_g,thx,thy,thz,phi,scalars)
    use mod_scal, only: scalar
    implicit none
    integer , intent(in), dimension(3) :: n,ng
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(0:) :: dzci_g,dzfi_g
    real(rp), intent(in), dimension(0:,0:,0:) :: thx,thy,thz,phi
    type(scalar), intent(in), dimension(:) :: scalars
    integer :: is,nmax
    !
    if(.not.any(scalars(:)%is_cht)) return
    do is=1,size(scalars)
      if(.not.scalars(is)%is_cht) cycle
      if(any(scalars(is)%cbc(0:1,3) == 'P').or.any(scalars(is)%cbc(0:1,1:2) /= 'P')) then
        call abort_cht('CHT scalars require periodic boundary conditions along x and y, and non-periodic along z.')
      end if
    end do
    if(.not.is_impdiff) return
    if(any(ng(1:2) < 3)) call abort_cht('CHT scalars with implicit diffusion require at least 3 grid points along x and y.')
    ngl(:) = ng(:)
    dxi2 = (ng(1)/l(1))**2
    dyi2 = (ng(2)/l(2))**2
    allocate(dzci_c(0:ng(3)+1),dzfi_c(0:ng(3)+1))
    dzci_c(:) = dzci_g(0:ng(3)+1)
    dzfi_c(:) = dzfi_g(0:ng(3)+1)
    !$acc enter data copyin(dzci_c,dzfi_c)
#if !defined(_OPENACC)
    nx_p(:) = xsize(:); ny_p(:) = ysize(:); nz_p(:) = zsize(:)
    iy_axis = 2
#else
    if(ipencil_axis /= 3) call abort_cht('on GPUs, CHT scalars with implicit diffusion require `ipencil_axis = 3`.')
    nx_p(:) = ap_x%shape(:); ny_p(:) = ap_y%shape(:); nz_p(:) = ap_z%shape(:)
    iy_axis = 1 ! axis-contiguous layout of the y pencil: (y,z,x)
#endif
    nmax = maxval([product(nx_p),product(ny_p),product(nz_p)])
#if !defined(_OPENACC)
    allocate(thx_x(nx_p(1),nx_p(2),nx_p(3)),phi_x(nx_p(1),nx_p(2),nx_p(3)), &
             thy_y(ny_p(1),ny_p(2),ny_p(3)),phi_y(ny_p(1),ny_p(2),ny_p(3)))
    allocate(pwork(1,1,1))
#else
    allocate(thx_x(nx_p(2),nx_p(1),nx_p(3)),phi_x(nx_p(2),nx_p(1),nx_p(3)), & ! permuted layouts
             thy_y(ny_p(2),ny_p(1),ny_p(3)),phi_y(ny_p(2),ny_p(1),ny_p(3)))
    allocate(pwork(nmax,1,1))
#endif
    allocate(thz_z(nz_p(1),nz_p(2),nz_p(3)),phi_z(nz_p(1),nz_p(2),nz_p(3)))
    allocate(cwork(nmax,1,1),uwork(nmax,1,1))
    !$acc enter data create(thx_x,thy_y,thz_z,phi_x,phi_y,phi_z,cwork,uwork,pwork)
    call scatter(n,thx,'thx')
    call scatter(n,thy,'thy')
    call scatter(n,thz,'thz')
    call scatter(n,phi,'phi')
  end subroutine cht_init
  !
  subroutine scatter(n,q,task)
    !
    ! stores a cell-centred field in the pencil(s) selected by `task`
    !
    implicit none
    integer , intent(in) :: n(3)
    real(rp), intent(in), dimension(0:,0:,0:) :: q
    character(len=*), intent(in) :: task
    real(rp), allocatable, dimension(:,:,:) :: p
    allocate(p(n(1),n(2),n(3)))
    !$acc enter data create(p)
#if defined(_OPENACC)
    !$acc kernels default(present) async(1)
    solver_buf_0(:) = 0. ! the padded parts of the pencils (x and y) are then zero
    solver_buf_1(:) = 0.
    !$acc end kernels
#endif
    !$acc kernels default(present) async(1)
    p(:,:,:) = q(1:n(1),1:n(2),1:n(3))
    !$acc end kernels
    call transposes(n,p,task)
    !$acc wait(1)
    !$acc exit data delete(p)
  end subroutine scatter
  !
  subroutine cht_solve(n,f,alphaf,alphas,cs_cf,cbcz,p)
    !
    ! solves (1 - f*L_x)(1 - f*L_y)(1 - f*L_z) dT = p for one CHT scalar; p holds the right-hand
    ! side on input and dT on output; dT satisfies homogeneous boundary conditions along z
    !
    implicit none
    integer , intent(in) :: n(3)
    real(rp), intent(in) :: f,alphaf,alphas,cs_cf
    character(len=1), intent(in), dimension(0:1) :: cbcz
    real(rp), intent(inout), dimension(:,:,:) :: p
    call transposes(n,p,'solve',f,alphaf,alphas,cs_cf,cbcz(0) == 'D',cbcz(1) == 'D')
  end subroutine cht_solve
  !
  subroutine transposes(n,p,task,f,alphaf,alphas,cs_cf,is_d0,is_d1)
    !
    ! visits the pencils aligned with x, y and z (in an order that depends on the pencil of `p`);
    ! task = 'thx'/'thy'/'thz'/'phi': stores the field (in p) in the pencil(s) it belongs to
    ! task = 'solve'                : solves the tridiagonal systems along each direction, result back in p
    !
    implicit none
    integer , intent(in) :: n(3)
    real(rp), intent(inout), dimension(:,:,:) :: p
    character(len=*), intent(in) :: task
    real(rp), intent(in), optional :: f,alphaf,alphas,cs_cf
    logical , intent(in), optional :: is_d0,is_d1
#if !defined(_OPENACC)
    real(rp), allocatable, dimension(:,:,:) :: px,py,pz
    allocate(px(nx_p(1),nx_p(2),nx_p(3)),py(ny_p(1),ny_p(2),ny_p(3)),pz(nz_p(1),nz_p(2),nz_p(3)))
    select case(ipencil_axis)
    case(1)
      !$OMP PARALLEL WORKSHARE
      px(:,:,:) = p(:,:,:)
      !$OMP END PARALLEL WORKSHARE
      call visit('x',px)
      call transpose_x_to_y(px,py); call visit('y',py)
      call transpose_y_to_z(py,pz); call visit('z',pz)
      if(task == 'solve') then
        call transpose_z_to_y(pz,py); call transpose_y_to_x(py,px)
        !$OMP PARALLEL WORKSHARE
        p(:,:,:) = px(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      end if
    case(2)
      !$OMP PARALLEL WORKSHARE
      py(:,:,:) = p(:,:,:)
      !$OMP END PARALLEL WORKSHARE
      call visit('y',py)
      call transpose_y_to_z(py,pz); call visit('z',pz)
      call transpose_z_to_y(pz,py)
      call transpose_y_to_x(py,px); call visit('x',px)
      if(task == 'solve') then
        call transpose_x_to_y(px,py)
        !$OMP PARALLEL WORKSHARE
        p(:,:,:) = py(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      end if
    case(3)
      !$OMP PARALLEL WORKSHARE
      pz(:,:,:) = p(:,:,:)
      !$OMP END PARALLEL WORKSHARE
      call visit('z',pz)
      call transpose_z_to_y(pz,py); call visit('y',py)
      call transpose_y_to_x(py,px); call visit('x',px)
      if(task == 'solve') then
        call transpose_x_to_y(px,py); call transpose_y_to_z(py,pz)
        !$OMP PARALLEL WORKSHARE
        p(:,:,:) = pz(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      end if
    end select
#else
    !
    ! GPU: z pencil (ipencil_axis = 3) -> y -> x with the cuDecomp transposes of the Poisson grid
    ! descriptor (x pencil padded along x, y pencil axis-contiguous); the buffers are those of the
    ! Poisson solver
    !
    real(rp), pointer, contiguous, dimension(:,:,:) :: px,py,pz
    integer :: istat
    pz(1:nz_p(1),1:nz_p(2),1:nz_p(3)) => solver_buf_0(1:product(nz_p(:)))
    px(1:nx_p(1),1:nx_p(2),1:nx_p(3)) => solver_buf_0(1:product(nx_p(:)))
    if(cudecomp_is_t_in_place) then
      py(1:ny_p(1),1:ny_p(2),1:ny_p(3)) => solver_buf_0(1:product(ny_p(:)))
    else
      py(1:ny_p(1),1:ny_p(2),1:ny_p(3)) => solver_buf_1(1:product(ny_p(:)))
    end if
    !
    ! the pencils of the Poisson grid descriptor are padded along x and y: clear the buffers first,
    ! so that the line solves in the padded parts never operate on stale (possibly non-finite) data
    !
    !$acc kernels default(present) async(1)
    solver_buf_0(:) = 0.
    solver_buf_1(:) = 0.
    !$acc end kernels
    !$acc kernels default(present) async(1)
    pz(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3))
    !$acc end kernels
    call visit('z',pz)
    !$acc host_data use_device(pz,py,work)
    istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
    !$acc end host_data
    call visit('y',py)
    !$acc host_data use_device(py,px,work)
    istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
    !$acc end host_data
    call visit('x',px)
    if(task == 'solve') then
      !$acc host_data use_device(px,py,pz,work)
      istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
      !$acc end host_data
      !$acc kernels default(present) async(1)
      p(1:n(1),1:n(2),1:n(3)) = pz(1:n(1),1:n(2),1:n(3))
      !$acc end kernels
    end if
#endif
    contains
    subroutine visit(dir,q)
      character(len=1), intent(in) :: dir
      real(rp), intent(inout), dimension(:,:,:) :: q
      select case(task)
      case('thx')
        if(dir == 'x') then
#if !defined(_OPENACC)
          thx_x(:,:,:) = q(:,:,:)
#else
          call swap12(nx_p(1),nx_p(2),nx_p(3),q,thx_x)
#endif
        end if
      case('thy')
        if(dir == 'y') then
#if !defined(_OPENACC)
          thy_y(:,:,:) = q(:,:,:)
#else
          call swap12(ny_p(1),ny_p(2),ny_p(3),q,thy_y)
#endif
        end if
      case('thz')
        if(dir == 'z') then
          !$acc kernels default(present) async(1)
          thz_z(:,:,:) = q(:,:,:)
          !$acc end kernels
        end if
      case('phi')
        select case(dir)
        case('x')
#if !defined(_OPENACC)
          phi_x(:,:,:) = q(:,:,:)
#else
          call swap12(nx_p(1),nx_p(2),nx_p(3),q,phi_x)
#endif
        case('y')
#if !defined(_OPENACC)
          phi_y(:,:,:) = q(:,:,:)
#else
          call swap12(ny_p(1),ny_p(2),ny_p(3),q,phi_y)
#endif
        case('z')
          !$acc kernels default(present) async(1)
          phi_z(:,:,:) = q(:,:,:)
          !$acc end kernels
        end select
      case('solve')
        select case(dir)
        case('x')
#if !defined(_OPENACC)
          call solve_periodic_axis1(nx_p(1),nx_p(2),nx_p(3),ngl(1),f*dxi2,alphaf,alphas,cs_cf,thx_x,phi_x,cwork,uwork,q)
#else
          call solve_swapped(nx_p,ngl(1),f*dxi2,alphaf,alphas,cs_cf,thx_x,phi_x,q)
#endif
        case('y')
#if !defined(_OPENACC)
          if(iy_axis == 1) then
            call solve_periodic_axis1(ny_p(1),ny_p(2),ny_p(3),ngl(2),f*dyi2,alphaf,alphas,cs_cf,thy_y,phi_y,cwork,uwork,q)
          else
            call solve_periodic_axis2(ny_p(1),ny_p(2),ny_p(3),ngl(2),f*dyi2,alphaf,alphas,cs_cf,thy_y,phi_y,cwork,uwork,q)
          end if
#else
          call solve_swapped(ny_p,ngl(2),f*dyi2,alphaf,alphas,cs_cf,thy_y,phi_y,q)
#endif
        case('z')
          call solve_z(nz_p(1),nz_p(2),nz_p(3),f,alphaf,alphas,cs_cf,is_d0,is_d1,dzci_c,dzfi_c,thz_z,phi_z,cwork,q)
        end select
      end select
    end subroutine visit
  end subroutine transposes
  !
  pure function kface(rf,rs,th) result(k)
    !
    ! face conductivity from the fluid fraction `th` of the centre-to-centre segment
    ! (rf, rs: inverse fluid/solid conductivity)
    !
    !$acc routine seq
    real(rp), intent(in) :: rf,rs,th
    real(rp) :: k
    k = 1./(th*rf + (1.-th)*rs)
  end function kface
  !
  pure function capinv(cs_cf,phi) result(capi)
    !
    ! inverse relative heat capacity of a cell with fluid volume fraction `phi`
    !
    !$acc routine seq
    real(rp), intent(in) :: cs_cf,phi
    real(rp) :: capi
    capi = 1./(1. + (cs_cf-1.)*(1.-phi))
  end function capinv
  !
  subroutine swap12(n1,n2,n3,a,b)
    !
    ! local permutation of the first two array axes, b(j,i,k) = a(i,j,k) (tiled, so that both the
    ! reads and the writes of a tile stay in cache)
    !
    implicit none
    integer , intent(in) :: n1,n2,n3
    real(rp), intent(in ), dimension(n1,n2,n3) :: a
    real(rp), intent(out), dimension(n2,n1,n3) :: b
    integer :: i,j,k
    !$acc parallel loop tile(32,32,1) default(present) async(1)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
    do k=1,n3
      do i=1,n1
        do j=1,n2
          b(j,i,k) = a(i,j,k)
        end do
      end do
    end do
  end subroutine swap12
  !
  subroutine solve_swapped(np,nl,fd,alphaf,alphas,cs_cf,th,phi,q)
    !
    ! cyclic tridiagonal systems along the first axis of the pencil array q (shape np), solved in the
    ! permuted layout (np(2),np(1),np(3)) with one thread per line along the contiguous axis
    ! (coalesced accesses); th and phi are stored in the permuted layout
    !
    implicit none
    integer , intent(in) :: np(3),nl
    real(rp), intent(in) :: fd,alphaf,alphas,cs_cf
    real(rp), intent(in   ), dimension(np(2),np(1),np(3)) :: th,phi
    real(rp), intent(inout), dimension(np(1),np(2),np(3)) :: q
    call swap12(np(1),np(2),np(3),q,pwork)
    call solve_periodic_axis2(np(2),np(1),np(3),nl,fd,alphaf,alphas,cs_cf,th,phi,cwork,uwork,pwork)
    call swap12(np(2),np(1),np(3),pwork,q)
  end subroutine solve_swapped
  !
  subroutine solve_periodic_axis1(n1,n2,n3,nl,fd,alphaf,alphas,cs_cf,th,phi,cc,uu,p)
    !
    ! cyclic tridiagonal systems along the first array axis (uniform spacing d, fd = f/d**2),
    ! nl valid points per line (the array may be padded), Sherman-Morrison formula;
    ! th(i) is the face between points i and i+1 (th(nl): between nl and 1)
    !
    implicit none
    integer , intent(in) :: n1,n2,n3,nl
    real(rp), intent(in) :: fd,alphaf,alphas,cs_cf
    real(rp), intent(in   ), dimension(n1,n2,n3) :: th,phi
    real(rp), intent(inout), dimension(n1,n2,n3) :: cc,uu,p
    real(rp) :: a,b,c,m,a1,cn,gam,fact,rf,rs,kp,km,capi
    integer :: i,j,k,im
    rf = 1./alphaf
    rs = 1./alphas
    !$acc parallel loop collapse(2) default(present) private(a,b,c,m,a1,cn,gam,fact,kp,km,capi,im) async(1)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(i,a,b,c,m,a1,cn,gam,fact,kp,km,capi,im)
    do k=1,n3
      do j=1,n2
        !
        ! corner coefficients: a1 couples point 1 to nl, cn couples point nl to 1 (same face)
        !
        km = kface(rf,rs,th(nl,j,k))
        a1 = -fd*km*capinv(cs_cf,phi(1 ,j,k))
        cn = -fd*km*capinv(cs_cf,phi(nl,j,k))
        gam = 0.
        !
        ! forward sweep for the modified system A'x = p and A'z = u, u = (gam,0,...,0,cn)
        !
        !$acc loop seq
        do i=1,nl
          im = i-1; if(i == 1) im = nl
          kp = kface(rf,rs,th(i ,j,k))
          km = kface(rf,rs,th(im,j,k))
          capi = capinv(cs_cf,phi(i,j,k))
          a = -fd*capi*km
          c = -fd*capi*kp
          b = 1. - a - c
          if(i == 1) then
            gam = -b
            b = b - gam
            uu(i,j,k) = gam
            m = 1./b
          else
            if(i == nl) then
              b = b - a1*cn/gam
              uu(i,j,k) = cn
              c = 0.
            else
              uu(i,j,k) = 0.
            end if
            m = 1./(b - a*cc(i-1,j,k))
            p(i,j,k)  = p(i,j,k)  - a*p(i-1,j,k)
            uu(i,j,k) = uu(i,j,k) - a*uu(i-1,j,k)
          end if
          cc(i,j,k) = c*m
          p(i,j,k)  = p(i,j,k)*m
          uu(i,j,k) = uu(i,j,k)*m
        end do
        !$acc loop seq
        do i=nl-1,1,-1
          p(i,j,k)  = p(i,j,k)  - cc(i,j,k)*p(i+1,j,k)
          uu(i,j,k) = uu(i,j,k) - cc(i,j,k)*uu(i+1,j,k)
        end do
        fact = (p(1,j,k) + a1*p(nl,j,k)/gam)/(1. + uu(1,j,k) + a1*uu(nl,j,k)/gam)
        !$acc loop seq
        do i=1,nl
          p(i,j,k) = p(i,j,k) - fact*uu(i,j,k)
        end do
      end do
    end do
  end subroutine solve_periodic_axis1
  !
  subroutine solve_periodic_axis2(n1,n2,n3,nl,fd,alphaf,alphas,cs_cf,th,phi,cc,uu,p)
    !
    ! same as `solve_periodic_axis1`, along the second array axis
    !
    implicit none
    integer , intent(in) :: n1,n2,n3,nl
    real(rp), intent(in) :: fd,alphaf,alphas,cs_cf
    real(rp), intent(in   ), dimension(n1,n2,n3) :: th,phi
    real(rp), intent(inout), dimension(n1,n2,n3) :: cc,uu,p
    real(rp) :: a,b,c,m,a1,cn,gam,fact,rf,rs,kp,km,capi
    integer :: i,j,k,jm
    rf = 1./alphaf
    rs = 1./alphas
    !$acc parallel loop collapse(2) default(present) private(a,b,c,m,a1,cn,gam,fact,kp,km,capi,jm) async(1)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(j,a,b,c,m,a1,cn,gam,fact,kp,km,capi,jm)
    do k=1,n3
      do i=1,n1
        km = kface(rf,rs,th(i,nl,k))
        a1 = -fd*km*capinv(cs_cf,phi(i,1 ,k))
        cn = -fd*km*capinv(cs_cf,phi(i,nl,k))
        gam = 0.
        !$acc loop seq
        do j=1,nl
          jm = j-1; if(j == 1) jm = nl
          kp = kface(rf,rs,th(i,j ,k))
          km = kface(rf,rs,th(i,jm,k))
          capi = capinv(cs_cf,phi(i,j,k))
          a = -fd*capi*km
          c = -fd*capi*kp
          b = 1. - a - c
          if(j == 1) then
            gam = -b
            b = b - gam
            uu(i,j,k) = gam
            m = 1./b
          else
            if(j == nl) then
              b = b - a1*cn/gam
              uu(i,j,k) = cn
              c = 0.
            else
              uu(i,j,k) = 0.
            end if
            m = 1./(b - a*cc(i,j-1,k))
            p(i,j,k)  = p(i,j,k)  - a*p(i,j-1,k)
            uu(i,j,k) = uu(i,j,k) - a*uu(i,j-1,k)
          end if
          cc(i,j,k) = c*m
          p(i,j,k)  = p(i,j,k)*m
          uu(i,j,k) = uu(i,j,k)*m
        end do
        !$acc loop seq
        do j=nl-1,1,-1
          p(i,j,k)  = p(i,j,k)  - cc(i,j,k)*p(i,j+1,k)
          uu(i,j,k) = uu(i,j,k) - cc(i,j,k)*uu(i,j+1,k)
        end do
        fact = (p(i,1,k) + a1*p(i,nl,k)/gam)/(1. + uu(i,1,k) + a1*uu(i,nl,k)/gam)
        !$acc loop seq
        do j=1,nl
          p(i,j,k) = p(i,j,k) - fact*uu(i,j,k)
        end do
      end do
    end do
  end subroutine solve_periodic_axis2
  !
  subroutine solve_z(n1,n2,nz,f,alphaf,alphas,cs_cf,is_d0,is_d1,dzci,dzfi,th,phi,cc,p)
    !
    ! tridiagonal systems along z (third array axis, full extent, non-uniform grid), with
    ! homogeneous Dirichlet (dT(0) = -dT(1)) or Neumann (dT(0) = dT(1)) conditions at both ends;
    ! th(k) is the face between cells k and k+1; the outer boundary faces (below cell 1, above
    ! cell nz) lie in the material of the adjacent cell
    !
    implicit none
    integer , intent(in) :: n1,n2,nz
    real(rp), intent(in) :: f,alphaf,alphas,cs_cf
    logical , intent(in) :: is_d0,is_d1
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(n1,n2,nz) :: th,phi
    real(rp), intent(inout), dimension(n1,n2,nz) :: cc,p
    real(rp) :: a,b,c,m,rf,rs,kp,km,capi
    integer :: i,j,k
    rf = 1./alphaf
    rs = 1./alphas
    !$acc parallel loop collapse(2) default(present) private(a,b,c,m,kp,km,capi) async(1)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(k,a,b,c,m,kp,km,capi)
    do j=1,n2
      do i=1,n1
        !$acc loop seq
        do k=1,nz
          if(k == 1) then
            km = kface(rf,rs,merge(1._rp,0._rp,phi(i,j,1) >= 0.5_rp))
          else
            km = kface(rf,rs,th(i,j,k-1))
          end if
          if(k == nz) then
            kp = kface(rf,rs,merge(1._rp,0._rp,phi(i,j,nz) >= 0.5_rp))
          else
            kp = kface(rf,rs,th(i,j,k))
          end if
          capi = capinv(cs_cf,phi(i,j,k))
          a = -f*capi*km*dzci(k-1)*dzfi(k)
          c = -f*capi*kp*dzci(k  )*dzfi(k)
          b = 1. - a - c
          if(k == 1) then
            if(is_d0) b = b - a
            a = 0.
          end if
          if(k == nz) then
            if(is_d1) b = b - c
            c = 0.
          end if
          if(k == 1) then
            m = 1./b
            p(i,j,k) = p(i,j,k)*m
          else
            m = 1./(b - a*cc(i,j,k-1))
            p(i,j,k) = (p(i,j,k) - a*p(i,j,k-1))*m
          end if
          cc(i,j,k) = c*m
        end do
        !$acc loop seq
        do k=nz-1,1,-1
          p(i,j,k) = p(i,j,k) - cc(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do
  end subroutine solve_z
  !
  subroutine abort_cht(msg)
    character(len=*), intent(in) :: msg
    if(myid == 0) print*, 'ERROR: '//msg
    if(myid == 0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abort_cht
end module mod_cht
