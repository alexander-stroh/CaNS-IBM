! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_solve_helmholtz
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_types
  use mod_bound, only: updt_rhs_b
  use mod_param, only: is_impdiff_1d,eps,ibm_pen
#if !defined(_OPENACC)
  use mod_solver    , only: solver
  use mod_solver    , only: solver_gaussel_z
#else
  use mod_solver_gpu, only: solver => solver_gpu
  use mod_solver_gpu, only: solver_gaussel_z => solver_gaussel_z_gpu
#endif
  implicit none
  private
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  public solve_helmholtz,solve_helmholtz_ibm,rhs_bound
  contains
  subroutine solve_helmholtz(n,ng,hi,arrplan,normfft,alpha,lambdaxy,a,b,c,rhsbx,rhsby,rhsbz,is_bound,cbc,c_or_f,p)
    !
    ! this is a wrapper subroutine to solve 1D/3D helmholtz problems: p/alpha + lap(p) = rhs
    !
    integer ,    intent(in   ), dimension(3)                :: n,ng,hi
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), intent(in   ), dimension(2,2),    optional :: arrplan
#else
    integer    , intent(in   ), dimension(2,2),    optional :: arrplan
#endif
    real(rp),    intent(in   ),                    optional :: normfft
    real(rp),    intent(in   )                              :: alpha
    real(rp),    intent(in   ), dimension(:,:),    optional :: lambdaxy
    real(rp),    intent(in   ), dimension(:)                :: a,b,c
    real(rp),    intent(in   ), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    logical ,    intent(in   ), dimension(0:1,3)            :: is_bound
    character(len=1), intent(in), dimension(0:1,3)          :: cbc
    character(len=1), intent(in), dimension(3)              :: c_or_f
    real(rp),    intent(inout), dimension(:,:,:)            :: p
    real(rp), allocatable, dimension(:), save :: bb
    real(rp) :: alphai
    integer :: k
    !
    logical, save :: is_first = .true.
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(bb,mold=b)
      !$acc enter data create(bb) async(1)
    end if
    !
    call updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p,alpha)
    !
    alphai = alpha**(-1)
    !$acc parallel loop default(present) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)
    do k=1,size(b)
      bb(k) = b(k) + alphai
    end do
    !
    if(.not.is_impdiff_1d) then
      call solver(n,ng,arrplan,normfft*alphai,lambdaxy,a,bb,c,cbc,c_or_f,p)
    else
      call solver_gaussel_z(n,ng,hi,a,bb,c,cbc(:,3),c_or_f,alphai,p)
    end if
  end subroutine solve_helmholtz
  !
  subroutine solve_helmholtz_ibm(n,alpha,a,b,c,rhsbx,rhsby,rhsbz,is_bound,cbc,c_or_f,mask,am,sw,is_set_solid,p)
    !
    ! z-implicit (1D) Helmholtz solve with the immersed boundary inside the operator;
    ! each z column is solved directly (Thomas algorithm) with a point-dependent diagonal:
    !   - non-solid points: (am - alpha*lap_z(.))p_new = am*p, with am = A of the exact
    !     integration of the IB term (`calc_a_b`, A = 1 away from the wall); the caller passes
    !     p = rhs/A, so the IB term keeps its exponential integration and the z diffusion is
    !     Crank-Nicolson
    !   - solid points: am is replaced by the penalty ibm_pen, which pins p to the wall value sw
    !     to O(1/ibm_pen), so the band points see the staircase wall, as in the explicit operator;
    !     ibm_pen*(p - sw) is the wall impulse of the solid point (see `ibm_impulse`)
    ! the right-hand side, am*p (ibm_pen*sw at the solid points) plus the boundary terms of
    ! `updt_rhs_b`, is formed in the forward sweep; with is_set_solid the solid points are set
    ! to sw in the back substitution, as `apply_ibm_staircase` does
    ! requires z columns local to each task (n(3) = ng(3)) and a non-periodic z direction
    !
    integer ,    intent(in   ), dimension(3)                :: n
    real(rp),    intent(in   )                              :: alpha
    real(rp),    intent(in   ), dimension(:)                :: a,b,c
    real(rp),    intent(in   ), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    logical ,    intent(in   ), dimension(0:1,3)            :: is_bound
    character(len=1), intent(in), dimension(0:1,3)          :: cbc
    character(len=1), intent(in), dimension(3)              :: c_or_f
    logical ,    intent(in   ), dimension(0:,0:,0:)         :: mask
    real(rp),    intent(in   ), dimension(0:,0:,0:)         :: am
    real(rp),    intent(in   ), dimension(0:)               :: sw ! wall value per level k
    logical ,    intent(in   )                              :: is_set_solid
    real(rp),    intent(inout), dimension(0:,0:,0:)         :: p
    real(rp), allocatable, dimension(:,:,:), save :: d
    real(rp) :: alphai,dk,z,rhs,dm,pm,pk
    integer  :: i,j,k,nn,nx,ny
    logical  :: is_bx0,is_bx1,is_by0,is_by1,is_bz0,is_bz1
    !
    if(.not.allocated(d)) then
      allocate(d(n(1),n(2),n(3)))
      !$acc enter data create(d) async(1)
    end if
    alphai = alpha**(-1)
    !
    ! last point of each direction (where the upper boundary term goes; the last face is the
    ! boundary for a Dirichlet face-centred variable), and the boundary terms to add
    !
    nx = n(1); if(c_or_f(1) == 'f' .and. cbc(1,1) == 'D') nx = n(1)-1
    ny = n(2); if(c_or_f(2) == 'f' .and. cbc(1,2) == 'D') ny = n(2)-1
    nn = n(3); if(c_or_f(3) == 'f' .and. cbc(1,3) == 'D') nn = n(3)-1
    is_bx0 = present(rhsbx) .and. is_bound(0,1); is_bx1 = present(rhsbx) .and. is_bound(1,1)
    is_by0 = present(rhsby) .and. is_bound(0,2); is_by1 = present(rhsby) .and. is_bound(1,2)
    is_bz0 = present(rhsbz) .and. is_bound(0,3); is_bz1 = present(rhsbz) .and. is_bound(1,3)
    !
    !$acc parallel loop gang vector collapse(2) default(present) private(dk,z,rhs,dm,pm,pk) async(1)
    !$OMP PARALLEL DO   COLLAPSE(2) DEFAULT(shared)  PRIVATE(dk,z,rhs,dm,pm,pk,k)
    do j=1,n(2)
      do i=1,n(1)
        dm = 0.
        pm = 0.
        !$acc loop seq
        do k=1,n(3)
          if(mask(i,j,k)) then
            rhs = ibm_pen*sw(k)
          else
            rhs = am(i,j,k)*p(i,j,k)
          end if
          if(is_bx0 .and. i == 1 ) rhs = rhs + rhsbx(j,k,0)*alpha
          if(is_bx1 .and. i == nx) rhs = rhs + rhsbx(j,k,1)*alpha
          if(is_by0 .and. j == 1 ) rhs = rhs + rhsby(i,k,0)*alpha
          if(is_by1 .and. j == ny) rhs = rhs + rhsby(i,k,1)*alpha
          if(is_bz0 .and. k == 1 ) rhs = rhs + rhsbz(i,j,0)*alpha
          if(is_bz1 .and. k == nn) rhs = rhs + rhsbz(i,j,1)*alpha
          if(k <= nn) then
            dk = merge(ibm_pen,am(i,j,k),mask(i,j,k))
            z = 1./(b(k)+alphai*dk-a(k)*dm+eps)
            dm = c(k)*z
            pm = (rhs*alphai-a(k)*pm)*z
            d(i,j,k) = dm
            p(i,j,k) = pm
          else
            p(i,j,k) = rhs ! boundary face, not solved
            if(is_set_solid .and. mask(i,j,k)) p(i,j,k) = sw(k)
          end if
        end do
        pk = p(i,j,nn)
        if(is_set_solid .and. mask(i,j,nn)) p(i,j,nn) = sw(nn)
        !$acc loop seq
        do k=nn-1,1,-1
          pk = p(i,j,k) - d(i,j,k)*pk
          p(i,j,k) = pk
          if(is_set_solid .and. mask(i,j,k)) p(i,j,k) = sw(k)
        end do
      end do
    end do
  end subroutine solve_helmholtz_ibm
end module mod_solve_helmholtz
