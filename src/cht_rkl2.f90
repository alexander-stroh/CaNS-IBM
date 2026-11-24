! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_cht_rkl2
  !
  ! RKL2 super time-stepping of the conjugate conduction of a CHT scalar (`is_cht_rkl2`; Meyer, Balsara & Aslam,
  ! J. Comput. Phys. 257 (2014) 594): an explicit s-stage method for dT/dt = L(T) whose stability polynomial
  ! is a shifted Legendre polynomial, second order, stable for dt*lambda_max <= (s**2+s-2)/2 on the negative
  ! real axis. L is the whole conduction operator over the cell heat capacity (series face flux, the tangential
  ! face flux with `is_cht_tang`, and the constant sources). Each stage is one application of L and one halo
  ! update, a conservative flux update; no linear system. The main loop splits it from the Runge-Kutta step of
  ! the convection (Strang): half a step of RKL2, the three Runge-Kutta stages, half a step of RKL2.
  !
  !   Y0 = T, M0 = L(Y0), Y1 = Y0 + mu1~ h M0, mu1~ = w1/3, w1 = 4/(s**2+s-2)
  !   Yj = mu_j Y(j-1) + nu_j Y(j-2) + (1-mu_j-nu_j) Y0 + mu_j~ h L(Y(j-1)) + gamma_j~ h M0, j = 2..s
  !   b_0 = b_1 = b_2 = 1/3, b_j = (j**2+j-2)/(2j(j+1)); mu_j = (2j-1)/j b_j/b_(j-1), nu_j = -(j-1)/j b_j/b_(j-2),
  !   mu_j~ = mu_j w1, gamma_j~ = -(1-b_(j-1)) mu_j~; T(t+h) = Ys
  !
  use mod_types
  use mod_bound   , only: boundp
  use mod_scal    , only: scal_cht_diff
  use mod_cht_tang, only: cht_tang_t,cht_tang_apply
  implicit none
  private
  public cht_rkl2_step
  real(rp), allocatable, dimension(:,:,:,:), save :: yb  ! Y0 and three rotating stages, with halos
  real(rp), allocatable, dimension(:,:,:)  , save :: m0,lw
  contains
  subroutine cht_rkl2_step(n,dli,dl,dzci,dzfi,dzc,cbc,bc,nb,is_bound,h,ns,alphaf,alphas,cs_cf,src_f,src_s,qw,adens, &
                           thx,thy,thz,phi,is_tang,tt,s)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli,dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzc
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp), intent(in), dimension(0:1,3) :: bc
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in) :: h,alphaf,alphas,cs_cf,src_f,src_s,qw
    integer , intent(in) :: ns
    real(rp), intent(in), dimension(0:,0:,0:) :: adens,thx,thy,thz,phi
    logical , intent(in) :: is_tang
    type(cht_tang_t), intent(in) :: tt
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    integer  :: i,j,k,jst,im1,im2,inew
    real(rp) :: w1,mu,nu,mut,gat,c0,b0,bm1,bm2
    !
    if(.not.allocated(yb)) then
      allocate(yb(0:n(1)+1,0:n(2)+1,0:n(3)+1,4),m0(n(1),n(2),n(3)),lw(n(1),n(2),n(3)))
      !$acc enter data create(yb,m0,lw)
    end if
    w1 = 4._rp/(ns*ns+ns-2)
    !
    ! Y0 = T (halos up to date), M0 = L(Y0), Y1 = Y0 + w1/3 h M0
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          yb(i,j,k,1) = s(i,j,k)
        end do
      end do
    end do
    call cht_op(yb(:,:,:,1),m0)
    mut = w1/3._rp
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          yb(i,j,k,2) = yb(i,j,k,1) + mut*h*m0(i,j,k)
        end do
      end do
    end do
    call boundp(cbc,n,bc,nb,is_bound,dl,dzc,yb(:,:,:,2))
    !
    ! stages 2..s: slots rotate among 2, 3, 4 (im2 = Y(j-2), im1 = Y(j-1), inew = Yj); Y(-1) = Y0 at j = 2
    !
    im2 = 1; im1 = 2
    do jst=2,ns
      inew = 2 + mod(jst-1,3) ! 3, 4, 2, 3, ...: never Y(j-1), Y(j-2) or Y0
      b0  = bcoef(jst); bm1 = bcoef(jst-1); bm2 = bcoef(jst-2)
      mu  = (2._rp*jst-1._rp)/jst*b0/bm1
      nu  = -(jst-1._rp)/jst*b0/bm2
      mut = mu*w1
      gat = -(1._rp-bm1)*mut
      c0  = 1._rp-mu-nu
      call cht_op(yb(:,:,:,im1),lw)
      call stage_update(yb(:,:,:,inew),yb(:,:,:,im1),yb(:,:,:,im2),yb(:,:,:,1),lw,m0,mu,nu,c0,mut*h,gat*h)
      call boundp(cbc,n,bc,nb,is_bound,dl,dzc,yb(:,:,:,inew))
      im2 = im1; im1 = inew
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          s(i,j,k) = yb(i,j,k,im1)
        end do
      end do
    end do
    contains
    !
    pure real(rp) function bcoef(j)
      integer, intent(in) :: j
      if(j < 2) then
        bcoef = 1._rp/3._rp
      else
        bcoef = (j*j+j-2._rp)/(2._rp*j*(j+1._rp))
      end if
    end function bcoef
    !
    subroutine cht_op(y,ly)
      !
      ! L(y): series conduction and sources over the cell capacity, plus the tangential term
      !
      real(rp), intent(in) , dimension(0:,0:,0:) :: y
      real(rp), intent(out), dimension(:,:,:) :: ly
      call scal_cht_diff(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,alphaf,alphas,cs_cf,thx,thy,thz,phi, &
                         src_f,src_s,qw,adens,y,ly)
      if(is_tang .and. tt%nc > 0) call cht_tang_apply(tt%nc,tt%cijk,tt%cf,tt%fijk,tt%fc,dli,dzfi,cs_cf,phi,y,ly)
    end subroutine cht_op
    !
    subroutine stage_update(yn,y1,y2,y0,ly,lm0,mu_,nu_,c0_,ah,bh)
      real(rp), intent(out), dimension(0:,0:,0:) :: yn
      real(rp), intent(in) , dimension(0:,0:,0:) :: y1,y2,y0
      real(rp), intent(in) , dimension(:,:,:) :: ly,lm0
      real(rp), intent(in) :: mu_,nu_,c0_,ah,bh
      integer :: i,j,k
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            yn(i,j,k) = mu_*y1(i,j,k) + nu_*y2(i,j,k) + c0_*y0(i,j,k) + ah*ly(i,j,k) + bh*lm0(i,j,k)
          end do
        end do
      end do
    end subroutine stage_update
  end subroutine cht_rkl2_step
end module mod_cht_rkl2
