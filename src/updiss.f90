! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_updiss
  !
  ! upwind dissipation in a band above the rough z-walls (`ibm_updiss > 0`):
  ! the dissipative part of third-order upwind advection in x and y,
  !   dq/dt += -d/dx( c*|U|*dx**3/12*d3q/dx3 ) - d/dy( c*|V|*dy**3/12*d3q/dy3 ),
  ! written in flux form (momentum conservative) and weighted by a band function phi(z);
  ! it damps the grid-scale (2*dx) wiggles that the non-dissipative central advection
  ! develops upstream of a sharp immersed body
  !
  use mod_bound, only: boundp
  use mod_param, only: ibm_updiss,ibm_updiss_band
  use mod_types
  implicit none
  private
  public updiss_init,mom_updiss
  character(len=1), dimension(0:1,3) :: cbc_g
  real(rp)        , dimension(0:1,3) :: bc_g
  integer         , dimension(0:1,3) :: nb_g
  logical         , dimension(0:1,3) :: is_bound_g
  real(rp)        , dimension(3)     :: dl_g
  real(rp), allocatable, dimension(:)     :: dzc_g,phic,phif
  real(rp), allocatable, dimension(:,:,:) :: gx,gy
  contains
  subroutine updiss_init(n,l,zc,zf,dzc,hmax,is_rough,cbcpre,bcpre,nb,is_bound,dl)
    !
    ! band weights at the cell centres (phic) and z faces (phif): 1 up to a distance
    ! ibm_updiss_band*hmax from a rough z-wall, then a cosine taper to 0 over hmax;
    ! the halos of the second differences use the pressure boundary conditions
    !
    implicit none
    integer , intent(in), dimension(3)     :: n
    real(rp), intent(in), dimension(3)     :: l
    real(rp), intent(in), dimension(0:)    :: zc,zf,dzc
    real(rp), intent(in)                   :: hmax
    logical , intent(in), dimension(0:1)   :: is_rough
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    real(rp), intent(in), dimension(0:1,3) :: bcpre
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3)     :: dl
    integer :: k
    !
    cbc_g(:,:) = cbcpre(:,:); bc_g(:,:) = bcpre(:,:); nb_g(:,:) = nb(:,:)
    is_bound_g(:,:) = is_bound(:,:); dl_g(:) = dl(:)
    allocate(dzc_g(0:n(3)+1),phic(0:n(3)+1),phif(0:n(3)+1))
    dzc_g(:) = dzc(0:n(3)+1)
    do k=0,n(3)+1
      phic(k) = band(zc(k))
      phif(k) = band(zf(k))
    end do
    allocate(gx(0:n(1)+1,0:n(2)+1,0:n(3)+1),gy(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    gx(:,:,:) = 0._rp; gy(:,:,:) = 0._rp
    !$acc enter data copyin(phic,phif,gx,gy,dzc_g)
  contains
    real(rp) function band(z)
      real(rp), intent(in) :: z
      real(rp) :: h,d
      band = 0._rp
      h = ibm_updiss_band*hmax
      d = huge(1._rp)
      if(is_rough(0)) d = min(d,z)
      if(is_rough(1)) d = min(d,l(3)-z)
      if(d <= h) then
        band = 1._rp
      else if(d < h+hmax) then
        band = 0.5_rp*(1._rp+cos(acos(-1._rp)*(d-h)/hmax))
      end if
    end function band
  end subroutine updiss_init
  !
  subroutine mom_updiss(n,dli,u,v,w,mask_u,mask_v,mask_w,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    logical , intent(in   ), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w
    real(rp), intent(inout), dimension(:,:,:)    :: dudt,dvdt,dwdt
    real(rp) :: fx,fy,uxp,uxm,vyp,vym
    integer  :: i,j,k
    !
    fx = ibm_updiss*dli(1)/12._rp
    fy = ibm_updiss*dli(2)/12._rp
    !
    ! u at (i+1/2,j,k): x fluxes at the cell centres, y fluxes at (i+1/2,j+1/2)
    !
    call diff2(n,phic,u)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(phic(k) > 0._rp .and. .not.mask_u(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i+1,j  ,k)))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i  ,j  ,k)))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i+1,j  ,k)))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i+1,j-1,k)))
            dudt(i,j,k) = dudt(i,j,k) - phic(k)*( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
    !
    ! v at (i,j+1/2,k): x fluxes at (i+1/2,j+1/2), y fluxes at the cell centres
    !
    call diff2(n,phic,v)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(phic(k) > 0._rp .and. .not.mask_v(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j+1,k)))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j+1,k)))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j+1,k)))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j  ,k)))
            dvdt(i,j,k) = dvdt(i,j,k) - phic(k)*( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
    !
    ! w at (i,j,k+1/2): x fluxes at (i+1/2,j,k+1/2), y fluxes at (i,j+1/2,k+1/2)
    !
    call diff2(n,phif,w)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(phif(k) > 0._rp .and. .not.mask_w(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j  ,k+1)))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j  ,k+1)))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j  ,k+1)))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j-1,k+1)))
            dwdt(i,j,k) = dwdt(i,j,k) - phif(k)*( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
  end subroutine mom_updiss
  !
  subroutine diff2(n,phi,q)
    !
    ! undivided second differences of q in x (gx) and y (gy) inside the band,
    ! with the halos exchanged/imposed for the flux stencil of the next pass
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: phi
    real(rp), intent(in), dimension(0:,0:,0:) :: q
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(phi(k) > 0._rp) then
            gx(i,j,k) = q(i+1,j,k)-2._rp*q(i,j,k)+q(i-1,j,k)
            gy(i,j,k) = q(i,j+1,k)-2._rp*q(i,j,k)+q(i,j-1,k)
          end if
        end do
      end do
    end do
    call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,gx)
    call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,gy)
  end subroutine diff2
end module mod_updiss
