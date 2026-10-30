! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_updiss
  !
  ! upwind dissipation near the immersed body (`ibm_updiss_mode /= 'none'`):
  ! the dissipative part of third-order upwind advection in x and y,
  !   dq/dt += -d/dx( c*phi_x*|U|*dx**3/12*d3q/dx3 ) - d/dy( c*phi_y*|V|*dy**3/12*d3q/dy3 ),
  ! written in flux form (momentum conservative, phi_x and phi_y averaged to the flux faces);
  ! it damps the grid-scale (2*dx) wiggles that the non-dissipative central advection
  ! develops upstream of a sharp immersed body. The band weights phi_x, phi_y in [0,1] are
  !  - 'slab': a function of z only, 1 up to ibm_updiss_band*hmax from a rough z-wall,
  !            then a cosine taper to 0 over a further hmax (phi_x = phi_y);
  !  - 'body': a function of the distance d to the body along the x (phi_x) or y (phi_y)
  !            grid line through the point, 1 up to d = ibm_updiss_ncell cells, then a cosine
  !            taper to 0 over a further ibm_updiss_ntaper cells; the wall is located with
  !            `isInbody`, as for the IB coefficients, and the solid points have phi = 1
  !
  use mpi
  use mod_bound   , only: boundp
  use mod_common_mpi, only: myid,ierr
  use mod_ibm     , only: isInbody
  use mod_param   , only: ibm_updiss,ibm_updiss_mode,ibm_updiss_band,ibm_updiss_ncell,ibm_updiss_ntaper, &
                          ibm_direction,amp_l,n_wave,l_0,phase_l
  use mod_thakkar , only: hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha
  use mod_types
  implicit none
  private
  public updiss_init,mom_updiss
  character(len=1), dimension(0:1,3) :: cbc_g
  real(rp)        , dimension(0:1,3) :: bc_g
  integer         , dimension(0:1,3) :: nb_g
  logical         , dimension(0:1,3) :: is_bound_g
  real(rp)        , dimension(3)     :: dl_g
  real(rp), allocatable, dimension(:)         :: dzc_g
  logical , allocatable, dimension(:,:)       :: kact      ! (0:n3+1,c): planes where component c is damped
  real(rp), allocatable, dimension(:,:,:,:,:) :: phi       ! (0:n1+1,0:n2+1,0:n3+1,x/y,c)
  real(rp), allocatable, dimension(:,:,:)     :: gx,gy
  contains
  subroutine updiss_init(n,lo,l,dl,zc,zf,dzc,hmax,cbcpre,bcpre,nb,is_bound)
    !
    ! band weights at the u, v and w points (halos included, so that no exchange is needed);
    ! the halos of the second differences use the pressure boundary conditions
    !
    implicit none
    integer , intent(in), dimension(3)     :: n,lo
    real(rp), intent(in), dimension(3)     :: l,dl
    real(rp), intent(in), dimension(0:)    :: zc,zf,dzc
    real(rp), intent(in)                   :: hmax
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    real(rp), intent(in), dimension(0:1,3) :: bcpre
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), parameter :: pi = acos(-1._rp)
    real(rp) :: xyz(3),d,frac(5)
    integer  :: i,j,k,c,a
    !
    cbc_g(:,:) = cbcpre(:,:); bc_g(:,:) = bcpre(:,:); nb_g(:,:) = nb(:,:)
    is_bound_g(:,:) = is_bound(:,:); dl_g(:) = dl(:)
    allocate(dzc_g(0:n(3)+1),kact(0:n(3)+1,3))
    allocate(phi(0:n(1)+1,0:n(2)+1,0:n(3)+1,2,3))
    dzc_g(:) = dzc(0:n(3)+1)
    phi(:,:,:,:,:) = 0._rp
    do c=1,3
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            xyz(1) = (lo(1)+i-1-0.5_rp+merge(0.5_rp,0._rp,c == 1))*dl(1)
            xyz(2) = (lo(2)+j-1-0.5_rp+merge(0.5_rp,0._rp,c == 2))*dl(2)
            xyz(3) = merge(zf(k),zc(k),c == 3)
            select case(trim(ibm_updiss_mode))
            case('slab')
              d = huge(1._rp)
              if(ibm_direction(0,3)) d = min(d,xyz(3))
              if(ibm_direction(1,3)) d = min(d,l(3)-xyz(3))
              phi(i,j,k,:,c) = taper(d,ibm_updiss_band*hmax,hmax)
            case('body')
              if(near_body(xyz(3))) then
                if(inbody(xyz)) then
                  phi(i,j,k,:,c) = 1._rp
                else
                  do a=1,2
                    phi(i,j,k,a,c) = taper(line_dist(xyz,a)/dl(a),ibm_updiss_ncell,ibm_updiss_ntaper)
                  end do
                end if
              end if
            end select
          end do
        end do
        kact(k,c) = any(phi(:,:,k,:,c) > 0._rp)
      end do
    end do
    !
    ! share of the domain where the u term acts (phi > 0) and the mean of phi, for x and y
    ! (the solid points, where phi = 1 but the term is skipped, are included)
    !
    frac(1:2) = [count(phi(1:n(1),1:n(2),1:n(3),1,1) > 0._rp),count(phi(1:n(1),1:n(2),1:n(3),2,1) > 0._rp)]
    frac(3:4) = [sum(phi(1:n(1),1:n(2),1:n(3),1,1)),sum(phi(1:n(1),1:n(2),1:n(3),2,1))]
    frac(5)   = product(1._rp*n(:))
    call MPI_ALLREDUCE(MPI_IN_PLACE,frac,5,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) print*,'upwind dissipation (',trim(ibm_updiss_mode),'): volume share, mean weight (x, y) = ',frac(1:4)/frac(5)
    allocate(gx(0:n(1)+1,0:n(2)+1,0:n(3)+1),gy(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    gx(:,:,:) = 0._rp; gy(:,:,:) = 0._rp
    !$acc enter data copyin(kact,phi,gx,gy,dzc_g)
  contains
    real(rp) function taper(d,h,w)
      !
      ! 1 for d <= h, cosine taper to 0 at d = h+w
      !
      real(rp), intent(in) :: d,h,w
      taper = 0._rp
      if(d <= h) then
        taper = 1._rp
      else if(d < h+w) then
        taper = 0.5_rp*(1._rp+cos(pi*(d-h)/w))
      end if
    end function taper
    logical function near_body(z)
      !
      ! the x and y grid lines at height z can meet the body (z-walls only: within hmax)
      !
      real(rp), intent(in) :: z
      near_body = any(ibm_direction(:,1:2)) .or. &
                  (ibm_direction(0,3) .and. z <= hmax) .or. (ibm_direction(1,3) .and. z >= l(3)-hmax)
    end function near_body
    logical function inbody(p)
      real(rp), intent(in) :: p(3)
      inbody = isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,p(1),p(2),p(3),n,l, &
                        hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha)
    end function inbody
    real(rp) function line_dist(p,a)
      !
      ! distance from the fluid point p to the body along direction a (both senses), found by
      ! stepping one cell at a time and bisecting the first fluid/solid interval; huge if the
      ! body is farther than the band
      !
      real(rp), intent(in) :: p(3)
      integer , intent(in) :: a
      integer  :: s,m,it
      real(rp) :: q(3),df,ds,dm
      line_dist = huge(1._rp)
      do s=-1,1,2
        q(:) = p(:)
        do m=1,ceiling(ibm_updiss_ncell+ibm_updiss_ntaper)+1
          q(a) = p(a)+s*m*dl(a)
          if(inbody(q)) then
            df = (m-1)*dl(a); ds = m*dl(a)
            do it=1,40
              dm = 0.5_rp*(df+ds)
              q(a) = p(a)+s*dm
              if(inbody(q)) then
                ds = dm
              else
                df = dm
              end if
            end do
            line_dist = min(line_dist,0.5_rp*(df+ds))
            exit
          end if
        end do
      end do
    end function line_dist
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
    fx = 0.5_rp*ibm_updiss*dli(1)/12._rp ! 0.5: face average of phi
    fy = 0.5_rp*ibm_updiss*dli(2)/12._rp
    !
    ! u at (i+1/2,j,k): x fluxes at the cell centres, y fluxes at (i+1/2,j+1/2)
    !
    call diff2(n,kact(:,1),u)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(kact(k,1) .and. .not.mask_u(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i+1,j  ,k)))*(phi(i,j,k,1,1)+phi(i+1,j,k,1,1))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i  ,j  ,k)))*(phi(i,j,k,1,1)+phi(i-1,j,k,1,1))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i+1,j  ,k)))*(phi(i,j,k,2,1)+phi(i,j+1,k,2,1))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i+1,j-1,k)))*(phi(i,j,k,2,1)+phi(i,j-1,k,2,1))
            dudt(i,j,k) = dudt(i,j,k) - ( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
    !
    ! v at (i,j+1/2,k): x fluxes at (i+1/2,j+1/2), y fluxes at the cell centres
    !
    call diff2(n,kact(:,2),v)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(kact(k,2) .and. .not.mask_v(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j+1,k)))*(phi(i,j,k,1,2)+phi(i+1,j,k,1,2))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j+1,k)))*(phi(i,j,k,1,2)+phi(i-1,j,k,1,2))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j+1,k)))*(phi(i,j,k,2,2)+phi(i,j+1,k,2,2))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j  ,k)))*(phi(i,j,k,2,2)+phi(i,j-1,k,2,2))
            dvdt(i,j,k) = dvdt(i,j,k) - ( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
    !
    ! w at (i,j,k+1/2): x fluxes at (i+1/2,j,k+1/2), y fluxes at (i,j+1/2,k+1/2)
    !
    call diff2(n,kact(:,3),w)
    !$acc parallel loop collapse(3) default(present) private(uxp,uxm,vyp,vym) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uxp,uxm,vyp,vym)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(kact(k,3) .and. .not.mask_w(i,j,k)) then
            uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j  ,k+1)))*(phi(i,j,k,1,3)+phi(i+1,j,k,1,3))
            uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j  ,k+1)))*(phi(i,j,k,1,3)+phi(i-1,j,k,1,3))
            vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j  ,k+1)))*(phi(i,j,k,2,3)+phi(i,j+1,k,2,3))
            vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j-1,k+1)))*(phi(i,j,k,2,3)+phi(i,j-1,k,2,3))
            dwdt(i,j,k) = dwdt(i,j,k) - ( &
                          fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                          fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
          end if
        end do
      end do
    end do
  end subroutine mom_updiss
  !
  subroutine diff2(n,act,q)
    !
    ! undivided second differences of q in x (gx) and y (gy) on the damped planes,
    ! with the halos exchanged/imposed for the flux stencil of the next pass
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:) :: act
    real(rp), intent(in), dimension(0:,0:,0:) :: q
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(act(k)) then
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
