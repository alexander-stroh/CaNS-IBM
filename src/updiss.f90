! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_updiss
  !
  ! upwind dissipation near the immersed body (`ibm_updiss_mode /= 'none'`):
  ! the dissipative part of third-order upwind advection (for the velocity components and,
  ! with `ibm_updiss_scal`, for the scalars),
  !   dq/dt += -d/dx( c*phi*|U|*dx**3/12*d3q/dx3 ) - d/dy( c*phi*|V|*dy**3/12*d3q/dy3 )
  !            -d/dz( c*phi*|W|*dz**3/12*d3q/dz3 ),
  ! written in flux form (momentum conservative, phi averaged to the flux faces); in z with
  ! undivided differences on the stretched grid. It damps the grid-scale (2*dx) wiggles that
  ! the non-dissipative central advection develops upstream of a sharp immersed body.
  ! The weight phi in [0,1] ('body') is a layer around the immersed wall, grown from the IB
  ! mask of each velocity component (and of the scalar): level 1 are the fluid points with a solid neighbour (the
  ! points with IB coefficients), level m the fluid points next to level m-1 (neighbours in
  ! all 26 directions, so that the layer is isotropic in grid cells). phi = 1 up to level
  ! ibm_updiss_ncell, then a cosine taper over a further ibm_updiss_ntaper levels; phi = 1 in
  ! the solid. It is independent of the geometry and of the wall orientation
  !
  ! the term is evaluated only on compact lists of the points where it acts and, for the
  ! second differences, of their neighbours
  !
  use mpi
  use mod_bound   , only: boundp
  use mod_common_mpi, only: myid,ierr
  use mod_param   , only: ibm_updiss,ibm_updiss_ncell,ibm_updiss_ntaper,ibm_updiss_scal
  use mod_types
  implicit none
  private
  public updiss_init,mom_updiss,scal_updiss
  character(len=1), dimension(0:1,3) :: cbc_g
  real(rp)        , dimension(0:1,3) :: bc_g
  integer         , dimension(0:1,3) :: nb_g
  logical         , dimension(0:1,3) :: is_bound_g
  real(rp)        , dimension(3)     :: dl_g
  real(rp), allocatable, dimension(:)       :: dzc_g,dzci_g,dzfi_g
  real(rp), allocatable, dimension(:,:,:,:) :: phi       ! (0:n1+1,0:n2+1,0:n3+1,c), c = 4: scalar
  real(rp), allocatable, dimension(:,:,:)   :: gx,gy,gz
  integer , allocatable, dimension(:,:,:)   :: ipt,ipg   ! (3,np,c): points of the term, of the differences
  integer , dimension(4) :: npt,npg
  contains
  subroutine updiss_init(n,dl,dzc,dzf,mask_u,mask_v,mask_w,mask_s,cbcpre,bcpre,nb,is_bound)
    !
    ! weights at the u, v and w points (and at the cell centres, with `ibm_updiss_scal`;
    ! halos included) and the point lists;
    ! the halos of the second differences use the pressure boundary conditions
    !
    implicit none
    integer , intent(in), dimension(3)     :: n
    real(rp), intent(in), dimension(3)     :: dl
    real(rp), intent(in), dimension(0:)    :: dzc,dzf
    logical , intent(in), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w,mask_s
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    real(rp), intent(in), dimension(0:1,3) :: bcpre
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), parameter :: pi = acos(-1._rp)
    real(rp), allocatable, dimension(:,:,:,:) :: lev
    logical , allocatable, dimension(:,:,:,:) :: ist,isg
    real(rp) :: frac(3)
    integer  :: i,j,k,c,m,nl,nc
    !
    cbc_g(:,:) = cbcpre(:,:); bc_g(:,:) = bcpre(:,:); nb_g(:,:) = nb(:,:)
    is_bound_g(:,:) = is_bound(:,:); dl_g(:) = dl(:)
    allocate(dzc_g(0:n(3)+1),dzci_g(0:n(3)+1),dzfi_g(0:n(3)+1))
    dzc_g(:) = dzc(0:n(3)+1); dzci_g(:) = 1._rp/dzc(0:n(3)+1); dzfi_g(:) = 1._rp/dzf(0:n(3)+1)
    nc = merge(4,3,ibm_updiss_scal)
    allocate(phi(0:n(1)+1,0:n(2)+1,0:n(3)+1,nc))
    phi(:,:,:,:) = 0._rp
    !
    ! level of each point: 0 in the solid, m in the m-th fluid layer, nl+3 beyond;
    ! levels nl+1 and nl+2 (phi = 0) are needed by the flux stencil
    !
    nl = ibm_updiss_ncell+ibm_updiss_ntaper
    allocate(lev(0:n(1)+1,0:n(2)+1,0:n(3)+1,nc))
    lev(:,:,:,1) = merge(0._rp,nl+3._rp,mask_u(:,:,:))
    lev(:,:,:,2) = merge(0._rp,nl+3._rp,mask_v(:,:,:))
    lev(:,:,:,3) = merge(0._rp,nl+3._rp,mask_w(:,:,:))
    if(nc == 4) lev(:,:,:,4) = merge(0._rp,nl+3._rp,mask_s(:,:,:))
    do c=1,nc
      do m=1,nl+2
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if(lev(i,j,k,c) > m .and. any(lev(i-1:i+1,j-1:j+1,k-1:k+1,c) == m-1)) lev(i,j,k,c) = m
            end do
          end do
        end do
        call halo(lev(:,:,:,c))
      end do
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            m = nint(lev(i,j,k,c))
            if(m <= ibm_updiss_ncell) then
              phi(i,j,k,c) = 1._rp
            else if(m <= nl) then
              phi(i,j,k,c) = 0.5_rp*(1._rp+cos(pi*(m-ibm_updiss_ncell)/(ibm_updiss_ntaper+1._rp)))
            end if
          end do
        end do
      end do
    end do
    !
    ! point lists: the fluid points with a nonzero face weight (the term) and the points
    ! next to them (the second differences)
    !
    allocate(ist(n(1),n(2),n(3),nc),isg(n(1),n(2),n(3),nc))
    isg(:,:,:,:) = lev(1:n(1),1:n(2),1:n(3),:) <= nl+2
    ist(:,:,:,:) = lev(1:n(1),1:n(2),1:n(3),:) <= nl+1
    deallocate(lev)
    ist(:,:,:,1) = ist(:,:,:,1) .and. .not.mask_u(1:n(1),1:n(2),1:n(3))
    ist(:,:,:,2) = ist(:,:,:,2) .and. .not.mask_v(1:n(1),1:n(2),1:n(3))
    ist(:,:,:,3) = ist(:,:,:,3) .and. .not.mask_w(1:n(1),1:n(2),1:n(3))
    if(nc == 4) ist(:,:,:,4) = ist(:,:,:,4) .and. .not.mask_s(1:n(1),1:n(2),1:n(3))
    npt(:) = 0; npg(:) = 0
    do c=1,nc
      npt(c) = count(ist(:,:,:,c)); npg(c) = count(isg(:,:,:,c))
    end do
    allocate(ipt(3,max(maxval(npt),1),nc),ipg(3,max(maxval(npg),1),nc))
    ipt(:,:,:) = 1; ipg(:,:,:) = 1
    do c=1,nc
      call fill(ist(:,:,:,c),ipt(:,:,c))
      call fill(isg(:,:,:,c),ipg(:,:,c))
    end do
    !
    ! share of the domain where the u term acts (phi > 0, fluid) and the mean of phi there
    !
    frac(1) = count(phi(1:n(1),1:n(2),1:n(3),1) > 0._rp .and. .not.mask_u(1:n(1),1:n(2),1:n(3)))
    frac(2) = sum(phi(1:n(1),1:n(2),1:n(3),1),mask=.not.mask_u(1:n(1),1:n(2),1:n(3)))
    frac(3) = product(1._rp*n(:))
    call MPI_ALLREDUCE(MPI_IN_PLACE,frac,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) print*,'upwind dissipation near the body: share of grid points (u), mean weight there = ', &
                         frac(1)/frac(3),frac(2)/max(frac(1),1._rp)
    deallocate(ist,isg)
    allocate(gx(0:n(1)+1,0:n(2)+1,0:n(3)+1),gy(0:n(1)+1,0:n(2)+1,0:n(3)+1),gz(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    gx(:,:,:) = 0._rp; gy(:,:,:) = 0._rp; gz(:,:,:) = 0._rp
    !$acc enter data copyin(phi,ipt,ipg,gx,gy,gz,dzc_g,dzci_g,dzfi_g)
  contains
    subroutine halo(p)
      !
      ! exchanges the halos of the levels; at the domain boundaries the halos keep the
      ! values from the mask (0 in the solid, nl+3 in the fluid)
      !
      real(rp), intent(inout), dimension(0:,0:,0:) :: p
      real(rp), allocatable, dimension(:,:,:) :: p0
      allocate(p0,source=p)
      !$acc data copy(p)
      call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,p)
      !$acc end data
      if(is_bound(0,1)) p(0     ,:,:) = p0(0     ,:,:)
      if(is_bound(1,1)) p(n(1)+1,:,:) = p0(n(1)+1,:,:)
      if(is_bound(0,2)) p(:,0     ,:) = p0(:,0     ,:)
      if(is_bound(1,2)) p(:,n(2)+1,:) = p0(:,n(2)+1,:)
      if(is_bound(0,3)) p(:,:,0     ) = p0(:,:,0     )
      if(is_bound(1,3)) p(:,:,n(3)+1) = p0(:,:,n(3)+1)
    end subroutine halo
    subroutine fill(is,ip)
      logical, intent(in ), dimension(:,:,:) :: is
      integer, intent(out), dimension(:,:)   :: ip
      integer :: i,j,k,np
      np = 0
      do k=1,size(is,3)
        do j=1,size(is,2)
          do i=1,size(is,1)
            if(is(i,j,k)) then
              np = np+1
              ip(:,np) = [i,j,k]
            end if
          end do
        end do
      end do
    end subroutine fill
  end subroutine updiss_init
  !
  subroutine mom_updiss(n,dli,u,v,w,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:)    :: dudt,dvdt,dwdt
    real(rp) :: fx,fy,fz,uxp,uxm,vyp,vym,wzp,wzm,ph
    integer  :: i,j,k,ip,np
    !
    fx = 0.5_rp*ibm_updiss*dli(1)/12._rp ! 0.5: face average of phi
    fy = 0.5_rp*ibm_updiss*dli(2)/12._rp
    fz = 0.5_rp*ibm_updiss/12._rp        ! times 1/dz of the control volume
    !
    ! u at (i+1/2,j,k): x fluxes at the cell centres, y fluxes at (i+1/2,j+1/2), z fluxes at (i+1/2,j,k+1/2)
    !
    call diff2(n,1,u)
    np = npt(1)
    !$acc parallel loop gang vector default(present) private(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm)
    do ip=1,np
      i = ipt(1,ip,1); j = ipt(2,ip,1); k = ipt(3,ip,1)
      ph  = phi(i,j,k,1)
      uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i+1,j  ,k)))*(ph+phi(i+1,j,k,1))
      uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i  ,j  ,k)))*(ph+phi(i-1,j,k,1))
      vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i+1,j  ,k)))*(ph+phi(i,j+1,k,1))
      vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i+1,j-1,k)))*(ph+phi(i,j-1,k,1))
      dudt(i,j,k) = dudt(i,j,k) - ( &
                    fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                    fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
      wzp = abs(0.5_rp*(w(i  ,j  ,k  )+w(i+1,j  ,k  )))*(ph+phi(i,j,k+1,1))
      wzm = abs(0.5_rp*(w(i  ,j  ,k-1)+w(i+1,j  ,k-1)))*(ph+phi(i,j,k-1,1))
      dudt(i,j,k) = dudt(i,j,k) - fz*dzfi_g(k)*(wzp*(gz(i,j,k+1)-gz(i,j,k))-wzm*(gz(i,j,k)-gz(i,j,k-1)))
    end do
    !
    ! v at (i,j+1/2,k): x fluxes at (i+1/2,j+1/2), y fluxes at the cell centres, z fluxes at (i,j+1/2,k+1/2)
    !
    call diff2(n,2,v)
    np = npt(2)
    !$acc parallel loop gang vector default(present) private(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm)
    do ip=1,np
      i = ipt(1,ip,2); j = ipt(2,ip,2); k = ipt(3,ip,2)
      ph  = phi(i,j,k,2)
      uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j+1,k)))*(ph+phi(i+1,j,k,2))
      uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j+1,k)))*(ph+phi(i-1,j,k,2))
      vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j+1,k)))*(ph+phi(i,j+1,k,2))
      vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j  ,k)))*(ph+phi(i,j-1,k,2))
      dvdt(i,j,k) = dvdt(i,j,k) - ( &
                    fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                    fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
      wzp = abs(0.5_rp*(w(i  ,j  ,k  )+w(i  ,j+1,k  )))*(ph+phi(i,j,k+1,2))
      wzm = abs(0.5_rp*(w(i  ,j  ,k-1)+w(i  ,j+1,k-1)))*(ph+phi(i,j,k-1,2))
      dvdt(i,j,k) = dvdt(i,j,k) - fz*dzfi_g(k)*(wzp*(gz(i,j,k+1)-gz(i,j,k))-wzm*(gz(i,j,k)-gz(i,j,k-1)))
    end do
    !
    ! w at (i,j,k+1/2): x fluxes at (i+1/2,j,k+1/2), y fluxes at (i,j+1/2,k+1/2), z fluxes at the cell centres
    !
    call diff2(n,3,w)
    np = npt(3)
    !$acc parallel loop gang vector default(present) private(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm)
    do ip=1,np
      i = ipt(1,ip,3); j = ipt(2,ip,3); k = ipt(3,ip,3)
      ph  = phi(i,j,k,3)
      uxp = abs(0.5_rp*(u(i  ,j  ,k)+u(i  ,j  ,k+1)))*(ph+phi(i+1,j,k,3))
      uxm = abs(0.5_rp*(u(i-1,j  ,k)+u(i-1,j  ,k+1)))*(ph+phi(i-1,j,k,3))
      vyp = abs(0.5_rp*(v(i  ,j  ,k)+v(i  ,j  ,k+1)))*(ph+phi(i,j+1,k,3))
      vym = abs(0.5_rp*(v(i  ,j-1,k)+v(i  ,j-1,k+1)))*(ph+phi(i,j-1,k,3))
      dwdt(i,j,k) = dwdt(i,j,k) - ( &
                    fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                    fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) )
      wzp = abs(0.5_rp*(w(i  ,j  ,k  )+w(i  ,j  ,k+1)))*(ph+phi(i,j,k+1,3))
      wzm = abs(0.5_rp*(w(i  ,j  ,k-1)+w(i  ,j  ,k  )))*(ph+phi(i,j,k-1,3))
      dwdt(i,j,k) = dwdt(i,j,k) - fz*dzci_g(k)*(wzp*(gz(i,j,k+1)-gz(i,j,k))-wzm*(gz(i,j,k)-gz(i,j,k-1)))
    end do
  end subroutine mom_updiss
  !
  subroutine scal_updiss(n,dli,u,v,w,s,dsdt)
    !
    ! the same term for a scalar at the cell centres: x fluxes at (i+1/2,j,k), y fluxes at
    ! (i,j+1/2,k), z fluxes at (i,j,k+1/2), i.e. at the velocity points
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,s
    real(rp), intent(inout), dimension(:,:,:)    :: dsdt
    real(rp) :: fx,fy,fz,uxp,uxm,vyp,vym,wzp,wzm,ph
    integer  :: i,j,k,ip,np
    !
    fx = 0.5_rp*ibm_updiss*dli(1)/12._rp
    fy = 0.5_rp*ibm_updiss*dli(2)/12._rp
    fz = 0.5_rp*ibm_updiss/12._rp
    call diff2(n,4,s)
    np = npt(4)
    !$acc parallel loop gang vector default(present) private(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(i,j,k,ph,uxp,uxm,vyp,vym,wzp,wzm)
    do ip=1,np
      i = ipt(1,ip,4); j = ipt(2,ip,4); k = ipt(3,ip,4)
      ph  = phi(i,j,k,4)
      uxp = abs(u(i  ,j,k))*(ph+phi(i+1,j,k,4))
      uxm = abs(u(i-1,j,k))*(ph+phi(i-1,j,k,4))
      vyp = abs(v(i,j  ,k))*(ph+phi(i,j+1,k,4))
      vym = abs(v(i,j-1,k))*(ph+phi(i,j-1,k,4))
      wzp = abs(w(i,j,k  ))*(ph+phi(i,j,k+1,4))
      wzm = abs(w(i,j,k-1))*(ph+phi(i,j,k-1,4))
      dsdt(i,j,k) = dsdt(i,j,k) - ( &
                    fx*(uxp*(gx(i+1,j,k)-gx(i,j,k))-uxm*(gx(i,j,k)-gx(i-1,j,k))) + &
                    fy*(vyp*(gy(i,j+1,k)-gy(i,j,k))-vym*(gy(i,j,k)-gy(i,j-1,k))) + &
                    fz*dzfi_g(k)*(wzp*(gz(i,j,k+1)-gz(i,j,k))-wzm*(gz(i,j,k)-gz(i,j,k-1))) )
    end do
  end subroutine scal_updiss
  !
  subroutine diff2(n,c,q)
    !
    ! undivided second differences of q in x (gx), y (gy) and z (gz) at the points of
    ! component c that the flux stencil needs, with the halos exchanged/imposed
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: c
    real(rp), intent(in), dimension(0:,0:,0:) :: q
    integer :: i,j,k,ip,np
    !
    np = npg(c)
    !$acc parallel loop gang vector default(present) private(i,j,k) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(i,j,k)
    do ip=1,np
      i = ipg(1,ip,c); j = ipg(2,ip,c); k = ipg(3,ip,c)
      gx(i,j,k) = q(i+1,j,k)-2._rp*q(i,j,k)+q(i-1,j,k)
      gy(i,j,k) = q(i,j+1,k)-2._rp*q(i,j,k)+q(i,j-1,k)
      gz(i,j,k) = q(i,j,k+1)-2._rp*q(i,j,k)+q(i,j,k-1)
    end do
    call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,gx)
    call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,gy)
    call boundp(cbc_g,n,bc_g,nb_g,is_bound_g,dl_g,dzc_g,gz)
  end subroutine diff2
end module mod_updiss
