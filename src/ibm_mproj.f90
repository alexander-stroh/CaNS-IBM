! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_ibm_mproj
  !
  ! IB-consistent ("masked") pressure projection (`ibm_proj = 'masked'`), cf. Seo & Mittal, J. Comput. Phys. 230
  ! (2011) 7347, with binary face apertures: with M = 0 at the solid velocity points (IBM masks) and 1 elsewhere,
  !
  !   L_M pp = D(M G pp) = D(u*)/dt_rk,   u = M (u* - dt_rk G pp),
  !
  ! so that the solid velocity points stay at zero and the divergence vanishes in every cell, including the
  ! cut cells (the plain full-grid projection leaves a transpiration velocity -dt_rk G pp on the solid faces,
  ! zeroed again at the next stage: a mass source/sink at the interface cells); L_M is the Laplacian with
  ! homogeneous Neumann conditions on the staircase surface; it is solved by conjugate gradients
  ! preconditioned with the constant-coefficient FFT solver of the plain Laplacian L0 (both negative
  ! semi-definite); cells with all six faces closed (inside the solid) are decoupled and kept at zero
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
  implicit none
  private
  public mproj_init,mproj_apply,mproj_dot,mproj_zero_closed,mproj_axpy,mproj_copy,mproj_zero_solid,mproj_set, &
         mp_r,mp_z,mp_d,mp_q,mp_x,mp_it,mp_res,mp_it_sum,mp_ncall
  logical , allocatable, dimension(:,:,:) :: closed
  real(rp), allocatable, dimension(:,:,:) :: mp_r,mp_z,mp_d,mp_q,mp_x
  integer  :: mp_it = 0, mp_it_sum = 0, mp_ncall = 0 ! iterations of the last solve, sums for the average
  real(rp) :: mp_res = 0._rp                         ! relative residual of the last solve
  contains
  subroutine mproj_init(n,mask_u,mask_v,mask_w,myid)
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w ! .true. = solid
    integer , intent(in) :: myid
    integer :: i,j,k,nc
    allocate(closed(n(1),n(2),n(3)))
    allocate(mp_r(0:n(1)+1,0:n(2)+1,0:n(3)+1),mp_z(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             mp_d(0:n(1)+1,0:n(2)+1,0:n(3)+1),mp_q(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             mp_x(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    nc = 0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ! the domain walls (k = 0, n(3) faces) are closed for the pressure (homogeneous Neumann)
          closed(i,j,k) = mask_u(i,j,k).and.mask_u(i-1,j,k).and.mask_v(i,j,k).and.mask_v(i,j-1,k).and. &
                          (mask_w(i,j,k).or.k == n(3)).and.(mask_w(i,j,k-1).or.k == 1)
          if(closed(i,j,k)) nc = nc + 1
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,nc,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) print*, '*** IBM masked projection: pressure cells enclosed by the solid: ', nc, ' ***'
    mp_r(:,:,:) = 0._rp; mp_z(:,:,:) = 0._rp; mp_d(:,:,:) = 0._rp; mp_q(:,:,:) = 0._rp; mp_x(:,:,:) = 0._rp
    !$acc enter data copyin(closed,mp_r,mp_z,mp_d,mp_q,mp_x)
  end subroutine mproj_init
  !
  subroutine mproj_apply(n,dli,dzci,dzfi,mask_u,mask_v,mask_w,x,y)
    !
    ! y = D(M G x) in the interior (x with valid halos), 0 in the enclosed cells
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    logical , intent(in), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w
    real(rp), intent(in ), dimension(0:,0:,0:) :: x
    real(rp), intent(out), dimension(0:,0:,0:) :: y
    real(rp) :: dxi2,dyi2,fe,fw,fn,fs,ft,fb
    integer :: i,j,k
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    !$acc parallel loop collapse(3) default(present) private(fe,fw,fn,fs,ft,fb) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) PRIVATE(fe,fw,fn,fs,ft,fb)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          fe = merge(0._rp,x(i+1,j,k)-x(i  ,j,k),mask_u(i  ,j,k))
          fw = merge(0._rp,x(i  ,j,k)-x(i-1,j,k),mask_u(i-1,j,k))
          fn = merge(0._rp,x(i,j+1,k)-x(i,j  ,k),mask_v(i,j  ,k))
          fs = merge(0._rp,x(i,j  ,k)-x(i,j-1,k),mask_v(i,j-1,k))
          ft = merge(0._rp,(x(i,j,k+1)-x(i,j,k  ))*dzci(k  ),mask_w(i,j,k  ))
          fb = merge(0._rp,(x(i,j,k  )-x(i,j,k-1))*dzci(k-1),mask_w(i,j,k-1))
          y(i,j,k) = (fe-fw)*dxi2 + (fn-fs)*dyi2 + (ft-fb)*dzfi(k)
          if(closed(i,j,k)) y(i,j,k) = 0._rp
        end do
      end do
    end do
  end subroutine mproj_apply
  !
  function mproj_dot(n,a,b) result(s)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:,0:,0:) :: a,b
    real(rp) :: s
    integer :: i,j,k
    s = 0._rp
    !$acc wait(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:s)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) REDUCTION(+:s)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s = s + a(i,j,k)*b(i,j,k)
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,s,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end function mproj_dot
  !
  subroutine mproj_zero_closed(n,a)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(inout), dimension(0:,0:,0:) :: a
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(closed(i,j,k)) a(i,j,k) = 0._rp
        end do
      end do
    end do
  end subroutine mproj_zero_closed
  !
  subroutine mproj_axpy(n,is_halo,al,x,y)
    !
    ! y = y + al*x (with the halo layers if is_halo)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in) :: is_halo
    real(rp), intent(in) :: al
    real(rp), intent(in   ), dimension(0:,0:,0:) :: x
    real(rp), intent(inout), dimension(0:,0:,0:) :: y
    integer :: i,j,k,h
    h = merge(1,0,is_halo)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1-h,n(3)+h
      do j=1-h,n(2)+h
        do i=1-h,n(1)+h
          y(i,j,k) = y(i,j,k) + al*x(i,j,k)
        end do
      end do
    end do
  end subroutine mproj_axpy
  !
  subroutine mproj_copy(n,is_halo,be,x,y)
    !
    ! y = x + be*y (with the halo layers if is_halo)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in) :: is_halo
    real(rp), intent(in) :: be
    real(rp), intent(in   ), dimension(0:,0:,0:) :: x
    real(rp), intent(inout), dimension(0:,0:,0:) :: y
    integer :: i,j,k,h
    h = merge(1,0,is_halo)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1-h,n(3)+h
      do j=1-h,n(2)+h
        do i=1-h,n(1)+h
          y(i,j,k) = x(i,j,k) + be*y(i,j,k)
        end do
      end do
    end do
  end subroutine mproj_copy
  !
  subroutine mproj_set(n,val,a)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: val
    real(rp), intent(out), dimension(0:,0:,0:) :: a
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          a(i,j,k) = val
        end do
      end do
    end do
  end subroutine mproj_set
  !
  subroutine mproj_zero_solid(n,mask,u)
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:,0:,0:) :: mask
    real(rp), intent(inout), dimension(0:,0:,0:) :: u
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          if(mask(i,j,k)) u(i,j,k) = 0._rp
        end do
      end do
    end do
  end subroutine mproj_zero_solid
end module mod_ibm_mproj
