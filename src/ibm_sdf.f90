! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_ibm_sdf
  !
  ! IBM geometry from a signed-distance level set phi (> 0 in the fluid, <= 0 in the solid, |grad(phi)| = 1),
  ! as an alternative to the analytic geometry (`isInbody` and bisection, see `mod_ibm`), selected with
  ! `geometry = 'sdf'` in the `ibm` namelist; the level set is given at the cell centres in a raw binary
  ! file (`sdf_file`, double precision, x fastest, ng(1)*ng(2)*ng(3) values, the format of `sdfp.bin`)
  !
  ! from phi, all geometric quantities of the IBM and of the conjugate heat transfer are local:
  !   - masks: phi <= 0 at the respective (staggered) position
  !   - distance to the wall along a grid line, between a fluid node (phi_f > 0) and a solid neighbour
  !     (phi_s <= 0): root of the quadratic interpolant of phi through the two nodes and the next node
  !     on the fluid side (linear interpolation, theta = phi_f/(phi_f-phi_s), if that node is missing or
  !     on the other side of another interface); the linear estimate is exact for planar walls, but not
  !     consistent for grid lines that graze a curved wall (error ~ curvature*spacing/slope), the
  !     quadratic one is second order also there
  !   - normal: n = grad(phi)/|grad(phi)| (central differences); closest wall point x - phi*n
  !   - fluid volume fraction of a cell cut by the wall: cut of the cell with the zero level of the
  !     local quadratic (Taylor) expansion of phi, phi_c + g.x' + x'.H.x'/2 (gradient g and Hessian H
  !     by central differences, i.e. including the curvature of the wall), evaluated on sub-columns
  !     along z
  ! the level set at the velocity positions is interpolated linearly from the cell centres
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: myid,ierr
  implicit none
  private
  public sdf_read,sdf_staggered,sdf_mask,set_ibm_2nd_sdf,calc_grad_dist_sdf,set_cht_geometry_sdf
  public frac_q,sdf_stag,sdf_point,sdf_taylor,sdf_taylor_eval
  contains
  subroutine sdf_read(fname,ng,lo,n,zc,ps)
    !
    ! reads the cell-centred level set into ps(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2) (two halo layers):
    ! periodic along x and y, linear extrapolation beyond the domain along z
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(out), dimension(-1:,-1:,-1:) :: ps
    real(rp), allocatable, dimension(:,:,:) :: buf
    integer :: fh,kmin,kmax,nk,i,j,k,ig,jg,kg
    integer(kind=MPI_OFFSET_KIND) :: disp
    real(rp) :: w
    logical :: is_file
    !
    inquire(file=trim(fname),exist=is_file)
    if(.not.is_file) then
      if(myid == 0) print*, 'ERROR: level set file `'//trim(fname)//'` not found.'
      call MPI_FINALIZE(ierr)
      error stop
    end if
    if(n(3) < 2) then
      if(myid == 0) print*, 'ERROR: `geometry = sdf` requires at least two cells per rank along z.'
      call MPI_FINALIZE(ierr)
      error stop
    end if
    kmin = max(1    ,lo(3)-2)
    kmax = min(ng(3),lo(3)+n(3)+1)
    nk = kmax-kmin+1
    allocate(buf(ng(1),ng(2),kmin:kmax))
    call MPI_FILE_OPEN(MPI_COMM_SELF,trim(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    disp = int(kmin-1,MPI_OFFSET_KIND)*ng(1)*ng(2)*storage_size(1._rp)/8
    call MPI_FILE_READ_AT(fh,disp,buf,ng(1)*ng(2)*nk,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    do k=-1,n(3)+2
      kg = lo(3)+k-1
      if(kg < kmin .or. kg > kmax) cycle
      do j=-1,n(2)+2
        jg = modulo(lo(2)+j-2,ng(2))+1
        do i=-1,n(1)+2
          ig = modulo(lo(1)+i-2,ng(1))+1
          ps(i,j,k) = buf(ig,jg,kg)
        end do
      end do
    end do
    !
    ! beyond the domain along z: linear extrapolation (uniform spacing assumed for the second layer)
    !
    if(lo(3) == 1) then
      w = (zc(1)-zc(0))/(zc(2)-zc(1))
      ps(:,:, 0) = ps(:,:,1) + (ps(:,:,1)-ps(:,:,2))*w
      ps(:,:,-1) = 2.*ps(:,:,0) - ps(:,:,1)
    end if
    if(lo(3)+n(3)-1 == ng(3)) then
      w = (zc(n(3)+1)-zc(n(3)))/(zc(n(3))-zc(n(3)-1))
      ps(:,:,n(3)+1) = ps(:,:,n(3)) + (ps(:,:,n(3))-ps(:,:,n(3)-1))*w
      ps(:,:,n(3)+2) = 2.*ps(:,:,n(3)+1) - ps(:,:,n(3))
    end if
  end subroutine sdf_read
  !
  subroutine sdf_staggered(n,dzc,dzf,ps,phi_s,phi_u,phi_v,phi_w)
    !
    ! level set at the cell centres (phi_s) and at the u, v, w positions (linear interpolation),
    ! on the arrays with one halo layer (0:n+1)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    real(rp), intent(out), dimension(0:,0:,0:) :: phi_s,phi_u,phi_v,phi_w
    integer :: i,j,k
    real(rp) :: wz
    do k=0,n(3)+1
      wz = 0.5_rp
      if(k <= n(3)) wz = 0.5_rp*dzf(k)/dzc(k) ! zf(k) between zc(k) and zc(k+1)
      do j=0,n(2)+1
        do i=0,n(1)+1
          phi_s(i,j,k) = ps(i,j,k)
          phi_u(i,j,k) = 0.5_rp*(ps(i,j,k)+ps(i+1,j,k))
          phi_v(i,j,k) = 0.5_rp*(ps(i,j,k)+ps(i,j+1,k))
          phi_w(i,j,k) = (1._rp-wz)*ps(i,j,k) + wz*ps(i,j,k+1)
        end do
      end do
    end do
  end subroutine sdf_staggered
  !
  subroutine sdf_mask(phi,mask)
    implicit none
    real(rp), intent(in ), dimension(0:,0:,0:) :: phi
    logical , intent(out), dimension(0:,0:,0:) :: mask
    mask(:,:,:) = phi(:,:,:) <= 0._rp
  end subroutine sdf_mask
  !
  pure function frac(pf,ps) result(th)
    !
    ! fraction of a segment from a fluid node (pf > 0) towards a solid node (ps <= 0) that lies in the
    ! fluid, linear interpolation
    !
    real(rp), intent(in) :: pf,ps
    real(rp) :: th
    th = pf/(pf-ps)
  end function frac
  !
  pure function frac_q(pf,ps,pm,d,dm) result(th)
    !
    ! as `frac`, from the quadratic interpolant through the fluid node (t = 0, pf > 0), the solid node
    ! (t = d, ps <= 0) and the next node on the fluid side (t = -dm, pm); falls back to the linear
    ! estimate if that node is not in the fluid or the quadratic has no root in (0,d]
    !
    real(rp), intent(in) :: pf,ps,pm,d,dm
    real(rp) :: th,d1,d2,a,b,c,disc,q,r1,r2,t
    th = pf/(pf-ps)
    if(pm <= 0._rp) return
    d1 = (ps-pf)/d                       ! Newton form: p(t) = pf + d1*t + d2*t*(t-d)
    d2 = ((pm-pf)/(-dm) - d1)/(-dm-d)
    a = d2; b = d1-d2*d; c = pf
    if(abs(a)*d*d < 1.e-12_rp*max(abs(b)*d,abs(c))) return
    disc = b*b-4._rp*a*c
    if(disc < 0._rp) return
    q = -0.5_rp*(b+sign(sqrt(disc),b))
    r1 = q/a
    r2 = huge(1._rp)
    if(q /= 0._rp) r2 = c/q
    t = huge(1._rp)
    if(r1 > 0._rp .and. r1 <= d) t = r1
    if(r2 > 0._rp .and. r2 <= d) t = min(t,r2)
    if(t <= d) th = t/d
  end function frac_q
  !
  subroutine set_ibm_2nd_sdf(n,dl,dzc,dzf,diz,phi,mask_id,laplacian_id,band_id)
    !
    ! the wall-distance correction of the Laplacian of the 2nd-order IBM (as `set_ibm_2nd`, with the
    ! distance from the level set): for a fluid node with a solid neighbour along a grid line,
    !   lambda = (1/(dl_a*dl_b))*(dl_n/d - 1),  d = theta*dl_n,
    ! dl_n the spacing to the neighbour and dl_a*dl_b the factor of the discrete Laplacian
    ! (dl**2 along x and y; dzc*dzf along z, as in `set_ibm_2nd`/`calc_lambda`)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    integer , intent(in) :: diz
    real(rp), intent(in), dimension(0:,0:,0:) :: phi
    logical , intent(in), dimension(0:,0:,0:) :: mask_id
    real(rp), intent(inout), dimension(0:,0:,0:) :: laplacian_id
    logical , intent(inout), dimension(0:,0:,0:) :: band_id
    real(rp), parameter :: eps = 1.e-10_rp
    integer :: i,j,k
    real(rp) :: th,dzc_l,dzf_l,sp,sm
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(mask_id(i,j,k)) cycle
          if(mask_id(i+1,j,k)) then
            th = max(frac_q(phi(i,j,k),phi(i+1,j,k),phi(i-1,j,k),dl(1),dl(1)),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/dl(1)**2)*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
          if(mask_id(i-1,j,k)) then
            th = max(frac_q(phi(i,j,k),phi(i-1,j,k),phi(i+1,j,k),dl(1),dl(1)),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/dl(1)**2)*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
          if(mask_id(i,j+1,k)) then
            th = max(frac_q(phi(i,j,k),phi(i,j+1,k),phi(i,j-1,k),dl(2),dl(2)),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/dl(2)**2)*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
          if(mask_id(i,j-1,k)) then
            th = max(frac_q(phi(i,j,k),phi(i,j-1,k),phi(i,j+1,k),dl(2),dl(2)),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/dl(2)**2)*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
          !
          ! along z: node spacings to the neighbours above (sp) and below (sm), zc (diz = 0) or zf (diz = 1)
          !
          if(diz == 0) then
            sp = dzc(k); sm = dzc(k-1)
          else
            sp = dzf(k+1); sm = dzf(k)
          end if
          if(mask_id(i,j,k+1)) then
            if(diz == 0) then
              dzc_l = dzc(k); dzf_l = dzf(k)
            else
              dzf_l = dzf(k+1); dzc_l = dzc(k)
            end if
            th = max(frac_q(phi(i,j,k),phi(i,j,k+1),phi(i,j,k-1),sp,sm),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/(dzc_l*dzf_l))*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
          if(mask_id(i,j,k-1)) then
            if(diz == 0) then
              dzc_l = dzc(k-1); dzf_l = dzf(k)
            else
              dzf_l = dzf(k); dzc_l = dzc(k)
            end if
            th = max(frac_q(phi(i,j,k),phi(i,j,k-1),phi(i,j,k+1),sm,sp),eps)
            laplacian_id(i,j,k) = laplacian_id(i,j,k) + (1._rp/(dzc_l*dzf_l))*(1._rp/th-1._rp); band_id(i,j,k) = .true.
          end if
        end do
      end do
    end do
  end subroutine set_ibm_2nd_sdf
  !
  subroutine grad_n(phi,i,j,k,dl,zm,zp,nv)
    !
    ! unit normal from central differences of the level set (pointing into the fluid)
    !
    implicit none
    real(rp), intent(in), dimension(0:,0:,0:) :: phi
    integer , intent(in) :: i,j,k
    real(rp), intent(in) :: dl(3),zm,zp
    real(rp), intent(out) :: nv(3)
    real(rp) :: nn
    nv(1) = (phi(i+1,j,k)-phi(i-1,j,k))/(2._rp*dl(1))
    nv(2) = (phi(i,j+1,k)-phi(i,j-1,k))/(2._rp*dl(2))
    nv(3) = (phi(i,j,k+1)-phi(i,j,k-1))/(zp-zm)
    nn = sqrt(sum(nv(:)**2))
    if(nn > 0._rp) then
      nv(:) = nv(:)/nn
    else
      nv(:) = [0._rp,0._rp,1._rp]
    end if
  end subroutine grad_n
  !
  subroutine calc_grad_dist_sdf(grad_dist_id,lo,n,dl,zc,zf,band_id,dix,diy,diz,phi)
    !
    ! as `calc_grad_dist`, from the level set: for each band point (interior, in the loop order of
    ! `calc_grad_dist`) the closest wall point x - phi*n (columns 1-3), the normal n (4-6), the
    ! indices (7-9) and the wall distance phi (10)
    !
    implicit none
    real(rp), intent(inout), dimension(:,:) :: grad_dist_id
    integer , intent(in), dimension(3) :: lo,n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: zc,zf
    logical , intent(in), dimension(0:,0:,0:) :: band_id
    integer , intent(in) :: dix,diy,diz
    real(rp), intent(in), dimension(0:,0:,0:) :: phi
    integer :: i,j,k,c
    real(rp) :: x,y,z,zm,zp,nv(3)
    grad_dist_id(:,:) = 0._rp
    c = 0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(.not.band_id(i,j,k)) cycle
          c = c+1
          x = (real(lo(1)+i-1,rp)-0.5_rp+0.5_rp*dix)*dl(1)
          y = (real(lo(2)+j-1,rp)-0.5_rp+0.5_rp*diy)*dl(2)
          if(diz /= 1) then
            z = zc(k); zm = zc(k-1); zp = zc(k+1)
          else
            z = zf(k); zm = zf(k-1); zp = zf(k+1)
          end if
          call grad_n(phi,i,j,k,dl,zm,zp,nv)
          grad_dist_id(c,1:3) = [x,y,z] - phi(i,j,k)*nv(:)
          grad_dist_id(c,4:6) = nv(:)
          grad_dist_id(c,7:9) = [real(i,rp),real(j,rp),real(k,rp)]
          grad_dist_id(c,10)  = phi(i,j,k)
        end do
      end do
    end do
  end subroutine calc_grad_dist_sdf
  !
  subroutine set_cht_geometry_sdf(n,dl,zc,dzc,dzf,ps,is_2nd,thx,thy,thz,phi)
    !
    ! as `set_cht_geometry` (see `mod_ibm`), from the level set at the cell centres (ps, two halo layers):
    !   thx(i,j,k): fluid fraction of the segment between the centres of cells i and i+1 (thy, thz alike):
    !               quadratic interpolation of the level set (`frac_q`, with `is_2nd`) or interface on
    !               the cell face (staircase)
    !   phi(i,j,k): fluid volume fraction: cut of the cell with the zero level of the quadratic Taylor
    !               expansion of the level set about the cell centre (`is_2nd`), or 0/1 (staircase)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: zc,dzc,dzf
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    logical , intent(in) :: is_2nd
    real(rp), intent(out), dimension(0:,0:,0:) :: thx,thy,thz,phi
    integer, parameter :: ns = 20 ! sub-columns per direction for the volume fraction
    real(rp), allocatable, dimension(:) :: ze
    integer  :: i,j,k,a,b,ncut
    real(rp) :: pc,pn,dz,hz,rdiag,fsum,xo,yo,g(3),hs(3,3),qa,qb,qc
    real(rp) :: zm,zp,dm,dp
    !
    ! cell-centre z coordinates with two halo layers (linear extrapolation of the spacing)
    !
    allocate(ze(-1:n(3)+2))
    ze(0:n(3)+1) = zc(0:n(3)+1)
    ze(-1) = 2._rp*ze(0)-ze(1)
    ze(n(3)+2) = 2._rp*ze(n(3)+1)-ze(n(3))
    ncut = 0
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          pc = ps(i,j,k)
          thx(i,j,k) = merge(1._rp,0._rp,pc > 0._rp)
          thy(i,j,k) = thx(i,j,k)
          thz(i,j,k) = thx(i,j,k)
          if(i <= n(1)) then
            pn = ps(i+1,j,k)
            if((pc > 0._rp).neqv.(pn > 0._rp)) then
              if(.not.is_2nd) then
                thx(i,j,k) = 0.5_rp
              else if(pc > 0._rp) then
                thx(i,j,k) = frac_q(pc,pn,ps(i-1,j,k),dl(1),dl(1))
              else
                thx(i,j,k) = frac_q(pn,pc,ps(i+2,j,k),dl(1),dl(1))
              end if
            end if
          end if
          if(j <= n(2)) then
            pn = ps(i,j+1,k)
            if((pc > 0._rp).neqv.(pn > 0._rp)) then
              if(.not.is_2nd) then
                thy(i,j,k) = 0.5_rp
              else if(pc > 0._rp) then
                thy(i,j,k) = frac_q(pc,pn,ps(i,j-1,k),dl(2),dl(2))
              else
                thy(i,j,k) = frac_q(pn,pc,ps(i,j+2,k),dl(2),dl(2))
              end if
            end if
          end if
          if(k <= n(3)) then
            pn = ps(i,j,k+1)
            if((pc > 0._rp).neqv.(pn > 0._rp)) then
              if(.not.is_2nd) then
                if(pc > 0._rp) then
                  thz(i,j,k) = 0.5_rp*dzf(k  )/dzc(k)
                else
                  thz(i,j,k) = 0.5_rp*dzf(k+1)/dzc(k)
                end if
              else if(pc > 0._rp) then
                thz(i,j,k) = frac_q(pc,pn,ps(i,j,k-1),ze(k+1)-ze(k),ze(k)-ze(k-1))
              else
                thz(i,j,k) = frac_q(pn,pc,ps(i,j,k+2),ze(k+1)-ze(k),ze(k+2)-ze(k+1))
              end if
            end if
          end if
          !
          ! fluid volume fraction, for cells closer to the wall than half their diagonal
          !
          phi(i,j,k) = merge(1._rp,0._rp,pc > 0._rp)
          if(is_2nd .and. k >= 1 .and. k <= n(3) .and. i >= 1 .and. i <= n(1) .and. j >= 1 .and. j <= n(2)) then
            dz = dzf(k); hz = 0.5_rp*dz
            rdiag = 0.5_rp*sqrt(dl(1)**2+dl(2)**2+dz**2)
            if(abs(pc) < rdiag) then
              !
              ! gradient and Hessian of the level set at the cell centre (central differences)
              !
              zm = ze(k-1); zp = ze(k+1); dm = ze(k)-zm; dp = zp-ze(k)
              g(1) = (ps(i+1,j,k)-ps(i-1,j,k))/(2._rp*dl(1))
              g(2) = (ps(i,j+1,k)-ps(i,j-1,k))/(2._rp*dl(2))
              g(3) = (ps(i,j,k+1)*dm**2 - ps(i,j,k-1)*dp**2 + pc*(dp**2-dm**2))/(dm*dp*(dm+dp))
              hs(1,1) = (ps(i+1,j,k)-2._rp*pc+ps(i-1,j,k))/dl(1)**2
              hs(2,2) = (ps(i,j+1,k)-2._rp*pc+ps(i,j-1,k))/dl(2)**2
              hs(3,3) = 2._rp*(ps(i,j,k+1)*dm + ps(i,j,k-1)*dp - pc*(dm+dp))/(dm*dp*(dm+dp))
              hs(1,2) = (ps(i+1,j+1,k)-ps(i+1,j-1,k)-ps(i-1,j+1,k)+ps(i-1,j-1,k))/(4._rp*dl(1)*dl(2))
              hs(1,3) = (ps(i+1,j,k+1)-ps(i+1,j,k-1)-ps(i-1,j,k+1)+ps(i-1,j,k-1))/(2._rp*dl(1)*(dm+dp))
              hs(2,3) = (ps(i,j+1,k+1)-ps(i,j+1,k-1)-ps(i,j-1,k+1)+ps(i,j-1,k-1))/(2._rp*dl(2)*(dm+dp))
              fsum = 0._rp
              do b=1,ns
                do a=1,ns
                  xo = ((a-0.5_rp)/ns-0.5_rp)*dl(1)
                  yo = ((b-0.5_rp)/ns-0.5_rp)*dl(2)
                  !
                  ! level set along the sub-column: qa*z'**2 + qb*z' + qc, z' in [-hz,hz]
                  !
                  qa = 0.5_rp*hs(3,3)
                  qb = g(3) + hs(1,3)*xo + hs(2,3)*yo
                  qc = pc + g(1)*xo + g(2)*yo + 0.5_rp*(hs(1,1)*xo**2 + 2._rp*hs(1,2)*xo*yo + hs(2,2)*yo**2)
                  fsum = fsum + fluid_length(qa,qb,qc,hz)/dz
                end do
              end do
              phi(i,j,k) = fsum/real(ns**2,rp)
              ncut = ncut+1
            end if
          end if
        end do
      end do
    end do
    if(myid == 0 .and. is_2nd) print*, 'CHT geometry (level set): cells with a sub-cell volume fraction (rank 0): ', ncut
  end subroutine set_cht_geometry_sdf
  !
  function sdf_stag(ps,n,dzc,dzf,c,i,j,k,val) result(ok)
    !
    ! level set at the node (i,j,k) of the grid of variable c (0: cell centres, 1, 2, 3: u, v, w), as in
    ! `sdf_staggered`, from the cell-centred level set with two halo layers; .false. if the node is out of
    ! the range covered by ps
    !
    implicit none
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    integer , intent(in) :: c,i,j,k
    real(rp), intent(out) :: val
    logical :: ok
    integer :: ie(3)
    real(rp) :: wz
    ie(:) = 0
    if(c >= 1) ie(c) = 1
    val = 0._rp
    ok = i >= -1 .and. j >= -1 .and. k >= -1 .and. i+ie(1) <= n(1)+2 .and. j+ie(2) <= n(2)+2 .and. k+ie(3) <= n(3)+2
    if(.not.ok) return
    select case(c)
    case(0)
      val = ps(i,j,k)
    case(1)
      val = 0.5_rp*(ps(i,j,k)+ps(i+1,j,k))
    case(2)
      val = 0.5_rp*(ps(i,j,k)+ps(i,j+1,k))
    case default
      wz = 0.5_rp
      if(k >= 0 .and. k <= n(3)) wz = 0.5_rp*dzf(k)/dzc(k)
      val = (1._rp-wz)*ps(i,j,k) + wz*ps(i,j,k+1)
    end select
  end function sdf_stag
  !
  function sdf_point(ps,lo,n,dl,ze,x,y,z) result(val)
    !
    ! level set at the point (x,y,z) (global coordinates, inside the local block extended by one cell):
    ! trilinear interpolation of the cell-centred level set; ze(-1:n(3)+2) the cell-centre z coordinates
    !
    implicit none
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    integer , intent(in), dimension(3) :: lo,n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(-1:) :: ze
    real(rp), intent(in) :: x,y,z
    real(rp) :: val,fx,fy,wx,wy,wz
    integer :: i0,j0,k0
    fx = x/dl(1) - (lo(1)-1) + 0.5_rp; i0 = max(-1,min(n(1)+1,floor(fx))); wx = fx-i0
    fy = y/dl(2) - (lo(2)-1) + 0.5_rp; j0 = max(-1,min(n(2)+1,floor(fy))); wy = fy-j0
    k0 = -1
    do while(k0 < n(3)+1)
      if(ze(k0+1) > z) exit
      k0 = k0+1
    end do
    wz = (z-ze(k0))/(ze(k0+1)-ze(k0))
    val = (1._rp-wz)*((1._rp-wy)*((1._rp-wx)*ps(i0,j0  ,k0  )+wx*ps(i0+1,j0  ,k0  ))  + &
                            wy *((1._rp-wx)*ps(i0,j0+1,k0  )+wx*ps(i0+1,j0+1,k0  ))) + &
                 wz *((1._rp-wy)*((1._rp-wx)*ps(i0,j0  ,k0+1)+wx*ps(i0+1,j0  ,k0+1))  + &
                            wy *((1._rp-wx)*ps(i0,j0+1,k0+1)+wx*ps(i0+1,j0+1,k0+1)))
  end function sdf_point
  !
  subroutine sdf_taylor(ps,i,j,k,dl,ze,coef)
    !
    ! coefficients of the local quadratic (Taylor) expansion of the level set about the centre of cell
    ! (i,j,k), as in `set_cht_geometry_sdf`: coef = [phi_c, g(1:3), H11, H22, H33, H12, H13, H23]
    !
    implicit none
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    integer , intent(in) :: i,j,k
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(-1:) :: ze
    real(rp), intent(out) :: coef(10)
    real(rp) :: pc,dm,dp
    pc = ps(i,j,k)
    dm = ze(k)-ze(k-1); dp = ze(k+1)-ze(k)
    coef(1)  = pc
    coef(2)  = (ps(i+1,j,k)-ps(i-1,j,k))/(2._rp*dl(1))
    coef(3)  = (ps(i,j+1,k)-ps(i,j-1,k))/(2._rp*dl(2))
    coef(4)  = (ps(i,j,k+1)*dm**2 - ps(i,j,k-1)*dp**2 + pc*(dp**2-dm**2))/(dm*dp*(dm+dp))
    coef(5)  = (ps(i+1,j,k)-2._rp*pc+ps(i-1,j,k))/dl(1)**2
    coef(6)  = (ps(i,j+1,k)-2._rp*pc+ps(i,j-1,k))/dl(2)**2
    coef(7)  = 2._rp*(ps(i,j,k+1)*dm + ps(i,j,k-1)*dp - pc*(dm+dp))/(dm*dp*(dm+dp))
    coef(8)  = (ps(i+1,j+1,k)-ps(i+1,j-1,k)-ps(i-1,j+1,k)+ps(i-1,j-1,k))/(4._rp*dl(1)*dl(2))
    coef(9)  = (ps(i+1,j,k+1)-ps(i+1,j,k-1)-ps(i-1,j,k+1)+ps(i-1,j,k-1))/(2._rp*dl(1)*(dm+dp))
    coef(10) = (ps(i,j+1,k+1)-ps(i,j+1,k-1)-ps(i,j-1,k+1)+ps(i,j-1,k-1))/(2._rp*dl(2)*(dm+dp))
  end subroutine sdf_taylor
  !
  pure function sdf_taylor_eval(coef,xo,yo,zo) result(val)
    !
    ! the expansion of `sdf_taylor` at the offset (xo,yo,zo) from the cell centre
    !
    real(rp), intent(in) :: coef(10),xo,yo,zo
    real(rp) :: val
    val = coef(1) + coef(2)*xo + coef(3)*yo + coef(4)*zo + &
          0.5_rp*(coef(5)*xo**2 + coef(6)*yo**2 + coef(7)*zo**2) + coef(8)*xo*yo + coef(9)*xo*zo + coef(10)*yo*zo
  end function sdf_taylor_eval
  !
  pure function fluid_length(qa,qb,qc,hz) result(len)
    !
    ! length of the part of [-hz,hz] where qa*z**2 + qb*z + qc > 0
    !
    real(rp), intent(in) :: qa,qb,qc,hz
    real(rp) :: len,r(4),zz,disc,q,t1,t2,tmp
    integer  :: nr,m
    r(1) = -hz; nr = 1
    if(abs(qa)*hz*hz > 1.e-12_rp*(abs(qb)*hz+abs(qc))) then
      disc = qb*qb-4._rp*qa*qc
      if(disc > 0._rp) then
        q = -0.5_rp*(qb+sign(sqrt(disc),qb))
        t1 = q/qa; t2 = qc/q
        if(t1 > t2) then; tmp = t1; t1 = t2; t2 = tmp; end if
        if(t1 > -hz .and. t1 < hz) then; nr = nr+1; r(nr) = t1; end if
        if(t2 > -hz .and. t2 < hz) then; nr = nr+1; r(nr) = t2; end if
      end if
    else if(qb /= 0._rp) then
      t1 = -qc/qb
      if(t1 > -hz .and. t1 < hz) then; nr = nr+1; r(nr) = t1; end if
    end if
    nr = nr+1; r(nr) = hz
    len = 0._rp
    do m=1,nr-1
      zz = 0.5_rp*(r(m)+r(m+1))
      if(qa*zz*zz+qb*zz+qc > 0._rp) len = len + (r(m+1)-r(m))
    end do
  end function fluid_length
end module mod_ibm_sdf
