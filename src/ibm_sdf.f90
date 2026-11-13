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
  public sdf_interface_density,set_cht_face_tensor,cht_th_tensor,cht_kdiff,set_cht_face_o2
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
  subroutine set_cht_face_tensor(n,dl,zc,dzf,ps,ax,ay,az,wx,wy,wz,pn,fn)
    !
    ! geometry of the anisotropic face conductivity of the CHT scalars (`cht_kface = 'tensor'`), from the
    ! level set at the cell centres (ps, two halo layers); for the face between cells i and i+1 (ax, wx;
    ! y and z alike, same indexing as thx, thy, thz):
    !   ax: fluid area fraction of the face (the face is cut with the zero level of the quadratic Taylor
    !       expansions of the level set about the two adjacent cell centres, averaged; ns sub-lines, exact
    !       length along each from the quadratic, see `fluid_length`; the sub-lines run along the in-plane
    !       direction with the larger variation of the level set over the face, so that the fraction is
    !       exact for a planar interface parallel to either in-plane direction)
    !   wx: n_x**2, the squared x component of the interface normal (gradient of the same expansions at
    !       the face centre)
    ! with these, the face conductivity is the diagonal element of the effective conductivity tensor of a
    ! locally planar interface, k = w*k_series(th) + (1-w)*k_parallel(a) (`cht_th_tensor`): series
    ! (harmonic) across the interface, area-weighted (arithmetic) along it, both exact for a planar
    ! interface normal or parallel to the face direction. Faces farther from the interface than half the
    ! diagonal of the face get w = 1 (series, i.e. unchanged), a = 0 or 1
    ! pn (optional, `cht_kface = 'tensor_full'`): products of the normal components at the faces, for the
    ! off-diagonal elements (k_series - k_parallel)*n_d*n_e of the tensor: pn(:,:,:,1:2) = n_x*n_y, n_x*n_z
    ! (x faces), pn(:,:,:,3:4) = n_y*n_x, n_y*n_z (y faces), pn(:,:,:,5:6) = n_z*n_x, n_z*n_y (z faces);
    ! zero away from the interface
    ! fn (optional, `cht_kface = 'tensor_o2'`): unit interface normal at the faces, fn(:,:,:,1:3,d) at the
    ! faces along d (zero away from the interface); the faces are then taken up to a full cell diagonal
    ! from the interface, to cover all faces whose flux stencil (see `scal_cht_cross`) is cut
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: zc,dzf
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    real(rp), intent(out), dimension(0:,0:,0:) :: ax,ay,az,wx,wy,wz
    real(rp), intent(out), dimension(0:,0:,0:,:), optional :: pn
    real(rp), intent(out), dimension(0:,0:,0:,:,:), optional :: fn
    integer, parameter :: ns = 20 ! sub-lines per face
    real(rp), allocatable, dimension(:) :: ze
    real(rp) :: c1(10),c2(10),g(3),pf,r,o1,o2,o3,qa,qb,qc,fsum,dz
    integer  :: i,j,k,m,nface
    allocate(ze(-1:n(3)+2))
    ze(0:n(3)+1) = zc(0:n(3)+1)
    ze(-1) = 2._rp*ze(0)-ze(1)
    ze(n(3)+2) = 2._rp*ze(n(3)+1)-ze(n(3))
    nface = 0
    if(present(pn)) pn(:,:,:,:) = 0._rp
    if(present(fn)) fn(:,:,:,:,:) = 0._rp
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          ax(i,j,k) = merge(1._rp,0._rp,ps(i,j,k) > 0._rp); ay(i,j,k) = ax(i,j,k); az(i,j,k) = ax(i,j,k)
          wx(i,j,k) = 1._rp; wy(i,j,k) = 1._rp; wz(i,j,k) = 1._rp
          dz = dzf(k)
          !
          ! x face (i+1/2): sub-lines along z at ns positions in y (or along y at ns positions in z, if the
          ! level set varies more along y over the face); offsets from the centres of i and i+1
          !
          if(i <= n(1)) then
            pf = 0.5_rp*(ps(i,j,k)+ps(i+1,j,k)); r = 0.5_rp*sqrt(dl(2)**2+dz**2)
            if(present(fn)) r = sqrt(dl(1)**2+dl(2)**2+dz**2)
            if(abs(pf) < r .or. ((ps(i,j,k) > 0._rp).neqv.(ps(i+1,j,k) > 0._rp))) then
              call sdf_taylor(ps,i,j,k,dl,ze,c1); call sdf_taylor(ps,i+1,j,k,dl,ze,c2)
              o1 = 0.5_rp*dl(1); fsum = 0._rp
              g = 0.5_rp*(grad_at(c1,o1,0._rp,0._rp) + grad_at(c2,-o1,0._rp,0._rp))
              do m=1,ns
                if(abs(g(2))*dl(2) <= abs(g(3))*dz) then
                  o2 = ((m-0.5_rp)/ns-0.5_rp)*dl(2)
                  ! level set along z' in [-dz/2,dz/2] at (x_f, y_j+o2): average of both expansions
                  qa = 0.25_rp*(c1(7)+c2(7))
                  qb = 0.5_rp*((c1(4)+c1(9)*o1+c1(10)*o2) + (c2(4)-c2(9)*o1+c2(10)*o2))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,o1,o2,0._rp) + sdf_taylor_eval(c2,-o1,o2,0._rp))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dz)/dz
                else
                  o3 = ((m-0.5_rp)/ns-0.5_rp)*dz
                  ! level set along y' in [-dy/2,dy/2] at (x_f, z_k+o3)
                  qa = 0.25_rp*(c1(6)+c2(6))
                  qb = 0.5_rp*((c1(3)+c1(8)*o1+c1(10)*o3) + (c2(3)-c2(8)*o1+c2(10)*o3))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,o1,0._rp,o3) + sdf_taylor_eval(c2,-o1,0._rp,o3))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dl(2))/dl(2)
                end if
              end do
              ax(i,j,k) = fsum/ns
              if(sum(g**2) > 0._rp) wx(i,j,k) = g(1)**2/sum(g**2)
              if(present(pn) .and. sum(g**2) > 0._rp) pn(i,j,k,1:2) = g(1)*[g(2),g(3)]/sum(g**2)
              if(present(fn) .and. sum(g**2) > 0._rp) fn(i,j,k,:,1) = g(:)/sqrt(sum(g**2))
              nface = nface+1
            end if
          end if
          !
          ! y face (j+1/2): sub-lines along z at ns positions in x (or along x at ns positions in z)
          !
          if(j <= n(2)) then
            pf = 0.5_rp*(ps(i,j,k)+ps(i,j+1,k)); r = 0.5_rp*sqrt(dl(1)**2+dz**2)
            if(present(fn)) r = sqrt(dl(1)**2+dl(2)**2+dz**2)
            if(abs(pf) < r .or. ((ps(i,j,k) > 0._rp).neqv.(ps(i,j+1,k) > 0._rp))) then
              call sdf_taylor(ps,i,j,k,dl,ze,c1); call sdf_taylor(ps,i,j+1,k,dl,ze,c2)
              o2 = 0.5_rp*dl(2); fsum = 0._rp
              g = 0.5_rp*(grad_at(c1,0._rp,o2,0._rp) + grad_at(c2,0._rp,-o2,0._rp))
              do m=1,ns
                if(abs(g(1))*dl(1) <= abs(g(3))*dz) then
                  o1 = ((m-0.5_rp)/ns-0.5_rp)*dl(1)
                  qa = 0.25_rp*(c1(7)+c2(7))
                  qb = 0.5_rp*((c1(4)+c1(9)*o1+c1(10)*o2) + (c2(4)+c2(9)*o1-c2(10)*o2))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,o1,o2,0._rp) + sdf_taylor_eval(c2,o1,-o2,0._rp))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dz)/dz
                else
                  o3 = ((m-0.5_rp)/ns-0.5_rp)*dz
                  ! level set along x' in [-dx/2,dx/2] at (y_f, z_k+o3)
                  qa = 0.25_rp*(c1(5)+c2(5))
                  qb = 0.5_rp*((c1(2)+c1(8)*o2+c1(9)*o3) + (c2(2)-c2(8)*o2+c2(9)*o3))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,0._rp,o2,o3) + sdf_taylor_eval(c2,0._rp,-o2,o3))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dl(1))/dl(1)
                end if
              end do
              ay(i,j,k) = fsum/ns
              if(sum(g**2) > 0._rp) wy(i,j,k) = g(2)**2/sum(g**2)
              if(present(pn) .and. sum(g**2) > 0._rp) pn(i,j,k,3:4) = g(2)*[g(1),g(3)]/sum(g**2)
              if(present(fn) .and. sum(g**2) > 0._rp) fn(i,j,k,:,2) = g(:)/sqrt(sum(g**2))
              nface = nface+1
            end if
          end if
          !
          ! z face (k+1/2, at zf(k)): sub-lines along x at ns positions in y (or along y at ns positions in x)
          !
          if(k <= n(3)) then
            pf = 0.5_rp*(ps(i,j,k)+ps(i,j,k+1)); r = 0.5_rp*sqrt(dl(1)**2+dl(2)**2)
            if(present(fn)) r = sqrt(dl(1)**2+dl(2)**2+dz**2)
            if(abs(pf) < r .or. ((ps(i,j,k) > 0._rp).neqv.(ps(i,j,k+1) > 0._rp))) then
              call sdf_taylor(ps,i,j,k,dl,ze,c1); call sdf_taylor(ps,i,j,k+1,dl,ze,c2)
              o1 = 0.5_rp*dzf(k); o3 = 0.5_rp*dzf(k+1) ! face offsets from the centres of k (up) and k+1 (down)
              fsum = 0._rp
              g = 0.5_rp*(grad_at(c1,0._rp,0._rp,o1) + grad_at(c2,0._rp,0._rp,-o3))
              do m=1,ns
                o2 = ((m-0.5_rp)/ns-0.5_rp)*merge(dl(2),dl(1),abs(g(2))*dl(2) <= abs(g(1))*dl(1))
                if(abs(g(2))*dl(2) <= abs(g(1))*dl(1)) then
                  ! level set along x' in [-dx/2,dx/2] at (y_j+o2, z_f)
                  qa = 0.25_rp*(c1(5)+c2(5))
                  qb = 0.5_rp*((c1(2)+c1(8)*o2+c1(9)*o1) + (c2(2)+c2(8)*o2-c2(9)*o3))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,0._rp,o2,o1) + sdf_taylor_eval(c2,0._rp,o2,-o3))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dl(1))/dl(1)
                else
                  ! level set along y' in [-dy/2,dy/2] at (x_i+o2, z_f)
                  qa = 0.25_rp*(c1(6)+c2(6))
                  qb = 0.5_rp*((c1(3)+c1(8)*o2+c1(10)*o1) + (c2(3)+c2(8)*o2-c2(10)*o3))
                  qc = 0.5_rp*(sdf_taylor_eval(c1,o2,0._rp,o1) + sdf_taylor_eval(c2,o2,0._rp,-o3))
                  fsum = fsum + fluid_length(qa,qb,qc,0.5_rp*dl(2))/dl(2)
                end if
              end do
              az(i,j,k) = fsum/ns
              if(sum(g**2) > 0._rp) wz(i,j,k) = g(3)**2/sum(g**2)
              if(present(pn) .and. sum(g**2) > 0._rp) pn(i,j,k,5:6) = g(3)*[g(1),g(2)]/sum(g**2)
              if(present(fn) .and. sum(g**2) > 0._rp) fn(i,j,k,:,3) = g(:)/sqrt(sum(g**2))
              nface = nface+1
            end if
          end if
        end do
      end do
    end do
    if(myid == 0) print*, 'CHT tensor face conductivity: faces near the interface (rank 0): ', nface
  contains
    pure function grad_at(c,xo,yo,zo) result(gg)
      real(rp), intent(in) :: c(10),xo,yo,zo
      real(rp) :: gg(3)
      gg(1) = c(2) + c(5)*xo + c(8)*yo + c(9)*zo
      gg(2) = c(3) + c(8)*xo + c(6)*yo + c(10)*zo
      gg(3) = c(4) + c(9)*xo + c(10)*yo + c(7)*zo
    end function grad_at
  end subroutine set_cht_face_tensor
  !
  pure function cht_th_tensor(th,a,w,rf,rs) result(the)
    !
    ! equivalent segment fraction of the tensor face conductivity k = w*k_s(th) + (1-w)*k_p(a),
    ! k_s(th) = 1/(th*rf + (1-th)*rs), k_p(a) = a/rf + (1-a)/rs (rf, rs: inverse fluid/solid
    ! conductivity): the value with k_s(the) = k, so that all CHT kernels keep the series form;
    ! unchanged (the = th) for w = 1 or equal conductivities
    !
    real(rp), intent(in) :: th,a,w,rf,rs
    real(rp) :: the,k
    the = th
    if(w >= 1._rp .or. rf == rs) return
    k = w/(th*rf + (1._rp-th)*rs) + (1._rp-w)*(a/rf + (1._rp-a)/rs)
    the = min(1._rp,max(0._rp,(1._rp/k - rs)/(rf - rs)))
  end function cht_th_tensor
  !
  subroutine set_cht_face_o2(n,dzc,is_bz,thx,thy,thz,ax,ay,az,fn,rf,rs,thxo,thyo,thzo,kc)
    !
    ! second-order face conductivities of a CHT scalar (`cht_kface = 'tensor_o2'`). Near the interface the
    ! temperature is taken as piecewise linear with the jump conditions of a locally planar interface
    ! (continuous temperature and tangential gradient, continuous normal heat flux):
    !
    !   T = T0 + g.(x-x0) + (r-1)*(n.g)*min(psi,0),   r = k_f/k_s,
    !
    ! with g the gradient on the fluid side, n the unit normal (into the fluid) and psi the signed distance.
    ! The difference quotients of the face-flux stencil (normal: the two adjacent cells; tangential: the
    ! four surrounding cells, as in `scal_cht_cross`) of this field are G_d = g_d + (r-1)*(n.g)*n_d*s_d,
    ! s_d: solid fraction of the stencil segments along d (from thx, thy, thz); hence n.g = (n.G)/sigma,
    ! sigma = 1 + (r-1)*sum_d(n_d**2*s_d), and the exact flux through, e.g., an x face of fluid area
    ! fraction a, a*k_f*g_x + (1-a)*k_s*(g_x + (r-1)*(n.g)*n_x), becomes
    !
    !   F_x = K_xx*G_x + K_xy*G_y + K_xz*G_z,   K_xx = k_a + c*n_x, K_xy = c*n_y, K_xz = c*n_z,
    !   k_a = a*k_f + (1-a)*k_s,   c = n_x*((1-a)*(k_f-k_s) - k_a*(r-1)*s_x)/sigma.
    !
    ! The flux is exact for a piecewise linear temperature at a planar interface of any orientation (the
    ! series/tensor closures have an O(1) flux error there, i.e. first-order temperatures); it reduces to
    ! the series conductivity for a grid-aligned interface and to k_f or k_s where the stencil is not cut.
    ! Output, in the form used by the CHT kernels: the diagonal element as equivalent segment fraction
    ! (thxo, thyo, thzo: k = 1/(th*rf + (1-th)*rs), rf, rs: inverse fluid/solid conductivity) and the
    ! off-diagonal elements kc(:,:,:,1:6) (as `kcr` of `tensor_full`: 1,2 = K_xy, K_xz at the x faces,
    ! 3,4 = K_yx, K_yz, 5,6 = K_zx, K_zy). Faces on the domain boundaries along z: diagonal element only (no
    ! tangential flux through the boundary, as for `tensor_full`).
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dzc
    logical , intent(in), dimension(0:1) :: is_bz ! subdomain at the lower/upper domain boundary along z
    real(rp), intent(in), dimension(0:,0:,0:) :: thx,thy,thz,ax,ay,az
    real(rp), intent(in), dimension(0:,0:,0:,:,:) :: fn
    real(rp), intent(in) :: rf,rs
    real(rp), intent(out), dimension(0:,0:,0:) :: thxo,thyo,thzo
    real(rp), intent(out), dimension(0:,0:,0:,:) :: kc
    real(rp) :: s(3),ko(3),kd,w
    integer  :: i,j,k,nface,nclip
    thxo(:,:,:) = thx(:,:,:); thyo(:,:,:) = thy(:,:,:); thzo(:,:,:) = thz(:,:,:)
    kc(:,:,:,:) = 0._rp
    if(abs(rf-rs) <= 1.e-12_rp*rf) return
    nface = 0; nclip = 0
    do k=1,n(3)
      w = dzc(k-1)/(dzc(k-1)+dzc(k)) ! weight of the segment k-1 -> k in k-1 -> k+1
      do j=1,n(2)
        do i=0,n(1) ! x faces
          if(sum(fn(i,j,k,:,1)**2) == 0._rp) cycle
          s(1) = 1._rp - thx(i,j,k)
          s(2) = 1._rp - 0.25_rp*(thy(i,j-1,k)+thy(i+1,j-1,k)+thy(i,j,k)+thy(i+1,j,k))
          s(3) = 1._rp - 0.5_rp*(w*(thz(i,j,k-1)+thz(i+1,j,k-1)) + (1._rp-w)*(thz(i,j,k)+thz(i+1,j,k)))
          call k_o2(ax(i,j,k),s,fn(i,j,k,:,1),1,kd,ko)
          thxo(i,j,k) = (1._rp/kd - rs)/(rf - rs)
          kc(i,j,k,1:2) = ko(2:3)
        end do
      end do
      do j=0,n(2) ! y faces
        do i=1,n(1)
          if(sum(fn(i,j,k,:,2)**2) == 0._rp) cycle
          s(2) = 1._rp - thy(i,j,k)
          s(1) = 1._rp - 0.25_rp*(thx(i-1,j,k)+thx(i-1,j+1,k)+thx(i,j,k)+thx(i,j+1,k))
          s(3) = 1._rp - 0.5_rp*(w*(thz(i,j,k-1)+thz(i,j+1,k-1)) + (1._rp-w)*(thz(i,j,k)+thz(i,j+1,k)))
          call k_o2(ay(i,j,k),s,fn(i,j,k,:,2),2,kd,ko)
          thyo(i,j,k) = (1._rp/kd - rs)/(rf - rs)
          kc(i,j,k,3:4) = ko([1,3])
        end do
      end do
    end do
    do k=0,n(3) ! z faces
      do j=1,n(2)
        do i=1,n(1)
          if(sum(fn(i,j,k,:,3)**2) == 0._rp) cycle
          s(3) = 1._rp - thz(i,j,k)
          s(1) = 1._rp - 0.25_rp*(thx(i-1,j,k)+thx(i-1,j,k+1)+thx(i,j,k)+thx(i,j,k+1))
          s(2) = 1._rp - 0.25_rp*(thy(i,j-1,k)+thy(i,j-1,k+1)+thy(i,j,k)+thy(i,j,k+1))
          call k_o2(az(i,j,k),s,fn(i,j,k,:,3),3,kd,ko)
          thzo(i,j,k) = (1._rp/kd - rs)/(rf - rs)
          kc(i,j,k,5:6) = ko(1:2)
          if((k == 0 .and. is_bz(0)) .or. (k == n(3) .and. is_bz(1))) kc(i,j,k,5:6) = 0._rp
        end do
      end do
    end do
    if(myid == 0) print*, 'CHT second-order face conductivity: faces with a cut stencil, with a clipped ', &
                          'diagonal element (rank 0): ',nface,nclip
  contains
    subroutine k_o2(a,s,nn,id,kd,ko)
      !
      ! diagonal (kd) and off-diagonal (ko, entry id unused) conductivities of a face along id
      !
      real(rp), intent(in)  :: a,s(3),nn(3)
      integer , intent(in)  :: id
      real(rp), intent(out) :: kd,ko(3)
      real(rp) :: kf,ks,r,ka,sig,c
      kf = 1._rp/rf; ks = 1._rp/rs; r = kf/ks
      ka = a*kf + (1._rp-a)*ks
      sig = 1._rp + (r-1._rp)*sum(nn(:)**2*s(:))
      c = nn(id)*((1._rp-a)*(kf-ks) - ka*(r-1._rp)*s(id))/sig
      ko(:) = c*nn(:)
      kd = ka + ko(id)
      if(any(abs(ko(:)) > 1.e-12_rp*max(kf,ks))) nface = nface+1
      !
      ! the diagonal element stays between k_f and k_s for a consistent planar geometry; enforced
      !
      if(kd < min(kf,ks)*(1._rp-1.e-9_rp) .or. kd > max(kf,ks)*(1._rp+1.e-9_rp)) nclip = nclip+1
      kd = min(max(kd,min(kf,ks)),max(kf,ks))
    end subroutine k_o2
  end subroutine set_cht_face_o2
  !
  pure function cht_kdiff(th,a,rf,rs) result(dk)
    !
    ! k_series(th) - k_parallel(a), the factor of n_d*n_e in the off-diagonal elements of the face
    ! conductivity tensor (`cht_kface = 'tensor_full'`); <= 0 for a face cut in proportion
    !
    real(rp), intent(in) :: th,a,rf,rs
    real(rp) :: dk
    dk = 1._rp/(th*rf + (1._rp-th)*rs) - (a/rf + (1._rp-a)/rs)
  end function cht_kdiff
  !
  subroutine sdf_interface_density(n,dl,zc,dzf,ps,eps,adens,area)
    !
    ! interface area per unit volume at the cell centres from the level set (ps, two halo layers):
    !
    !   adens = delta_h(ps)*|grad ps|,   delta_h(d) = (1 + cos(pi*d/h))/(2*h) for |d| < h, 0 otherwise,
    !
    ! h = eps*max(dx,dy,dz(k)); the volume integral of adens is the interface area (exact for a plane with
    ! a smooth, resolved delta); `area` returns that integral over the whole domain
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: zc,dzf
    real(rp), intent(in), dimension(-1:,-1:,-1:) :: ps
    real(rp), intent(in) :: eps
    real(rp), intent(out), dimension(0:,0:,0:) :: adens
    real(rp), intent(out) :: area
    real(rp), parameter :: pi = acos(-1._rp)
    integer  :: i,j,k
    real(rp) :: d,h,gx,gy,gz
    adens(:,:,:) = 0._rp
    area = 0._rp
    do k=1,n(3)
      h = eps*max(dl(1),dl(2),dzf(k))
      do j=1,n(2)
        do i=1,n(1)
          d = ps(i,j,k)
          if(abs(d) < h) then
            gx = (ps(i+1,j,k)-ps(i-1,j,k))/(2._rp*dl(1))
            gy = (ps(i,j+1,k)-ps(i,j-1,k))/(2._rp*dl(2))
            gz = (ps(i,j,k+1)-ps(i,j,k-1))/(zc(k+1)-zc(k-1))
            adens(i,j,k) = (1._rp + cos(pi*d/h))/(2._rp*h)*sqrt(gx**2 + gy**2 + gz**2)
            area = area + adens(i,j,k)*dl(1)*dl(2)*dzf(k)
          end if
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,area,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine sdf_interface_density
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
