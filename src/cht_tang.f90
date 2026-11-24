! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_cht_tang
  !
  ! Conjugate heat transfer (CHT): the tangential face flux (`is_cht_tang`) and the super time-stepping of the
  ! conjugate conduction (`is_cht_rkl2`). Port of the AMPHIBIOUS scheme (CHT_TANGENTIAL with CHT_3D,
  ! CHT_TANG_GUARD, CHT_RKL2).
  !
  ! The tangential face flux. The series face flux of `scal_cht`, k_face*(T_N-T_P)/h, prices a face that the
  ! interface cuts or clips obliquely at the wrong conductivity, which leaves the scheme first order at a curved
  ! interface. The missing part is added face by face:
  !
  !   F = k_face*(T_N-T_P)/h + (k_area-k_face)*s_t
  !
  ! - k_area = f*k_f + (1-f)*k_s, f the fluid area fraction of the true face (the two adjacent cells' quadratic
  !   Taylor models of the level set, averaged; exact roots along one in-plane direction, midpoint columns across);
  ! - s_t = e_d.grad_G(T), the tangential gradient along the face normal e_d, from a least-squares fit of the
  !   conjugate-linear model T = a + sigma1*s1 + sigma2*s2 + q*n/k on the 2x3x3 block of cells that lie in the
  !   3x3x3 stencils of both cells of the face. (n, s1, s2) are coordinates in a local frame of the interface:
  !   the total-least-squares plane through the interface crossings of the grid segments within two local
  !   spacings of the face centre; k is the conductivity of the block cell, from the solver's own mask.
  !
  ! The fit weights depend on the geometry and on k_s/k_f only: s_t is a fixed linear combination of the 18
  ! block values, stored per face as c_m = (k_area-k_face)*w_m. Each face is computed once, canonically from its
  ! lower cell and its axis, and applied to both of its cells with opposite signs (conservative). A face keeps
  ! the series flux where the fit is refused (fewer than 3 crossings, rank-deficient fit, block reaching a domain
  ! boundary along z) and, with `is_cht_tang_guard`, where the term would exceed `cht_tang_cap` times either
  ! cell's series diagonal. The term needs explicit time integration: `is_cht_tang` requires `is_cht_rkl2`.
  !
  ! Storage: the band only. Per face its lower cell, axis and 18 coefficients; per band cell the indices of its
  ! (up to six) faces with a term. The application is a gather loop over the band cells (one writer per cell,
  ! no atomics); both cells of a face evaluate the identical sum.
  !
  ! RKL2. The whole conjugate conduction (series flux, tangential flux and the constant sources) is advanced by
  ! RKL2 super time-stepping (Meyer, Balsara & Aslam, J. Comput. Phys. 257, 2014), Strang-split around the
  ! Runge-Kutta step of the convection: half a step of RKL2, the three Runge-Kutta stages, half a step of RKL2.
  ! The stage count s of each half step follows from dt and a Gershgorin bound G of the operator:
  ! (s**2+s-2)/2 >= G*dt/2. Stable for any dt, second order in time, explicit and conservative at every stage.
  !
  ! References: theory/docs/cht_method_results.tex (the tangential face flux), timestep_cht.tex (RKL2);
  ! AMPHIBIOUS src/modules/preproc-computeCHT3D.cpl (tangFace3, tangFrame3), doc/IMPLEMENTATION.md sec. 9-12.
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: myid,ierr
  implicit none
  private
  public cht_tang_t,cht_tang_setup,cht_tang_apply,cht_tang_dump,cht_rkl2_stages,cht_sdf_read_halo,cht_halo_width,cht_zw
  !
  type cht_tang_t
    integer :: nf = 0                                  ! faces with a tangential term (this rank)
    integer :: nc = 0                                  ! band cells with at least one such face (this rank)
    integer , allocatable, dimension(:,:) :: fijk      ! (4,nf): lower cell i,j,k (local) and axis 1..3
    real(rp), allocatable, dimension(:,:) :: fc        ! (18,nf): (k_area-k_face)*(weight of s_t), block order
    integer , allocatable, dimension(:,:) :: cijk      ! (3,nc): the band cells
    integer , allocatable, dimension(:,:) :: cf        ! (6,nc): their faces -x,+x,-y,+y,-z,+z (0: no term)
    real(rp) :: gersh = 0._rp                          ! Gershgorin bound of the conduction operator (1/time)
  end type cht_tang_t
  !
  integer , parameter :: nfcol  = 20        ! midpoint columns across a face for its fluid fraction
  real(rp), parameter :: frame_rad = 2._rp  ! crossings within this many local spacings fix the frame
  real(rp), parameter :: rank_tol = 1.e-3_rp! refuse a fit whose equilibrated Cholesky diagonal spans more
  real(rp), parameter :: tang_tol = 1.e-8_rp! no term where |k_area-k_face| <= tang_tol*max(k_s/k_f,1)*k_f
  integer , parameter :: band_half = 2      ! band: cells whose (2*band_half+1)**3 patch holds both materials
  !
  contains
  !
  ! ------------------------------------------------------------------------------------------- geometry input
  !
  subroutine cht_halo_width(dl,dzc_g,ng,nh)
    !
    ! halo widths of the level set the set-up needs: the frame window (frame_rad times the largest local
    ! spacing, in cells of the smallest spacing along each axis, plus one), the link fraction of a segment
    ! (two more cells, `frac_q`) and one for the window's last segment
    !
    implicit none
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc_g
    integer , intent(in), dimension(3) :: ng
    integer , intent(out), dimension(3) :: nh
    real(rp) :: hmax,dzmin
    dzmin = minval(dzc_g(1:ng(3)-1))
    hmax  = max(dl(1),dl(2),maxval(dzc_g(1:ng(3)-1)))
    nh(1) = ceiling(frame_rad*hmax/dl(1)) + 4
    nh(2) = ceiling(frame_rad*hmax/dl(2)) + 4
    nh(3) = ceiling(frame_rad*hmax/dzmin) + 4
    nh(:) = max(nh(:),band_half+1)
  end subroutine cht_halo_width
  !
  subroutine cht_sdf_read_halo(fname,ng,lo,n,nh,zc_g,ps,zw)
    !
    ! as `sdf_read` (mod_ibm_sdf), with nh(:) halo layers: ps(1-nh(1):n(1)+nh(1), ...), periodic along x and y,
    ! linear extrapolation beyond the domain along z (the first two layers exactly as `sdf_read`, so that the
    ! link fractions recomputed from it equal the stored ones); zw(1-nh(3):n(3)+nh(3)) the cell-centre z
    ! coordinates, extrapolated linearly beyond the ghost layer
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,n,nh
    real(rp), intent(in), dimension(0:) :: zc_g
    real(rp), intent(out), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: ps
    real(rp), intent(out), dimension(1-nh(3):) :: zw
    real(rp), allocatable, dimension(:,:,:) :: buf
    integer :: fh,kmin,kmax,nk,i,j,k,ig,jg,kg
    integer(kind=MPI_OFFSET_KIND) :: disp
    real(rp) :: w
    !
    call cht_zw(ng,lo,n,nh,zc_g,zw)
    kmin = max(1    ,lo(3)-nh(3))
    kmax = min(ng(3),lo(3)+n(3)-1+nh(3))
    nk = kmax-kmin+1
    allocate(buf(ng(1),ng(2),kmin:kmax))
    call MPI_FILE_OPEN(MPI_COMM_SELF,trim(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    disp = int(kmin-1,MPI_OFFSET_KIND)*ng(1)*ng(2)*storage_size(1._rp)/8
    call MPI_FILE_READ_AT(fh,disp,buf,ng(1)*ng(2)*nk,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    ps(:,:,:) = 0._rp
    do k=1-nh(3),n(3)+nh(3)
      kg = lo(3)+k-1
      if(kg < kmin .or. kg > kmax) cycle
      do j=1-nh(2),n(2)+nh(2)
        jg = modulo(lo(2)+j-2,ng(2))+1
        do i=1-nh(1),n(1)+nh(1)
          ig = modulo(lo(1)+i-2,ng(1))+1
          ps(i,j,k) = buf(ig,jg,kg)
        end do
      end do
    end do
    !
    ! beyond the domain along z: the two layers of `sdf_read`, then linear continuation
    !
    if(lo(3)-nh(3) < 1) then
      k = 1-lo(3)+1                       ! local index of global k = 1
      w = (zc_g(1)-zc_g(0))/(zc_g(2)-zc_g(1))
      ps(:,:,k-1) = ps(:,:,k) + (ps(:,:,k)-ps(:,:,k+1))*w
      do kg=-1,lo(3)-nh(3),-1
        ps(:,:,kg-lo(3)+1) = 2.*ps(:,:,kg-lo(3)+2) - ps(:,:,kg-lo(3)+3)
      end do
    end if
    if(lo(3)+n(3)-1+nh(3) > ng(3)) then
      k = ng(3)-lo(3)+1                   ! local index of global k = ng(3)
      w = (zc_g(ng(3)+1)-zc_g(ng(3)))/(zc_g(ng(3))-zc_g(ng(3)-1))
      ps(:,:,k+1) = ps(:,:,k) + (ps(:,:,k)-ps(:,:,k-1))*w
      do kg=ng(3)+2,lo(3)+n(3)-1+nh(3)
        ps(:,:,kg-lo(3)+1) = 2.*ps(:,:,kg-lo(3)) - ps(:,:,kg-lo(3)-1)
      end do
    end if
  end subroutine cht_sdf_read_halo
  !
  subroutine cht_zw(ng,lo,n,nh,zc_g,zw)
    !
    ! cell-centre z coordinates zw(1-nh(3):n(3)+nh(3)), extrapolated linearly beyond the ghost layers
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,n,nh
    real(rp), intent(in), dimension(0:) :: zc_g
    real(rp), intent(out), dimension(1-nh(3):) :: zw
    integer :: k,kg
    do k=1-nh(3),n(3)+nh(3)
      kg = lo(3)+k-1
      if(kg < 0) then
        zw(k) = zc_g(0) + kg*(zc_g(1)-zc_g(0))
      else if(kg > ng(3)+1) then
        zw(k) = zc_g(ng(3)+1) + (kg-ng(3)-1)*(zc_g(ng(3)+1)-zc_g(ng(3)))
      else
        zw(k) = zc_g(kg)
      end if
    end do
  end subroutine cht_zw
  !
  ! ------------------------------------------------------------------------------------------------- set-up
  !
  subroutine cht_tang_setup(n,lo,ng,nh,dl,dzf,ps,zw,thx,thy,thz,phi,alphaf,alphas,cs_cf, &
                            is_tang,is_guard,cap,iscal,tt)
    !
    ! builds the face list of one CHT scalar (`is_tang`) and the Gershgorin bound of its conduction operator
    ! (always; with is_tang = F only the series flux enters it)
    !
    use mod_ibm_sdf, only: frac_q,fluid_length
    implicit none
    integer , intent(in), dimension(3) :: n,lo,ng,nh
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), intent(in), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: ps
    real(rp), intent(in), dimension(1-nh(3):) :: zw
    real(rp), intent(in), dimension(0:,0:,0:) :: thx,thy,thz,phi
    real(rp), intent(in) :: alphaf,alphas,cs_cf,cap
    logical , intent(in) :: is_tang,is_guard
    integer , intent(in) :: iscal
    type(cht_tang_t), intent(inout) :: tt
    !
    logical , allocatable, dimension(:,:,:) :: band
    integer , allocatable, dimension(:,:,:) :: cidx
    integer , allocatable, dimension(:,:)   :: fijk,cf
    real(rp), allocatable, dimension(:,:)   :: fc
    real(rp), allocatable, dimension(:,:,:) :: rowt
    integer  :: i,j,k,d,iL,jL,kL,iU,jU,kU,nfmax,nf,nc,m,l,ii,jj,kk
    integer  :: nrej_cross,nrej_rank,nrej_bnd,nrej_gain,norphan
    integer  :: cnt(6)
    real(rp) :: fac,w(18),wl,wu,dg_l,dg_u,g,gmax_rej,worst_gain,worst_rs,rs,ra,capi,dia,gersh
    real(rp) :: rv(3)
    logical  :: ok
    character(len=8) :: reason
    !
    ! band (including the halo cells, for the orphan test): the cells whose 5x5x5 patch holds both materials
    !
    allocate(band(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    band(:,:,:) = .false.
    if(is_tang) then
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          band(i,j,k) = .false.
          do kk=k-band_half,k+band_half
            do jj=j-band_half,j+band_half
              do ii=i-band_half,i+band_half
                if((ps(ii,jj,kk) > 0._rp).neqv.(ps(i,j,k) > 0._rp)) band(i,j,k) = .true.
              end do
            end do
          end do
        end do
      end do
    end do
    end if
    !
    nrej_cross = 0; nrej_rank = 0; nrej_bnd = 0; nrej_gain = 0; norphan = 0
    gmax_rej = 0._rp
    nf = 0
    if(is_tang) then
      nfmax = 3*count(band(1:n(1),1:n(2),1:n(3))) + 3*product(n(1:2)) + 3*n(1)*n(3) + 3*n(2)*n(3) + 16
      allocate(fijk(4,nfmax),fc(18,nfmax))
      do kL=0,n(3)
        do jL=0,n(2)
          do iL=0,n(1)
            do d=1,3
              if(d == 1 .and. (jL < 1 .or. kL < 1)) cycle
              if(d == 2 .and. (iL < 1 .or. kL < 1)) cycle
              if(d == 3 .and. (iL < 1 .or. jL < 1)) cycle
              iU = iL; jU = jL; kU = kL
              if(d == 1) iU = iL+1
              if(d == 2) jU = jL+1
              if(d == 3) kU = kL+1
              !
              ! only faces of a local band cell
              !
              if(.not.((is_local(iL,jL,kL) .and. band(iL,jL,kL)) .or. (is_local(iU,jU,kU) .and. band(iU,jU,kU)))) cycle
              call tang_face(iL,jL,kL,d,fac,w,ok,reason)
              if(.not.ok) then
                select case(trim(reason))
                case('cross'); nrej_cross = nrej_cross+1
                case('rank') ; nrej_rank  = nrej_rank +1
                case('bound'); nrej_bnd   = nrej_bnd  +1
                end select
                cycle
              end if
              if(is_guard) then
                wl = cell_width(iL,jL,kL,d); wu = cell_width(iU,jU,kU,d)
                dg_l = step1_diag(iL,jL,kL); dg_u = step1_diag(iU,jU,kU)
                g = abs(fac)*sum(abs(w(:)))*max(1._rp/(wl*max(dg_l,tiny(1._rp))),1._rp/(wu*max(dg_u,tiny(1._rp))))
                if(g > cap) then
                  nrej_gain = nrej_gain+1
                  gmax_rej = max(gmax_rej,g)
                  cycle
                end if
              end if
              !
              ! a face whose coefficients all vanish exactly (a face lying on a flat interface: no tangential
              ! component along its normal) is dropped, so that a grid-aligned wall keeps the series operator as is
              !
              if(all(fac*w(:) == 0._rp)) cycle
              if(.not.(band(iL,jL,kL) .and. band(iU,jU,kU))) norphan = norphan+1
              nf = nf+1
              fijk(:,nf) = [iL,jL,kL,d]
              fc(:,nf) = fac*w(:)
            end do
          end do
        end do
      end do
    end if
    !
    ! band-cell lists: every local cell that owns a face with a term
    !
    allocate(cidx(n(1),n(2),n(3))); cidx(:,:,:) = 0
    nc = 0
    do m=1,nf
      do l=0,1
        i = fijk(1,m); j = fijk(2,m); k = fijk(3,m); d = fijk(4,m)
        if(l == 1) then
          if(d == 1) i = i+1
          if(d == 2) j = j+1
          if(d == 3) k = k+1
        end if
        if(.not.is_local(i,j,k)) cycle
        if(cidx(i,j,k) == 0) then
          nc = nc+1
          cidx(i,j,k) = nc
        end if
      end do
    end do
    allocate(cf(6,max(nc,1))); cf(:,:) = 0
    if(allocated(tt%cijk)) deallocate(tt%cijk,tt%cf,tt%fijk,tt%fc)
    allocate(tt%cijk(3,max(nc,1)),tt%cf(6,max(nc,1)),tt%fijk(4,max(nf,1)),tt%fc(18,max(nf,1)))
    tt%cijk(:,:) = 1; tt%fijk(:,:) = 1; tt%fc(:,:) = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(cidx(i,j,k) > 0) tt%cijk(:,cidx(i,j,k)) = [i,j,k]
        end do
      end do
    end do
    do m=1,nf
      do l=0,1
        i = fijk(1,m); j = fijk(2,m); k = fijk(3,m); d = fijk(4,m)
        if(l == 1) then
          if(d == 1) i = i+1
          if(d == 2) j = j+1
          if(d == 3) k = k+1
        end if
        if(.not.is_local(i,j,k)) cycle
        cf(2*d-1+(1-l),cidx(i,j,k)) = m ! the cell is the lower one (l = 0): its + face; the upper: its - face
      end do
    end do
    tt%nf = nf; tt%nc = nc
    tt%cf(:,:) = cf(:,:)
    if(nf > 0) then
      tt%fijk(:,1:nf) = fijk(:,1:nf)
      tt%fc(:,1:nf) = fc(:,1:nf)
    end if
    !
    ! diagnostics: the tangential row against the series diagonal, the row sums (the fit reproduces a uniform
    ! field: every face's weights sum to zero), and the Gershgorin bound of the conduction operator
    !
    allocate(rowt(n(1),n(2),n(3))); rowt(:,:,:) = 0._rp
    worst_rs = 0._rp
    do m=1,nf
      rs = sum(tt%fc(:,m)); ra = sum(abs(tt%fc(:,m)))
      if(ra > 0._rp) worst_rs = max(worst_rs,abs(rs)/ra)
      do l=0,1
        i = tt%fijk(1,m); j = tt%fijk(2,m); k = tt%fijk(3,m); d = tt%fijk(4,m)
        if(l == 1) then
          if(d == 1) i = i+1
          if(d == 2) j = j+1
          if(d == 3) k = k+1
        end if
        if(.not.is_local(i,j,k)) cycle
        rowt(i,j,k) = rowt(i,j,k) + ra/cell_width(i,j,k,d)
      end do
    end do
    worst_gain = 0._rp
    gersh = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dia  = step1_diag(i,j,k)
          capi = 1._rp/(1._rp + (cs_cf-1._rp)*(1._rp-phi(i,j,k)))
          if(rowt(i,j,k) > 0._rp) worst_gain = max(worst_gain,rowt(i,j,k)/dia)
          gersh = max(gersh,capi*(2._rp*dia + rowt(i,j,k)))
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,gersh,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    tt%gersh = gersh
    if(is_tang) then
      cnt(:) = [nf,nrej_cross,nrej_rank,nrej_bnd,nrej_gain,norphan]
      call MPI_ALLREDUCE(MPI_IN_PLACE,cnt,6,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      rv(:) = [worst_gain,worst_rs,gmax_rej]
      call MPI_ALLREDUCE(MPI_IN_PLACE,rv,3,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        print'(a,i3,a,i10,a,3i8)', ' CHT tangential term, scalar',iscal,': faces with a term',cnt(1), &
              '; refused for: < 3 crossings, rank, block at a z boundary',cnt(2:4)
        if(is_guard) print'(a,f6.2,a,i8,a,es10.3)', '   row-gain guard (cap',cap,'): faces kept on the series flux', &
                                                    cnt(5),', worst refused gain',rv(3)
        print'(a,es10.3,a,es10.3,a,i6)', '   worst |tangential row|/|series diagonal|',rv(1), &
              '; worst relative face row sum',rv(2),'; faces with a term shared with a cell outside the band',cnt(6)
        print'(a)', '   (face counts are summed over the ranks: faces on a rank boundary count twice)'
      end if
      if(cnt(6) > 0) then
        if(myid == 0) print*, 'ERROR: CHT tangential term: faces with a term touch a cell outside the band.'
        call MPI_FINALIZE(ierr)
        error stop
      end if
      if(rv(2) > 1.e-10_rp) then
        if(myid == 0) print*, 'ERROR: CHT tangential term: the face weights do not annihilate constants.'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    end if
    if(myid == 0) print'(a,i3,a,es12.5)', ' CHT scalar',iscal,': Gershgorin bound of the conduction operator',tt%gersh
    deallocate(band,cidx,cf,rowt)
    if(allocated(fijk)) deallocate(fijk,fc)
    !
    contains
    !
    logical function is_local(i,j,k)
      integer, intent(in) :: i,j,k
      is_local = i >= 1 .and. i <= n(1) .and. j >= 1 .and. j <= n(2) .and. k >= 1 .and. k <= n(3)
    end function is_local
    !
    ! coordinates: x and y relative to the face's lower cell (integer offsets times the spacing, so that every
    ! rank and every periodic image computes the same numbers), z global
    !
    real(rp) function zcell(k)
      integer, intent(in) :: k
      zcell = zw(k)
    end function zcell
    real(rp) function zface(k) ! the face between cells k and k+1
      integer, intent(in) :: k
      if(k >= 0 .and. k <= n(3)) then
        zface = zcell(k) + 0.5_rp*dzf(k)
      else
        zface = 0.5_rp*(zcell(k)+zcell(k+1))
      end if
    end function zface
    real(rp) function cell_width(i,j,k,d)
      integer, intent(in) :: i,j,k,d
      if(d == 1) then
        cell_width = dl(1)
      else if(d == 2) then
        cell_width = dl(2)
      else
        cell_width = dzf(k)
      end if
    end function cell_width
    !
    ! fluid fraction of the segment from cell (i,j,k) to its neighbour along d, as `set_cht_geometry_sdf`
    !
    real(rp) function link_th(i,j,k,d)
      integer, intent(in) :: i,j,k,d
      real(rp) :: pc,pn
      pc = ps(i,j,k)
      if(d == 1) then
        pn = ps(i+1,j,k)
      else if(d == 2) then
        pn = ps(i,j+1,k)
      else
        pn = ps(i,j,k+1)
      end if
      link_th = merge(1._rp,0._rp,pc > 0._rp)
      if((pc > 0._rp).eqv.(pn > 0._rp)) return
      if(d == 1) then
        if(pc > 0._rp) then
          link_th = frac_q(pc,pn,ps(i-1,j,k),dl(1),dl(1))
        else
          link_th = frac_q(pn,pc,ps(i+2,j,k),dl(1),dl(1))
        end if
      else if(d == 2) then
        if(pc > 0._rp) then
          link_th = frac_q(pc,pn,ps(i,j-1,k),dl(2),dl(2))
        else
          link_th = frac_q(pn,pc,ps(i,j+2,k),dl(2),dl(2))
        end if
      else
        if(pc > 0._rp) then
          link_th = frac_q(pc,pn,ps(i,j,k-1),zw(k+1)-zw(k),zw(k)-zw(k-1))
        else
          link_th = frac_q(pn,pc,ps(i,j,k+2),zw(k+1)-zw(k),zw(k+2)-zw(k+1))
        end if
      end if
    end function link_th
    !
    ! series face conductivity of the segment from (i,j,k) along d: the stored value where it exists
    ! (identical to `scal_cht`), else recomputed from the same inputs
    !
    real(rp) function kface(i,j,k,d)
      integer, intent(in) :: i,j,k,d
      real(rp) :: th
      if(i >= 0 .and. i <= n(1)+1 .and. j >= 0 .and. j <= n(2)+1 .and. k >= 0 .and. k <= n(3)+1 .and. &
         ((d == 1 .and. i <= n(1)) .or. (d == 2 .and. j <= n(2)) .or. (d == 3 .and. k <= n(3)))) then
        if(d == 1) th = thx(i,j,k)
        if(d == 2) th = thy(i,j,k)
        if(d == 3) th = thz(i,j,k)
      else
        th = link_th(i,j,k,d)
      end if
      kface = 1._rp/(th/alphaf + (1._rp-th)/alphas)
    end function kface
    !
    ! series diagonal of cell (i,j,k): sum over its six arms of k_face/(h*W), as `scal_cht`
    !
    real(rp) function step1_diag(i,j,k)
      integer, intent(in) :: i,j,k
      step1_diag = (kface(i-1,j,k,1)+kface(i,j,k,1))/dl(1)**2 + &
                   (kface(i,j-1,k,2)+kface(i,j,k,2))/dl(2)**2 + &
                   (kface(i,j,k,3)/(zw(k+1)-zw(k)) + kface(i,j,k-1,3)/(zw(k)-zw(k-1)))/dzf(k)
    end function step1_diag
    !
    ! quadratic Taylor model of the level set about the centre of cell (i,j,k) (`sdf_taylor`):
    ! value c, gradient gr(3), Hessian hs(3,3)
    !
    subroutine taylor(i,j,k,c,gr,hs)
      integer , intent(in) :: i,j,k
      real(rp), intent(out) :: c,gr(3),hs(3,3)
      real(rp) :: dm,dp
      c  = ps(i,j,k)
      dm = zw(k)-zw(k-1); dp = zw(k+1)-zw(k)
      gr(1) = (ps(i+1,j,k)-ps(i-1,j,k))/(2._rp*dl(1))
      gr(2) = (ps(i,j+1,k)-ps(i,j-1,k))/(2._rp*dl(2))
      gr(3) = (ps(i,j,k+1)*dm**2 - ps(i,j,k-1)*dp**2 + c*(dp**2-dm**2))/(dm*dp*(dm+dp))
      hs(1,1) = (ps(i+1,j,k)-2._rp*c+ps(i-1,j,k))/dl(1)**2
      hs(2,2) = (ps(i,j+1,k)-2._rp*c+ps(i,j-1,k))/dl(2)**2
      hs(3,3) = 2._rp*(ps(i,j,k+1)*dm + ps(i,j,k-1)*dp - c*(dm+dp))/(dm*dp*(dm+dp))
      hs(1,2) = (ps(i+1,j+1,k)-ps(i+1,j-1,k)-ps(i-1,j+1,k)+ps(i-1,j-1,k))/(4._rp*dl(1)*dl(2))
      hs(1,3) = (ps(i+1,j,k+1)-ps(i+1,j,k-1)-ps(i-1,j,k+1)+ps(i-1,j,k-1))/(2._rp*dl(1)*(dm+dp))
      hs(2,3) = (ps(i,j+1,k+1)-ps(i,j+1,k-1)-ps(i,j-1,k+1)+ps(i,j-1,k-1))/(2._rp*dl(2)*(dm+dp))
      hs(2,1) = hs(1,2); hs(3,1) = hs(1,3); hs(3,2) = hs(2,3)
    end subroutine taylor
    !
    ! fluid area fraction of the face between cell (iL,jL,kL) and its neighbour along d: the average of the two
    ! cells' Taylor models, exact roots along the in-plane direction where the model varies more, nfcol
    ! midpoint columns across it
    !
    real(rp) function face_fluid(iL,jL,kL,d)
      integer, intent(in) :: iL,jL,kL,d
      integer  :: iU,jU,kU,p,q,dex,dmid,ic
      real(rp) :: cL,gL(3),hL(3,3),cU,gU(3),hU(3,3),oL(3),oU(3),half(3),gav(3)
      real(rp) :: xL(3),xU(3),qa,qb,qc,a,f
      iU = iL; jU = jL; kU = kL
      if(d == 1) iU = iL+1
      if(d == 2) jU = jL+1
      if(d == 3) kU = kL+1
      call taylor(iL,jL,kL,cL,gL,hL)
      call taylor(iU,jU,kU,cU,gU,hU)
      ! face centre relative to the two cell centres, and the half-extents of the face
      oL(:) = 0._rp; oU(:) = 0._rp
      half(1) = 0.5_rp*dl(1); half(2) = 0.5_rp*dl(2); half(3) = 0.5_rp*dzf(kL)
      if(d == 1) then
        oL(1) =  0.5_rp*dl(1); oU(1) = -0.5_rp*dl(1)
      else if(d == 2) then
        oL(2) =  0.5_rp*dl(2); oU(2) = -0.5_rp*dl(2)
      else
        oL(3) = zface(kL)-zcell(kL); oU(3) = zface(kL)-zcell(kU)
      end if
      ! the two in-plane axes p < q; exact along the one where the averaged gradient is larger
      p = merge(2,1,d == 1); q = merge(2,3,d == 3)
      gav(:) = 0.5_rp*(gL(:) + matmul(hL,oL) + gU(:) + matmul(hU,oU))
      dex = q; dmid = p
      if(abs(gav(p)) > abs(gav(q))) then
        dex = p; dmid = q
      end if
      f = 0._rp
      do ic=1,nfcol
        a = ((ic-0.5_rp)/nfcol-0.5_rp)*2._rp*half(dmid)
        xL(:) = oL(:); xU(:) = oU(:)
        xL(dmid) = xL(dmid)+a; xU(dmid) = xU(dmid)+a
        ! P(x0 + t e) = P(x0) + t (g + H x0).e + t**2/2 e.H.e, averaged over the two models
        qa = 0.25_rp*(hL(dex,dex)+hU(dex,dex))
        qb = 0.5_rp*(gL(dex) + dot_product(hL(dex,:),xL) + gU(dex) + dot_product(hU(dex,:),xU))
        qc = 0.5_rp*(cL + dot_product(gL,xL) + 0.5_rp*dot_product(xL,matmul(hL,xL)) + &
                     cU + dot_product(gU,xU) + 0.5_rp*dot_product(xU,matmul(hU,xU)))
        f = f + fluid_length(qa,qb,qc,half(dex))/(2._rp*half(dex))
      end do
      face_fluid = f/nfcol
      if(face_fluid < 1.e-9_rp) face_fluid = 0._rp
      if(face_fluid > 1._rp-1.e-9_rp) face_fluid = 1._rp
    end function face_fluid
    !
    ! the interface frame near the face centre: total least squares on the crossings of every grid segment
    ! within frame_rad local spacings (normal: eigenvector of the smallest eigenvalue of the scatter matrix),
    ! anchored at the foot of the perpendicular from the face centre; positions relative to the lower cell
    !
    subroutine frame(iL,jL,kL,fcen,pf,nf_,ok)
      integer , intent(in) :: iL,jL,kL
      real(rp), intent(in) :: fcen(3)
      real(rp), intent(out) :: pf(3),nf_(3)
      logical , intent(out) :: ok
      real(rp) :: hx,hy,hz,rad,s0,s1(3),s2(3,3),bar(3),cv(3,3),vv(3,3),xa(3),xb(3),dq(3),t,n0,nrm
      integer  :: mx,my,mz,iw,jw,kw,dd,imin,a
      ok = .false.
      hx = dl(1); hy = dl(2); hz = zcell(kL+1)-zcell(kL)
      rad = frame_rad*max(hx,hy,hz)
      mx = ceiling(rad/hx)+1; my = ceiling(rad/hy)+1; mz = ceiling(rad/hz)+1
      s0 = 0._rp; s1(:) = 0._rp; s2(:,:) = 0._rp
      do kw=kL-mz,kL+mz
        do jw=jL-my,jL+my
          do iw=iL-mx,iL+mx
            xa(:) = [(iw-iL)*dl(1),(jw-jL)*dl(2),zcell(kw)-zcell(kL)]
            do dd=1,3
              xb(:) = xa(:)
              if(dd == 1) xb(1) = (iw+1-iL)*dl(1)
              if(dd == 2) xb(2) = (jw+1-jL)*dl(2)
              if(dd == 3) xb(3) = zcell(kw+1)-zcell(kL)
              if(dd == 1) then
                if((ps(iw,jw,kw) > 0._rp).eqv.(ps(iw+1,jw,kw) > 0._rp)) cycle
              else if(dd == 2) then
                if((ps(iw,jw,kw) > 0._rp).eqv.(ps(iw,jw+1,kw) > 0._rp)) cycle
              else
                if((ps(iw,jw,kw) > 0._rp).eqv.(ps(iw,jw,kw+1) > 0._rp)) cycle
              end if
              t = link_th(iw,jw,kw,dd)
              if(ps(iw,jw,kw) > 0._rp) then
                dq(:) = xa(:) + t*(xb(:)-xa(:)) - fcen(:)
              else
                dq(:) = xb(:) + t*(xa(:)-xb(:)) - fcen(:)
              end if
              if(sum(dq(:)**2) <= rad*rad) then
                s0 = s0+1._rp
                s1(:) = s1(:)+dq(:)
                do a=1,3
                  s2(:,a) = s2(:,a) + dq(:)*dq(a)
                end do
              end if
            end do
          end do
        end do
      end do
      if(s0 < 3._rp) return
      bar(:) = s1(:)/s0
      do a=1,3
        cv(:,a) = s2(:,a)/s0 - bar(:)*bar(a)
      end do
      call jacobi3(cv,vv)
      imin = 1
      do a=2,3
        if(cv(a,a) < cv(imin,imin)) imin = a
      end do
      nrm = sqrt(sum(vv(:,imin)**2))
      if(nrm <= 0._rp) return
      nf_(:) = vv(:,imin)/nrm
      n0 = -dot_product(bar,nf_)
      pf(:) = fcen(:) - n0*nf_(:)
      ok = .true.
    end subroutine frame
    !
    ! one face, canonically from its lower cell: fac = k_area-k_face and the weights w of s_t on the block
    !
    subroutine tang_face(iL,jL,kL,d,fac,w,ok,reason)
      integer , intent(in) :: iL,jL,kL,d
      real(rp), intent(out) :: fac,w(18)
      logical , intent(out) :: ok
      character(len=8), intent(out) :: reason
      integer  :: m,r,ib,jb,kb,p,q,iref,kg
      real(rp) :: kf,ka,f,fcen(3),pf(3),nv(3),t1(3),t2(3),nt,bx(3),bk,dq(3)
      real(rp) :: bm(18,4),colw(4),mm(4,4),rhs(4),wa(18),wb(18),rk
      logical  :: okf
      ok = .false.; fac = 0._rp; w(:) = 0._rp; reason = 'none'
      kf = kface(iL,jL,kL,d)
      f = face_fluid(iL,jL,kL,d)
      ka = f*alphaf + (1._rp-f)*alphas
      if(abs(ka-kf) <= tang_tol*max(alphas/alphaf,1._rp)*alphaf) return
      ! the block must lie inside the domain along z
      do m=1,18
        call block_offset(m,d,ib,jb,kb)
        kg = lo(3)-1+kL+kb
        if(kg < 1 .or. kg > ng(3)) then
          reason = 'bound'; return
        end if
      end do
      ! face centre, relative to the lower cell centre
      fcen(:) = 0._rp
      if(d == 1) fcen(1) = 0.5_rp*dl(1)
      if(d == 2) fcen(2) = 0.5_rp*dl(2)
      if(d == 3) fcen(3) = zface(kL)-zcell(kL)
      call frame(iL,jL,kL,fcen,pf,nv,okf)
      if(.not.okf) then
        reason = 'cross'; return
      end if
      ! tangent basis: the axis least aligned with n, projected on the plane, and n x t1
      iref = 1
      do m=2,3
        if(abs(nv(m)) < abs(nv(iref))) iref = m
      end do
      do m=1,3
        t1(m) = merge(1._rp,0._rp,m == iref) - nv(iref)*nv(m)
      end do
      nt = sqrt(sum(t1(:)**2)); t1(:) = t1(:)/nt
      t2(1) = nv(2)*t1(3)-nv(3)*t1(2); t2(2) = nv(3)*t1(1)-nv(1)*t1(3); t2(3) = nv(1)*t1(2)-nv(2)*t1(1)
      ! conjugate-linear least squares on the block, columns equilibrated
      do m=1,18
        call block_offset(m,d,ib,jb,kb)
        bx(:) = [ib*dl(1),jb*dl(2),zcell(kL+kb)-zcell(kL)]
        bk = merge(alphaf,alphas,ps(iL+ib,jL+jb,kL+kb) > 0._rp)
        dq(:) = bx(:)-pf(:)
        bm(m,1) = 1._rp
        bm(m,2) = dot_product(dq,t1)
        bm(m,3) = dot_product(dq,t2)
        bm(m,4) = dot_product(dq,nv)/bk
      end do
      do p=1,4
        nt = sum(bm(:,p)**2)
        if(nt <= 0._rp) then
          reason = 'rank'; return
        end if
        colw(p) = 1._rp/sqrt(nt)
        bm(:,p) = bm(:,p)*colw(p)
      end do
      ! the two tangential amplitudes: rows 2 and 3 of (B'B)^-1 B'
      do r=2,3
        do q=1,4
          do p=1,4
            mm(p,q) = sum(bm(:,p)*bm(:,q))
          end do
        end do
        rhs(:) = 0._rp; rhs(r) = 1._rp
        call spd_solve(mm,rhs,okf,rk)
        if(.not.okf .or. rk < rank_tol) then
          reason = 'rank'; return
        end if
        do m=1,18
          if(r == 2) wa(m) = colw(2)*sum(bm(m,:)*rhs(:))
          if(r == 3) wb(m) = colw(3)*sum(bm(m,:)*rhs(:))
        end do
      end do
      ! s_t = (t1.e_d) sigma1 + (t2.e_d) sigma2
      w(:) = t1(d)*wa(:) + t2(d)*wb(:)
      fac = ka-kf
      ok = .true.
    end subroutine tang_face
  end subroutine cht_tang_setup
  !
  pure subroutine block_offset(m,d,ib,jb,kb)
    !
    ! the m-th cell (1..18) of the block of a face along d, relative to the face's lower cell: two cells along d
    ! (a = 0, 1) times the 3x3 of the two other axes (b, c = -1..1, in the order x, y, z with d removed)
    !
    integer, intent(in)  :: m,d
    integer, intent(out) :: ib,jb,kb
    integer :: a,b,c,r
    a = (m-1)/9; r = mod(m-1,9); b = r/3-1; c = mod(r,3)-1
    if(d == 1) then
      ib = a; jb = b; kb = c
    else if(d == 2) then
      ib = b; jb = a; kb = c
    else
      ib = b; jb = c; kb = a
    end if
  end subroutine block_offset
  !
  subroutine spd_solve(a,b,ok,rk)
    !
    ! Cholesky solve of a 4x4 symmetric positive definite system in place; rk = smallest over largest diagonal
    ! of the factor (the rank test of AMPHIBIOUS `spdSolve`)
    !
    real(rp), intent(inout) :: a(4,4),b(4)
    logical , intent(out) :: ok
    real(rp), intent(out) :: rk
    integer  :: i,j,m
    real(rp) :: s,dmin,dmax
    ok = .false.; rk = 0._rp
    dmin = huge(1._rp); dmax = 0._rp
    do i=1,4
      s = a(i,i)
      do m=1,i-1
        s = s-a(i,m)**2
      end do
      if(s <= 0._rp) return
      a(i,i) = sqrt(s)
      dmin = min(dmin,a(i,i)); dmax = max(dmax,a(i,i))
      do j=i+1,4
        s = a(j,i)
        do m=1,i-1
          s = s-a(j,m)*a(i,m)
        end do
        a(j,i) = s/a(i,i)
      end do
    end do
    rk = dmin/dmax
    do i=1,4
      s = b(i)
      do m=1,i-1
        s = s-a(i,m)*b(m)
      end do
      b(i) = s/a(i,i)
    end do
    do i=4,1,-1
      s = b(i)
      do m=i+1,4
        s = s-a(m,i)*b(m)
      end do
      b(i) = s/a(i,i)
    end do
    ok = .true.
  end subroutine spd_solve
  !
  subroutine jacobi3(a,v)
    !
    ! eigenvectors of a symmetric 3x3 matrix by cyclic Jacobi rotations (AMPHIBIOUS `jacobi3`): a is
    ! diagonalised in place, the columns of v are the eigenvectors
    !
    real(rp), intent(inout) :: a(3,3)
    real(rp), intent(out)   :: v(3,3)
    integer , parameter :: jp(3) = [1,1,2], jq(3) = [2,3,3]
    integer  :: sw,ir,p,q,kj
    real(rp) :: off,dia,theta,t,c,s,akp,akq,apq,vkp,vkq
    v(:,:) = 0._rp
    v(1,1) = 1._rp; v(2,2) = 1._rp; v(3,3) = 1._rp
    do sw=1,50
      off = abs(a(1,2))+abs(a(1,3))+abs(a(2,3))
      dia = abs(a(1,1))+abs(a(2,2))+abs(a(3,3))
      if(off <= 1.e-15_rp*dia .or. off == 0._rp) return
      do ir=1,3
        p = jp(ir); q = jq(ir)
        if(a(p,q) /= 0._rp) then
          theta = (a(q,q)-a(p,p))/(2._rp*a(p,q))
          t = sign(1._rp,theta)/(abs(theta)+sqrt(theta*theta+1._rp))
          if(theta == 0._rp) t = 1._rp/(abs(theta)+sqrt(theta*theta+1._rp))
          c = 1._rp/sqrt(t*t+1._rp); s = t*c
          do kj=1,3
            if(kj /= p .and. kj /= q) then
              akp = a(kj,p); akq = a(kj,q)
              a(kj,p) = c*akp-s*akq; a(p,kj) = a(kj,p)
              a(kj,q) = s*akp+c*akq; a(q,kj) = a(kj,q)
            end if
          end do
          apq = a(p,q)
          a(p,p) = a(p,p)-t*apq; a(q,q) = a(q,q)+t*apq; a(p,q) = 0._rp; a(q,p) = 0._rp
          do kj=1,3
            vkp = v(kj,p); vkq = v(kj,q)
            v(kj,p) = c*vkp-s*vkq; v(kj,q) = s*vkp+c*vkq
          end do
        end if
      end do
    end do
  end subroutine jacobi3
  !
  ! ---------------------------------------------------------------------------------------------- application
  !
  subroutine cht_tang_apply(nc,cijk,cf,fijk,fc,dli,dzfi,cs_cf,phi,s,dsdt)
    !
    ! adds the divergence of the tangential face fluxes, over the cell heat capacity, to dsdt:
    ! dsdt(c) += capi * sum over the faces of c of (+-) (sum_m fc(m,f) s(block_m(f)))/W
    ! (+ where c is the face's lower cell, - where it is the upper one; W the cell width along the face normal)
    !
    implicit none
    integer , intent(in) :: nc
    integer , intent(in), dimension(:,:) :: cijk,cf,fijk
    real(rp), intent(in), dimension(:,:) :: fc
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), intent(in) :: cs_cf
    real(rp), intent(in), dimension(0:,0:,0:) :: phi,s
    real(rp), intent(inout), dimension(:,:,:) :: dsdt
    integer  :: ic,i,j,k,l,f,iL,jL,kL,d,m,a,b,c,r,ib,jb,kb
    real(rp) :: acc,st,wi,capi,cs
    cs = cs_cf
    !$acc parallel loop gang vector default(present) private(i,j,k,l,f,iL,jL,kL,d,m,a,b,c,r,ib,jb,kb,acc,st,wi,capi) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared) PRIVATE(i,j,k,l,f,iL,jL,kL,d,m,a,b,c,r,ib,jb,kb,acc,st,wi,capi)
    do ic=1,nc
      i = cijk(1,ic); j = cijk(2,ic); k = cijk(3,ic)
      acc = 0._rp
      !$acc loop seq
      do l=1,6
        f = cf(l,ic)
        if(f == 0) cycle
        iL = fijk(1,f); jL = fijk(2,f); kL = fijk(3,f); d = fijk(4,f)
        st = 0._rp
        !$acc loop seq
        do m=1,18
          a = (m-1)/9; r = mod(m-1,9); b = r/3-1; c = mod(r,3)-1
          if(d == 1) then
            ib = a; jb = b; kb = c
          else if(d == 2) then
            ib = b; jb = a; kb = c
          else
            ib = b; jb = c; kb = a
          end if
          st = st + fc(m,f)*s(iL+ib,jL+jb,kL+kb)
        end do
        if(d == 1) then
          wi = dli(1)
        else if(d == 2) then
          wi = dli(2)
        else
          wi = dzfi(k)
        end if
        if(mod(l,2) == 0) then ! the + face: the cell is the face's lower cell
          acc = acc + st*wi
        else
          acc = acc - st*wi
        end if
      end do
      capi = 1._rp/(1._rp + (cs-1._rp)*(1._rp-phi(i,j,k)))
      dsdt(i,j,k) = dsdt(i,j,k) + acc*capi
    end do
  end subroutine cht_tang_apply
  !
  integer function cht_rkl2_stages(mu) result(s)
    !
    ! the fewest RKL2 stages stable for dt*lambda_max = mu on the negative real axis: (s**2+s-2)/2 >= mu
    !
    implicit none
    real(rp), intent(in) :: mu
    s = 2
    do while(0.5_rp*(s*s+s-2) < mu)
      s = s+1
    end do
  end function cht_rkl2_stages
  !
  subroutine cht_tang_dump(fname,n,tt)
    !
    ! the face list of one scalar, for the order and conservation tests (serial runs): nf, then per face the
    ! lower cell (i,j,k), the axis and the 18 coefficients
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    type(cht_tang_t), intent(in) :: tt
    integer :: iunit,m
    open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='replace')
    write(iunit) int(tt%nf,4),int(n(:),4)
    do m=1,tt%nf
      write(iunit) int(tt%fijk(:,m),4),tt%fc(:,m)
    end do
    close(iunit)
  end subroutine cht_tang_dump
end module mod_cht_tang
