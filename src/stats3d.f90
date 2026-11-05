! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_stats3d
  !
  ! online 3D time averages (same quantities and files as cans_world):
  ! <u>, <v>, <w> interpolated to the cell centres, <uu>, <vv>, <ww>, <uv>, <uw>, <vw>, <p>, and per scalar <s>, <s^2>,
  ! <u s>, <w s> (u, w at the cell centres; files us_<nnn>, ws_<nnn>);
  ! mean_n = mean_{n-1} + (x_n - mean_{n-1})/n, every `iaverage` steps for time >= `average_start_time`.
  ! output at every checkpoint: data/<name>_mean.bin (visualization), data/<name>_mean_res.bin (restart),
  ! data/count_avg.out (number of samples)
  !
  use mpi
  use mod_types
  use mod_load  , only: load_one
  use mod_output, only: write_visu_3d
  implicit none
  private
  public stats3d_init,stats3d_update,stats3d_save
  integer, parameter :: nvel = 10
  character(len=3), parameter, dimension(nvel) :: cvel = ['vex','vey','vez','u2 ','v2 ','w2 ','uv ','uw ','vw ','p  ']
  character(len=17), parameter, dimension(nvel) :: dvel = ['Velocity_Mean_X  ','Velocity_Mean_Y  ','Velocity_Mean_Z  ', &
                                                          'Velocity_Mean_U2 ','Velocity_Mean_V2 ','Velocity_Mean_W2 ', &
                                                          'Velocity_Mean_UV ','Velocity_Mean_UW ','Velocity_Mean_VW ', &
                                                          'Pressure_Mean    ']
  real(rp), allocatable, dimension(:,:,:,:) :: vm,sm,s2m,usm,wsm
  integer :: count_avg = 0
  contains
  subroutine stats3d_init(n,nscal,is_restart,datadir,ng,lo,hi,myid)
    !
    ! allocates the averages; on restart, continues from data/<name>_mean_res.bin and data/count_avg.out if present
    ! (a missing <s^2> file restarts it from zero, consistent only if count_avg is reset to 0)
    !
    integer , intent(in) :: n(3),nscal,ng(3),lo(3),hi(3),myid
    logical , intent(in) :: is_restart
    character(len=*), intent(in) :: datadir
    character(len=3) :: scalnum
    logical :: is_file
    integer :: m,is,iunit,ierr
    allocate(vm (0:n(1)+1,0:n(2)+1,0:n(3)+1,nvel), &
             sm (0:n(1)+1,0:n(2)+1,0:n(3)+1,max(nscal,1)), &
             s2m(0:n(1)+1,0:n(2)+1,0:n(3)+1,max(nscal,1)), &
             usm(0:n(1)+1,0:n(2)+1,0:n(3)+1,max(nscal,1)), &
             wsm(0:n(1)+1,0:n(2)+1,0:n(3)+1,max(nscal,1)))
    vm(:,:,:,:) = 0._rp; sm(:,:,:,:) = 0._rp; s2m(:,:,:,:) = 0._rp; usm(:,:,:,:) = 0._rp; wsm(:,:,:,:) = 0._rp
    count_avg = 0
    inquire(file=trim(datadir)//'count_avg.out',exist=is_file)
    if(is_restart.and.is_file) then
      if(myid == 0) then
        open(newunit=iunit,file=trim(datadir)//'count_avg.out',status='old',action='read')
        read(iunit,*) count_avg
        close(iunit)
      end if
      call MPI_BCAST(count_avg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      do m=1,nvel
        call load_res(trim(datadir)//trim(cvel(m))//'_mean_res.bin',vm(:,:,:,m))
      end do
      do is=1,nscal
        write(scalnum,'(i3.3)') is
        call load_res(trim(datadir)//'s_' //scalnum//'_mean_res.bin',sm (:,:,:,is))
        call load_res(trim(datadir)//'s2_'//scalnum//'_mean_res.bin',s2m(:,:,:,is))
        call load_res(trim(datadir)//'us_'//scalnum//'_mean_res.bin',usm(:,:,:,is))
        call load_res(trim(datadir)//'ws_'//scalnum//'_mean_res.bin',wsm(:,:,:,is))
      end do
    end if
    if(myid == 0) print*, '*** 3D time averages: starting at sample count = ', count_avg, ' ***'
    !$acc enter data copyin(vm,sm,s2m,usm,wsm)
    contains
    subroutine load_res(fname,a)
      character(len=*), intent(in) :: fname
      real(rp), intent(inout), dimension(0:,0:,0:) :: a
      real(rp) :: tdum
      integer  :: idum
      inquire(file=fname,exist=is_file)
      if(is_file) then
        call load_one('r',fname,MPI_COMM_WORLD,ng,[1,1,1],lo,hi,a,'mean',tdum,idum)
      else
        if(myid == 0) print*, 'Warning: ', fname, ' not found, this average restarts from zero.'
      end if
    end subroutine load_res
  end subroutine stats3d_init
  !
  subroutine stats3d_update(n,nscal,u,v,w,p,s)
    use mod_scal, only: scalar
    integer , intent(in) :: n(3),nscal
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w,p
    type(scalar), intent(in), dimension(:) :: s
    real(rp) :: inv
    integer :: is
    count_avg = count_avg + 1
    inv = 1._rp/count_avg
    call avg_vel(n,inv,u,v,w,p,vm)
    do is=1,nscal
      call avg_scal(n,inv,u,w,s(is)%val,sm(:,:,:,is),s2m(:,:,:,is),usm(:,:,:,is),wsm(:,:,:,is))
    end do
  end subroutine stats3d_update
  !
  subroutine avg_vel(n,inv,u,v,w,p,a)
    integer , intent(in) :: n(3)
    real(rp), intent(in) :: inv
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,p
    real(rp), intent(inout), dimension(0:,0:,0:,:) :: a
    real(rp) :: ui,vi,wi
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) private(ui,vi,wi) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) PRIVATE(ui,vi,wi)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ui = 0.5_rp*(u(i-1,j,k)+u(i,j,k))
          vi = 0.5_rp*(v(i,j-1,k)+v(i,j,k))
          wi = 0.5_rp*(w(i,j,k-1)+w(i,j,k))
          a(i,j,k, 1) = a(i,j,k, 1) + (ui         - a(i,j,k, 1))*inv
          a(i,j,k, 2) = a(i,j,k, 2) + (vi         - a(i,j,k, 2))*inv
          a(i,j,k, 3) = a(i,j,k, 3) + (wi         - a(i,j,k, 3))*inv
          a(i,j,k, 4) = a(i,j,k, 4) + (ui*ui      - a(i,j,k, 4))*inv
          a(i,j,k, 5) = a(i,j,k, 5) + (vi*vi      - a(i,j,k, 5))*inv
          a(i,j,k, 6) = a(i,j,k, 6) + (wi*wi      - a(i,j,k, 6))*inv
          a(i,j,k, 7) = a(i,j,k, 7) + (ui*vi      - a(i,j,k, 7))*inv
          a(i,j,k, 8) = a(i,j,k, 8) + (ui*wi      - a(i,j,k, 8))*inv
          a(i,j,k, 9) = a(i,j,k, 9) + (vi*wi      - a(i,j,k, 9))*inv
          a(i,j,k,10) = a(i,j,k,10) + (p(i,j,k)   - a(i,j,k,10))*inv
        end do
      end do
    end do
  end subroutine avg_vel
  !
  subroutine avg_scal(n,inv,u,w,s,a,a2,aus,aws)
    integer , intent(in) :: n(3)
    real(rp), intent(in) :: inv
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,w,s
    real(rp), intent(inout), dimension(0:,0:,0:) :: a,a2,aus,aws
    real(rp) :: ui,wi
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) private(ui,wi) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) PRIVATE(ui,wi)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          a (i,j,k) = a (i,j,k) + (s(i,j,k)   - a (i,j,k))*inv
          a2(i,j,k) = a2(i,j,k) + (s(i,j,k)**2 - a2(i,j,k))*inv
          ui = 0.5_rp*(u(i-1,j,k)+u(i,j,k))
          wi = 0.5_rp*(w(i,j,k-1)+w(i,j,k))
          aus(i,j,k) = aus(i,j,k) + (ui*s(i,j,k) - aus(i,j,k))*inv
          aws(i,j,k) = aws(i,j,k) + (wi*s(i,j,k) - aws(i,j,k))*inv
        end do
      end do
    end do
  end subroutine avg_scal
  !
  subroutine stats3d_save(datadir,ng,lo,hi,nscal,time,istep,xc_g,yc_g,zc_g,myid)
    integer , intent(in) :: ng(3),lo(3),hi(3),nscal,istep,myid
    real(rp), intent(in) :: time
    character(len=*), intent(in) :: datadir
    real(rp), intent(in), dimension(0:) :: xc_g,yc_g,zc_g
    real(rp) :: tt
    integer  :: it,m,is,iunit
    character(len=3) :: scalnum
    !$acc wait
    !$acc update self(vm,sm,s2m,usm,wsm)
    tt = time; it = istep
    do m=1,nvel
      call write_visu_3d(datadir,trim(cvel(m))//'_mean.bin','log_visu_3d.out',trim(dvel(m)), &
                         ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,vm(:,:,:,m),xc_g,yc_g,zc_g)
      call load_one('w',trim(datadir)//trim(cvel(m))//'_mean_res.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi, &
                    vm(:,:,:,m),'mean',tt,it)
    end do
    do is=1,nscal
      write(scalnum,'(i3.3)') is
      call write_visu_3d(datadir,'s_' //scalnum//'_mean.bin','log_visu_3d.out','Scalar_S'//scalnum, &
                         ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,sm(:,:,:,is),xc_g,yc_g,zc_g)
      call write_visu_3d(datadir,'s2_'//scalnum//'_mean.bin','log_visu_3d.out','Scalar_S2_'//scalnum, &
                         ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,s2m(:,:,:,is),xc_g,yc_g,zc_g)
      call load_one('w',trim(datadir)//'s_' //scalnum//'_mean_res.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi, &
                    sm(:,:,:,is),'mean',tt,it)
      call load_one('w',trim(datadir)//'s2_'//scalnum//'_mean_res.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi, &
                    s2m(:,:,:,is),'mean',tt,it)
      call write_visu_3d(datadir,'us_'//scalnum//'_mean.bin','log_visu_3d.out','Scalar_US'//scalnum, &
                         ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,usm(:,:,:,is),xc_g,yc_g,zc_g)
      call write_visu_3d(datadir,'ws_'//scalnum//'_mean.bin','log_visu_3d.out','Scalar_WS'//scalnum, &
                         ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,wsm(:,:,:,is),xc_g,yc_g,zc_g)
      call load_one('w',trim(datadir)//'us_'//scalnum//'_mean_res.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi, &
                    usm(:,:,:,is),'mean',tt,it)
      call load_one('w',trim(datadir)//'ws_'//scalnum//'_mean_res.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi, &
                    wsm(:,:,:,is),'mean',tt,it)
    end do
    if(myid == 0) then
      open(newunit=iunit,file=trim(datadir)//'count_avg.out',status='replace')
      write(iunit,*) count_avg
      close(iunit)
    end if
  end subroutine stats3d_save
end module mod_stats3d
