module mod_ibm
    use mod_types 
    use mod_param
    use mod_thakkar
    use decomp_2d, only: update_halo
    implicit none
    !
    ! exact wall force per cell (`FV-fric` output, see `init_wall_force`): for each velocity
    ! component c, the points whose IB impulse is assigned along the wall directions
    ! (band points: the lambda part; solid points next to the fluid: the local fluxes), with
    ! up to wf_ntmax target cells; wf_imp accumulates the impulse per cell and component
    !
    integer , parameter :: wf_ntmax = 12
    integer :: wf_np(3) = 0
    integer , allocatable, dimension(:,:,:)   :: wf_ijk       ! (3,np,c)
    integer , allocatable, dimension(:,:)     :: wf_nt        ! (np,c)
    integer , allocatable, dimension(:,:,:,:) :: wf_tc        ! (3,wf_ntmax,np,c)
    real(rp), allocatable, dimension(:,:,:)   :: wf_tw        ! (wf_ntmax,np,c)
    real(rp), allocatable, dimension(:,:,:,:) :: wf_imp       ! (0:n1+1,0:n2+1,0:n3+1,c)
    logical , allocatable, dimension(:,:,:)   :: wf_cell      ! cells that may hold wall force
    contains 
    ! we initalize the ibm coef. here
    ! at the main.f90
    ! check if the real location(x,y,z) is in the given body shape 
    logical function isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
        implicit none
        logical , intent(in), dimension(0:1,3)   :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)   :: amp_l
        integer , intent(in), dimension(0:1,3)   :: n_wave
        real(rp), intent(in), dimension(0:1,3)   :: l_0
        real(rp), intent(in), dimension(0:1,3)   :: phase_l
        integer , intent(in), dimension(3)       :: n
        real(rp), intent(in), dimension(3)       :: l
        real(rp),intent(in)                      :: x,y,z 
        integer                                  :: t
        real(rp)                                 :: xyz(3)
        integer                                  :: i,ii
        real(rp),intent(in),optional             :: hmap(0:,0:)
        real(rp),intent(in),optional             :: l1_hmap,l2_hmap 
        integer,intent(in),optional              :: n1_hmap,n2_hmap
        xyz = [x,y,z]
        isInBody=.false.
        do t = 1,3
            i=modulo(t,3)+1
            ii=modulo(t+1,3)+1
            if(ibm_direction(0,t))then
                if(xyz(t)<=wall_height(0,t,xyz(i),xyz(ii),amp_l,n_wave,l_0,phase_l,l,&
                                       hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                    isInBody=.true.
                endif
            endif
            if(ibm_direction(1,t))then
                if(xyz(t)>=wall_height(1,t,xyz(i),xyz(ii),amp_l,n_wave,l_0,phase_l,l,&
                                       hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                    isInBody=.true.
                endif
            endif   
        end do 
    end function isInbody
    real(rp) function wall_height(side,t,x1,x2,amp_l,n_wave,l_0,phase_l,l,&
                                  hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap) result(zw)
        !
        ! position, along direction t, of the IB wall on `side` (0: solid below, 1: solid above)
        ! at the in-plane coordinates x1 (direction modulo(t,3)+1) and x2 (direction modulo(t+1,3)+1)
        !
        implicit none
        integer , intent(in)                     :: side,t
        real(rp), intent(in)                     :: x1,x2
        real(rp), intent(in), dimension(0:1,3)   :: amp_l
        integer , intent(in), dimension(0:1,3)   :: n_wave
        real(rp), intent(in), dimension(0:1,3)   :: l_0
        real(rp), intent(in), dimension(0:1,3)   :: phase_l
        real(rp), intent(in), dimension(3)       :: l
        real(rp),intent(in),optional             :: hmap(0:,0:)
        real(rp),intent(in),optional             :: l1_hmap,l2_hmap 
        integer,intent(in),optional              :: n1_hmap,n2_hmap
        real(rp),parameter                       :: pi = 3.141592653589793_rp
        real(rp)                                 :: height,loc1,loc2
        integer                                  :: i,ii
        real(rp)                                 :: dl_hmap(2)
        integer                                  :: i1(2),i2(2)
        real(rp)                                 :: r1(2),r2(2)
        real(rp)                                 :: w1(2),w2(2)
        i=modulo(t,3)+1
        ii=modulo(t+1,3)+1
        if(.not. use_hmap)then
            select case(trim(ibm_shape))
            case('sinusoid')
                height=amp_l(side,i)*0.5_rp*(1._rp+sin(2._rp*pi*&
                       real(n_wave(side,i)*x1/l(i),rp)+phase_l(side,i)))+&
                       amp_l(side,ii)*0.5_rp*(1._rp+sin(2._rp*pi*&
                       real(n_wave(side,ii)*x2/l(ii),rp)+phase_l(side,ii)))
            case('cones')
                !
                ! staggered truncated cones (Theobald et al., WCCM-ECCOMAS 2020): height amp_l(side,i),
                ! base radius 1.25*height, flank slope 55 degrees; n_wave(side,i) x n_wave(side,ii) units
                ! of two cones each, centred at (1/4,1/4) and (3/4,3/4) of the unit
                !
                block
                  real(rp) :: lu1,lu2,d1,d2,r,rb,rt
                  lu1 = l(i )/n_wave(side,i )
                  lu2 = l(ii)/n_wave(side,ii)
                  d1 = x1-0.25_rp*lu1; d1 = d1-lu1*nint(d1/lu1)
                  d2 = x2-0.25_rp*lu2; d2 = d2-lu2*nint(d2/lu2)
                  r  = sqrt(d1**2+d2**2)
                  d1 = x1-0.75_rp*lu1; d1 = d1-lu1*nint(d1/lu1)
                  d2 = x2-0.75_rp*lu2; d2 = d2-lu2*nint(d2/lu2)
                  r  = min(r,sqrt(d1**2+d2**2))
                  rb = 1.25_rp*amp_l(side,i)
                  rt = rb-amp_l(side,i)/tan(55._rp*pi/180._rp)
                  height = amp_l(side,i)*min(1._rp,max(0._rp,(rb-r)/(rb-rt)))
                end block
            case default ! 'macdonald'
                if(side == 0)then
                    height=amp_l(side,i)*(1._rp + &
                           cos(2._rp*pi*real(n_wave(side,i)*x1/l(i),rp))*&
                           cos(2._rp*pi*real(n_wave(side,ii)*x2/l(ii),rp)))
                else
                    height=amp_l(side,i)*(1._rp - &
                           cos(2._rp*pi*real(n_wave(side,i)*x1/l(i),rp))*&
                           cos(2._rp*pi*real(n_wave(side,ii)*x2/l(ii),rp)))
                endif
            end select
        else
            loc1=x1; loc2=x2 ! get_hmap_loc shifts them for side 1
            if(trim(hmap_mode)=="fit")then
                call get_hmap_loc(side,n1_hmap,n2_hmap,loc1,loc2,l(i),l(ii),dl_hmap,i1,i2,r1,r2,w1,w2)
            elseif(trim(hmap_mode)=="normal")then
                ! here we need to change the shape of the hmap otherwise we will use different sized dl_hmap...
                call get_hmap_loc(side,n1_hmap,n2_hmap,loc1,loc2,l1_hmap,l2_hmap,dl_hmap,i1,i2,r1,r2,w1,w2)
            else
                print*,"unknown hmap mode entered..."
                stop
            endif
            ! Bilinear interpolation 
            height=(hmap(i1(1),i2(2))*(w1(1))+hmap(i1(2),i2(2))*(w1(2)))*(w2(2))+&
                   (hmap(i1(1),i2(1))*(w1(1))+hmap(i1(2),i2(1))*(w1(2)))*(w2(1))
        endif
        if(side == 0)then
            zw=l_0(0,t)+height
        else
            zw=l_0(1,t)-height
        endif
    end function wall_height
    subroutine wall_height_range(ng,l,ibm_direction,amp_l,n_wave,l_0,phase_l,hmin,hmax,&
                                 hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
        !
        ! range of the roughness height above the bottom wall (side 0) and below the top
        ! wall (side 1) of the z-walls, sampled on a grid 4x finer than the mesh in x and y
        !
        implicit none
        integer , intent(in), dimension(3)       :: ng
        real(rp), intent(in), dimension(3)       :: l
        logical , intent(in), dimension(0:1,3)   :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)   :: amp_l
        integer , intent(in), dimension(0:1,3)   :: n_wave
        real(rp), intent(in), dimension(0:1,3)   :: l_0
        real(rp), intent(in), dimension(0:1,3)   :: phase_l
        real(rp), intent(out)                    :: hmin,hmax
        real(rp),intent(in),optional             :: hmap(0:,0:)
        real(rp),intent(in),optional             :: l1_hmap,l2_hmap 
        integer,intent(in),optional              :: n1_hmap,n2_hmap
        integer , parameter                      :: nref = 4
        integer                                  :: side,i,j
        real(rp)                                 :: x,y,h
        hmin = huge(1._rp); hmax = -huge(1._rp)
        do side=0,1
            if(.not.ibm_direction(side,3)) cycle
            do j=1,nref*ng(2)
                y = (j-1)*l(2)/(nref*ng(2)) ! fine-grid nodes: include the mesh faces and centres
                do i=1,nref*ng(1)
                    x = (i-1)*l(1)/(nref*ng(1))
                    h = wall_height(side,3,x,y,amp_l,n_wave,l_0,phase_l,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                    if(side == 1) h = l(3)-h
                    hmin = min(hmin,h)
                    hmax = max(hmax,h)
                end do
            end do
        end do
        if(hmax < hmin) then ! no z-wall
            hmin = 0._rp; hmax = 0._rp
        end if
    end subroutine wall_height_range
    subroutine get_hmap_loc(side,n1_hmap,n2_hmap,loc_1,loc_2,&
                            l1_hmap,l2_hmap,dl_hmap,i1,i2,r1,r2,w1,w2)
        integer,intent(in)                              :: n1_hmap,n2_hmap
        real(rp),intent(inout)                          :: loc_1,loc_2
        real(rp),intent(in)                             :: l1_hmap,l2_hmap
        real(rp),intent(out)                            :: dl_hmap(2)
        integer,intent(out)                             :: i1(2),i2(2)
        real(rp),intent(out)                            :: r1(2),r2(2)
        real(rp),intent(out)                            :: w1(2),w2(2)
        integer,intent(in)                              :: side
        !     i1(1),i2(2)------------|------------i1(2),i2(2)
        !          |                                    |
        !          |        loc_1,loc_2--->hmap_val     |
        !          |                                    |
        !     i1(1),i2(1)------------|-------------i1(2),i2(1)
        ! additionally we need to shift upper wall lx_hmap/2,ly_hmap/2
        ! so we add a check
        i1(:)=0;i2(:)=0;
        r1(:)=0;r2(:)=0;
        dl_hmap(:)=0._rp;
        dl_hmap(1)=l1_hmap/n1_hmap
        dl_hmap(2)=l2_hmap/n2_hmap

        select case(side)
            case(0)
                loc_1=loc_1
                loc_2=loc_2
            case(1)
                loc_1=modulo(loc_1+(0.5*n1_hmap*dl_hmap(1)),(n1_hmap*dl_hmap(1)))
                loc_2=modulo(loc_2+(0.5*n2_hmap*dl_hmap(2)),(n2_hmap*dl_hmap(2)))
        end select
        
        i1(1)=floor(real(loc_1/dl_hmap(1),kind=rp))
        i1(2)=i1(1)+1

        r1(1)=i1(1)*dl_hmap(1)
        r1(2)=i1(2)*dl_hmap(1)

        i1(1)=modulo(i1(1),n1_hmap)
        i1(2)=modulo(i1(2),n1_hmap)

        i2(1)=floor(real(loc_2/dl_hmap(2),kind=rp))
        i2(2)=i2(1)+1

        r2(1)=i2(1)*dl_hmap(2)
        r2(2)=i2(2)*dl_hmap(2)

        i2(1)=modulo(i2(1),n2_hmap)
        i2(2)=modulo(i2(2),n2_hmap)

        w1(1)=abs(loc_1-r1(2))/dl_hmap(1)
        w1(2)=abs(loc_1-r1(1))/dl_hmap(1)
        w2(1)=abs(loc_2-r2(2))/dl_hmap(2)
        w2(2)=abs(loc_2-r2(1))/dl_hmap(2)

    end subroutine get_hmap_loc

    ! staircase masks: .true. at the solid points
    ! we change the diL depending the velocity mask we are handling 
    ! e.g. we need to apply 1,0,0 for mask_u and we need to apply 0,1,0 for mask_v 
    ! 0,0,1 for mask_w
    subroutine set_ibm_staircase(lo,mask_id,dix,diy,diz,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,zc,zf)
        implicit none
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        logical,intent(inout),dimension(0:,0:,0:)   :: mask_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k
        integer                                     :: ii,jj,kk
        real(rp)                                    :: x,y,z
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional  :: zc,zf

        do k = lbound(mask_id,3),ubound(mask_id,3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    if(diz/=1)then
                    ! this means we are looking for either u or v so their location is at z center
                    ! we gonna use the senter of zc
                        z = zc(k)
                        ! we use k inestead of kk since kk is the global and k is the local array index
                    else
                    ! else than we are looking for the w which is located on the z face
                        z = zf(k)
                    endif
                    if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap).eqv..true.)then
                        mask_id(i,j,k) = .true. ! this means we are in the solid
                    endif
                end do 
            end do 
        end do
    end subroutine set_ibm_staircase
    subroutine set_ibm_band(laplacian_id,idx,lap_b)
        !
        ! compact list of the points with a nonzero laplacian (the band points next to the wall):
        ! indices idx(1:3,:) and laplacian lap_b(:), where `calc_a_b` computes A and B
        ! (A = B = 1 elsewhere)
        !
        implicit none
        real(rp),intent(in),dimension(0:,0:,0:)         :: laplacian_id
        integer ,intent(out),allocatable,dimension(:,:) :: idx
        real(rp),intent(out),allocatable,dimension(:)   :: lap_b
        integer :: i,j,k,ib
        allocate(idx(3,count(laplacian_id /= 0._rp)),lap_b(count(laplacian_id /= 0._rp)))
        ib = 0
        do k=lbound(laplacian_id,3),ubound(laplacian_id,3)
            do j=lbound(laplacian_id,2),ubound(laplacian_id,2)
                do i=lbound(laplacian_id,1),ubound(laplacian_id,1)
                    if(laplacian_id(i,j,k) /= 0._rp) then
                        ib = ib+1
                        idx(:,ib) = [i,j,k]
                        lap_b(ib) = laplacian_id(i,j,k)
                    end if
                end do
            end do
        end do
    end subroutine set_ibm_band
    !2nd order scheme--laplacian settings
    subroutine set_ibm_2nd(lo,mask_id,laplacian_id,dix,diy,diz&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,&
        zc,zf,dzc,dzf,use_hmap,band_id)
        implicit none
        logical,intent(in)                          :: mask_id(0:,0:,0:)
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        real(rp),intent(inout),dimension(0:,0:,0:)  :: laplacian_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k,ip,im,jp,jm,kp,km
        integer                                     :: ii,jj,kk,c
        real(rp)                                    :: x,y,z,xp,xm,yp,ym,zp,zm
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp)                                    :: lambda
        integer                                     :: n_dir
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional  :: zc,zf
        real(rp),intent(in),dimension(0:),optional  :: dzc,dzf   
        real(rp)                                    :: dzf_l,dzc_l
        logical, intent(in)                         :: use_hmap
        real(rp)                                    :: hmax,hmin
        real(rp)                                    :: dl_int,l_int
        integer                                     :: n_hidden,ncand
        logical                                     :: calc_inBetween
        logical                                     :: ibm_diagnostic
        logical,dimension(0:,0:,0:),intent(inout)   :: band_id
        calc_inBetween=.false.
        ibm_diagnostic=.false.
        if(use_hmap)then
            if (.not.present(hmap)) then
                error stop "use_hmap=T but hmap not present"
            endif
            n_hidden=0;
            ncand=0;
            dl_int=0._rp;
            l_int=0._rp;
            hmax=0._rp;
            hmin=0._rp;
            hmin=minval(hmap);
            hmax=maxval(hmap);
            !hmax=hmax-hmin;
        else
            print*,"some problem between use_hmap and hmap"! doesnt mean problem
            n_hidden=0;
            ncand=0;
            dl_int=0._rp;
            l_int=0._rp;
            hmax=0._rp;
            hmin=0._rp;
        endif
        do k = 1,n(3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    dzf_l=0._rp;dzc_l=0._rp
                    xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);
                    call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                    do n_dir=1,6
                        select case(n_dir)
                            case(1)
                                ! xp
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xp,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band and fluid ! these are for friction calc.
                                endif
                                !mask_id(i,j,k) already has the velocity body information
                                if(calc_inBetween)then 
                                    ! this means we are close to body region and both cells that we check are fluid now
                                    ! we check if in between part has body!
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                xp,y,z,n,l,&
                                                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            ! we check 15 times
                                            dl_int=dl(1)*real(c,kind=rp)/15
                                            l_int=x+dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit ! we need to exit otherwise we gonna keep adding to laplacian !note#2 not sure what should be done
                                            endif
                                        end do
                                        
                                    endif
                                endif
                            case(2)
                                ! xm
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xm,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,xm,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                xm,y,z,n,l,&
                                                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            ! we check 4 times since 5th time is on the neigbour which it is in fluid
                                            dl_int=dl(1)*real(c,kind=rp)/15
                                            l_int=x-dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif
                            case(3)
                                ! yp
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,yp,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,yp,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                x,yp,z,n,l,&
                                                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            dl_int=dl(2)*real(c,kind=rp)/15
                                            l_int=y+dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif
                            case(4)
                                ! ym
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,ym,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,ym,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                x,ym,z,n,l,&
                                                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            dl_int=dl(2)*real(c,kind=rp)/15
                                            l_int=y-dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif       
                            case(5)
                                ! zp 
                                if(diz==0)then
                                    dzc_l=dzc(k)
                                    dzf_l=dzf(k)
                                else
                                    dzf_l=dzf(k+1)
                                    dzc_l=dzc(k)
                                endif
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zp,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,zp,3,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                            case(6)
                                ! zm
                                if(diz==0)then
                                    dzc_l=dzc(k-1)
                                    dzf_l=dzf(k)
                                else
                                    dzf_l=dzf(k)
                                    dzc_l=dzc(k)
                                endif 
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zm,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,zm,3,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                        end select
                    end do
                end do 
            end do 
        end do
        if(ibm_diagnostic)then
            c=0
            print*,"max height of the hmap: ",hmax,hmin
            print*, "possible problematic locations: ",n_hidden,ncand
            print*, "number of band locations: ",count(band_id)
            print*, "number of IBM locations: ",count(mask_id)
            do k=1,n(3)
                do j=1,n(2)
                    do i=1,n(1)
                        if(band_id(i,j,k).and.mask_id(i,j,k))then
                            print*,"some problems with ibm processing"
                            print*,"some points are in the band and solid mask"
                            c=c+1
                        endif
                    enddo
                enddo
            enddo
            print*,"overlapping cells: ", c
        endif
    end subroutine set_ibm_2nd

    subroutine calc_lambda(x,y,z,l_n,case_num,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf,dzc,diz,wall_loc)
        implicit none
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer,intent(in)                          ::case_num
        real(rp),intent(in)                         :: x,y,z,l_n
        real(rp),intent(out)                        :: lambda
        real(rp)                                    :: l_fluid,l_solid,l_int,l_diff
        real(rp)                                    :: eps
        integer                                     ::  n_iter
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),optional                :: dzf,dzc
        integer,intent(in),optional                 :: diz   
        real(rp)                                    :: dz
        real(rp),intent(out),optional,dimension(3)  :: wall_loc
        lambda=0._rp
        dz=0._rp
        if(case_num==3)then
            if(diz==0)then
                ! this means we are looking for u or v
                ! their location is at zc 
                dz=dzc
            else
                ! this means we are looking for w
                ! its location is at zf 
                dz=dzf
            endif
        endif
        select case(case_num)
            case(1)
                eps = 1.e-10_rp*dl(1)
            case(2)
                eps = 1.e-10_rp*dl(2)
            case(3)
                eps = 1.e-10_rp*dz
        end select
        select case(case_num)
            case(1)
                l_fluid = x;l_solid = l_n;l_int=0._rp;l_diff=0._rp
            case(2)
                l_fluid = y;l_solid = l_n;l_int=0._rp;l_diff=0._rp
            case(3)
                l_fluid = z;l_solid = l_n;l_int=0._rp;l_diff=0._rp
        end select
        do n_iter=1,60
            l_int=real((l_solid+l_fluid)/2._rp,kind=rp)
            select case(case_num)
                case(1)!x
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(2)!y
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(3)!z
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,l_int,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif
            end select
        enddo
        select case(case_num)
            case(1)
                l_diff=abs(x-l_int)
            case(2)
                l_diff=abs(y-l_int)
            case(3)
                l_diff=abs(z-l_int)
        end select
        if(l_diff<=eps)then
            l_diff=eps
        endif
        select case(case_num)
            
            case(1)
                lambda=real((1._rp/dl(1)**2)*(dl(1)/l_diff-1._rp),kind=rp)
            case(2)
                lambda=real((1._rp/dl(2)**2)*(dl(2)/l_diff-1._rp),kind=rp)
            case(3)
                ! imporant note here  we need both dzc and
                ! dzf since their gradinet either at cell center or face or vice versa
                !                      !
                lambda=real((1._rp/(dzc*dzf))*((dz)/l_diff-1._rp),kind=rp)
        end select 
        ! return the wall location
        if(present(wall_loc))then
            select case(case_num)
            case(1)
                wall_loc=[l_int,y,z]
            case(2)
                wall_loc=[x,l_int,z]
            case(3)
                wall_loc=[x,y,l_int]
            end select
        endif 
    end subroutine calc_lambda
    subroutine calc_grad_dist(grad_dist_id,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_id,dix,diy,diz,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf,dzc)
        ! routine creates the tangent and normal plane for 3D IBM
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id                      
        logical , intent(in), dimension(0:1,3)                  :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)                  :: amp_l
        integer , intent(in), dimension(0:1,3)                  :: n_wave
        real(rp), intent(in), dimension(0:1,3)                  :: l_0
        real(rp), intent(in), dimension(0:1,3)                  :: phase_l
        real(rp), intent(in   ), dimension(3)                   :: l
        integer , intent(in   ), dimension(3)                   :: n
        real(rp),intent(in),optional                            :: hmap(0:,0:)
        integer,intent(in),optional                             :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                            :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional              :: dzf,dzc

        logical, intent(in), dimension(0:,0:,0:)                :: band_id
        integer,intent(in)                                      :: dix,diy,diz
        integer                                                 :: i,j,k,c,wc,nc
        real(rp),dimension(0:2,6,0:26)                          :: wall_loc
        logical,dimension(6,0:26)                               :: wall_loc_log
        integer,intent(in),dimension(3)                         :: lo
        real(rp),intent(in),dimension(0:)                       :: zc,zf
        real(rp), intent(in),dimension(3)                       :: dl
        real(rp)                                                :: x,y,z
        real(rp)                                                :: xp,xm,yp,ym,zp,zm
        real(rp)                                                :: lambda!just to fullfill the calc_lamda
        integer                                                 :: di,dj,dk,in,jn,kn,count
        real(rp),dimension(3)                                   :: wall_loc_real
        real(rp)                                                :: dzf_l,dzc_l
        real(rp),dimension(3,2)                                 :: plane
        real(rp),dimension(3)                                   :: nv,cp,bp
        real(rp)                                                :: angle_plane,out_plane,grad

        ! init vars
        grad_dist_id(:,:) = 0._rp 
        count = 0
        angle_plane=0._rp;out_plane=0._rp;grad=0._rp
        wall_loc_log(:,:)=.false.
        wall_loc_real(:)=0._rp
        wall_loc(:,:,:)=0._rp
        x=0._rp;y=0._rp;z=0._rp;
        dzc_l=0._rp;dzf_l=0._rp
        xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
        ! calc part
        do k=(lbound(band_id,3)+1),(ubound(band_id,3)-1)
            do j=(lbound(band_id,2)+1),(ubound(band_id,2)-1)
                do i=(lbound(band_id,1)+1),(ubound(band_id,1)-1)
                    if(band_id(i,j,k))then ! this means we are in band
                        wall_loc(:,:,:)=0._rp
                        wall_loc_log(:,:)=.false.
                        wc=0
                       
                        do nc=0,26 ! time to check all possible near cells
                            
                            di=modulo(nc,3)-1
                            dj=modulo(nc/3,3)-1
                            dk=modulo(nc/9,3)-1
                            in=i+di;jn=j+dj;kn=k+dk;
                            if(band_id(in,jn,kn))then
                                x=0._rp;y=0._rp;z=0._rp;xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
                                call get_grid_loc(lo,in,jn,kn,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm) 
                                do c=1,6    
                                    select case(c)
                                        case(1)!xp
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xp,y,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 

                                                call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif   
                                        case(2)!i-1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xm,y,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 

                                                call calc_lambda(x,y,z,xm,1,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif   
                                        case(3)!j+1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,yp,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                            
                                                call calc_lambda(x,y,z,yp,2,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif    
                                        case(4)!j-1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,ym,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,ym,2,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp
                                            end if    
                                        case(5)!k+1
                                            if(diz==0)then
                                                dzc_l=dzc(kn)
                                                dzf_l=dzf(kn)
                                            else
                                                dzf_l=dzf(kn+1)
                                                dzc_l=dzc(kn)
                                            endif
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zp,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,zp,3,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp 
                                            endif  
                                        case(6)!k-1
                                            if(diz==0)then
                                                dzc_l=dzc(kn-1)
                                                dzf_l=dzf(kn)
                                            else
                                                dzf_l=dzf(kn)
                                                dzc_l=dzc(kn)
                                            endif
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zm,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,zm,3,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp
                                            endif
                                        case default
                                            print*,"something is off.... at ibm.f90"   
                                    end select
                                end do
                            endif
                        end do
                        if(wc<=3)cycle
                        count=count+1
                        call get_plane(wall_loc,wall_loc_log,plane)
                        !now we have the normal and center now lets get the distance
                        !lets get the location of the band point
                        x=0._rp;y=0._rp;z=0._rp;xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
                        call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                        !output of these are the location of the bandpoint
                        bp(1)=x;bp(2)=y;bp(3)=z;
                        cp(1)=bp(1)-plane(1,1);
                        cp(2)=bp(2)-plane(2,1);
                        cp(3)=bp(3)-plane(3,1);
                        !now we have exact distance from the center point of the plane
                        !to the band point! lets get them scalar product
                        nv=plane(:,2)
                        call comp_sca(nv,cp,angle_plane,out_plane)
                        if(out_plane<0._rp)nv=-nv ! we check if the output is neg.
                        ! this means the normal vector is facing inwards, we want outwards facing
                        call comp_sca(nv,cp,angle_plane,out_plane)
                        grad_dist_id(count,1)=plane(1,1)
                        grad_dist_id(count,2)=plane(2,1)
                        grad_dist_id(count,3)=plane(3,1)
                        grad_dist_id(count,4)=nv(1)
                        grad_dist_id(count,5)=nv(2)
                        grad_dist_id(count,6)=nv(3)
                        grad_dist_id(count,7)=i!these are the locations of the velocities
                        grad_dist_id(count,8)=j!these are the locations of the velocities
                        grad_dist_id(count,9)=k!these are the locations of the velocities
                        grad_dist_id(count,10)=out_plane
                        grad_dist_id(count,11)=0._rp
                        grad_dist_id(count,12)=0._rp
                        ! here we only compute the geometrical distances! 
                        ! since we have a stat. wall we dont need to compute the distances
                        ! in each time step!
                    endif
                end do 
            end do 
        end do
    end subroutine calc_grad_dist
    subroutine calc_shear_st(grad_dist_id,vel_id,visc_loc)
        ! we assume the walls arent moving therefore we do the distance calc.
        ! time and compute the gradient here at each wanted time step
        real(rp),dimension(0:,0:,0:),intent(in)                 :: vel_id
        real(rp),intent(in)                                     :: visc_loc ! not the module visc: device code
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id
        integer                                                 :: i,j,k,kk
        real(rp)                                                :: grad,dist_grad,shear_st
        grad=0._rp;i=0;j=0;k=0;
        !$acc parallel loop private(grad,dist_grad,shear_st,kk,i,j,k)
        do kk=lbound(grad_dist_id,1),ubound(grad_dist_id,1)
            grad=0._rp
            ! now we go through the data
            i=grad_dist_id(kk,7)
            j=grad_dist_id(kk,8)
            k=grad_dist_id(kk,9)
            dist_grad=grad_dist_id(kk,10)
            !we have the location of the i,j,k for the velocity
            call comp_grad(dist_grad,vel_id(i,j,k),grad)
            shear_st=visc_loc*grad
            grad_dist_id(kk,11)=grad
            grad_dist_id(kk,12)=shear_st
        end do 
        
    end subroutine calc_shear_st
    subroutine comp_grad(d,vel,grad)
        !$acc routine seq
        real(rp), intent(in)                      :: vel
        real(rp),intent(in)                       :: d
        real(rp),intent(out)                      :: grad
        if(abs(d)>1e-10)then
            grad=vel*d**(-1)
        else
            grad=0._rp
        endif
    end subroutine comp_grad
    subroutine get_plane(wall_loc,wall_loc_log,plane)
        real(rp),intent(in),dimension(0:2,6,0:26)           :: wall_loc
        logical,intent(in),dimension(6,0:26)                :: wall_loc_log
        integer                                             :: o,p,m,n,co,g,f   
        real(rp),dimension(3)                               :: center
        real(rp),dimension(0:2,6,0:26)                      :: vec
        real(rp),dimension(3)                               :: vec1,vec2,vec3,n1,n2,n3,n_av,n_ref,norm_nav
        real(rp),dimension(0:2,6,0:26)                      :: normal
        logical,dimension(6,0:26)                           :: normal_loc
        real(rp)                                            :: angle1,angle2,angle3
        real(rp),intent(out),dimension(3,2)                 :: plane
        real(rp)                                            :: out_dum,out_n,angle_dum
        real(rp),dimension(3,3)                             :: cov
        integer                                             :: c_n

        !first lets create a center
        angle1=90._rp;angle2=90._rp;angle3=90._rp;
        plane(:,:)=0._rp
        normal_loc(:,:)=.false.
        vec(:,:,:) = 0._rp
        vec1(:) = 0._rp;vec2(:) = 0._rp;vec3(:) = 0._rp;
        n1(:) = 0._rp;n2(:) = 0._rp;n3(:) = 0._rp;
        n_av(:) = 0._rp;n_ref(:) = 0._rp;norm_nav(:)=0._rp;
        center(:)=0._rp
        co=0;c_n=0
        do n=0,26
            do m=1,6
                if(wall_loc_log(m,n))then
                    center(:)=center(:)+wall_loc(:,m,n)
                    co=co+1
                endif
            enddo
        enddo
        if(co>1)then
            center(:)=center(:)/co
        endif
        !now lets get the vectors from the center
        do n=0,26
            do m=1,6
                if(wall_loc_log(m,n))then
                    vec(:,m,n)=wall_loc(:,m,n)-center(:)
                endif
            enddo
        enddo
        ! least-squares plane through the wall points: the normal is the eigenvector of
        ! the covariance matrix C = sum(vec vec^T) with the smallest eigenvalue
        cov(:,:)=0._rp
        do n=0,26
            do m=1,6
                if(wall_loc_log(m,n))then
                    do g=1,3
                        cov(:,g)=cov(:,g)+vec(:,m,n)*vec(g-1,m,n)
                    enddo
                endif
            enddo
        enddo
        call smallest_eigvec_sym3(cov,n_av)
        call nomrlalize_vec(n_av,norm_nav)
        plane(:,1)=center(:)
        plane(:,2)=norm_nav(:)
        ! now we have our normals lets check each normal w/ each vector 
    end subroutine get_plane
    subroutine smallest_eigvec_sym3(a,ev)
        ! eigenvector of the smallest eigenvalue of a symmetric 3x3 matrix
        ! (closed-form eigenvalues, eigenvector from the largest cross product of the
        ! rows of a - lambda*I); returns 0 for a zero matrix
        real(rp),intent(in),dimension(3,3)  :: a
        real(rp),intent(out),dimension(3)   :: ev
        real(rp),parameter                  :: pi = 3.141592653589793_rp
        real(rp)                            :: p1,p2,p,q,r,phi,lam,mag,mag_c
        real(rp),dimension(3,3)             :: b
        real(rp),dimension(3)               :: c
        integer                             :: i1,i2
        ev(:)=0._rp
        p1=a(1,2)**2+a(1,3)**2+a(2,3)**2
        q=(a(1,1)+a(2,2)+a(3,3))/3._rp
        p2=(a(1,1)-q)**2+(a(2,2)-q)**2+(a(3,3)-q)**2+2._rp*p1
        if(p2<=tiny(1._rp))return ! a = q*I: no preferred direction
        p=sqrt(p2/6._rp)
        b(:,:)=a(:,:)
        b(1,1)=b(1,1)-q;b(2,2)=b(2,2)-q;b(3,3)=b(3,3)-q
        b(:,:)=b(:,:)/p
        r=0.5_rp*(b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2)) &
                 -b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1)) &
                 +b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1)))
        r=min(max(r,-1._rp),1._rp)
        phi=acos(r)/3._rp
        lam=q+2._rp*p*cos(phi+2._rp*pi/3._rp) ! smallest eigenvalue
        b(:,:)=a(:,:)
        b(1,1)=b(1,1)-lam;b(2,2)=b(2,2)-lam;b(3,3)=b(3,3)-lam
        mag=0._rp
        do i1=1,2
            do i2=i1+1,3
                c(1)=b(i1,2)*b(i2,3)-b(i1,3)*b(i2,2) ! not comp_cross: it normalises
                c(2)=b(i1,3)*b(i2,1)-b(i1,1)*b(i2,3)
                c(3)=b(i1,1)*b(i2,2)-b(i1,2)*b(i2,1)
                mag_c=norm2(c)
                if(mag_c>mag)then
                    mag=mag_c
                    ev(:)=c(:)
                endif
            enddo
        enddo
    end subroutine smallest_eigvec_sym3
    subroutine comp_cross(vec1,vec2,vec3)
        real(rp),intent(in),dimension(3) :: vec1,vec2
        real(rp),intent(out),dimension(3) :: vec3
        real(rp)                          :: i,j,k
        real(rp)                          :: n
        k=(vec1(1)*vec2(2))-(vec1(2)*vec2(1))
        j=(vec1(3)*vec2(1))-(vec1(1)*vec2(3))
        i=(vec1(2)*vec2(3))-(vec1(3)*vec2(2))
        n=sqrt(i**2+j**2+k**2)
        if(abs(n)<1e-10)then
            vec3=0._rp
        else
            i=i*n**(-1);j=j*n**(-1);k=k*n**(-1)
            vec3(:)=[i,j,k]
        endif
        
    end subroutine comp_cross
    subroutine nomrlalize_vec(vec1,vec2)
        real(rp),intent(in),dimension(3)     :: vec1 
        real(rp),intent(out),dimension(3)    :: vec2
        real(rp)                             :: i,j,k,n
        n=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
        if(n>epsilon(1._rp))then
            vec2(1)=vec1(1)*n**(-1)
            vec2(2)=vec1(2)*n**(-1)
            vec2(3)=vec1(3)*n**(-1)
        else
            vec2(1)=0._rp;vec2(2)=0._rp;vec2(3)=0._rp
        endif
    end subroutine nomrlalize_vec
    subroutine get_magni(vec1,mag)
        real(rp),intent(in),dimension(3)     :: vec1 
        real(rp),intent(out) :: mag
        mag=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
    end subroutine get_magni
    subroutine comp_sca(vec1,vec2,angle,out)
        real(rp),intent(in),dimension(3) :: vec1,vec2
        real(rp),intent(out)             :: angle
        real(rp),intent(out)             :: out
        real(rp)                         :: b1,b2
        out=vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
        b1=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
        b2=sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)
        if(abs(b1)>1e-10.and.abs(b2)>1e-10)then
            angle=out*(b1*b2)**(-1)
            angle=acos(angle)*180/(3.1415) !radians to degree
        else
            !we give out as epsilon since we will use it to divide
            out=epsilon(0._rp)
        endif
    end subroutine comp_sca
    subroutine get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
        !$acc routine seq
        integer,intent(in)                      :: i,j,k
        real(rp),intent(out)                    :: x,y,z,xp,xm,yp,ym,zp,zm
        integer,intent(in)                      :: dix,diy,diz
        integer                                 :: ii,jj,kk
        integer,intent(in),dimension(3)         :: lo
        real(rp),intent(in),dimension(0:),optional  :: zc,zf
        real(rp), intent(in),dimension(3)       :: dl
        x=0;y=0;z=0;xp=0;xm=0;yp=0;ym=0;zp=0;zm=0;    
        ii = lo(1)+i-1
        jj = lo(2)+j-1
        kk = lo(3)+k-1
        x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
        y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
        xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);
        if(diz/=1)then
            ! this means we are looking for either u or v so their location is at z center
            ! we gonna use the senter of zc
            z = zc(k)
            zp = zc(k+1)
            zm = zc(k-1)
            ! we use k inestead of kk since kk is the global and k is the local array index
        else
            ! else than we are looking for the w which is located on the z face
            z = zf(k)
            zp = zf(k+1)
            zm = zf(k-1)
        endif                 
    end subroutine get_grid_loc
    ! here we do the wall pressure calc stuff

    subroutine get_wall_pres(p,p_grad,lo,dl,zc,zf,mask_s,band_s)
        real(rp),intent(in),dimension(0:,0:,0:)     :: p                        
        real(rp),dimension(:,:),intent(inout)       :: p_grad
        integer,intent(in),dimension(3)             :: lo
        real(rp),intent(in),dimension(0:)           :: zc,zf
        real(rp), intent(in),dimension(3)           :: dl
        logical,intent(in)                          :: mask_s(0:,0:,0:)
        logical, intent(in), dimension(0:,0:,0:)    :: band_s
        real(rp),dimension(3)                       :: c0,p_normalv
        integer                                     :: kk,dix,diy,diz 
        integer                                     :: i,j,k,m,in,jn,kn,di,dj,dk,count
        real(rp)                                    :: dist,x,y,z,xp,xm,yp,ym,zp,zm
        real(rp),dimension(3)                       :: vec                                    
        real(rp)                                    :: pw,b !b is the the pressure gradient
        real(rp)                                    :: s0,s1,s2,t0,t1
        kk=0;dix=0;diy=0;diz=0;m=0;count=0;
        dist=0._rp;vec(:)=0._rp
        ! now we have the pressures that are on the surface but arent in the solid
        ! lets get their planes etc. we can directly call the grad distance
        ! to create a normal and a plane for each pressure band point 
        ! like we did w/ the velocities
        ! this subroutine directly returns us an array where the pressure band locations that we located
        ! have their normal and a center point 

        !$acc parallel loop private(m,c0,p_normalv,i,j,k,s0,s1,s2,t0,t1) &
        !$acc& private(vec,pw,b,dist,x,y,z,xp,xm,yp,ym,zp,zm) &
        !$acc& private(di,dj,dk,in,jn,kn,count,dix,diy,diz)
        do kk=lbound(p_grad,1),ubound(p_grad,1)
            c0(1)=p_grad(kk,1)
            c0(2)=p_grad(kk,2)
            c0(3)=p_grad(kk,3)
            p_normalv(1)=p_grad(kk,4)
            p_normalv(2)=p_grad(kk,5)
            p_normalv(3)=p_grad(kk,6)
            i=p_grad(kk,7)
            j=p_grad(kk,8)
            k=p_grad(kk,9)
            ! so lets search the around the 
            s0=0._rp;s1=0._rp;s2=0._rp;t0=0._rp;t1=0._rp;count=0._rp;
            do m=0,26
                di=modulo(m,3)-1
                dj=modulo(m/3,3)-1
                dk=modulo(m/9,3)-1
                in=i+di;jn=j+dj;kn=k+dk
                !so we check if the neighbour of the band poin is in fluid
                if(kn==k.and.jn==j.and.in==i)cycle ! we avoid the band point we are on
                if(in<lbound(p,1).or.in>ubound(p,1)) cycle
                if(jn<lbound(p,2).or.jn>ubound(p,2)) cycle
                if(kn<lbound(p,3).or.kn>ubound(p,3)) cycle 
                if(mask_s(in,jn,kn))cycle!mask_s is only true if the point is in body
                if(band_s(in,jn,kn))cycle!if the location is also a band point cycle
                call get_grid_loc(lo,in,jn,kn,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                ! comp vec
                count=count+1
                vec(1)=x-c0(1);vec(2)=y-c0(2);vec(3)=z-c0(3)
                dist=dot_product(vec,p_normalv)
                !now we have the distance to the center
                s0=count
                s2=s2+dist**2
                s1=s1+dist
                t0=t0+p(in,jn,kn)
                t1=t1+dist*p(in,jn,kn)
            end do
            if(s0<3 .or. abs(s0*s2-s1**2) <= epsilon(1._rp)*s0*s2)cycle ! too few points, or all at the same distance
            ! we use direct mse from the gatti's source 
            pw=(s2*t0-s1*t1)*(s0*s2-s1**2)**(-1)
            b=(s0*t1-s1*t0)*(s0*s2-s1**2)**(-1)
            p_grad(kk,11)=pw
            p_grad(kk,12)=b
        end do 
    end subroutine get_wall_pres
    subroutine init_fric_cubes(band_s,mask_s,lo,dl,zc,zf,dzc,&
                                        dzf,xg,yg,dix,diy,diz,ibm_direction,amp_l,&
                                        n_wave,l_0,phase_l,n,l,hmap,l1_hmap,l2_hmap,&
                                        n1_hmap,n2_hmap,output)
        logical, intent(in),dimension(0:,0:,0:)             :: band_s
        logical,intent(in)                                  :: mask_s(0:,0:,0:)
        integer,intent(in),dimension(3)                     :: lo
        real(rp), intent(in),dimension(3)                   :: dl
        real(rp),intent(in),dimension(0:)                   :: zc,zf,xg,yg,dzf,dzc
        logical , intent(in), dimension(0:1,3)              :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)              :: amp_l
        integer , intent(in), dimension(0:1,3)              :: n_wave
        real(rp), intent(in), dimension(0:1,3)              :: l_0
        real(rp), intent(in), dimension(0:1,3)              :: phase_l
        integer , intent(in), dimension(3)                  :: n
        real(rp), intent(in), dimension(3)                  :: l
        real(rp),intent(in),optional                        :: hmap(0:,0:)
        real(rp),intent(in),optional                        :: l1_hmap,l2_hmap 
        integer,intent(in),optional                         :: n1_hmap,n2_hmap
        integer,intent(in)                                  :: dix,diy,diz
        real(rp),intent(inout),dimension(0:,0:,0:,:)        :: output
        integer                                             :: i,j,k
        output(:,:,:,:) = 0._rp;
        do k=lbound(band_s,3)+1,ubound(band_s,3)-1
            do j=lbound(band_s,2)+1,ubound(band_s,2)-1
                do i=lbound(band_s,1)+1,ubound(band_s,1)-1
                    if(band_s(i,j,k))then
                        call get_fric_cubes(band_s,i,j,k,mask_s,lo,dl,zc,zf,dzc,dzf,xg,yg,dix,diy,diz,&
                                            ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,hmap,l1_hmap,&
                                            l2_hmap,n1_hmap,n2_hmap,output)
                    endif
                end do 
            end do
        enddo


    end subroutine init_fric_cubes
    subroutine get_fric_cubes(band_s,i,j,k,mask_s,lo,dl,zc,zf,dzc,dzf,xg,yg,dix,diy,diz,ibm_direction,amp_l,n_wave,l_0,&
        phase_l,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,output)
        logical, intent(in),dimension(0:,0:,0:)             :: band_s
        logical,intent(in)                                  :: mask_s(0:,0:,0:)
        integer,intent(in),dimension(3)                     :: lo
        real(rp), intent(in),dimension(3)                   :: dl
        real(rp),intent(in),dimension(0:)                   :: zc,zf,xg,yg,dzf,dzc
        logical , intent(in), dimension(0:1,3)              :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)              :: amp_l
        integer , intent(in), dimension(0:1,3)              :: n_wave
        real(rp), intent(in), dimension(0:1,3)              :: l_0
        real(rp), intent(in), dimension(0:1,3)              :: phase_l
        integer , intent(in), dimension(3)                  :: n
        real(rp), intent(in), dimension(3)                  :: l
        real(rp),intent(in),optional                        :: hmap(0:,0:)
        real(rp),intent(in),optional                        :: l1_hmap,l2_hmap 
        integer,intent(in),optional                         :: n1_hmap,n2_hmap
        integer,intent(in)                                  :: i,j,k,dix,diy,diz
        integer                                             :: ic,jc,kc
        real(rp)                                            :: x,y,z,xp,xm,yp,ym,zp,zm,x_int,y_int,z_int
        real(rp)                                            :: dx_int,dy_int,dz_int
        real(rp)                                            :: zw_lo,zw_hi
        integer                                             :: co
        real(rp)                                            :: cell_vol,fluid_vol,solid_vol
        integer, parameter                                  :: nsub = 20 ! sub-cells per direction
        logical,dimension(nsub,nsub,nsub)                   :: mask_subcell
        integer                                             :: s_xp,s_xm,s_yp,s_ym,s_zp,s_zm 
        real(rp)                                            :: A_dx,A_dy,A_dz,A_x,A_y,A_z
        real(rp)                                            :: A_xsi,A_ysi,A_zsi
        real(rp)                                            :: x_in,x_out,y_in,y_out,z_in,z_out
        real(rp),dimension(3)                               :: normal
        real(rp),dimension(2,3)                             :: A_inlets
        real(rp)                                            :: A_tot 
        real(rp),intent(inout),dimension(0:,0:,0:,:)          :: output
        logical :: diag,inloop
        !options
        diag=.false.;inloop=.false.;
        s_xp=0;s_xm=0;s_yp=0;s_ym=0;s_zp=0;s_zm=0;
        A_dx=0._rp;A_dy=0._rp;A_dz=0._rp;
        A_x=0._rp;A_y=0._rp;A_z=0._rp;
        A_xsi=0;A_ysi=0;A_zsi=0;
        A_tot  = 0._rp;normal = 0._rp
        if(band_s(i,j,k))then ! we have a band point
            if(i==0.or.j==0.or.k==0)return !guard
            mask_subcell(:,:,:)=.false.
            x=0._rp;y=0._rp;z=0._rp;xp=0._rp;
            xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp;
            x_int=0._rp;y_int=0._rp;z_int=0._rp;
            co = 0;
            call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
            ! now lets check the inside of the band point with smaller cubes
            dz_int=dzf(k);dy_int=(yg(j)-yg(j-1));dx_int=(xg(i)-xg(i-1))!get dx,dy,dz from global
            A_dx=dy_int*dz_int/nsub**2;A_dy=dx_int*dz_int/nsub**2;A_dz=dx_int*dy_int/nsub**2
            if(.not.any(ibm_direction(:,1:2))) then
              ! z-walls only: the solid test is monotone in z, so the wall heights are evaluated
              ! once per sub-column (the same comparisons as in isInbody)
              do jc=1,nsub
                y_int=y-0.5_rp*dy_int+(real(jc,kind=rp)-0.5_rp)/nsub*dy_int
                do ic=1,nsub
                  x_int=x-0.5_rp*dx_int+(real(ic,kind=rp)-0.5_rp)/nsub*dx_int
                  zw_lo = -huge(1._rp); zw_hi = huge(1._rp)
                  if(ibm_direction(0,3)) zw_lo = wall_height(0,3,x_int,y_int,amp_l,n_wave,l_0,phase_l,l, &
                                                             hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                  if(ibm_direction(1,3)) zw_hi = wall_height(1,3,x_int,y_int,amp_l,n_wave,l_0,phase_l,l, &
                                                             hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                  do kc=1,nsub
                    z_int=z-0.5_rp*dz_int+(real(kc,kind=rp)-0.5_rp)/nsub*dz_int ! sub-cell centre
                    if(z_int<=zw_lo.or.z_int>=zw_hi)then
                      mask_subcell(ic,jc,kc)=.true.
                      co=co+1
                    endif
                  end do
                end do
              end do
            else
            do kc=1,nsub
                z_int=z-0.5_rp*dz_int+(real(kc,kind=rp)-0.5_rp)/nsub*dz_int ! sub-cell centre
                do jc=1,nsub
                    y_int=y-0.5_rp*dy_int+(real(jc,kind=rp)-0.5_rp)/nsub*dy_int
                    do ic=1,nsub
                        x_int=x-0.5_rp*dx_int+(real(ic,kind=rp)-0.5_rp)/nsub*dx_int
                        if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x_int,y_int,z_int,n,l,&
                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                            mask_subcell(ic,jc,kc)=.true. ! this subcell is in solid
                            co=co+1
                        endif
                    end do
                end do
            end do
            end if
            do kc=1,nsub
                do jc=1,nsub
                    do ic=1,nsub
                            if(mask_subcell(ic,jc,kc))then
                                if(ic<nsub)then
                                    if(.not.mask_subcell(ic+1,jc,kc))then ! neigbour isnt in solid
                                        s_xp=s_xp+1! thus its area is what we look
                                    endif
                                endif
                                if(ic>1)then
                                    if(.not.mask_subcell(ic-1,jc,kc))then
                                        s_xm=s_xm+1
                                    endif
                                endif
                                if(jc<nsub)then
                                    if(.not.mask_subcell(ic,jc+1,kc))then
                                        s_yp=s_yp+1
                                    endif
                                endif
                                if(jc>1)then
                                    if(.not.mask_subcell(ic,jc-1,kc))then
                                        s_ym=s_ym+1
                                    endif
                                endif
                                if(kc<nsub)then
                                    if(.not.mask_subcell(ic,jc,kc+1))then
                                        s_zp=s_zp+1
                                    endif
                                endif
                                if(kc>1)then
                                    if(.not.mask_subcell(ic,jc,kc-1))then
                                        s_zm=s_zm+1
                                    endif
                                endif
                            endif                       
                    end do 
                end do 
            end do
            !print*,co
            if(co==0)then ! this means they are band but they sit perfectly
                if(i>lbound(mask_s,1))then
                    if(mask_s(i-1,j,k))then
                        normal=[1._rp,0._rp,0._rp]
                        A_xsi=dy_int*dz_int
                        A_tot=A_xsi
                    endif
                endif
                if(i<ubound(mask_s,1))then
                    if(mask_s(i+1,j,k))then
                        normal=[-1._rp,0._rp,0._rp]
                        A_xsi=dy_int*dz_int
                        A_tot=A_xsi
                    endif
                endif
                if(j>lbound(mask_s,2))then
                    if(mask_s(i,j-1,k))then
                        normal=[0._rp,1._rp,0._rp]
                        A_ysi=dx_int*dz_int
                        A_tot=A_ysi
                    endif
                endif
                if(j<ubound(mask_s,2))then
                    if(mask_s(i,j+1,k))then
                        normal=[0._rp,-1._rp,0._rp]
                        A_ysi=dx_int*dz_int
                        A_tot=A_ysi
                    endif
                endif
                if(k>lbound(mask_s,3))then
                    if(mask_s(i,j,k-1))then
                        normal=[0._rp,0._rp,1._rp]
                        A_zsi=dx_int*dy_int
                        A_tot=A_zsi
                    endif
                endif
                if(k<ubound(mask_s,3))then
                    if(mask_s(i,j,k+1))then
                        normal=[0._rp,0._rp,-1._rp]
                        A_zsi=dx_int*dy_int
                        A_tot=A_zsi
                    endif
                endif
                ! disable this for flat wall... 
                ! otherwise overestimates for curvy structure
            endif
            !lets get the inlet values
            x_in=count(.not.mask_subcell(1,:,:))*A_dx
            x_out=count(.not.mask_subcell(nsub,:,:))*A_dx
            y_in=count(.not.mask_subcell(:,1,:))*A_dy
            y_out=count(.not.mask_subcell(:,nsub,:))*A_dy
            z_in=count(.not.mask_subcell(:,:,1))*A_dz
            z_out=count(.not.mask_subcell(:,:,nsub))*A_dz
            ! surface area
            if(co>0)then
                A_x=(s_xp+s_xm)*A_dx;A_y=(s_yp+s_ym)*A_dy;A_z=(s_zp+s_zm)*A_dz;
                A_xsi=(s_xp-s_xm)*A_dx;A_ysi=(s_yp-s_ym)*A_dy;A_zsi=(s_zp-s_zm)*A_dz;
                A_tot=sqrt(A_xsi**2+A_ysi**2+A_zsi**2)
            endif
            if(co>0)call nomrlalize_vec([A_xsi,A_ysi,A_zsi],normal) ! we already fill the otherones before
            ! now lets calc the vol.
            cell_vol=dx_int*dy_int*dz_int 
            solid_vol=(real(co,rp)/nsub**3)*cell_vol
            fluid_vol=(1*cell_vol)-solid_vol
            if(inloop)then
                A_inlets(1,1)=x_in;A_inlets(2,1)=x_out;
                A_inlets(1,2)=y_in;A_inlets(2,2)=y_out;
                A_inlets(1,3)=z_in;A_inlets(2,3)=z_out;
            endif
            if(co==0)then
            !    A_tot     = 0._rp
                solid_vol = 0._rp
                fluid_vol = cell_vol
            !    normal    = 0._rp
                if(mask_s(i-1,j,k)) x_in  = 0._rp
                if(mask_s(i+1,j,k)) x_out = 0._rp

                if(mask_s(i,j-1,k)) y_in  = 0._rp
                if(mask_s(i,j+1,k)) y_out = 0._rp

                if(mask_s(i,j,k-1)) z_in  = 0._rp
                if(mask_s(i,j,k+1)) z_out = 0._rp
                A_xsi = x_out - x_in
                A_ysi = y_out - y_in
                A_zsi = z_out - z_in

                A_tot = sqrt(A_xsi**2 + A_ysi**2 + A_zsi**2)
                normal = 0._rp
                if (A_tot > tiny(1._rp)) then
                    normal = [A_xsi, A_ysi, A_zsi] / A_tot
                endif
            endif
            if(diag)then 
                print*,"cell volume is:",cell_vol
                print*,"cell fluid volume is:",fluid_vol
                print*,"cell solid volume is:",solid_vol
                print*,"total surface area in x:",A_x
                print*,"total surface area in y:",A_y
                print*,"total surface area in z:",A_z
                print*,"signed surface area in x:",A_xsi
                print*,"signed surface area in y:",A_ysi
                print*,"signed surface area in z:",A_zsi
                print*,"approx. surface area:",A_tot
                print*,"inlet area in x:",x_in
                print*,"outlet area in x:",x_out
                print*,"inlet area in y:",y_in
                print*,"outlet area in y:",y_out
                print*,"inlet area in z:",z_in
                print*,"outlet area in z:",z_out
                print*,"normal is x,y,z:",normal
            endif
            output(i,j,k,1)=cell_vol
            output(i,j,k,2)=fluid_vol
            output(i,j,k,3)=solid_vol
            output(i,j,k,4)=A_x
            output(i,j,k,5)=A_y
            output(i,j,k,6)=A_z
            output(i,j,k,7)=A_xsi
            output(i,j,k,8)=A_ysi
            output(i,j,k,9)=A_zsi
            output(i,j,k,10)=A_tot
            output(i,j,k,11)=x_in
            output(i,j,k,12)=x_out
            output(i,j,k,13)=y_in
            output(i,j,k,14)=y_out
            output(i,j,k,15)=z_in
            output(i,j,k,16)=z_out
            output(i,j,k,17)=normal(1)
            output(i,j,k,18)=normal(2)
            output(i,j,k,19)=normal(3)            
        endif
    end subroutine get_fric_cubes
    subroutine init_wall_force(n,lo,l,dl,zc,zf,dzc,dzf,mask_u,mask_v,mask_w,mask_s,ibm_direction,amp_l,n_wave,l_0,phase_l,&
                               hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
        !
        ! set-up of the exact wall force per cell (`wf_gather`, `wf_scatter`, `calc_wall_force`):
        !  - wf_cell: the cells next to the fluid/solid interface (both kinds in their 3x3x3 neighbourhood)
        !  - per velocity component, the points next to the interface: band points (fluid, with a solid
        !    neighbour) and solid points with a fluid neighbour; their impulse is split over the wall
        !    directions with the weights lambda_d/lambda (band) or 1/h_d^2 (solid), and each direction
        !    goes to the cell that contains the wall point between the two nodes (`calc_lambda`);
        !    a wall point on a face normal to the staggered direction is split between the two cells
        !
        integer , intent(in), dimension(3)        :: n,lo
        real(rp), intent(in), dimension(3)        :: l,dl
        real(rp), intent(in), dimension(0:)       :: zc,zf,dzc,dzf
        logical , intent(in), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w,mask_s
        logical , intent(in), dimension(0:1,3)    :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)    :: amp_l
        integer , intent(in), dimension(0:1,3)    :: n_wave
        real(rp), intent(in), dimension(0:1,3)    :: l_0
        real(rp), intent(in), dimension(0:1,3)    :: phase_l
        real(rp), intent(in), optional            :: hmap(0:,0:)
        real(rp), intent(in), optional            :: l1_hmap,l2_hmap
        integer , intent(in), optional            :: n1_hmap,n2_hmap
        logical , allocatable, dimension(:,:,:)   :: mask
        integer  :: c,i,j,k,d,ax,a,ip,nt,npmax,own(3),idx,ie(3),dix,diy,diz
        real(rp) :: x,y,z,xp,xm,yp,ym,zp,zm,pos(3),xn(3),lam,wsum,wd(6),wl(3,6),dzf_l,dzc_l,wcoord
        logical  :: has(6),is_solid,is_split,in_cell(-1:1)
        integer  :: sgn
        real(rp) :: xs(3),xs2(3)
        allocate(wf_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        wf_cell(:,:,:) = .false.
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    wf_cell(i,j,k) = any(mask_s(i-1:i+1,j-1:j+1,k-1:k+1)).and..not.all(mask_s(i-1:i+1,j-1:j+1,k-1:k+1))
                end do
            end do
        end do
        npmax = 0
        do c=1,3
            select case(c)
            case(1); mask = mask_u
            case(2); mask = mask_v
            case(3); mask = mask_w
            end select
            npmax = max(npmax,count(is_iface(mask,n)))
        end do
        allocate(wf_ijk(3,npmax,3),wf_nt(npmax,3),wf_tc(3,wf_ntmax,npmax,3),wf_tw(wf_ntmax,npmax,3))
        wf_ijk(:,:,:) = 1; wf_nt(:,:) = 0; wf_tc(:,:,:,:) = 1; wf_tw(:,:,:) = 0._rp
        allocate(wf_imp(0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
        wf_imp(:,:,:,:) = 0._rp
        do c=1,3
            ie(:) = 0; ie(c) = 1
            dix = ie(1); diy = ie(2); diz = ie(3)
            select case(c)
            case(1); mask = mask_u
            case(2); mask = mask_v
            case(3); mask = mask_w
            end select
            ip = 0
            do k=1,n(3)
                do j=1,n(2)
                    do i=1,n(1)
                        is_solid = mask(i,j,k)
                        has(:) = .false.; wd(:) = 0._rp; wl(:,:) = 0._rp
                        call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                        pos(:) = [x,y,z]
                        do d=1,6
                            ax = (d+1)/2
                            xn(:) = pos(:)
                            select case(d)
                            case(1); xn(1) = xp; has(d) = mask(i+1,j,k).neqv.is_solid
                            case(2); xn(1) = xm; has(d) = mask(i-1,j,k).neqv.is_solid
                            case(3); xn(2) = yp; has(d) = mask(i,j+1,k).neqv.is_solid
                            case(4); xn(2) = ym; has(d) = mask(i,j-1,k).neqv.is_solid
                            case(5); xn(3) = zp; has(d) = mask(i,j,k+1).neqv.is_solid
                            case(6); xn(3) = zm; has(d) = mask(i,j,k-1).neqv.is_solid
                            end select
                            if(.not.has(d)) cycle
                            dzf_l = 0._rp; dzc_l = 0._rp ! as in `set_ibm_2nd`
                            if(d == 5) then
                                if(diz == 0) then
                                    dzc_l = dzc(k); dzf_l = dzf(k)
                                else
                                    dzf_l = dzf(k+1); dzc_l = dzc(k)
                                end if
                            else if(d == 6) then
                                if(diz == 0) then
                                    dzc_l = dzc(k-1); dzf_l = dzf(k)
                                else
                                    dzf_l = dzf(k); dzc_l = dzc(k)
                                end if
                            end if
                            if(.not.is_solid) then
                                call calc_lambda(pos(1),pos(2),pos(3),xn(ax),ax,lam,ibm_direction,amp_l,n_wave,l_0,phase_l, &
                                                 n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wl(:,d))
                                wd(d) = lam
                            else
                                call calc_lambda(xn(1),xn(2),xn(3),pos(ax),ax,lam,ibm_direction,amp_l,n_wave,l_0,phase_l, &
                                                 n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wl(:,d))
                                wd(d) = 1._rp/(xn(ax)-pos(ax))**2
                            end if
                        end do
                        if(.not.any(has)) cycle
                        ip = ip+1
                        wf_ijk(:,ip,c) = [i,j,k]
                        wsum = sum(wd(:))
                        if(wsum <= 0._rp) then ! no lambda (wall on the node): equal weights
                            wd(:) = merge(1._rp,0._rp,has(:)); wsum = sum(wd(:))
                        end if
                        nt = 0
                        do d=1,6
                            if(.not.has(d)) cycle
                            ax = (d+1)/2
                            own(:) = [i,j,k]
                            is_split = .false.
                            do a=1,3
                                if(a == c) then
                                    if(ax == c) then
                                        if(mod(d,2) == 1) own(a) = own(a)+1 ! staggered point between cells own and own+1
                                    else
                                        is_split = .true.
                                    end if
                                else if(a == ax) then
                                    wcoord = wl(a,d)
                                    if(a < 3) then
                                        idx = ceiling(wcoord/dl(a)) - lo(a) + 1
                                        own(a) = min(max(idx,own(a)-1),own(a)+1)
                                    else
                                        if(wcoord >  zf(own(3)  )) own(3) = own(3)+1
                                        if(wcoord <= zf(own(3)-1)) own(3) = own(3)-1
                                    end if
                                end if
                            end do
                            if(.not.is_split) then
                                nt = nt+1
                                wf_tc(:,nt,ip,c) = own(:)
                                wf_tw(nt,ip,c) = wd(d)/wsum
                            else
                                !
                                ! the wall point is on the face between cells own and own+1 (along c):
                                ! give it to the cell(s) in which the wall crosses the segment between
                                ! the two nodes, shifted by a quarter cell off the face
                                !
                                in_cell(:) = .false.
                                do sgn=-1,1,2
                                    xs(:) = pos(:); xs(c) = pos(c) + sgn*0.25_rp*h_c(c,k,sgn)
                                    xs2(:) = xn(:); xs2(c) = xs(c)
                                    in_cell(sgn) = isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xs(1),xs(2),xs(3),n,l, &
                                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap).neqv. &
                                                   isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xs2(1),xs2(2),xs2(3),n,l, &
                                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                end do
                                if(.not.any(in_cell)) in_cell(-1:1:2) = .true. ! wall along the face
                                do sgn=-1,1,2
                                    if(.not.in_cell(sgn)) cycle
                                    nt = nt+1
                                    wf_tc(:,nt,ip,c) = own(:)
                                    if(sgn == 1) wf_tc(c,nt,ip,c) = own(c)+1
                                    wf_tw(nt,ip,c) = wd(d)/wsum/count(in_cell)
                                end do
                            end if
                        end do
                        wf_nt(ip,c) = nt
                    end do
                end do
            end do
            wf_np(c) = ip
        end do
        !$acc enter data copyin(wf_cell,wf_ijk,wf_nt,wf_tc,wf_tw,wf_imp)
    contains
        real(rp) function h_c(c,k,sgn)
            ! width, along c, of the cell on side sgn of the staggered point (cells k, k+1 along z)
            integer, intent(in) :: c,k,sgn
            if(c < 3) then
                h_c = dl(c)
            else
                h_c = dzf(k+(sgn+1)/2)
            end if
        end function h_c
        function is_iface(mask,n) result(res)
            logical, intent(in), dimension(0:,0:,0:) :: mask
            integer, intent(in), dimension(3)        :: n
            logical, dimension(n(1),n(2),n(3))       :: res
            integer :: i,j,k
            do k=1,n(3)
                do j=1,n(2)
                    do i=1,n(1)
                        res(i,j,k) = any(mask(i-1:i+1:2,j,k).neqv.mask(i,j,k)) .or. &
                                     any(mask(i,j-1:j+1:2,k).neqv.mask(i,j,k)) .or. &
                                     any(mask(i,j,k-1:k+1:2).neqv.mask(i,j,k))
                    end do
                end do
            end do
        end function is_iface
    end subroutine init_wall_force
    !
    subroutine wf_gather(c,n,dli,dzfi,grid_vol_ratio_f,mask,fac,p,pref,fs)
        !
        ! cell part of the impulse the IB removes at the solid points of component c (all of them,
        ! also deep in the solid): the pressure part -fac*grad(p) of each solid point, telescoped onto
        ! the cells (a cell gets fac*(p-pref)*area from a solid face on its + side and minus that from
        ! one on its - side, so it cancels inside the solid and is left on the cut cells), and the
        ! forcing part fs (impulse per unit volume), split over the two half-cells of the point;
        ! units as `ibm_impulse` (per unit box volume)
        !
        integer , intent(in)                      :: c
        integer , intent(in), dimension(3)        :: n
        real(rp), intent(in), dimension(3)        :: dli
        real(rp), intent(in), dimension(0:)       :: dzfi,grid_vol_ratio_f
        logical , intent(in), dimension(0:,0:,0:) :: mask
        real(rp), intent(in)                      :: fac,pref,fs
        real(rp), intent(in), dimension(0:,0:,0:) :: p
        real(rp) :: sp,sm,area
        integer  :: i,j,k,di,dj,dk
        di = 0; dj = 0; dk = 0
        select case(c)
        case(1); di = 1
        case(2); dj = 1
        case default; dk = 1
        end select
        !$acc parallel loop collapse(3) default(present) private(sp,sm,area) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(sp,sm,area)
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    if(wf_cell(i,j,k)) then
                        sp = merge(1._rp,0._rp,mask(i,j,k))
                        sm = merge(1._rp,0._rp,mask(i-di,j-dj,k-dk))
                        area = (di*dli(1)+dj*dli(2)+dk*dzfi(k))*grid_vol_ratio_f(k) ! face area/V_box
                        wf_imp(i,j,k,c) = wf_imp(i,j,k,c) + fac*(p(i,j,k)-pref)*area*(sp-sm) &
                                                          + fs*.5_rp*grid_vol_ratio_f(k)*(sp+sm)
                    end if
                end do
            end do
        end do
    end subroutine wf_gather
    !
    subroutine wf_scatter(part,c,dli,dzci,grid_vol_ratio,mask,a,b,fac,bf,phi,vel)
        !
        ! point part of the impulse the IB removes from component c, assigned to the wall cells
        ! (targets from `init_wall_force`), same parts as `ibm_impulse` ('proj': damped correction):
        !  - band points: the lambda part, 'pre' (1-b)*u, 'post'/'impl' (a-1)*u, 'proj' -fac*grad(phi)*(1-1/a)
        !  - solid points next to the fluid: the local fluxes, i.e. the impulse minus its pressure and
        !    forcing parts (`wf_gather`): 'leak' u + fac*grad(phi) (phi, fac of the previous stage),
        !    'post'/'post_s' a*u - fac*(bf - grad(phi)) (phi = p), 'impl' (ibm_pen-1)*u
        ! deeper solid points have no flux part
        !
        character(len=*), intent(in)              :: part
        integer , intent(in)                      :: c
        real(rp), intent(in), dimension(3)        :: dli
        real(rp), intent(in), dimension(0:)       :: dzci,grid_vol_ratio
        logical , intent(in), dimension(0:,0:,0:) :: mask
        real(rp), intent(in), dimension(0:,0:,0:) :: a,b,phi,vel
        real(rp), intent(in)                      :: fac,bf
        real(rp) :: q,g
        integer  :: ip,it,i,j,k,ii,jj,kk,ipart,di,dj,dk,np
        np = wf_np(c)
        select case(part)
        case('leak');   ipart = 1
        case('pre');    ipart = 2
        case('post');   ipart = 3
        case('post_s'); ipart = 4
        case('impl');   ipart = 5
        case default;   ipart = 6 ! 'proj'
        end select
        di = 0; dj = 0; dk = 0
        select case(c)
        case(1); di = 1
        case(2); dj = 1
        case default; dk = 1
        end select
        !$acc parallel loop default(present) private(q,g,i,j,k,ii,jj,kk) async(1)
        !$OMP PARALLEL DO   DEFAULT(shared)  PRIVATE(q,g,i,j,k,ii,jj,kk,it)
        do ip=1,np
            i = wf_ijk(1,ip,c); j = wf_ijk(2,ip,c); k = wf_ijk(3,ip,c)
            g = 0._rp
            if(ipart /= 2 .and. ipart /= 5 .and. fac /= 0._rp) &
              g = (phi(i+di,j+dj,k+dk)-phi(i,j,k))*(di*dli(1)+dj*dli(2)+dk*dzci(k))
            q = 0._rp
            if(mask(i,j,k)) then
                select case(ipart)
                case(1);   q = vel(i,j,k) + fac*g
                case(3,4); q = a(i,j,k)*vel(i,j,k) - fac*(bf - g)
                case(5);   q = (ibm_pen-1._rp)*vel(i,j,k)
                end select
            else
                select case(ipart)
                case(2);   q = (1._rp-b(i,j,k))*vel(i,j,k)
                case(3,5); q = (a(i,j,k)-1._rp)*vel(i,j,k)
                case(6);   q = -fac*g*(1._rp-1._rp/a(i,j,k))
                end select
            end if
            q = q*grid_vol_ratio(k)
            !$acc loop seq
            do it=1,wf_nt(ip,c)
                ii = wf_tc(1,it,ip,c); jj = wf_tc(2,it,ip,c); kk = wf_tc(3,it,ip,c)
#if defined(_OPENACC)
                !$acc atomic update
#else
                !$OMP ATOMIC UPDATE
#endif
                wf_imp(ii,jj,kk,c) = wf_imp(ii,jj,kk,c) + wf_tw(it,ip,c)*q
            end do
        end do
    end subroutine wf_scatter
    !
    subroutine wf_halo_add(n)
        !
        ! adds the impulse accumulated in the halo cells to the cells of the neighbouring tasks
        ! (periodic wrap along the pencil axis), then clears the halos; host data
        !
        use mpi
        integer, intent(in), dimension(3) :: n
        real(rp), allocatable, dimension(:,:,:,:) :: sb,rb
        integer :: idir,ierr,nn,nbr(0:1)
        do idir=1,3
            select case(idir)
            case(1); allocate(sb(1,0:n(2)+1,0:n(3)+1,3),rb(1,0:n(2)+1,0:n(3)+1,3))
            case(2); allocate(sb(0:n(1)+1,1,0:n(3)+1,3),rb(0:n(1)+1,1,0:n(3)+1,3))
            case(3); allocate(sb(0:n(1)+1,0:n(2)+1,1,3),rb(0:n(1)+1,0:n(2)+1,1,3))
            end select
            nn = size(sb)
            if(idir == ipencil_axis) then
                if(cbcpre(0,idir)//cbcpre(1,idir) == 'PP') then
                    call wf_add_plane(idir,n(idir),wf_get_plane(idir,0))
                    call wf_add_plane(idir,1      ,wf_get_plane(idir,n(idir)+1))
                end if
            else
                nbr(:) = merge(nb(:,idir),MPI_PROC_NULL,nb(:,idir) >= 0) ! cuDecomp marks a missing neighbour by -1
                sb = wf_get_plane(idir,n(idir)+1); rb = 0._rp
                call MPI_SENDRECV(sb,nn,MPI_REAL_RP,nbr(1),0,rb,nn,MPI_REAL_RP,nbr(0),0, &
                                  MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                call wf_add_plane(idir,1,rb)
                sb = wf_get_plane(idir,0); rb = 0._rp
                call MPI_SENDRECV(sb,nn,MPI_REAL_RP,nbr(0),1,rb,nn,MPI_REAL_RP,nbr(1),1, &
                                  MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                call wf_add_plane(idir,n(idir),rb)
            end if
            select case(idir)
            case(1); wf_imp(0,:,:,:) = 0._rp; wf_imp(n(1)+1,:,:,:) = 0._rp
            case(2); wf_imp(:,0,:,:) = 0._rp; wf_imp(:,n(2)+1,:,:) = 0._rp
            case(3); wf_imp(:,:,0,:) = 0._rp; wf_imp(:,:,n(3)+1,:) = 0._rp
            end select
            deallocate(sb,rb)
        end do
    contains
        function wf_get_plane(idir,m) result(pl)
            integer, intent(in) :: idir,m
            real(rp), allocatable, dimension(:,:,:,:) :: pl
            select case(idir)
            case(1); pl = wf_imp(m:m,:,:,:)
            case(2); pl = wf_imp(:,m:m,:,:)
            case(3); pl = wf_imp(:,:,m:m,:)
            end select
        end function wf_get_plane
        subroutine wf_add_plane(idir,m,pl)
            integer, intent(in) :: idir,m
            real(rp), intent(in), dimension(:,:,:,:) :: pl
            select case(idir)
            case(1); wf_imp(m:m,:,:,:) = wf_imp(m:m,:,:,:) + pl
            case(2); wf_imp(:,m:m,:,:) = wf_imp(:,m:m,:,:) + pl
            case(3); wf_imp(:,:,m:m,:) = wf_imp(:,:,m:m,:) + pl
            end select
        end subroutine wf_add_plane
    end subroutine wf_halo_add
    !
    subroutine calc_wall_force(n,lo,l,dti,pi_drive,mask_s,cubeInput,output)
        !
        ! exact wall force per cell for one time step, from the IB impulse accumulated in wf_imp
        ! over its stages: F = imp*V_box/dt - pi_drive*V_solid (the driving force on the solid part
        ! of the cell is not a wall force; V_solid from the sub-cell sampling of `get_fric_cubes`).
        ! Summed over the cells it is the total IB force minus the driving force on the solid, i.e.
        ! the drag of `ibm_force.out`/`ibm_wall_stress_*.out`. F is the force of the fluid on the
        ! wall (positive = drag); its component along the normal n (out of the solid) is reported as
        ! the normal (pressure) part, the rest as the tangential (shear) part. A cell holds the force
        ! of the nodes next to the wall portion it contains, so single-cell values are local but
        ! not a pointwise traction. Clears wf_imp.
        !
        integer , intent(in), dimension(3)         :: n,lo
        real(rp), intent(in), dimension(3)         :: l,pi_drive
        real(rp), intent(in)                       :: dti
        logical , intent(in), dimension(0:,0:,0:)  :: mask_s
        real(rp), intent(in), dimension(0:,0:,0:,:):: cubeInput
        real(rp), intent(out), dimension(:,:)      :: output
        real(rp) :: A_tot,normal(3),F(3),f_n,f_nv(3),f_tv(3),f_t,r(3)
        logical  :: nbo(6)
        integer  :: i,j,k,nc
        !$acc wait(1)
        !$acc update self(wf_imp)
        call wf_halo_add(n)
        !
        ! a cell without wall area (by the sub-cell sampling, e.g. a solid cell whose wall lies
        ! within the sampling distance of a face) passes its force to the face neighbours of
        ! the other kind (solid/fluid centre), which hold that wall face
        !
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    if(.not.wf_cell(i,j,k) .or. cubeInput(i,j,k,10) > tiny(1._rp)) cycle
                    nbo(:) = [mask_s(i+1,j,k),mask_s(i-1,j,k),mask_s(i,j+1,k),mask_s(i,j-1,k), &
                              mask_s(i,j,k+1),mask_s(i,j,k-1)].neqv.mask_s(i,j,k)
                    if(.not.any(nbo)) cycle
                    r(:) = (wf_imp(i,j,k,:) - pi_drive(:)*cubeInput(i,j,k,3)/(product(l(:))*dti))/count(nbo)
                    if(nbo(1)) wf_imp(i+1,j,k,:) = wf_imp(i+1,j,k,:) + r(:)
                    if(nbo(2)) wf_imp(i-1,j,k,:) = wf_imp(i-1,j,k,:) + r(:)
                    if(nbo(3)) wf_imp(i,j+1,k,:) = wf_imp(i,j+1,k,:) + r(:)
                    if(nbo(4)) wf_imp(i,j-1,k,:) = wf_imp(i,j-1,k,:) + r(:)
                    if(nbo(5)) wf_imp(i,j,k+1,:) = wf_imp(i,j,k+1,:) + r(:)
                    if(nbo(6)) wf_imp(i,j,k-1,:) = wf_imp(i,j,k-1,:) + r(:)
                    wf_imp(i,j,k,:) = wf_imp(i,j,k,:) - count(nbo)*r(:)
                end do
            end do
        end do
        call wf_halo_add(n)
        nc = 0
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    A_tot = cubeInput(i,j,k,10)
                    if(A_tot <= tiny(1._rp)) cycle
                    nc = nc+1
                    normal(:) = cubeInput(i,j,k,17:19)
                    F(:) = wf_imp(i,j,k,:)*product(l(:))*dti - pi_drive(:)*cubeInput(i,j,k,3)
                    f_n = dot_product(normal,F)
                    f_nv(:) = f_n*normal(:)
                    f_tv(:) = F(:)-f_nv(:)
                    f_t = norm2(f_tv)
                    output(nc,1:3)   = real([lo(1)+i-1,lo(2)+j-1,lo(3)+k-1],rp)
                    output(nc,4)     = A_tot
                    output(nc,5:7)   = normal(:)
                    output(nc,8:10)  = F(:)
                    output(nc,11)    = f_n
                    output(nc,12:14) = f_nv(:)
                    output(nc,15)    = f_t
                    output(nc,16:18) = f_tv(:)
                    output(nc,19)    = f_n/A_tot
                    output(nc,20:22) = f_nv(:)/A_tot
                    output(nc,23)    = f_t/A_tot
                    output(nc,24:26) = f_tv(:)/A_tot
                end do
            end do
        end do
        wf_imp(:,:,:,:) = 0._rp
        !$acc update device(wf_imp)
    end subroutine calc_wall_force
    !
    subroutine out_ibm_wall_stress(fname,n,lo,ng,l,dl,pi_drive,imp_map)
        !
        ! local wall stress of the lower wall in the plan view, from the IB force of one
        ! time step integrated over each (i,j) column of the lower half (imp_map, per unit
        ! box volume and time, see `ibm_impulse`) minus the driving force on the solid
        ! part of the column: tau_x(i,j) = F_x(i,j)/(dx*dy) - pi_x*h_s(x,y), h_s being the
        ! wall height above z = 0 (not the staircase height: the band cells hold the cut
        ! cell's solid part); tau_x at (xf(i),yc(j)), tau_y at (xc(i),yf(j)); z-walls only
        !
        use mpi
        implicit none
        character(len=*), intent(in)                 :: fname
        integer , intent(in), dimension(3)           :: n,lo,ng
        real(rp), intent(in), dimension(3)           :: l,dl,pi_drive
        real(rp), intent(in), dimension(:,:,:)       :: imp_map
        real(rp), allocatable, dimension(:,:,:)      :: tau
        integer :: i,j,ig,jg,iunit,ierr,myid
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        allocate(tau(ng(1),ng(2),2))
        tau(:,:,:) = 0._rp
        if(lo(3) == 1) then ! one rank per (i,j) column subtracts the solid part
            do j=1,n(2)
                jg = lo(2)+j-1
                do i=1,n(1)
                    ig = lo(1)+i-1
                    tau(ig,jg,1) = -pi_drive(1)*wall_height(0,3,ig*dl(1),(jg-0.5_rp)*dl(2),amp_l,n_wave,l_0,phase_l,l, &
                                                            hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha)
                    tau(ig,jg,2) = -pi_drive(2)*wall_height(0,3,(ig-0.5_rp)*dl(1),jg*dl(2),amp_l,n_wave,l_0,phase_l,l, &
                                                            hmap_tha,lx_tha,ly_tha,nx_hmap_tha,ny_hmap_tha)
                end do
            end do
        end if
        do j=1,n(2)
            jg = lo(2)+j-1
            do i=1,n(1)
                ig = lo(1)+i-1
                tau(ig,jg,:) = tau(ig,jg,:) + imp_map(i,j,1:2)*product(l(:))/(dl(1)*dl(2))
            end do
        end do
        if(myid == 0) then
            call MPI_REDUCE(MPI_IN_PLACE,tau,size(tau),MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            open(newunit=iunit,file=fname,status='replace',action='write')
            do j=1,ng(2)
                do i=1,ng(1)
                    write(iunit,'(2i8,2(1x,es24.16e3))') i,j,tau(i,j,1),tau(i,j,2)
                end do
            end do
            close(iunit)
        else
            call MPI_REDUCE(tau,tau,size(tau),MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        end if
    end subroutine out_ibm_wall_stress
    subroutine writeFVfric(fname,output)
        implicit none
        real(rp),dimension(:,:),intent(in)       :: output
        character(len=*)                         :: fname   
        integer                                  :: iunit,wr
        character(len=*), parameter              :: fmt_dp = '(*(es24.16e3,1x))', &
                                                    fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
        character(len=*), parameter              :: fmt_rp = fmt_dp
#else
        character(len=*), parameter              :: fmt_rp = fmt_sp
#endif
        open(newunit=iunit,file=fname,status='replace',action='write')
            do wr=1,ubound(output,1)
                write(iunit,fmt_rp) output(wr,1:26)
            end do
        close(iunit)
    end subroutine writeFVfric
    !
    !CaNS in CaNS 

    !
    !
    !
    subroutine write_data(fname,myid,grad_id)
        real(rp),dimension(:,:),intent(in)       :: grad_id
        integer,intent(in)                       :: myid
        character(len=*)                         :: fname   
        integer                                  :: iunit,wr
        character(len=*), parameter              :: fmt_dp = '(*(es24.16e3,1x))', &
                                                    fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
        character(len=*), parameter              :: fmt_rp = fmt_dp
#else
        character(len=*), parameter              :: fmt_rp = fmt_sp
#endif
        if(myid == 0) then
        open(newunit=iunit,file=fname)
            do wr=1,ubound(grad_id,1)
                write(iunit,fmt_rp) grad_id(wr,1),grad_id(wr,2),grad_id(wr,3),&
                                    grad_id(wr,4),grad_id(wr,5),grad_id(wr,6),&
                                    grad_id(wr,7),grad_id(wr,8),grad_id(wr,9),&
                                    grad_id(wr,10),grad_id(wr,11),grad_id(wr,12)
            end do
        close(iunit)
        end if
    end subroutine write_data
    subroutine apply_ibm_staircase(field,mask_id,dt)
        implicit none
        real(rp),intent(inout),dimension(0:,0:,0:)  :: field
        logical,intent(in),dimension(0:,0:,0:)      :: mask_id
        real(rp),intent(in)                         :: dt
        integer :: i,j,k
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k = lbound(field,3)+1,ubound(field,3)-1
            do j = lbound(field,2)+1,ubound(field,2)-1
                do i = lbound(field,1)+1,ubound(field,1)-1
                    if (mask_id(i,j,k).eqv..true.)then
                        field(i,j,k) = 0._rp 
                    endif
                end do 
            end do 
        end do
    end subroutine apply_ibm_staircase
    subroutine apply_ibm_staircase_scalar(scalar_field,scalar_mask_id,scalar_bc)
        implicit none
        real(rp),intent(inout),dimension(0:,0:,0:)  :: scalar_field
        logical,intent(in),dimension(0:,0:,0:)      :: scalar_mask_id
        real(rp),intent(in)                         :: scalar_bc
        integer :: i,j,k
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k = lbound(scalar_field,3)+1,ubound(scalar_field,3)-1
            do j = lbound(scalar_field,2)+1,ubound(scalar_field,2)-1
                do i = lbound(scalar_field,1)+1,ubound(scalar_field,1)-1
                    if (scalar_mask_id(i,j,k).eqv..true.)then
                        scalar_field(i,j,k) = scalar_bc 
                    endif
                end do 
            end do 
        end do
    end subroutine apply_ibm_staircase_scalar
end module mod_ibm