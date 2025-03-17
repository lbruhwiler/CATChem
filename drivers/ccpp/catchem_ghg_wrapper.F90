!>\file catchem_ghg_wrapper.F90
!! This file is GSDChem ghg wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 05/2020
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov
!! 03/2025, Adapt for GHGs, lori.bruhwiler@noaa.gov

 module catchem_ghg_wrapper

   use mpi      
   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config   only : num_chem, num_moist, epsilc

   implicit none

   private

   public :: catchem_ghg_wrapper_init, catchem_ghg_wrapper_run, catchem_ghg_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_ghg_wrapper_init()
      end subroutine catchem_ghg_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_ghg_wrapper_finalize Argument Table
!!
      subroutine catchem_ghg_wrapper_finalize()
      end subroutine catchem_ghg_wrapper_finalize

!> \defgroup catchem_ghg_group CATChem ghg wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_ghg_wrapper CATChem ghg wrapper Module  
!> \ingroup catchem_ghg_group
!! This is the CATChem ghg wrapper Module
!! \section arg_table_catchem_ghg_wrapper_run Argument Table
!! \htmlinclude catchem_ghg_wrapper_run.html
!!
!>\section catchem_ghg_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_ghg_wrapper_run(im, kte, kme, ktau, dt, garea,rlat, rlon, julian, idat, pgr,  &
                   pr3d, ph3d, prl3d, tk3d, spechum, ghgi_in, ch4chm1, ch4chm2, ch4chm3,ch4loss,  &
                   ntrac, ntco2,ntco2_bgd,ntco2_land,ntco2_fossil,ntco2_fire,ntco2_ocn,ntch4, ntsf6, &
                   ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, gq0, qgrs, errmsg, errflg)


    implicit none

    integer,        intent(in) :: im,kte,kme,ktau,idat(8)
    integer,        intent(in) :: ntrac
    integer,        intent(in) :: ntco2,ntco2_bgd,ntco2_land,ntco2_fossil,ntco2_fire,ntco2_ocn
    integer,        intent(in) :: ntch4,ntsf6,ntqv,ntcw,ntiw,ntrw,ntsw,ntgl
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1
    integer, parameter :: num_ghg=3,p_co2=1,p_ch4=2,p_sf6=3
    real,    parameter :: mw_co2=44.01, mw_ch4=16.04, mw_sf6=146.06, mw_dry=28.97

    real(kind_phys), dimension(im,64, 3),     intent(in)    :: ghgi_in
    real(kind_phys), dimension(im,90, 2),     intent(in)    :: ch4chm1
    real(kind_phys), dimension(im,31, 2),     intent(in)    :: ch4chm3
    real(kind_phys), dimension(im,60, 2),     intent(in)    :: ch4chm2
    real(kind_phys), dimension(im,2),         intent(out)   :: ch4loss
    real(kind_phys), dimension(im),           intent(in)    :: garea, rlat,rlon, pgr
    real(kind_phys), dimension(im,kme),       intent(in)    :: ph3d, pr3d
    real(kind_phys), dimension(im,kte),       intent(in)    :: prl3d, tk3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
!    integer,                                  intent(in)    :: chem_in_opt
    character(len=*),                         intent(out)   :: errmsg
    integer,                                  intent(out)   :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme)         :: rri, t_phy, p_phy, dz8w, p8w, t8w, rho_phy
    real(kind_phys), dimension(ims:im, jms:jme)             :: xlat, xlong, dxy


!>- vapor & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

    integer :: ide, ime, ite, kde, julday

!>- ghg variables
    !integer, parameter :: chem_in_opt = 0  ! 0 for coldstart, 1 for restart
    logical, parameter :: readrestart = .false.
    integer, parameter :: nvl_gocart  = 64  ! number of input levels from gocart file

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: pm10, pm2_5_dry, pm2_5_dry_ec

!>- chemical background variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: backg_oh

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: oh_t
    real(kind_phys), dimension(ims:im, jms:jme) :: ttday, tcosz

    real(kind_phys) :: dtstep, gmt
    real(kind_phys), dimension(1:num_ghg) :: ppm2ugkg

    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25!, ebu_oc
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: oh_bg, h2o2_bg, no3_bg

    !-- column masses for checking mass conservation

    real(kind_phys) :: colghg


    integer :: current_month
    real(kind_phys) :: dtstep
    real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!>-- local variables
    logical :: call_ch4chem
    integer :: i, j, jp, k, kp, n
    integer :: numprocs,nproc,ierr,tag
    integer :: status(MPI_STATUS_SIZE)
    real    :: totmass,q_corr


    errmsg = ''
    errflg = 0

    gmt = real(idat(5))
    julday = real(julian)


    chem_opt          = chem_opt_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
!    ppm2ugkg(p_co2) = 1.e+03_kind_phys * mw_co2 / mwdry ! should it be ug/kg or kg/kg????
!    ppm2ugkg(p_ch4) = 1.e+03_kind_phys * mw_ch4 / mwdry
!    ppm2ugkg(p_sf6) = 1.e+03_kind_phys * mw_sf6 / mwdry
    ppm2ugkg(p_co2) = 1._kind_phys * mw_co2 / mw_dry
    ppm2ugkg(p_ch4) = 1._kind_phys * mw_ch4 / mw_dry
    ppm2ugkg(p_sf6) = 1._kind_phys * mw_sf6 / mw_dry


    ! frequency of call to ch4chem can be changed using mod(ktau,1) - e.g.
    ! change 1 to 2

!!!    call_ch4chem      = (mod(ktau, 1) == 0) .or. (ktau == 1)

! this code doesn't work because intel converts "T" to -1, take absolute val.

!!!    if (ktau > 1) then
!!!      dtstep = abs(call_ch4chem) * dt
!!!    else
      dtstep = dt
!!!    end if

   !print *,'ntco2:',ntco2
   !print *,'ntsf6:',ntsf6
   !print *,'ntch4:',ntch4
   !print *,'ntqv:',ntqv
   !stop

! qg0 is the tracer concentration to be updated by physics

!    if (chem_in_opt.eq.0) then
     if (ktau.eq.1) then
      do k=kts,kte
!       do i=its,ite
       do i = 1,im
!        gq0(i,k,ntco2  )=ppm2ugkg(p_co2   ) * max(epsilc,ghgi_in(i,k,1))
!        gq0(i,k,ntch4  )=ppm2ugkg(p_ch4   ) * max(epsilc,ghgi_in(i,k,2))
!        gq0(i,k,ntsf6  )=ppm2ugkg(p_sf6   ) * max(epsilc,ghgi_in(i,k,3))
!correction factor to retain constant field (incoming is wr to dry air, so m.r.s
!will get smaller if moisture or condensates are present)
        q_corr=1.-gq0(i,k,ntqv)-gq0(i,k,ntcw)-gq0(i,k,ntiw)-gq0(i,k,ntrw)-gq0(i,k,ntsw)-gq0(i,k,ntgl)
        gq0(i,k,ntco2  )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*1.e-6*q_corr
        gq0(i,k,ntco2_bgd  )= ppm2ugkg(p_co2) * 400.*1.e-6*q_corr
!        gq0(i,k,ntco2_bgd  )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*q_corr
!        gq0(i,k,ntco2_land  )= ppm2ugkg(p_co2) * 400.*1.e-6*qcorr
        gq0(i,k,ntco2_land  )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*1.e-6*q_corr
!        gq0(i,k,ntco2_fossil  )= ppm2ugkg(p_co2) * 400.*1.e-6*qcorr
        gq0(i,k,ntco2_fossil  )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*1.e-6*q_corr
!        gq0(i,k,ntco2_fire  )= ppm2ugkg(p_co2) * 400.*1.e-6*qcorr
        gq0(i,k,ntco2_fire  )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*1.e-6*q_corr
!        gq0(i,k,ntco2_ocn  )= ppm2ugkg(p_co2) * 400.*1.e-6*qcorr
        gq0(i,k,ntco2_ocn )=ppm2ugkg(p_co2   ) * ghgi_in(i,k,1)*1.e-6*q_corr
        gq0(i,k,ntch4  )=ppm2ugkg(p_ch4   ) * ghgi_in(i,k,2)*1.e-9*q_corr
!        gq0(i,k,ntch4  )=ghgi_in(i,k,2)
        gq0(i,k,ntsf6  )=ppm2ugkg(p_sf6   ) * ghgi_in(i,k,3)*q_corr
!        gq0(i,k,ntsf6  )=ghgi_in(i,k,3)
       enddo
      enddo
!      print *,'initializing co2:',gq0(1,1,ntco2),ktau,epsilc,q_corr
!      print *,'initializing co2 bg:',gq0(1,1,ntco2_bgd),ktau,epsilc,q_corr
!      print *,'initializing co2 land:',gq0(1,1,ntco2_land),ktau,epsilc,q_corr
!      print *,'initializing co2 fossil:',gq0(1,1,ntco2_fossil),ktau,epsilc,q_corr
!      print *,'initializing co2 fire:',gq0(1,1,ntco2_fire),ktau,epsilc,q_corr
!      print *,'initializing co2 ocn:',gq0(1,1,ntco2_ocn),ktau,epsilc,q_corr
!      print *,'initializing ch4:',gq0(1,1,ntch4),ktau,epsilc,q_corr
!      print *,'initializing sf6:',gq0(1,1,ntsf6),ktau,epsilc,q_corr
     endif
!    endif
!     print *,'Initialized GHGs'
!!    if (chem_in_opt.eq.0) then
!     if (ktau.le.1) then
!      do k=kts,kte
!       do i=its,ite
!!        gq0(i,k,ntco2  )=ppm2ugkg(p_co2   ) * max(epsilc,ghgi_in(i,k,1))
!!        gq0(i,k,ntch4  )=ppm2ugkg(p_ch4   ) * max(epsilc,ghgi_in(i,k,2))
!!        gq0(i,k,ntsf6  )=ppm2ugkg(p_sf6   ) * max(epsilc,ghgi_in(i,k,3))
!        gq0(i,k,ntco2  ) = 400e-6
!!        gq0(i,k,ntco2  )=ppm2ugkg(p_co2   ) * ghgi_in(i,1,1)  !ghgi_in(i,k,1)
!        gq0(i,k,ntch4  )=ppm2ugkg(p_ch4   ) * ghgi_in(i,1,2)  !ghgi_in(i,k,2)
!        gq0(i,k,ntsf6  )=ppm2ugkg(p_sf6   ) * ghgi_in(i,1,3)  !ghgi_in(i,k,3)
!       enddo
!      enddo
!      print *,'initializing ghg:',gq0(1,1,ntco2),ktau,chem_in_opt,epsilc
!     endif
!!    endif

!!!    call ghg_ch4_chem(im,kte,kme,ntrac,ntch4,ntqv,ntcw,ntiw,ntrw,ntsw,ntgl,ktau,dtstep,garea,rlat,rlon,&

!!!                      pgr,tk3d,prl3d,pr3d,gq0,ch4chm1,ch4chm2,ch4chm3,ch4loss)

!!!print *,'shaper(qgrs):',shape(qgrs)
    do k=kts,kte
!     do i=its,ite
       do i = 1, im
        qgrs(i,k,ntco2 )  = gq0(i,k,ntco2 )
        qgrs(i,k,ntco2_bgd  )= gq0(i,k,ntco2_bgd )
        qgrs(i,k,ntco2_land  )= gq0(i,k,ntco2_land )
        qgrs(i,k,ntco2_fossil  )= gq0(i,k,ntco2_fossil )
        qgrs(i,k,ntco2_fire  )= gq0(i,k,ntco2_fire )
        qgrs(i,k,ntco2_ocn  )= gq0(i,k,ntco2_ocn )
        qgrs(i,k,ntch4 )  = gq0(i,k,ntch4 )
        qgrs(i,k,ntsf6 )  = gq0(i,k,ntsf6 )
     enddo
    enddo
!    print *,'outgoing ch4',ppm2ugkg(p_ch4),gq0(1,1,ntch4)
!    print *,'outgoing co2',ppm2ugkg(p_ch4),gq0(1,1,ntco2)
!    stop

!   print *,'mpi info:',mpicomm,mpirank,mpiroot
! calculate column totals (kg)

!   call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
!   print *,'numprocs=',numprocs
!   call MPI_Comm_rank(MPI_COMM_WORLD, nproc, ierr)
!   print *,'proc=',nproc
!   call MPI_Comm_size(mpicomm, mpisize, ierr)
!   print *,'mpisize=',mpisize
!
   colghg=0.0
!   do k=kts,kte
    do i=its,ite
!      colghg=colghg+(1./g)*(pr3d(i,k)-pr3d(i,k+1))*garea(i)
!      colghg=colghg+(1./g)*(pgr(i))*garea(i)
      colghg=colghg+garea(i)
!      if (i.eq.1) then
!       print *,i,k,pr3d(i,k),pr3d(i,k+1),pr3d(i,k)-pr3d(i,k+1),pgr(i)
!      endif
    enddo
!   enddo
!   print *,nproc,' column=', colghg,pgr(1)
!   stop
!   if (mpirank.ne.mpiroot) then
!     tag=2001
!     call MPI_SEND(colghg,1,MPI_REAL,mpiroot,tag,mpicomm,ierr)
!     call MPI_SEND(colghg,1,MPI_REAL,mpiroot,0,mpicomm,ierr)
!   endif
!   if (mpirank.eq.mpiroot) then
!     totmass=colghg
!     do j=1,mpisize-1
!!       call MPI_RECV(colghg,1,MPI_REAL,j,tag,mpicomm,status,ierr)
!       call MPI_RECV(colghg,1,MPI_REAL,j,0,mpicomm,status,ierr)
!       totmass=totmass+colghg
!     enddo
!     print *,'totmass=',totmass
!   endif
!   stop

  end subroutine catchem_ghg_wrapper_run


!> @}

  subroutine ghg_ch4_chem(im,kte,kme,ntrac,ntch4,ntqv,ntcw,ntiw,ntrw,ntsw,ntgl,ktau,dtstep,garea,rlat,rlon,&
                          pgr,tk3d,prl3d,pr3d,gq0,ch4chm1,ch4chm2,ch4chm3,ch4loss)

   implicit none
   integer                              :: i,j,k,im,kte,kme,ktau,ntch4,ntqv,ntcw,ntiw,ntrw,ntsw,ntgl,ntrac
   real(kind_phys)                      :: dtstep
   real(kind_phys), dimension(im,90, 2),     intent(in)    :: ch4chm1
   real(kind_phys), dimension(im,31, 2),     intent(in)    :: ch4chm3
   real(kind_phys), dimension(im,60,2),      intent(in)    :: ch4chm2
   real(kind_phys), dimension(im),           intent(in)    :: garea, rlat,rlon, pgr
   real(kind_phys), dimension(im,kte),       intent(in)    :: prl3d, tk3d
   real(kind_phys), dimension(im,2),         intent(out)   :: ch4loss
   real(kind_phys), dimension(im,kme),       intent(in)    :: pr3d
   real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0

   real, allocatable, dimension(:,:)    :: xmdry
   real, allocatable, dimension(:,:)    :: krate,ch4_loss_rate
   real                                 :: trop_prs
!   type(budget)                         :: ghgbud
!   real, allocatable, dimension(:,:,:)  :: tempstrat,tempOH,tempOHp,tempCl,tempstratp,tempClp
   real, allocatable, dimension(:,:)    :: Sloss,OHloss,Clloss,OHprs,delt_ch4,tau
   real, parameter                      :: mwair=28.97
   real, parameter                      :: OH_scale=0.901
   real, dimension(61)                  :: at,bt

   at =[0., 0., 7.367743, 65.889244, 210.39389, 467.333588, 855.361755, 1385.912598,           &
        2063.779785, 2887.696533, 3850.91333, 4941.77832, 6144.314941, 7438.803223,            &
        8802.356445, 10209.500977, 11632.758789, 13043.21875, 14411.124023,                    &
        15706.447266, 16899.46875, 17961.357422, 18864.75, 19584.330078, 20097.402344,         &
        20384.480469, 20429.863281, 20222.205078, 19755.109375, 19027.695313,                  &
        18045.183594, 16819.474609, 15379.805664, 13775.325195, 12077.446289,                  &
        10376.126953, 8765.053711, 7306.631348, 6018.019531, 4906.708496, 3960.291504,         &
        3196.421631, 2579.888672, 2082.273926, 1680.640259, 1356.474609, 1094.834717,          &
        883.660522, 713.218079, 575.651001, 464.618134, 373.971924, 298.495789,                &
        234.779053, 180.584351, 134.483307, 95.636963, 63.647804, 38.425343, 20., 0.]
   bt= [1., 0.99763012, 0.99401945, 0.9882701, 0.97966272, 0.96764523, 0.95182151,             &
        0.93194032, 0.90788388, 0.87965691, 0.84737492, 0.81125343, 0.77159661,                &
        0.72878581, 0.68326861, 0.63554746, 0.58616841, 0.53570992, 0.48477158,                &
        0.43396294, 0.38389215, 0.33515489, 0.28832296, 0.24393314, 0.2024759,                 &
        0.16438432, 0.13002251, 0.09967469, 0.07353383, 0.05169041, 0.03412116,                &
        0.02067788, 0.01114291, 0.00508112, 0.00181516, 0.00046139, 7.582e-05, 0., 0.,         &
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]



! vertical interpolation of incoming fields

! stratloss (ch4chm1)
   allocate(Sloss(im,kte))
   call vert_interp(im,90,kte,prl3d,ch4chm1(:,:,2)*100.,ch4chm1(:,:,1),Sloss)
! OH (ch4chm2)
   allocate(OHprs(im,60))
   do i=1,im
    do k=1,60
      OHprs(i,k)=at(k)+pgr(i)*bt(k)
    enddo
   enddo
   allocate(OHloss(im,kte))
   call vert_interp(im,60,kte,prl3d,OHprs,ch4chm2(:,:,1),OHloss)
! Clloss (ch4chm3)
   allocate(Clloss(im,kte))
   call vert_interp(im,31,kte,prl3d,ch4chm3(:,:,2),ch4chm3(:,:,1),Clloss)

   allocate(krate(im,kte))            ! all incoming chemloss info defined on same local domain
   allocate(ch4_loss_rate(im,kte))
   allocate(delt_ch4(im,kte))
   allocate(tau(im,kte))

! stratospheric loss due to OH, Cl and O1D

   ch4_loss_rate=SLoss        !s^-1
!   print *,'ch4_loss_rate(strat)=',ch4_loss_rate(1,:)

! tropospheric loss due to OH and Cl

   do k=1,kte
    do i=1,im

       trop_prs=30000.0-21500.0*cos(rlat(i))**2      !tropopause pressure approximation, Lawrence et al., 2002, lat in deg.

       if (prl3d(i,k).ge.trop_prs) then
         krate(i,k)=2.45e-12*exp(-1775./(tk3d(i,k)))  !cm^3 molec^-1 s^-1
         ch4_loss_rate(i,k)=OH_scale*krate(i,k)*OHloss(i,k)   !cm^3 molec^-1 s^-1 * molec * cm^-3 = s^-1
         krate(i,k)=2.36e-12*((tk3d(i,k)/298.)**1.37)*exp(-939./tk3d(i,k)) !cm^3 molec^-1 s^-1
         ch4_loss_rate(i,k)=ch4_loss_rate(i,k)+krate(i,k)*Clloss(i,k) !cm^3 molec^-1 s^-1 * molec * cm^-3 = s^-1
       endif

    enddo
   enddo

!  print *,'ch4_loss_rate(strat+trop)=',ch4_loss_rate(1,:),rlat(1)*180./3.141596,rlon(1)*180./3.141596

! apply chemical loss to CH4

   ch4loss=0.0
   allocate(xmdry(im,kte))
   do k=1,kte
    do i=1,im

     trop_prs=30000.0-21500.0*cos(rlat(i))**2      !tropopause pressure approximation, Lawrence et al., 2002, lat in deg.

! calculate total kg of dry air in each gbox
     xmdry(i,k) = ((pr3d(i,k)-pr3d(i,k+1))/g)*garea(i)* &
           (1-gq0(i,1,ntqv)-gq0(i,1,ntcw)-gq0(i,1,ntiw)-gq0(i,1,ntrw)-&
           gq0(i,1,ntsw)-gq0(i,1,ntgl))

     delt_ch4(i,k)=gq0(i,k,ntch4)*(exp(-ch4_loss_rate(i,k)*dtstep)-1.)
     if (delt_ch4(i,k).gt.0.) then
        print *,'delt_ch4 has the wrong sign',delt_ch4(i,k),gq0(i,k,ntch4),ch4_loss_rate(i,k),dtstep,i,k
        stop
     endif
!     if (i.eq.1.and.k.eq.1) print *,'in ch4chem:',gq0(1,1,ntch4),ch4_loss_rate(1,1),delt_ch4(1,1),prl3d(1,1)
     tau(i,k)=-dtstep/(alog(delt_ch4(i,k)/gq0(i,k,ntch4)+1))
     tau(i,k)=tau(i,k)/(86400.*365.)
!     gq0(i,k,ntch4)=gq0(i,k,ntch4)+delt_ch4(i,k)
!     if (prl3d(i,k).ge.trop_prs) then
    if (prl3d(i,k).ge.0.) then
       ch4loss(i,1)=ch4loss(i,1)+delt_ch4(i,k)*xmdry(i,k)
!       print *,'trop data',i,k,prl3d(i,k),trop_prs,delt_ch4(i,k),xmdry(i,k),ch4loss(i,1)
     else
!       ch4loss(i,2)=ch4loss(i,2)+delt_ch4(i,k)*xmdry(i,k)
       ch4loss(i,2)=garea(i)
     endif
!     ghgbud%chem(2)=ghgbud%chem(2)+delt_ch4(i,j,k)*xmdry(i,j,k)*mwair*1.e-3

    enddo
   enddo
   print *,'trop loss:',maxval(ch4loss(:,1)),minval(ch4loss(:,1)),rlat(1)*180./3.141596,rlon(1)*180./3.141596
   print *,'strat loss:',maxval(ch4loss(:,2)),minval(ch4loss(:,2))

!   print *,'xmdry',xmdry(:,5)
!   print *,'delt_ch4',delt_ch4(:,5)
!   print *,'new value:',gq0(1,1,ntch4)
  end subroutine ghg_ch4_chem
  subroutine vert_interp( im, nlev, clev, prs, fieldcsprs, fieldcs, field_out )

! nlev - number of ll levs
! im - horizontal grid of points 
! clev - number of FV3 levs 
! prs - FV3 pressure on local domain
! fieldcsprs - incoming pressure array on cs horiz grid but with ll levs
! fieldcs  - incoming mixing ratio array on cs horiz grid but with ll levs
! vertically interpolated array on cs grid

      integer, intent(in)                :: nlev,clev,im
      real                               :: fieldcsprs(:,:),fieldcs(:,:)
      real*8                             :: prs(:,:)
      real, intent(out)                  :: field_out(:,:)
      integer                            :: i,j,k,l
      real                               :: dz

!      print *,'incoming pressure levels',prs_cs(is,js,:)
!      print *,'log of incoming pressure levels',alog(prs_cs(is,js,:))
!      print *,'incoming init pressure levels',fieldcsprs(is,js,:)
!      print *,'log of incoming init pressure levels',alog(fieldcsprs(is,js,:))

      do i=1,im
       do k=1,clev

! check for Zero pressures at top of fv3 atm (why does this happen???)

          if (prs(i,k).le.0.) then
                prs(i,k)=prs(i,k-1)
          endif
          if (prs(i,k).gt.fieldcsprs(i,1)) then
            field_out(i,k)=fieldcs(i,1)
          else if (prs(i,k).le.fieldcsprs(i,nlev)) then
            field_out(i,k)=fieldcs(i,nlev)
          else
            do l=1,nlev-1
              if (prs(i,k).le.fieldcsprs(i,l).and.prs(i,k).gt.fieldcsprs(i,l+1)) then
                dz=(alog(prs(i,k))-alog(fieldcsprs(i,l)))/(alog(fieldcsprs(i,l+1))-alog(fieldcsprs(i,l)))
                field_out(i,k)=fieldcs(i,l)+dz*(fieldcs(i,l+1)-fieldcs(i,l))
              endif
            enddo
          endif
       enddo
      enddo

  end subroutine vert_interp
  

!> @}
  end module catchem_dust_wrapper
