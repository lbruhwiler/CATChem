!>\file catchem_ghg_emi_wrapper.F90
!! This file is the catchem-ghg emission wrapper with CCPP coupling to FV3
!! Lori.Bruhwiler@noaa.gov 02/2024

 module catchem_ghg_emi_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_ghg_config,     only : num_chem,num_moist,epsilc

   implicit none

   private

   public :: catchem_ghg_emi_wrapper_init, catchem_ghg_emi_wrapper_run, catchem_ghg_emi_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_ghg_emi_wrapper_init()
      end subroutine catchem_ghg_emi_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_ghg_emi_wrapper_finalize Argument Table
!!
      subroutine catchem_ghg_emi_wrapper_finalize()
      end subroutine catchem_ghg_emi_wrapper_finalize

!> \defgroup gsd_chem_anthropogenic_group GSD Chem seas wrapper Module
!! This is the gsd chemistry
!>\defgroup catchem_ghg_emi_wrapper GML EMI driver Module  
!> \ingroup catchem_ghg_emi_group
!! This is the GML EMI wrapper Module
!! \section arg_table_catchem_ghg_emi_wrapper_run Argument Table
!! \htmlinclude catchem_ghg_emi_wrapper_run.html
!!
!>\section catchem_ghg_emi_wrapper GML Emi Scheme General Algorithm
!> @{
    subroutine catchem_ghg_emi_wrapper_run(im, kte, kme, ktau, dt, garea, &
                 rlat,rlon,julian, pr3d, ph3d, prl3d, tk3d, spechum,idat,    &
                 ntrac,ntco2,ntco2_bgd,ntco2_land,ntco2_fossil,ntco2_fire,ntco2_ocn,ntch4,ntsf6,ntqv,ntcw,ntiw,ntrw,ntsw,  &
                 ntgl,ghg_emi,ghgem,gq0,qgrs,tile_num,mpirank,mpiroot,errmsg,errflg,naux2d,naux3d,aux2d,aux3d)  !,mpirank,mpiroot)


    implicit none

    integer,        intent(in) :: im,kte,kme,ktau,tile_num,idat(8)
    integer,        intent(in) :: ntrac
    integer,        intent(in) :: ntco2,ntco2_bgd,ntco2_land,ntco2_fossil,ntco2_fire,ntco2_ocn,ntch4,ntsf6
    integer,        intent(in) :: ntqv,ntcw,ntiw,ntrw,ntsw,ntgl
    real(kind_phys),intent(in) :: dt,julian

    ! MPI information
    integer,                   intent(in   ) :: mpirank
    integer,                   intent(in   ) :: mpiroot

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1
    integer, parameter :: num_ghg=3,p_co2=1,p_ch4=2,p_sf6=3
    real,    parameter :: mw_co2=44.01, mw_ch4=16.04, mw_sf6=146.06, mw_dry=28.97

    integer, intent(in) :: naux2d,naux3d
    real(kind_phys), intent(inout) :: aux2d(:,:)
    real(kind_phys), intent(inout) :: aux3d(:,:,:)

    real(kind_phys), dimension(im), intent(in)        :: garea,rlat,rlon
    real(kind_phys), dimension(im, 24, 6), intent(in) :: ghg_emi
    real(kind_phys), dimension(im,6), intent(inout) :: ghgem
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: prl3d, tk3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,       &
                     p_phy, z_at_w, dz8w, p8w, rho_phy

    real(kind_phys), dimension(im)         :: xmdry, xmwet
    integer :: ide, ime, ite, kde
    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
    real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang

    integer :: i, j, jp, k, kp, n, ihr
  
    errmsg = ''
    errflg = 0

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg        = 1._kind_phys
    ppm2ugkg(p_co2) = 1._kind_phys * mw_co2 / mw_dry
!    ppm2ugkg(p_ch4) = 1._kind_phys * mw_ch4 / mw_dry
    ppm2ugkg(p_ch4) = 1._kind_phys * mw_co2 / mw_dry
    ppm2ugkg(p_sf6) = 1._kind_phys * mw_sf6 / mw_dry


    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
    !*!  dtstep = call_chemistry * dt
    dtstep = dt
    else
      dtstep = dt
    end if

!   using modulo and just one diurnal cycle for now
    ihr=int((ktau*dt)/(3600.))+1
    ihr = mod(ihr,24)+1
    if(mpirank == mpiroot) print *,'time vars:',ihr,ktau,dt

! calculate total kg of dry air in each gbox

    do i=1,im
        xmdry(i) = ((pr3d(i,1)-pr3d(i,2))/g)*garea(i)* &
           (1-gq0(i,1,ntqv)-gq0(i,1,ntcw)-gq0(i,1,ntiw)-gq0(i,1,ntrw)-&
           gq0(i,1,ntsw)-gq0(i,1,ntgl))
        xmwet(i) = ((pr3d(i,1)-pr3d(i,2))/g)*garea(i)
    enddo

    do i=1,im
!!!       ghgem(i,1)=1.e-3*mw_co2*ghg_emi(i,ihr,1)*garea(i)*dt   !in kg
!!!       ghgem(i,2)=1.e-3*mw_ch4*ghg_emi(i,ihr,2)*garea(i)*dt
!!!       ghgem(i,3)=1.e-3*mw_sf6*ghg_emi(i,ihr,3)*garea(i)*dt
!!! like abem:

       ghgem(i,1)=1.e-3*mw_co2*ghg_emi(i,ihr,1)   !  kg m^-2 s^-1
       ghgem(i,2)=1.e-3*mw_co2*ghg_emi(i,ihr,2)
       ghgem(i,3)=1.e-3*mw_co2*ghg_emi(i,ihr,3)
       ghgem(i,4)=1.e-3*mw_co2*ghg_emi(i,ihr,4)
       ghgem(i,5)=1.e-3*mw_ch4*ghg_emi(i,ihr,5)
       ghgem(i,6)=1.e-3*mw_sf6*ghg_emi(i,ihr,6)

!       ghgem(i,1)=1.e-3*mw_co2*ghg_emi(i,1)
!       ghgem(i,2)=1.e-3*mw_co2*ghg_emi(i,2)
!       ghgem(i,3)=1.e-3*mw_co2*ghg_emi(i,3)
!       ghgem(i,4)=1.e-3*mw_co2*ghg_emi(i,4)

!       aux2d(i,1) = ghgem(i,1)
!       aux2d(i,2) = ghgem(i,2)
!       aux2d(i,3) = ghgem(i,3)
!       aux2d(i,4) = ghgem(i,4)

!       ghgem(i,1)=0.0
!        ghgem(i,2)=1.e-3*mw_ch4*ghg_emi(i,ihr,2)
!       ghgem(i,2)=0.0
!        ghgem(i,3)=1.e-3*mw_sf6*ghg_emi(i,ihr,3)
!       ghgem(i,3)=0.0

!       if(mpirank == mpiroot) print *,'adding:',ghgem(i,1)*garea(i)*dt/xmwet(i),' to ',gq0(i,1,ntco2_land)

!        print *,'ntco2,ntch4,ntsf6:',ntco2,':',ntch4,':',ntsf6

       gq0(i,1,ntco2)=gq0(i,1,ntco2_land)+(ghgem(i,1)+ghgem(i,2)+ghgem(i,3)+ghgem(i,4))*garea(i)*dt/xmwet(i)   !kg/kg
       gq0(i,1,ntco2_land)=gq0(i,1,ntco2_land) +ghgem(i,1)*garea(i)*dt/xmwet(i)   !kg/kg
       gq0(i,1,ntco2_fire)=gq0(i,1,ntco2_fire)+ghgem(i,2)*garea(i)*dt/xmwet(i)   !kg/kg
       gq0(i,1,ntco2_fossil)=gq0(i,1,ntco2_fossil)+ghgem(i,3)*garea(i)*dt/xmwet(i)   !kg/kg
       gq0(i,1,ntco2_ocn)=gq0(i,1,ntco2_ocn)+ghgem(i,4)*garea(i)*dt/xmwet(i)   !kg/kg
       gq0(i,1,ntch4)=gq0(i,1,ntch4)+ghgem(i,5)*garea(i)*dt/xmwet(i)
       gq0(i,1,ntsf6)=gq0(i,1,ntsf6)+ghgem(i,6)*garea(i)*dt/xmwet(i)
    enddo
    if(mpirank == mpiroot) print *,mpirank,'max/min of co2_land emi:',maxval(ghgem(:,1)),minval(ghgem(:,1)),ntco2_land
    if(mpirank == mpiroot) print *,'max/min of co2_fire emi:',maxval(ghgem(:,2)),minval(ghgem(:,2)),ntco2_fire
    if(mpirank == mpiroot) print *,'max/min of co2_fossil emi:',maxval(ghgem(:,3)),minval(ghgem(:,3)),ntco2_fossil
    if(mpirank == mpiroot) print *,'max/min of co2_ocn emi:',maxval(ghgem(:,4)),minval(ghgem(:,4)),ntco2_ocn
    if(mpirank == mpiroot) print *,'max/min of ch4 emi:',maxval(ghgem(:,5)),minval(ghgem(:,5)),ntch4
    if(mpirank == mpiroot) print *,'max/min of sf6 emi:',maxval(ghgem(:,6)),minval(ghgem(:,6)),ntsf6


!    ! -- put chem stuff back into tracer array
!    THERE IS NO REASON TO ALTER UNIT, ALREADY KG/KG
!     do i=its,ite
!      do i=1,im
!       gq0(i,1,ntco2  )=ppm2ugkg(p_co2   ) * max(epsilc,gq0(i,1,p_co2))
!       gq0(i,1,ntco2_land  )=ppm2ugkg(p_co2   ) * max(epsilc,gq0(i,1,ntco2_land))
!       gq0(i,1,ntco2_fossil  )=ppm2ugkg(p_co2   ) * max(epsilc,gq0(i,1,ntco2_fossil))
!       gq0(i,1,ntco2_fire  )=ppm2ugkg(p_co2   ) * max(epsilc,gq0(i,1,ntco2_fire))
!       gq0(i,1,ntco2_ocn  )=ppm2ugkg(p_co2   ) * max(epsilc,gq0(i,1,ntco2_ocn))
!       gq0(i,1,ntch4  )=ppm2ugkg(p_ch4  ) * max(epsilc,gq0(i,1,p_ch4))
!       gq0(i,1,ntsf6  )=ppm2ugkg(p_sf6   ) * max(epsilc,gq0(i,1,p_sf6))
!     enddo

!     print *,'shaper(qgrs):',shape(qgrs)
     if(mpirank == mpiroot) print *,'before:',i,qgrs(i,k,ntco2_land  )*1e6,' afterr:',gq0(i,1,ntco2_land)*1e6
     if(mpirank == mpiroot) print *,'before:',i,qgrs(i,k,ntco2_ocn  )*1e6,' afterr:',gq0(i,1,ntco2_ocn)*1e6

!    do k=kts,kte
     do k=1,kte
!     do i=its,ite
      do i=1,im
       qgrs(i,k,ntco2 )=gq0(i,k,ntco2 )
       qgrs(i,k,ntco2_land  )=gq0(i,k,ntco2_land  )
       qgrs(i,k,ntco2_fossil  )=gq0(i,k,ntco2_fossil  )
       qgrs(i,k,ntco2_fire  )=gq0(i,k,ntco2_fire  )
       qgrs(i,k,ntco2_ocn  )=gq0(i,k,ntco2_ocn  )
       qgrs(i,k,ntch4  )=gq0(i,k,ntch4 )
       qgrs(i,k,ntsf6  )=gq0(i,k,ntsf6 )
     enddo
    enddo

   if(mpirank == mpiroot) print *,'DONE DONE DONE with catchem ghg emi....'
!
   end subroutine catchem_ghg_emi_wrapper_run
!> @}

  end module catchem_ghg_emi_wrapper
