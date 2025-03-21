module ghg_mod
!
!  This module developed by Lori Bruhwiler (NOAA GML) and Andrew Schuh (CIRA)
!  For serious  and unserious questions contact lori.bruhwiler@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li
!  06/2023, Restructure for CATChem, Jian.He@noaa.gov

  use catchem_constants, only : kind_chem, g=>con_g, pi=>con_pi
  use catchem_ghg_config, only: num_chem,num_moist,epsilc

  implicit none

  private

  public :: ghg_driver

contains

  subroutine ghg_driver(dt,           &
          chem_arr,rho_phy,dz8w,smois,         &
          delp,ssm,isltyp,vegfra,snowh,area,  &
          emis_dust,ust,znt,clay,sand, &
          rdrag,uthr,num_soil_layers,random_factor)

  IMPLICIT NONE

     INTEGER,      INTENT(IN   ) :: isltyp,  &
                                    num_soil_layers
     REAL(kind=kind_chem), INTENT(IN   ) :: dt,&
                                            rho_phy, &
                                            dz8w, &
                                            delp,&
                                            ssm,vegfra,snowh,&
                                            area, &
                                            ust,znt, &
                                            clay,sand, & 
                                            rdrag,uthr,&
                                            random_factor

     REAL(kind=kind_chem), DIMENSION( num_chem ),                 &
           INTENT(INOUT ) ::                                   chem_arr
     REAL(kind=kind_chem), DIMENSION(num_emis_dust),&
           OPTIONAL, INTENT(INOUT ) ::                  emis_dust
     REAL(kind=kind_chem), DIMENSION( num_soil_layers ) , &
        INTENT(INOUT) ::                               smois

    ! Local variables

    integer :: nmx,smx
    integer :: ilwi
    real(kind_chem) :: erodtot, gravsm, drylimit
    real(kind_chem), DIMENSION (5)   :: tc,bems
    real(kind_chem) :: airden,airmas,ustar,dxy
    real(kind_chem), dimension (3) :: massfrac
    real(kind_chem) :: conver,converi
    real(kind_chem) :: R

    ! threshold values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    nmx=ndust
    smx=nsalt


    ! Don't do dust over water!!!

    ilwi=1

    ! Total concentration at lowest model level. This is still hardcoded for 5 bins.
    tc(1)=chem_arr(p_dust_1)*conver
    tc(2)=chem_arr(p_dust_2)*conver
    tc(3)=chem_arr(p_dust_3)*conver
    tc(4)=chem_arr(p_dust_4)*conver
    tc(5)=chem_arr(p_dust_5)*conver

    ! Air mass and density at lowest model level.
    airmas=area * delp / g
    airden=rho_phy
    ustar=ust
    dxy=area

    ! Mass fractions of clay, silt, and sand.
    massfrac(1)=clay
    massfrac(2)=1-(clay+sand)
    massfrac(3)=sand


    ! Total erodibility.
    erodtot = ssm       ! SUM(erod(i,j,:))
             
    ! Don't allow roughness lengths greater than 20 cm to be lofted.
    ! This kludge accounts for land use types like urban areas and
    ! forests which would otherwise show up as high dust emitters.
    ! This is a placeholder for a more widely accepted kludge
    ! factor in the literature, which reduces lofting for rough areas.
    ! Forthcoming...

    IF (znt .gt. 0.2) then
      ilwi=0
    ENDIF

    ! limit where there is lots of vegetation
    if (vegfra .gt. .17) then
       ilwi = 0
    endif

    ! limit where there is snow on the ground
    if (snowh .gt. 0) then
       ilwi = 0
    endif

    ! Do not allow areas with bedrock, lava, or land-ice to loft

    IF (isltyp.eq. 15 .or. isltyp .eq. 16. .or. &
          isltyp .eq. 18) then
      ilwi=0
    ENDIF
    IF (isltyp .eq. 0)then
      ilwi=0
    endif

    if(ilwi == 1 ) return

    ! Calculate gravimetric soil moisture and drylimit.
    gravsm=100.*smois(1)/((1.-maxsmc(isltyp))*(2.65*(1.-clay)+2.50*clay))
    drylimit=14.0*clay*clay+17.0*clay

    ! get drag partition
    ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
    !     doi:10.1016/j.geomorph.2004.03.009
    if (dust_calcdrag .ne. 1) then
       call fengsha_drag(znt,R)
    else
      ! use the precalculated version derived from ASCAT; Prigent et al. (2012,2015)
      ! doi:10.1109/TGRS.2014.2338913 & doi:10.5194/amt-5-2703-2012
      ! pick only valid values
      if (rdrag > 0.) then
        R = real(rdrag, kind=kind_chem)
      else
        return 
      endif
    endif  

    ! Call dust emission routine.
    call source_dust(nmx, smx, dt, tc, ustar, massfrac, & 
                  erodtot, dxy, gravsm, airden, airmas, &
                  bems, g, drylimit, dust_alpha, dust_gamma, R, uthr, random_factor)

    chem_arr(p_dust_1)=tc(1)*converi
    chem_arr(p_dust_2)=tc(2)*converi
    chem_arr(p_dust_3)=tc(3)*converi
    chem_arr(p_dust_4)=tc(4)*converi
    chem_arr(p_dust_5)=tc(5)*converi

    ! -- for output diagnostics
    emis_dust(p_edust1)=bems(1)
    emis_dust(p_edust2)=bems(2)
    emis_dust(p_edust3)=bems(3)
    emis_dust(p_edust4)=bems(4)
    emis_dust(p_edust5)=bems(5)

    !

  end subroutine gocart_dust_fengsha_driver


  SUBROUTINE source_dust(nmx, smx, dt1, tc, ustar, massfrac, &
       erod, dxy, gravsm, airden, airmas, bems, g0, drylimit, alpha,  &
       gamma, R, uthres, random_factor)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
    ! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         SMX       Number of saltation bins                      (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
    ! *         SPOINT    Pointer to 3 soil classes                     (-)
    ! *         DEN_DUST  Dust density                                  (kg/m3)
    ! *         DEN_SALT  Saltation particle density                    (kg/m3)
    ! *         REFF_SALT Reference saltation particle diameter         (m)
    ! *         REFF_DUST Reference dust particle diameter              (m)
    ! *         LO_DUST   Lower diameter limits for dust bins           (m)
    ! *         UP_DUST   Upper diameter limits for dust bins           (m)
    ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
    ! *
    ! *  Parameters:
    ! *         CMB       Constant of proportionality                   (-)
    ! *         MMD_DUST  Mass median diameter of dust                  (m)
    ! *         GSD_DUST  Geometric standard deviation of dust          (-)
    ! *         LAMBDA    Side crack propagation length                 (m)
    ! *         CV        Normalization constant                        (-)
    ! *         G0        Gravitational acceleration                    (m/s2)
    ! *         G         Gravitational acceleration in cgs             (cm/s2)
    ! *
    ! *  Working:
    ! *         U_TS0     "Dry" threshold friction velocity             (m/s)
    ! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DEN       Dust density in cgs                           (g/cm3)
    ! *         DIAM      Dust diameter in cgs                          (cm)
    ! *         DMASS     Saltation mass distribution                   (-)
    ! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         SALT      Saltation flux                                (kg/m/s)
    ! *         DLNDP     Dust bin width                                (-)
    ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
    ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
    ! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
    ! *
    ! *  Output:
    ! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
    ! *         BEMS      Source of each dust type           (kg/timestep/cell)
    ! *
    ! ****************************************************************************

    INTEGER,            INTENT(IN)    :: nmx,smx
    REAL(kind_chem), INTENT(IN)    :: dt1
    REAL(kind_chem), INTENT(INOUT) :: tc(nmx)
    REAL(kind_chem), INTENT(IN)    :: ustar
    REAL(kind_chem), INTENT(IN)    :: massfrac(3)
    REAL(kind_chem), INTENT(IN)    :: erod
    REAL(kind_chem), INTENT(IN)    :: dxy
    REAL(kind_chem), INTENT(IN)    :: gravsm
    REAL(kind_chem), INTENT(IN)    :: random_factor
    REAL(kind_chem), INTENT(IN)    :: airden
    REAL(kind_chem), INTENT(IN)    :: airmas
    REAL(kind_chem), INTENT(OUT)   :: bems(nmx)
    REAL(kind_chem), INTENT(IN)    :: g0
    REAL(kind_chem), INTENT(IN)    :: drylimit
    !! Sandblasting mass efficiency, aka "fudge factor" (based on Tegen et al,
    !! 2006 and Hemold et al, 2007)
    !
    !  REAL, PARAMETER :: alpha=1.8E-8  ! (m^-1)
    REAL(kind_chem), INTENT(IN)    :: alpha
    ! Experimental optional exponential tuning constant for erodibility.
    ! 0 < gamma < 1 -> more relative impact by low erodibility regions.
    REAL(kind_chem), INTENT(IN)    :: gamma
    REAL(kind_chem), INTENT(IN)    :: R
    REAL(kind_chem), INTENT(IN)    :: uthres

    REAL(kind_chem)    :: den(smx), diam(smx)
    REAL(kind_chem)    :: dvol(nmx), distr_dust(nmx), dlndp(nmx)
    REAL(kind_chem)    :: dsurface(smx), ds_rel(smx)
    REAL(kind_chem)    :: u_ts0, u_ts, dsrc, dmass, dvol_tot
    REAL(kind_chem)    :: salt,emit, emit_vol, stotal
    REAL(kind_chem)    :: rhoa, g
    INTEGER   :: i, j, n

    ! Sandblasting mass efficiency, beta.
    ! Beta maxes out for clay fractions above 0.2 = betamax.

    REAL(kind_chem), PARAMETER :: betamax=5.25E-4
    REAL(kind_chem) :: beta
    integer :: styp

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(kind_chem), PARAMETER :: cmb=1.0
    ! REAL, PARAMETER :: cmb=2.61   ! from White,1979

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(kind_chem), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_chem), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_chem), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(kind_chem), PARAMETER :: cv=12.62D-6      ! normalization constant

    ! Calculate saltation surface area distribution from sand, silt, and clay
    ! mass fractions and saltation bin fraction. This will later become a
    ! modifier to the total saltation flux.  The reasoning here is that the
    ! size and availability of saltators affects saltation efficiency. Based
    ! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

    DO n=1,smx
       dmass=massfrac(spoint(n))*frac_salt(n)
       dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
    ENDDO

    ! The following equation yields relative surface area fraction.  It will only
    ! work if you are representing the "full range" of all three soil classes.
    ! For this reason alone, we have incorporated particle sizes that encompass
    ! the clay class, to account for the its relative area over the basal
    ! surface, even though these smaller bins would be unlikely to play any large
    ! role in the actual saltation process. - GC

    stotal=SUM(dsurface(:))
    DO n=1,smx
       ds_rel(n)=dsurface(n)/stotal
    ENDDO

    ! Calculate total dust emission due to saltation of sand sized particles.
    ! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    ! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    ! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally,
    ! calculate total dust emission (tot_emit), taking into account erodibility.

    ! Set DRY threshold friction velocity to input value
    u_ts0 = uthres

    g = g0*1.0E2
    emit=0.0

    DO n = 1, smx
       den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
       diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
       rhoa = airden*1.0D-3                ! (g cm^-3)

             ! FENGSHA uses the 13 category soil type from the USDA
             ! call calc_fengsha_styp(massfrac(1),massfrac(3),massfrac(2),styp)
             ! Fengsha uses threshold velocities based on dale gilletes data
             ! call fengsha_utst(styp,uthres,u_ts0)

             ! Friction velocity threshold correction function based on physical
             ! properties related to moisture tension. Soil moisture greater than
             ! dry limit serves to increase threshold friction velocity (making
             ! it more difficult to loft dust). When soil moisture has not reached
             ! dry limit, treat as dry

             IF (gravsm > drylimit) THEN
                u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm-drylimit)**0.68)) / R)
             ELSE
                u_ts = u_ts0 / R
             END IF

             ! Calculate total vertical mass flux (note beta has units of m^-1)
             ! Beta acts to tone down dust in areas with so few dust-sized particles that the
             ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
             ! producers, which is generally not the case.  Equation derived from wind-tunnel
             ! experiments (see MB95).

             beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
             if (massfrac(1) <= 0.2) then
                beta=10**(13.4*massfrac(1)-6.0)
             else
                beta = 2.E-4
             endif

             !---------------------------------------------------------------------
             ! formula of Draxler & Gillette (2001) Atmos. Environ.
             ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
             !
             ! where:
             !     F   = vertical emission flux  [g/m**2-s]
             !     K   = constant 2.0E-04                      [1/m]
             !     A   = 0~3.5  mean = 2.8  (fudge factor)
             !     U*  = friction velocity                     [m/s]
             !     Ut* = threshold friction velocity           [m/s]
             !
             !--------------------------------------------------------------------

             IF (ustar .gt. u_ts) then
                call fengsha_hflux(ustar,u_ts,beta, salt)
                salt = alpha * cmb * ds_rel(n) * airden / g0 * salt * (erod**gamma) * beta
             else
                salt = 0.
             endif
             ! EROD is taken into account above
             emit = emit + salt 
    END DO

    ! Now that we have the total dust emission, distribute into dust bins using
    ! lognormal distribution (Dr. Jasper Kok, in press), and
    ! calculate total mass emitted over the grid box over the timestep.
    !
    ! In calculating the Kok distribution, we assume upper and lower limits to each bin.
    ! For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
    ! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
    ! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
    ! These may be changed within module_data_gocart_dust.F, but make sure it is
    ! consistent with reff_dust values.  These values were taken from the original
    ! GOCART bin configuration. We use them here to calculate dust bin width, dlndp.
    ! dVol is the volume distribution. You know...if you were wondering. GC

    dvol_tot=0.
    DO n=1,nmx
       dlndp(n)=LOG(up_dust(n)/lo_dust(n))
       dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
            EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
       dvol_tot=dvol_tot+dvol(n)
       ! Convert mass flux to volume flux
       !emit_vol=emit/den_dust(n) ! (m s^-1)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
       !print *,"distr_dust(",n,")=",distr_dust(n)
    END DO

    ! Now distribute total vertical emission into dust bins and update concentration.

    DO n=1,nmx
             ! Calculate total mass emitted
             dsrc = emit*distr_dust(n)*dxy*dt1 *random_factor  ! (kg)
             IF (dsrc < 0.0) dsrc = 0.0

             ! Update dust mixing ratio at first model level.
             tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
             !   bems(i,j,n) = dsrc  ! diagnostic
             !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
             bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
    END DO

  END SUBROUTINE source_dust

  subroutine fengsha_utst(styp,uth, ut)
    integer,                            intent(in)  :: styp
    real(kind_chem), dimension(fengsha_maxstypes), intent(in)  :: uth
    real(kind_chem),                 intent(out) :: ut
    ut = uth(styp)
!     real (kind_chem) :: uth(13) = &
!          (/ 0.08,   & ! Sand          - 1
!          0.20,    & ! Loamy Sand      - 2
!          0.30,    & ! Sandy Loam      - 3
!          0.30,    & ! Silt Loam       - 4
!          0.35,    & ! Silt            - 5
!          0.60,    & ! Loam            - 6
!          0.30,    & ! Sandy Clay Loam - 7
!          0.35,    & ! Silty Clay Loam - 8
!          0.45,    & ! Clay Loam       - 9
!          0.45,    & ! Sandy Clay      - 10
!          0.45,    & ! Silty Clay      - 11
!          0.60,    & ! Clay            - 12
!          9.999 /)   ! Other           - 13
    return
  end subroutine fengsha_utst

  subroutine calc_fengsha_styp(clay, sand, silt, type)

    !---------------------------------------------------------------
    ! Function: calculate soil type based on USDA definition.
    ! Source: USDA soil texture calculator
    !
    ! Defintion of soil types:
    !
    !
    ! NOAH 1      2             3           4           5      6      7                 8                9           10           11           12
    ! PX   1      2             3           4           -      5      6                 7                8           9            10           11
    ! Soil "Sand" "Loamy Sand" "Sandy Loam" "Silt Loam" "Silt" "Loam" "Sandy Clay Loam" "Silt Clay Loam" "Clay Loam" "Sandy Clay" "Silty Clay" "Clay"
    !---------------------------------------------------------------
    REAL(kind_chem), intent(in) ::  clay, sand, silt
    integer, intent(out) ::  type
    real(kind_chem) :: cly, snd, slt

    type = 0

    snd = sand * 100.
    cly = clay * 100.
    slt = silt * 100.
    if (slt+1.5*cly .lt. 15)                                                                type = 1      ! snd
    if (slt+1.5*cly .ge. 15 .and.slt+1.5*cly .lt. 30)                                       type = 2      ! loamy snd
    if (cly .ge. 7 .and. cly .lt. 20 .and. snd .gt. 52 .and. slt+2*cly .ge. 30)             type = 3      ! sndy loam (cond 1)
    if (cly .lt. 7 .and. slt .lt. 50 .and. slt+2*cly .ge. 30)                               type = 3      ! sndy loam (cond 2)
    if (slt .ge. 50 .and. cly .ge. 12 .and.cly .lt. 27 )                                    type = 4      ! slt loam (cond 1)
    if (slt .ge. 50 .and. slt .lt. 80 .and.cly .lt. 12)                                     type = 4      ! slt loam (cond 2)
    if (slt .ge. 80 .and. cly .lt. 12)                                                      type = 5      ! slt
    if (cly .ge. 7  .and. cly .lt. 27 .and.slt .ge. 28 .and. slt .lt. 50 .and.snd .le. 52)  type = 6      ! loam
    if (cly .ge. 20 .and. cly .lt. 35 .and.slt .lt. 28 .and. snd .gt. 45)                   type = 7      ! sndy cly loam
    if (cly .ge. 27 .and. cly .lt. 40 .and.snd .lt. 20)                                     type = 8      ! slt cly loam
    if (cly .ge. 27 .and. cly .lt. 40 .and.snd .ge. 20 .and. snd .le. 45)                   type = 9      ! cly loam
    if (cly .ge. 35 .and. snd .gt. 45)                                                      type = 10     ! sndy cly
    if (cly .ge. 40 .and. slt .ge. 40)                                                      type = 11     ! slty cly
    if (cly .ge. 40 .and. snd .le. 45 .and.slt .lt. 40)                                     type = 12     ! clay
    return
  end subroutine calc_fengsha_styp

  subroutine fengsha_drag(z0,R)
    real(kind_chem), intent(in) :: z0
    real(kind_chem), intent(out) :: R
    real(kind_chem), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine fengsha_hflux(ust,utst, kvh, salt)
    !---------------------------------------------------------------------
    ! Function: Calculates the Horizontal Saltation Flux, Q, and then
    !           calculates the vertical flux.
    !
    ! formula of Draxler & Gillette (2001) Atmos. Environ.
    ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
    !
    ! where:
    !     F   = vertical emission flux  [g/m**2-s]
    !     K   = constant 2.0E-04                      [1/m]
    !     A   = 0~3.5  mean = 2.8  (fudge factor)
    !     U*  = friction velocity                     [m/s]
    !     Ut* = threshold friction velocity           [m/s]
    !
    !--------------------------------------------------------------------
    real(kind_chem), intent(in) :: ust, & ! friction velocity
                                     utst, & ! threshold friction velocity
                                      kvh    ! vertical to horizontal mass flux ratio

    real(kind_chem), intent(out) :: salt
    real(kind_chem) :: Q
    Q = ust * (ust * ust - utst * utst)
    salt = Q ! sdep * kvh * Q

    return
  end subroutine fengsha_hflux


end module dust_fengsha_mod
