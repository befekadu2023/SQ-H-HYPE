!> \file soil_proc.f90
!> Contains module soil_processes.

!>Water processes in soil in HYPE and some more
MODULE SOIL_PROCESSES

  !Copyright 2012-2016 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !-----------------------------------------------------------------------------------------

  USE STATETYPE_MODULE, ONLY :soilstatetype, snowicestatetype
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water, &
											remove_substance_only, &
											remove_substance_only_simult, &
                                          error_remove_water, &
                                          add_water,&
										  add_substance_only
  USE NPC_SOIL_PROCESSES, ONLY : atmdep_in_loss, &
                                 doc_percolation_reduction, &
                                 onpp_percolation_reduction
  USE ATMOSPHERIC_PROCESSES, ONLY : deltasaturationpressure_function
  !Uses also modvar, hypevariables, hype_indata

  IMPLICIT NONE
  PRIVATE 
  PUBLIC initiate_soil_water_state, &
         initiate_soil_water, &
         calculate_snow, &
         latentheat_tempfunction, &
		 set_evaporation_concentrations, &
         calculate_potential_evaporation, &
         calculate_actual_soil_evapotranspiration, &
         calculate_tile_drainage, &
         calculate_soil_runoff, &
         infiltration, &
         infiltrationfrozen, &      !for CryoProcesses
         percolation, &
         add_macropore_flow, &
         calculate_groundwater_table, &
         calculate_snowdepth, &
         calculate_soiltemp , &
         calculate_weighted_temperature, &
         calculate_frostdepth, &
         calculate_soil_moisture_deficit, &
         snowalbedo_function, &
         calculate_fractional_snowcover, &
         calculate_interception,	&
	evaporation_percolation_soilrunoff_linear_reservoir
         
  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(19)  !error message for location of remove_water call, 1-6 not used
  PARAMETER (errstring = (/'regional groundwater flow, soillayer 3        ',    & !1
                           'regional groundwater flow, soillayer 2        ',    &
                           'regional groundwater flow, soillayer 1        ',    &
                           'evapotranspiration lake, less than lake volume',    &
                           'evapotranspiration lake, more than lake volume',    &
                           'evapotranspiration lake, slowlake part used   ',       &
                           'evapotranspiration soil, soillayer 1          ',      & !7
                           'evapotranspiration soil, soillayer 2          ',      & !8
                           'runoff from soillayer 1, stream in layer 1    ',      &
                           'runoff from soillayer 2, stream in layer 2    ',      &
                           'runoff from soillayer 1, stream in layer 2    ',      &
                           'runoff from soillayer 3, stream in layer 3    ',      &
                           'runoff from soillayer 1, stream in layer 3    ',      &
                           'runoff from soillayer 2, stream in layer 3    ',      &   !14
                           'tile runoff, drainage pipe in soillayer 1     ',      &
                           'tile runoff, drainage pipe in soillayer 2     ',      &
                           'tile runoff, drainage pipe in soillayer 3     ',      &  !17
                           'percolation from soillayer 1                  ',      &  
                           'percolation from soillayer 2                  '/))

CONTAINS

  !>\brief Initiate soil water state variables when no saved state exist.
  !!
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions)
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_soil_water_state(soilstate)

    USE HYPEVARIABLES, ONLY : m_wcfc,m_wcwp,m_wcep,  &
                              m_wcfc1,m_wcfc2,m_wcfc3, &
                              m_wcwp1,m_wcwp2,m_wcwp3, &
                              m_wcep1,m_wcep2,m_wcep3, &
                              m_fccorr,m_wpcorr,m_epcorr
    USE MODVAR, ONLY : classdata,     &
                       nsub,          &
                       nclass,        &
                       maxsoillayers, &
                       soiliniwet,    &
                       soilthick,     &
                       soilpar,       &
                       genpar
    
    !Argument declaration
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate   !<Soil states
    
    !Local variables
    INTEGER i,j           !loop-variables (subbasin,class)
    REAL :: epmm(maxsoillayers,nclass),wpmm(maxsoillayers,nclass),fcmm(maxsoillayers,nclass)

    !>\b Algoritm \n
    !>Calculate size of water storage in soil (wp,fc and ep) in mm
    DO j = 1,nclass
      IF(soilpar(m_wcfc1,classdata(j)%soil) > 0.)THEN                               !Field capacity (mm)
        fcmm(1,j)=soilpar(m_wcfc1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_fccorr)       !First layer
        fcmm(2,j)=soilpar(m_wcfc2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_fccorr)       !Second layer
        fcmm(3,j)=soilpar(m_wcfc3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_fccorr)       !Third layer
      ELSE
        fcmm(:,j)=soilpar(m_wcfc,classdata(j)%soil) * soilthick(:,j) * 1000.*genpar(m_fccorr)        !All layers
      ENDIF
      IF(soilpar(m_wcwp1,classdata(j)%soil) > 0.)THEN                               !Wilting point (mm)
        wpmm(1,j)=soilpar(m_wcwp1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_wpcorr)       !First layer
        wpmm(2,j)=soilpar(m_wcwp2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_wpcorr)       !Second layer
        wpmm(3,j)=soilpar(m_wcwp3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_wpcorr)       !Third layer
      ELSE
        wpmm(:,j)=soilpar(m_wcwp,classdata(j)%soil) * soilthick(:,j) * 1000.*genpar(m_wpcorr)        !All layers
      ENDIF
      IF(soiliniwet)THEN
        IF(soilpar(m_wcep1,classdata(j)%soil) > 0.)THEN                             !Effective porosity (mm)
          epmm(1,j)=soilpar(m_wcep1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_epcorr)     !First layer
          epmm(2,j)=soilpar(m_wcep2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_epcorr)     !Second layer
          epmm(3,j)=soilpar(m_wcep3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_epcorr)     !Third layer
        ELSE
          epmm(:,j)=soilpar(m_wcep,classdata(j)%soil) * soilthick(:,j) * 1000.*genpar(m_epcorr)      !All layers
        ENDIF
      ELSE
        epmm = 0.
      ENDIF
    ENDDO
    !>Initiate soil water to filed capacity (fc+wp)
    DO i = 1,nsub
      soilstate%water(:,:,i) = fcmm(:,:) + wpmm(:,:) + epmm(:,:) 
    ENDDO

  END SUBROUTINE initiate_soil_water_state

  !>\brief Initiate soil water help parameters.
  !!
  !>\b Consequences Module hypevariables epotdist, soilrc, basinrrcscorr, basincevpam, 
  !> basincevpph, basinlp, pwmm, wpmm, fcmm and epmm is set.
  !!
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions, Soil water -
  !! Groundwater runoff) and Processes above ground (Evaporation)
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_soil_water()

    USE HYPEVARIABLES, ONLY : epotdist, &         !OUT
                              soilrc,   &         !OUT
                              wpmm,fcmm,epmm, &   !OUT
                              pwmm,  & !OUT
                              basinrrcscorr,   & !OUT
                              basincevpam,   & !OUT
                              basincevpph,   & !OUT
                              basinlp,   & !OUT
                              m_epotdist,    &
                              m_cevpam,m_cevpph,m_lp,  &
                              m_wcfc,m_wcwp,m_wcep,  &
                              m_wcfc1,m_wcfc2,m_wcfc3, &
                              m_wcwp1,m_wcwp2,m_wcwp3, &
                              m_wcep1,m_wcep2,m_wcep3, &
                              m_rrcs1,m_rrcs2,m_rrcs3, &
                              m_rrcscorr,  &
                              n_rrcsc,n_rrcs3,n_cevpa,n_cevpp,n_lp, &
                              m_fccorr,m_wpcorr,m_epcorr,m_rrccorr, &
                              m_rrccorr
    USE MODVAR, ONLY : classdata,         &
                       basin,             &
                       nsub, nclass,      &
                       maxsoillayers,     &
                       soildepth,         &
                       soilthick,         &
                       genpar,soilpar,regpar, &
                       conductregest
    USE HYPE_INDATA, ONLY : clust_group,catdes, &
                            ndes_reg, idx_catdes, wght_catdes, indx_par
         
    !Local variables
    INTEGER i,j,k,isb           !loop-variables (subbasin,class)
    REAL    coeff
    REAL    sums
    REAL    rc0,rc1,rc2,b   !help variables for recession coefficient calculation

    !>\b Algoritm \n
    !>Calculate distribution of potential evaporation between soil layers
    IF(.NOT.ALLOCATED(epotdist)) ALLOCATE(epotdist(2,nclass))
    coeff=genpar(m_epotdist)
    DO j = 1, nclass
      sums = soilthick(1,j)*EXP(-coeff*soildepth(1,j)/2.) + soilthick(2,j)*EXP(-coeff*(soildepth(1,j)+(soildepth(2,j)-soildepth(1,j))/2.))
      epotdist(1,j) = soilthick(1,j)*EXP(-coeff*soildepth(1,j)/2.) / sums
    ENDDO
    epotdist(2,:) = 1 - epotdist(1,:)

    !>Initiate soil water content parameters
    IF(.NOT.ALLOCATED(wpmm)) ALLOCATE(wpmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(fcmm)) ALLOCATE(fcmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(epmm)) ALLOCATE(epmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(pwmm)) ALLOCATE(pwmm(maxsoillayers,nclass))
    DO j = 1,nclass
      IF(soilpar(m_wcfc1,classdata(j)%soil) > 0)THEN                               !Field capacity (mm)
        fcmm(1,j)=soilpar(m_wcfc1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_fccorr)       !First layer
        fcmm(2,j)=soilpar(m_wcfc2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_fccorr)       !Second layer
        fcmm(3,j)=soilpar(m_wcfc3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_fccorr)       !Third layer
      ELSE
        fcmm(:,j)=soilpar(m_wcfc,classdata(j)%soil) * soilthick(:,j) * 1000.*genpar(m_fccorr)        !All layers
      ENDIF
      IF(soilpar(m_wcwp1,classdata(j)%soil) > 0)THEN                               !Wilting point (mm)
        wpmm(1,j)=soilpar(m_wcwp1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_wpcorr)       !First layer
        wpmm(2,j)=soilpar(m_wcwp2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_wpcorr)      !Second layer
        wpmm(3,j)=soilpar(m_wcwp3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_wpcorr)       !Third layer
      ELSE
        wpmm(:,j)=soilpar(m_wcwp,classdata(j)%soil) * soilthick(:,j) * 1000.        !All layers
      ENDIF
      IF(soilpar(m_wcep1,classdata(j)%soil) > 0)THEN                               !Effectiv porosity (mm)
        epmm(1,j)=soilpar(m_wcep1,classdata(j)%soil) * soilthick(1,j) * 1000.*genpar(m_epcorr)       !First layer
        epmm(2,j)=soilpar(m_wcep2,classdata(j)%soil) * soilthick(2,j) * 1000.*genpar(m_epcorr)       !Second layer
        epmm(3,j)=soilpar(m_wcep3,classdata(j)%soil) * soilthick(3,j) * 1000.*genpar(m_epcorr)       !Third layer
      ELSE
        epmm(:,j)=soilpar(m_wcep,classdata(j)%soil) * soilthick(:,j) * 1000.        !All layers
      ENDIF
    ENDDO
    pwmm = wpmm + fcmm + epmm

    !Set soil runoff recession correction
    IF(.NOT.ALLOCATED(basinrrcscorr)) ALLOCATE(basinrrcscorr(nsub))
    DO i = 1,nsub
      IF(basin(i)%parregion>0)THEN
        basinrrcscorr(i) = 1. + regpar(m_rrcscorr,basin(i)%parregion)   !Correction of recession coefficients
      ELSE
        basinrrcscorr(i)  = 1.
      ENDIF
        
      !Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        IF(indx_par(n_rrcsc).gt.0) THEN
          basinrrcscorr(i) = 1.
          DO k=1,ndes_reg(clust_group(i),indx_par(n_rrcsc))
            basinrrcscorr(i) = basinrrcscorr(i) + wght_catdes(clust_group(i),indx_par(n_rrcsc),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rrcsc),k))
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    !Initiate soil runoff recession coeffcients
    IF(.NOT.ALLOCATED(soilrc)) ALLOCATE(soilrc(maxsoillayers,nclass,nsub))
    !>Calculate adjustment factors
    rc2 = genpar(m_rrcs3)*genpar(m_rrccorr)                         !Runoff coefficient slope dependence   
    DO i = 1,nsub
      !>Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        IF(indx_par(n_rrcs3).gt.0)THEN
          rc2=0.0
          DO j=1,ndes_reg(clust_group(i),indx_par(n_rrcs3))
            rc2 = rc2 + wght_catdes(clust_group(i),indx_par(n_rrcs3),j)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rrcs3),j))
          ENDDO
        ENDIF
      ENDIF          

      !>Calculate soil runoff recession coeffcients for each soil layer, class and subbasin
      DO j=1,nclass
        rc0 = soilpar(m_rrcs1,classdata(j)%soil)*basinrrcscorr(i)*genpar(m_rrccorr)       !Runoff coefficient in surface layer 
        IF(rc0>1.) rc0 = 1.
        rc0 = rc0+rc2*basin(i)%slope    !runoff coefficient in upper soil layer (slope dependent)
        IF(rc0>1.) rc0 = 1.
        rc1 = soilpar(m_rrcs2,classdata(j)%soil)*basinrrcscorr(i)*genpar(m_rrccorr)       !Runoff coefficient in bottom layer 
        IF(rc1>1.) rc1 = 1.
        IF(rc1==0) rc1 = rc0
        b = LOG(rc0/rc1)/((soildepth(3,j)-soilthick(3,j)/ 2.) - soilthick(1,j)/ 2.)
        soilrc(1,j,i) = rc0
        soilrc(3,j,i) = rc1
        soilrc(2,j,i) = rc0 * EXP (-b* (soildepth(2,j) - soilthick(2,j)/2. - soilthick(1,j)/ 2.))
      ENDDO
    ENDDO
    
    !>Set evaporation sinus corrections
    IF(.NOT.ALLOCATED(basincevpam))THEN
      ALLOCATE(basincevpam(nsub))
      ALLOCATE(basincevpph(nsub))
    ENDIF
    basincevpam = genpar(m_cevpam)
    basincevpph = genpar(m_cevpph)
    !Replace parameter values with regional parameter estimates
    IF(conductregest)THEN
      IF(indx_par(n_cevpa).gt.0)THEN
        basincevpam = 0.
        DO isb = 1,nsub
          DO k=1,ndes_reg(clust_group(isb),indx_par(n_cevpa))
            basincevpam(isb) = basincevpam(isb) + wght_catdes(clust_group(isb),indx_par(n_cevpa),k)*catdes(isb,idx_catdes(clust_group(isb),indx_par(n_cevpa),k))
          ENDDO
        ENDDO
      ENDIF
      IF(indx_par(n_cevpp).gt.0)THEN
        basincevpph = 0.
        DO isb = 1,nsub
          DO k=1,ndes_reg(clust_group(isb),indx_par(n_cevpp))
            basincevpph(isb) = basincevpph(isb) + wght_catdes(clust_group(isb),indx_par(n_cevpp),k)*catdes(isb,idx_catdes(clust_group(isb),indx_par(n_cevpp),k))
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    !>Set evaporation subbasin parameter
    IF(.NOT.ALLOCATED(basinlp)) ALLOCATE(basinlp(nsub))
    basinlp = genpar(m_lp)
    IF(conductregest)THEN   !Replace parameter value with regional parameter estimates
      IF(indx_par(n_lp).gt.0)THEN
        basinlp=0.0
        DO isb = 1,nsub
          DO k=1,ndes_reg(clust_group(isb),indx_par(n_lp))
            basinlp(isb) = basinlp(isb) + wght_catdes(clust_group(isb),indx_par(n_lp),k)*catdes(isb,idx_catdes(clust_group(isb),indx_par(n_lp),k))
          ENDDO
        ENDDO
      ENDIF
    ENDIF      

  END SUBROUTINE initiate_soil_water

!>Function to calculate snow albedo depending on the snow age        
!>
!> \b Reference ModelDescription Chapter Land routines (Snow routines)
!------------------------------------------------------------------------
  FUNCTION snowalbedo_function(snowage,albmin,albmax,kexp) RESULT(albedo)
  
  !Argument declarations
  REAL,INTENT(IN)  :: snowage !<snow age (timesteps)
  REAL,INTENT(IN)  :: albmin  !<minimum albedo (typical value 0.4)
  REAL,INTENT(IN)  :: albmax  !<maximum albedo (typical value 0.9)
  REAL,INTENT(IN)  :: kexp    !<exponential factor (1/time step) (typical value 0.1 for daily timesteps)
  REAL             :: albedo  ! albedo, fractional reflection of shortwave radiation (-)
  !David Gustafsson, 2013-02-05
  
  !Calculate albedo with a simplified exponential function
  albedo = albmin+(albmax-albmin)*EXP(-kexp*snowage)
  
  END FUNCTION snowalbedo_function

  !>\brief Subroutine for calculation of fractional snow cover area
  !>
  !Based on Lindström&Gardelin(1999;2000) following the implementation in 
  !the Rossby centre RCA-model (Samuelsson et al, 2006)
  !
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Snow cover)
  !------------------------------------------------------------------------
  SUBROUTINE calculate_fractional_snowcover(iluse,elevstd,snow,snowmax,fsc)
  
    USE MODVAR, ONLY : landpar, genpar, seconds_per_timestep
    USE HYPEVARIABLES, ONLY :  m_fscmax,m_fscmin,m_fsclim, & 
                               m_fscdistmax, m_fscdist0,m_fscdist1, &
                               m_fsck1, m_fsckexp
    !Argument declarations
    INTEGER, INTENT(IN) :: iluse       !<index of landuse
    REAL, INTENT(IN)    :: elevstd     !<standard deviation of elevation (m)
    REAL, INTENT(IN)    :: snow        !<snow pack (mm)
    REAL, INTENT(INOUT) :: snowmax     !<maximum snow pack during winter (mm)
    REAL, INTENT(OUT)   :: fsc         !<fractional snowcover area (-)

    !Local variables
    REAL timestep_seconds,fscdist
    timestep_seconds = REAL(seconds_per_timestep)
    
    !Check snow status
    IF(snow.gt.0.)THEN   !Snow present
      !Check snowcover model
      IF(genpar(m_fscmax)==0)THEN
        fsc = 1.
      ELSE
        !Check snowpack development phase, and select corresponding FSC function
        IF(snowmax.le.0.)THEN
          !1) Accumulation phase, snowmax = 0
          !1.1) fsc = tangens-hyperbolic function, Eq 28 (Samuelsson 2006)
          fsc = MAX(genpar(m_fscmin),genpar(m_fscmax) * TANH(0.1 * snow))
          !1.2) Set snowmax = snow, if fsc >= fscmax - fsclim
          IF(fsc.ge.(genpar(m_fscmax)-genpar(m_fsclim)))THEN
            snowmax = snow
          ENDIF
        ELSE
          !2) Melting phase, snowmax>0 (onset in previous timesteps)
          !2.1) update snowmax
          IF(snow.GT.snowmax)THEN
            !update snowmax to new maximum snow value
            snowmax = snow
          ELSE
            !decrease snowmax towards end of melt season, eq. 31 (Samuelsson 2006)
            IF(snow.LT.genpar(m_fsck1)*snowmax)THEN
              snowmax = snowmax - (genpar(m_fsck1) * snowmax - snow)*(1.-EXP(-genpar(m_fsckexp) * timestep_seconds)) / genpar(m_fsck1)
            ENDIF 
          ENDIF
          !2.2) calculate snow distribution factor, Eq 30 (Samuelsson 2006)
          fscdist = MIN(landpar(m_fscdistmax,iluse),landpar(m_fscdist0,iluse) + landpar(m_fscdist1,iluse) * elevstd)
          !2.3) fsc=linear function, Eq 29 (Samuelsson 2006)
          fsc = MAX(genpar(m_fscmin),MIN(genpar(m_fscmax),snow / (snowmax * fscdist)))
        ENDIF
      ENDIF 
    ELSE  !No snow
      snowmax = 0.
      fsc = 0.
    ENDIF
    
  END SUBROUTINE calculate_fractional_snowcover
  
  !>Subroutine for calculation of interception and evaporation of intercepted snow and rain
  !!
  !> \b Reference ModelDescription Chapter Land routines (Interception routines)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_interception(iluse,snowfall,rainfall,cprec,snowfallthrough,rainfallthrough,csnowfallthrough,crainfallthrough,evapintprec,cevapintprec)
    USE MODVAR, ONLY : landpar, &
         numsubstances
    USE HYPEVARIABLES, ONLY : m_pcluse,m_rncluse,m_sncluse
    
    !ARGUMENTS
    INTEGER, INTENT(IN) :: iluse
    REAL, INTENT(IN)    :: snowfall
    REAL, INTENT(IN)    :: rainfall
    REAL, INTENT(IN)    :: cprec(numsubstances)
    REAL, INTENT(OUT)   :: snowfallthrough
    REAL, INTENT(OUT)   :: rainfallthrough
    REAL, INTENT(OUT)   :: csnowfallthrough(numsubstances)
    REAL, INTENT(OUT)   :: crainfallthrough(numsubstances)
    REAL, INTENT(OUT)   :: evapintprec
    REAL, INTENT(OUT)   :: cevapintprec(numsubstances)
    
    !LOCAL VARIABLES
    REAL evapsnowfrac, evaprainfrac, evapintsnow, evapintrain
      
    !First VERY simple function - just take a fraction of precipitation (sncluse and rncluse (or pcluse)) as evaporation of intercepted snow or rain
    
    !Parameters (PCLUSE have priority over SNCLUSE and RNCLUSE)
    IF(landpar(m_pcluse,iluse).GT.0.)THEN 
      evapsnowfrac = landpar(m_pcluse,iluse)
      evaprainfrac = landpar(m_pcluse,iluse)
    ELSE
      evapsnowfrac = landpar(m_sncluse,iluse)
      evaprainfrac = landpar(m_rncluse,iluse)
    ENDIF
    
    !Evaporation of intercepted rainfall and snowfall
    evapintsnow=snowfall * evapsnowfrac
    evapintrain=rainfall * evaprainfrac
    evapintprec=evapintsnow+evapintrain
    
    !Throughfall
    snowfallthrough = snowfall - evapintsnow
    rainfallthrough = rainfall - evapintrain
    
    !Concentration of elements in evaporation and throughfall
    IF(numsubstances.GT.0)THEN
      csnowfallthrough(numsubstances)=0.
      crainfallthrough(numsubstances)=0.
      cevapintprec(numsubstances)=0.
      !Concentration in throughfall rain
      IF(rainfallthrough.GT.0.)THEN
        !adjust concentration of non-evaporating substances
        crainfallthrough(:)=cprec(:)*rainfall/rainfallthrough
        !keep concentrations of evaporating substances (1)
        crainfallthrough(1)=cprec(1)
      ELSE
        crainfallthrough(:)=0.
      ENDIF
      !concentration in throughfall snow
      IF(snowfallthrough.GT.0.)THEN
        !adjust concentration of non-sublimating substances
        csnowfallthrough(:)=cprec(:)*snowfall/snowfallthrough
        !keep concentrations of sublimating substances 1
        csnowfallthrough(1)=cprec(1)
      ELSE
        csnowfallthrough(:)=0.
      ENDIF
      !Concentration in the intercepted precipitation evaporation
      IF(evapintprec.GT.0.)THEN
        cevapintprec(1) = cprec(1) 
      ENDIF
    ENDIF
  
  END SUBROUTINE calculate_interception
  
  !>Subroutine for calculation of changes in snow pack; snowfall addition, 
  !snow pack melting and snow age
  !!
  !> \b Reference ModelDescription Chapter Land routines (Snow routines)
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snow(iluse,snowfall,csnowfall,snow,csnow,temp,melt,cmelt,swrad,snowage,snowcover,epot,evap,cevap,effcov)
  
    USE MODVAR, ONLY : landpar,   &
         genpar,          &
         numsubstances,   &
         missing_value,   &
         modeloption,     &
         p_snowmelt,      &
         i_t2,            &
         i_t1,            &
         p_snowevap,      &
         i_sm
    
    USE HYPEVARIABLES, ONLY : m_ttmp,m_cmlt,m_snalbmin,m_snalbmax,m_snalbkexp,m_cmrad,m_fsceff,m_cmrefr,m_cmltcorr
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iluse      !<index of landuse
    REAL, INTENT(IN)    :: snowfall   !<precipitation as snow (mm/timestep)
    REAL, INTENT(IN)    :: csnowfall(numsubstances) !<concentration of precipitation as snow 
    REAL, INTENT(INOUT) :: snow       !<snow pack (mm)
    REAL, INTENT(INOUT) :: csnow(numsubstances) !<concentration of snow 
    REAL, INTENT(IN)    :: temp       !<air temperature (C)
    REAL, INTENT(OUT)   :: melt       !<snow melt (mm/timestep)
    REAL, INTENT(OUT)   :: cmelt(numsubstances)     !<substances of snow melt
    REAL, INTENT(IN)    :: swrad      !<shortwave radiation (MJ/m2/day?)
    REAL, INTENT(INOUT) :: snowage    !<age of snow (timesteps)
    REAL, INTENT(IN)    :: snowcover  !<snowcover fraction
    REAL, INTENT(IN)    :: epot       !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap       !<snow sublimation (mm/timestep)
    REAL, INTENT(OUT)   :: effcov     !<effective snowcover used for scaling snow and soil evaporation (0 if snowevap is switched off)
    REAL, INTENT(OUT)   :: cevap(numsubstances)   !<concentrations in snow sublimation
    
    !Local variables
    REAL tt       !threshold temperature for snow melt (and evaporation) (C)
    REAL cm       !coefficient for snow melt (mm/C/timestep)
    REAL newsnow
    REAL snowalbedo
    REAL snalbmax, snalbmin, snalbkexp
    REAL cmrad     !radiation index snow melt factor (mm/MJ/timestep)
    REAL fsceff,abla,abla0,cmrefr  

    !Set parameter values
    tt = landpar(m_ttmp,iluse)  !threshold temperature for snow melt
    cm = landpar(m_cmlt,iluse)*genpar(m_cmltcorr)  !Coefficient for snow melt
    fsceff = genpar(m_fsceff)   !efficiency of fractional snow cover to reduce melt and evap
    effcov = 1.-fsceff*(1.-snowcover) !effective snow cover used to scale melt and evap
    
    !Melting and Sublimation, select model
    ! -> To reduce the number of snowmelt options, snow cover melt scaling 
    !    is now included in all snowmelt models. Instead of options for each combination
    !    of melt and sublimation, the fraction of snow cover reduction is controlled
    !    by a new parameter ffscred (fraction of fsc reduction)
    ! -> For backward compitability, the previous snowmelt options values are still used 
    !    (0,1 temp index, 2 temp+rad index) - but note that p_snowmelt = 0 may now also 
    !    imply snowcover scaling and sublimation, depending on parameter ffscred and fepotsnow.
    ! -> Sublimation is calculated separately after the melt section, and is also controlled by ffscred.
    ! -> Ablation = melt + sublimation is introduced
    ! -> Minimization of ablation to current snow is made after calculation of (potential) 
    !    melt and sublimation. The reduction from potential to actual ablation is finally
    !    distributed on melt and sublimation.
    SELECT CASE(modeloption(p_snowmelt))
    CASE(0,1) ! Original temperature index model, with/without snowcover scaling
      IF(snow>0 .AND. temp >= tt) THEN
        melt = cm   * (temp - tt)  !potential melt
        melt = melt * effcov       !snowcover melt scaling (no reduction of snowcover=1 and/or fsceff=0)
      ELSE
        melt = 0.
      ENDIF
    CASE(2) ! Temperature AND Radiation index model, with/without snowcover scaling and sublimation
      !Set parameter values
      snalbmin  = landpar(m_snalbmin,iluse)
      snalbmax  = landpar(m_snalbmax,iluse)
      snalbkexp = landpar(m_snalbkexp,iluse)
      cmrad     = landpar(m_cmrad,iluse)
      cmrefr    = genpar(m_cmrefr)
    
      !Radiation melt component
      snowalbedo = snowalbedo_function(snowage,snalbmin,snalbmax,snalbkexp)
      melt = cmrad * swrad * (1.-snowalbedo)
      
      !Add Temperature component
      IF(snow>0. .AND. temp >= tt)THEN
        melt = melt + cm * (temp - tt)
      ENDIF
      
      !Refreezing component when temperatures below tt, as a fraction cmrefr of cm 
      IF(snow>0. .AND. temp < tt .AND. melt > 0.)THEN
        melt = melt - cmrefr * cm * (tt - temp)
        IF(melt<0.) melt = 0.
      ENDIF
      
      !Scale melt with fractional snow cover
      melt = melt * effcov
    CASE DEFAULT ! Original temperature index model, with snowcover scaling
      IF(snow>0 .AND. temp >= tt) THEN
        melt = cm   * (temp - tt)  !potential melt
        melt = melt * effcov       !snowcover melt scaling
      ELSE
        melt = 0.
      ENDIF
    END SELECT
    
    !Evaporation (sublimation), with/without snowcover scaling
    IF(modeloption(p_snowevap).GE.1)THEN
      IF(snow>0)THEN
        evap = epot * effcov
      ELSE
        effcov = 0. !make sure effcov = 0 if there is no snow, otherwise there will be no soil evaporation
        evap   = 0.
      ENDIF
    ELSE
      evap   = 0.
      effcov = 0. !make sure effcov = 0 if snowevap is switched off, otherwise there will be no soil evaporation
    ENDIF

    !Ablation = Melt + Sublimation
    abla0 = melt + evap                     !potential ablation (melt+evap)
    abla = MIN(abla0, snow)                 !minimize ablation to available snow
    IF(abla0.GT.0.)THEN
      IF(abla<abla0)THEN
        melt = melt * abla/abla0              !distribute ablation on melt and evap
        evap = evap * abla/abla0              !distribute ablation on melt and evap
      ENDIF
    ELSE
      melt = 0.
      evap = 0.
    ENDIF

    !New snow water equivalent, after snowfall, melting, and sublimation
    newsnow = snow + snowfall  - melt - evap
    IF(newsnow<0.) newsnow = 0.
    
    !Update concentrations of substances in snow 
    IF(numsubstances.GT.0)THEN
      !Meltwater concentrations
      IF(snow > 0.) THEN
        cmelt(:) = csnow(:)
      ELSE
        cmelt(:) = 0.
      ENDIF
      IF(i_t2>0) cmelt(i_t2)=0. !temp.conc. in meltwater = 0
      IF(i_sm>0) cmelt(i_sm)=1. !snowmelt concentration in meltwater is 1!
      !Evaporation concentrations
      IF(evap.GT.0)THEN
          cevap(:)=0.
          IF(i_t1>0)THEN
            cevap(i_t1)= csnow(i_t1) !018
          ENDIF
      ENDIF
      !Concentrations in snow, after snowfall, melting and sublimation
      IF(newsnow==0.)THEN
        csnow(:) = missing_value             
      ELSE
        !csnow(:) = (csnow(:)*(snow-melt) + csnowfall(:)*snowfall) / newsnow
        csnow(:) = (csnow(:)*snow - cmelt(:)*melt -cevap(:)*evap + csnowfall(:)*snowfall) / newsnow
      ENDIF
    ENDIF
    !Update snow water equivalent to new value
    snow = newsnow

  END SUBROUTINE calculate_snow
  
  !>Calculate air pressure (kPa) as a function of elevation, FAO(7)
  !-------------------------------------------------------------------------------
  REAL FUNCTION airpressure_elevationfunction(elev)
     
    !Argument decalaration
    REAL, INTENT(IN) :: elev   !<elevation
     
    airpressure_elevationfunction = 101.3 * ((293. - 0.0065 * elev)/293. ) ** 5.26
     
  END FUNCTION airpressure_elevationfunction
  
  !>Calculate latent heat of vaporization (MJ kg-1) as a function of temperature
  !-------------------------------------------------------------------------------
  REAL FUNCTION latentheat_tempfunction(temp)

    !Argument decalaration
    REAL, INTENT(IN) :: temp !<temperature (C)
     
    latentheat_tempfunction = 2.501 - 0.002361 * temp  !MJ/kg
     
  END FUNCTION latentheat_tempfunction
  
  !> Calculate psychrometric constant (kPa C^-1) as a function of
  !! pressure(elevation) and lambda (temperature), FAO
  !-------------------------------------------------------------------------------
  REAL FUNCTION psychrometric_constant(pa,lambda)
  
    !Argument decalarations
    REAL, INTENT(IN) :: pa       !<air pressure [kPa]
    REAL, INTENT(IN) :: lambda   !<latent heat of vaporaization
    
    !Parameter declaration
    REAL, PARAMETER :: cp = 0.001013  ! specific heat of moist air at constant pressure (MJ kg^-1 C^-1)
 
    psychrometric_constant = cp * pa / (0.622 * lambda)
 
  END FUNCTION psychrometric_constant

  !>Calculates potential evaporation or uses value supplied as input
  !
  !> \b Reference ModelDescription Processes above ground (Evaporation)
  !--------------------------------------------------------------
  SUBROUTINE calculate_potential_evaporation(i,j,temp,epot,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
  
    USE MODVAR, ONLY : basin,classdata, &
                       landpar,         &
                       genpar,          &
                       dayno,           &
                       pi,              &
                       xobsi,xobsindex, &
                       modeloption,p_petmodel,  &
                       classbasin
    USE HYPEVARIABLES, ONLY : o_reepot,   &
                              m_ttmp,m_cevp,  &
                              basincevpam,  &
                              basincevpph,  &
                              m_krs,m_kc,m_jhtadd,m_jhtscale,m_alfapt,m_fepotsnow, &
                              m_kccorr, m_fpsnocorr

    !Argument declarations
    INTEGER, INTENT(IN) :: i      !<index of current subbasin
    INTEGER, INTENT(IN) :: j      !<index of current class 
    REAL, INTENT(IN)    :: temp   !<air temperature
    REAL, INTENT(OUT)   :: epot   !<potential evapotranspiration mm/timestep
    REAL, INTENT(IN)    :: radext !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: swrad  !<downward shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: netrad !<net downward radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: actvap !<actual vapor pressure [kPa]
    REAL, INTENT(IN)    :: satvap !<saturated vapour pressure [kPa]
    REAL, INTENT(IN)    :: wind   !<wind speed [m/s]
    REAL, INTENT(OUT)   :: epotsnow !<potential evapotranspiration for snow mm/timestep
    
    !Local variables
    REAL tt       !threshold temperature for melting (C)
    REAL ce       !coefficient for evaporation (mm/C/timestep)
    REAL dsatvap ! Slope of saturation pressure curve [kPa/C]
    REAL gamma   ! psychrometric constant
    REAL lambda  ! latent heat of evaporation [MJ/kg]
    REAL pa      ! atmospheric pressure [kPa]
    REAL kc      ! crop coefficient used for the new PET functions
    REAL elev    ! elevation
    REAL turbidity ! atmospheric turbidity
    REAL fepotsnow ! fraction of potential evaporation used for snow
    
    !>\b Algorithm \n
    !>Set local parameters and corrections
    tt=landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation
    ce=landpar(m_cevp,classdata(j)%luse)       !Coefficient for potential evaporation
    fepotsnow = landpar(m_fepotsnow,classdata(j)%luse)*genpar(m_fpsnocorr)       !Coefficient for potential evaporation for snow
    ce = ce * (1 + basincevpam(i)*SIN(2.*pi*(dayno-basincevpph(i))/365.))
!    epotcorr = 1 + genpar(m_cevpam)*SIN(2.*pi*(dayno-genpar(m_cevpph))/365.) 
    
    !>Calculate additional input variables for the optional PET functions
    IF(modeloption(p_petmodel).GT.1)THEN
      !Slope of saturated vapour pressure curve, using mean temperature
      dsatvap = deltasaturationpressure_function(temp)
      !Latent heat of vaporization
      lambda = latentheat_tempfunction(temp)
      !Elevation
      elev = basin(i)%elev+classbasin(i,j)%deltah
      !Air pressure, assuming normal pressure at sea level
      pa = airpressure_elevationfunction(elev)
      !Psychrometric constant
      gamma = psychrometric_constant(pa,lambda)
      !Landuse scaling parameter, "crop coefficient"
      kc = landpar(m_kc,classdata(j)%luse)*genpar(m_kccorr)
      !Turbidity
      turbidity = swrad / radext
    ENDIF
      
    !>Calculate potential evaporation 
    SELECT CASE(modeloption(p_petmodel))
      CASE(0) !HYPE original model (with Xobs replacement, if available)
        IF(xobsindex(o_reepot,i)>0)THEN       
          epot = xobsi(xobsindex(o_reepot,i))
        ELSEIF(temp>tt)THEN
          epot = ce*(temp-tt)
!          epot = ce*(temp-tt)*epotcorr
        ELSE
          epot = 0.  
        ENDIF
      CASE(1) !HYPE original model (without Xobs replacement)
        IF(temp>tt)THEN
          epot = ce*(temp-tt)
!          epot = ce*(temp-tt)*epotcorr
        ELSE
          epot = 0.  
        ENDIF
      CASE(2) !Modified Jensen-Haise/McGuinness following Oudin et al (2005)
        !parameters suggested by Oudin et al, jhtadd = 5, jhtscale = 100
        epot = kc * max(0.,radext / (lambda) * (temp + genpar(m_jhtadd)) / genpar(m_jhtscale))
      CASE(3) !Hargreaves-Samani (known to overpredict in humid areas)
        ! The function is modified by DG to limit the "turbidity-factor" with the Ångström formula:
        ! 
        !   The original Hargreaves function is:
        !     epot = 0.0023 * radext / (lambda*rho) * (Tmax-Tmin)^0.5 * (temp + 17.8)
        !   and the Hargreaves turbidity for estimating swrad = krs * (Tmax-Tmin)^0.5
        !
        !   Thus, by replacing (Tmax-Tmin)^2 with turbidity/krs, we get a reasonable limitation of the Hargreaves (tmax-Tmin) impact
        !  (furthermore, if Tmax-min was missing, we actually use the clearsky turbidity at this point)
        !
        ! also note that rho = 1 and excluded in equations below...
        epot = max(0.,kc * 0.0023 * radext /(lambda) * turbidity / genpar(m_krs) * (temp+17.8))
      CASE(4) ! Priestly Taylor (known to underpredict in arid and semi-arid areas)
        epot = max(0.,kc * genpar(m_alfapt) * dsatvap * netrad / (lambda * (dsatvap+gamma)))
      CASE(5) ! FAO Penman Monteith reference crop evapotranspiration
        epot = max(0., kc * ((0.408 * dsatvap * netrad + gamma*900./(temp+273.)*wind*(satvap-actvap))/(dsatvap+gamma*(1.+0.34*wind))))
      CASE DEFAULT !TODO: print out a warning that PETMODEL don't have a useful number 0-5 and that original HYPE is used
        !HYPE original model (with Xobs replacement, if available)
        IF(xobsindex(o_reepot,i)>0)THEN       
          epot = xobsi(xobsindex(o_reepot,i))
        ELSEIF(temp>tt)THEN
          epot = ce*(temp-tt)
!          epot = ce*(temp-tt)*epotcorr
        ELSE
          epot = 0.  
        ENDIF  
      END SELECT
      !Potential evaporation used for snow evaporation (sublimation)
      epotsnow = fepotsnow * epot

  END SUBROUTINE calculate_potential_evaporation
  
  !>\brief Calculate and set concentration of evaporating water
  !
  !> \b Reference ModelDescription Chapter Processes above ground (Evaporation)
  !--------------------------------------------------------------------
  SUBROUTINE set_evaporation_concentrations(conc,cevap)
  
    USE MODVAR, ONLY : numsubstances

    !Argument declarations
    REAL, INTENT(IN)  :: conc(numsubstances)     !<concentration in water body (?)
    REAL, INTENT(OUT) :: cevap(numsubstances)    !<concentration in evapotranspiration (?)

    cevap = 0.
    IF(numsubstances == 0) RETURN
    !IF(simulate%substance(i_t1)) cevap(i_t1) = genpar(m_T1evap) * conc(i_t1) !!<concentration in water body (?)  !This version of HYPE doesn't make concentration adjustment to the evaporating substance
    !IF(simulate%substance(i_t2)) cevap(i_t2) = conc(i_t2)
    cevap(:) = conc(:)


  END SUBROUTINE set_evaporation_concentrations

  !>\brief Calculate and remove evapotranspiration from the soil upper
  !>two layers
  !
  !> \b Reference ModelDescription Chapter Processes above ground (Evaporation)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_actual_soil_evapotranspiration(i,j,temp,epot,wp,fc,  &
                                      epotfrac,soilstate,evap,evapflows,cevap, &
                                       barefrac,soiltemp)

    USE MODVAR, ONLY : basin,classdata, &
                       landpar,  &
                       numsubstances,   &
                       maxsoillayers,   &
                       i_t1,            &
                       soilthick,       &
                       i_sm,i_gm,i_rn
    USE HYPEVARIABLES, ONLY : basinlp, &
                              m_ttmp, & 
                              m_ttrig, &
                              m_tredA, &
                              m_tredB

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<index of current subbasin
    INTEGER, INTENT(IN) :: j                        !<index of current class
    REAL, INTENT(IN)    :: temp                     !<air temperature
    REAL, INTENT(IN)    :: epot                     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<field capacity (mm)
    REAL, INTENT(IN)    :: epotfrac(2)              !<relative distribution of potential evaporation between upper two soil layers (-)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: evap                     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evapflows(2)             !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubstances)     !<concentration in evapotranspiration (?)
    REAL, INTENT(IN)    :: barefrac                 !<fraction of soil that has evapotranspiration (-)
    REAL, INTENT(IN)    :: soiltemp(maxsoillayers)  !<soil temperature (deg)

    !Local variables
    INTEGER k   !loop-variable
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for melting (C)
    REAL evap1,evap2  !evapotranspiration of soillayer 1 and 2 (mm/timestep)
    REAL cevap1(numsubstances),cevap2(numsubstances)  !concentration of evapotranspiration
    REAL soiltemp_reduction(maxsoillayers)
    
    !>\b Algorithm \n
    !Default values output variables
    evap = 0.
    cevap = 0.
    evapflows = 0.

    !>Set local parameters
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    !>Set Soil Temperature reduction (only if tredA>0)
    soiltemp_reduction=1.
    IF(landpar(m_tredA,classdata(j)%luse).GT.0.)THEN
      DO k=1,MIN(2,maxsoillayers)
        IF(k.EQ.1 .OR. soilthick(MIN(2,maxsoillayers),j)>0)THEN
          IF(soiltemp(k).GT.landpar(m_ttrig,classdata(j)%luse))THEN
            soiltemp_reduction(k) = 1.-EXP(-landpar(m_tredA,classdata(j)%luse) * (soiltemp(k)-landpar(m_ttrig,classdata(j)%luse))**landpar(m_tredB,classdata(j)%luse))
          ELSE
            soiltemp_reduction(k) = 0.
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    
    !>If temperature above threshold:
    IF(temp>tt)THEN

       !>\li calculate actual evapotranspiration in the uppermost layer 
      IF(soilstate%water(1,j,i) - wp(1)> basinlp(i) * fc(1)) THEN
        evap1 = epot*epotfrac(1)*soiltemp_reduction(1)
      ELSEIF(soilstate%water(1,j,i)-wp(1) <= 0.0) THEN
        evap1 = 0.0
      ELSE
        evap1 = epot*epotfrac(1)*((soilstate%water(1,j,i)-wp(1))/(basinlp(i) * fc(1)))*soiltemp_reduction(1)
      ENDIF
      IF(evap1>soilstate%water(1,j,i)-wp(1)) evap1 = soilstate%water(1,j,i)-wp(1)
      DO k=1,numsubstances
        cevap1(k) = 0.
        IF(k==i_t1) cevap1(k) = soilstate%conc(k,1,j,i)  !t1 == O18 i den här modellen
        IF(k==i_sm) cevap1(k) = soilstate%conc(k,1,j,i)  !sm == snowmelt, water origin trace element model
        IF(k==i_gm) cevap1(k) = soilstate%conc(k,1,j,i)  !gm == glaciermelt, water origin trace element model
        IF(k==i_rn) cevap1(k) = soilstate%conc(k,1,j,i)  !rn == rainfall, water origin trace element model
      ENDDO
      evap1 = evap1*MIN(1.,barefrac)  !Scale evapotranspiration with fraction of bare soil

      !>\li Remove evapotranspiration of soillayer 1
      CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),evap1,cevap1,status)
      IF(status.NE.0) CALL error_remove_water(errstring(7),basin(i)%subid,i,j)

      !Second soillayer:
      IF(soilthick(2,j)>0)THEN

        !>\li Calculate actual evapotranspiration in the second soillayer
        IF(soilstate%water(2,j,i)-wp(2) > basinlp(i) * fc(2)) THEN
          evap2 = epot*epotfrac(2)*soiltemp_reduction(2)
        ELSEIF(soilstate%water(2,j,i)-wp(2) <= 0.0) THEN
          evap2 = 0.0
        ELSE
          evap2 = epot*epotfrac(2)*((soilstate%water(2,j,i)-wp(2))/(basinlp(i) * fc(2)))*soiltemp_reduction(2)
        ENDIF
        IF(evap2>soilstate%water(2,j,i)-wp(2)) evap2 = soilstate%water(2,j,i) - wp(2)
        DO k=1,numsubstances
          cevap2(k) = 0.
          IF(k==i_t1) cevap2(k) = soilstate%conc(k,2,j,i)  !t1 == O18 i den här modellen
          IF(k==i_sm) cevap2(k) = soilstate%conc(k,2,j,i)  !sm == snowmelt, water origin trace element model
          IF(k==i_gm) cevap2(k) = soilstate%conc(k,2,j,i)  !gm == glaciermelt, water origin trace element model
          IF(k==i_rn) cevap2(k) = soilstate%conc(k,2,j,i)  !rn == rainfall, water origin trace element model
        ENDDO
        evap2 = evap2*MIN(1.,barefrac)  !Scale evapotranspiration with fraction of bare soil

        !>\li Remove evapotranspiration of soillayer 2
        CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),evap2,cevap2,status)
        IF(status.NE.0) CALL error_remove_water(errstring(8),basin(i)%subid,i,j)
      ELSE
        evap2 = 0.
        cevap2 = 0.
      ENDIF

      !>Set output variables
      evap = evap1 + evap2
      evapflows(1) = evap1
      evapflows(2) = evap2
      IF(i_t1>0 .AND. evap>0) cevap(i_t1) = (cevap1(i_t1)*evap1 + cevap2(i_t1)*evap2)/evap
      IF(i_sm>0 .AND. evap>0) cevap(i_sm) = (cevap1(i_sm)*evap1 + cevap2(i_sm)*evap2)/evap
      IF(i_gm>0 .AND. evap>0) cevap(i_gm) = (cevap1(i_gm)*evap1 + cevap2(i_gm)*evap2)/evap
      IF(i_rn>0 .AND. evap>0) cevap(i_rn) = (cevap1(i_rn)*evap1 + cevap2(i_rn)*evap2)/evap
    ENDIF

  END SUBROUTINE calculate_actual_soil_evapotranspiration
  
  
  
  
  
  
  
  
  SUBROUTINE evaporation_percolation_soilrunoff_linear_reservoir(i,j,temp,epot,wp,fc,epotfrac,soilstate,evap,evapflows,cevap,barefrac, &
																	isoil,subid,ep,sthick,infilt,cinfilt,percflow,cpercflow,soilrunoff,csoilrunoff)							!!modified to simultaneously combine evaporation, percolation and soilrunoff calculation by BTW

    USE MODVAR, ONLY : basin,soilthick,realzero, &
	                   numsubstances,   &
                       maxsoillayers,   &
                       genpar,soilpar,landpar, &            
                       classdata,		&
					   dayno,			&														!!Added by BTW JUly11_2020
					   currentdate,		&														!!Added by BTW JUly12_2020
						timesteps_per_day														!!Added by BTW Feb24_2021

   USE HYPEVARIABLES, ONLY : basinlp, &
                              m_ttmp,  &
                              m_ttrig,m_tredA,m_tredB, &
							  m_perc1,m_perc2,    &
                              m_crate5,   &
                              m_onpercred, m_pppercred, &
							  soilrc

	USE gear_GlobVARs							!Added by BTW
	USE GEAR_IMPLICIT							!Added by BTW
	USE uawp, ONLY : dist_max, norm_max, awp	!Added by BTW
							 
	USE t_dgls, ONLY : 	dgl						!Added by BTW
	USE fgauss, ONLY : gauss					!Added by BTW

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    INTEGER, INTENT(IN) :: subid              !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<effective porosity volume (mm) (water avaliable for runoff)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<thickness of soil layers (m)
    REAL, INTENT(INOUT) :: percflow(2)        !<percolation (mm/time step)
    REAL, INTENT(INOUT) :: cpercflow(2,numsubstances) !<concentration of percolation (mm/time step)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
	
	!!Variables copied by BTW from calculate_infiltration subroutine
	REAL, INTENT(IN)    :: infilt            !<gross infiltration; rain+snowmelt (mm/timestep)
	REAL, INTENT(IN)    :: cinfilt(numsubstances)       !<concentration of infiltration
	
	!!Variables copied by BTW from calculate_soil_runoff subroutine
	!______________________________________________________________________
	
	!REAL, INTENT(IN)    :: ddepth         !<Depth of stream, drainagedepth (m)
    !TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(OUT)   :: soilrunoff(maxsoillayers)  !<runoff
    REAL, INTENT(OUT)   :: csoilrunoff(numsubstances,maxsoillayers)    !<concentration of runoff
	!________________________________________________________________________
    
	
	!!Variables copied by BTW from calculate_actual_soil_evaporation
	!####################################################################

    !Argument declarations
   
    REAL, INTENT(IN)    :: temp                     !<air temperature
    REAL, INTENT(IN)    :: epot                     !<potential evapotranspiration (mm/timestep)
    
    REAL, INTENT(IN)    :: epotfrac(2)              !<relative distribution of potential evaporation between upper two soil layers (-)
    
    REAL, INTENT(OUT)   :: evap                     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evapflows(2)             !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubstances)     !<concentration in evapotranspiration (?)
    REAL, INTENT(IN)    :: barefrac                 !<fraction of soil that has evapotranspiration (-)

    !Local variables
    INTEGER k   !loop-variable
    INTEGER status  !error status of subroutine
    REAL evap1,evap2  !evapotranspiration of soillayer 1 and 2 (mm/timestep)
    REAL cevap1(numsubstances),cevap2(numsubstances)  !concentration of evapotranspiration
    REAL soiltemp_reduction(maxsoillayers),soiltemp

   
	
	!####################################################################
    !Local variables
    INTEGER isl    !soil layer (1-3)
    
    REAL delta(maxsoillayers)    !groundwater level above drainage level (streamdepth) (m)
    REAL deltah    !groundwater level above drainage level (streamdepth) (m)
    REAL avail(maxsoillayers)    !water available for runoff and percolation (mm)
    REAL runoff1,runoff2,runoff3
    REAL crunoff1(numsubstances),crunoff2(numsubstances),crunoff3(numsubstances)
    
    !Local variables
    
    REAL perc1,perc2          !percolation soillayer 1 and 2 (mm/timestep)
    REAL maxperc1,maxperc2    !maximum percolation from soillayer 1 and maximum percolation to soillayer 2 (mm/timestep)
    REAL cperc(numsubstances) !concentration of percolation
    REAL firstpercflow(2),cfirstpercflow(2,numsubstances)
	
	!!LOcal variables  defined by BTW
	INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)

	real :: recession_perc1, recession_perc2
	real :: recession_evap1, recession_evap2
	REAL :: perc1_day1,perc2_day1         !percolation until sunrise (mm)
	REAL :: perc1_day2,perc2_day2         !percolation from surise to suset (mm)
	REAL :: perc1_day3,perc2_day3        !percolation from surise to suset (mm)
    REAL :: runoff1_day1,runoff2_day1
	REAL :: runoff1_day2,runoff2_day2
	REAL :: runoff1_day3,runoff2_day3
	REAL :: runoff3_day1,runoff3_day2
	REAL :: evap1_day,evap2_day
	REAL :: evap1_lower
	REAL :: evap2_lower
	REAL :: evap1_day0,evap2_day0
	REAL :: avail1_day, avail2_day
	REAL :: evap1_default, evap2_default
    REAL :: extra_perc1, extra_perc2
	REAL :: soil_water_layer1_initial
	REAL :: soil_water_layer2_initial
	REAL :: soil_water_layer3_initial
	REAL :: recession_soilrunoff1
	REAL :: recession_soilrunoff2
	REAL :: totStorage_removed_1half1,totStorage_removed_1half2 
	REAL :: totStorage_removed_2half1,totStorage_removed_2half2
	REAL :: totalfraction1,totalfraction2
	REAL :: Q_integrated,avail_1_2,avail_2_1,avail_2_2
	REAL :: avail_1_3,avail_2_3
	REAL :: dayfrac
	REAL :: dayfrac1,dayfrac2,dayfrac3
    REAL :: coef
	REAL :: k_tot_inv
	REAL(dp) :: sunrise						!!https://en.wikipedia.org/wiki/Sunrise_equation
	REAL(dp) :: sunset
	REAL(dp) :: Mean_solar_noon				!! is an approximation of mean solar time at {\displaystyle n}n expressed as a Julian day with the day fraction
	REAL(dp) :: long_west					!! is the longitude west (west is negative, east is positive) of the observer on the Earth
	REAL(dp) :: Solar_mean_anomaly			!! is the solar mean anomaly
	REAL(dp) :: center_value				!! is the Equation of the center value needed to calculate Ecliptic_longitude
	REAL(dp) :: Ecliptic_longitude
	REAL(dp) :: Solar_transit				!!  is the Julian date for the local true solar transit (or solar noon).
	REAL(dp) :: Declination_of_the_Sun
	REAL(dp) :: hour_angle					!! is the hour angle from the observer's zenith
	REAL(dp) :: latitude_i
	REAL(dp) :: current_Julian_day			!!  number of days since Jan 1st, 2000 12:00.
	REAL(dp) :: YEAR_current
	REAL(dp) :: MONTH_current
	REAL(dp) :: DAY_current
	REAL(dp) :: J_date
	REAL :: time_diff_UT
	REAL :: eq_time
	REAL :: decl 
	REAL :: year_fractional
	REAL :: var1
	REAL :: var2
	REAL :: var3 
	REAL :: var4
	REAL :: storage_sunrise
	REAL :: storage_sunset
	REAL :: storage2_sunrise
	REAL :: storage2_sunset
	REAL :: S1,S2,S3
	Real(dp) :: Arguments_as_real(15)
	Integer :: hh
	REAL ::	crunoff(numsubstances)
	REAL :: correctK,ktot
	REAL :: evapx,Strx
	
	
	
	
	
	
	
	
	Mean_solar_noon = 0.
	long_west = basin(i)%longitude					!!longitude of the basin centroid
	latitude_i = basin(i)%latitude					!!latitude of the basin centroid
	center_value = 0.
	Solar_mean_anomaly = 0.
	Ecliptic_longitude = 0.
	Solar_transit = 0.
	Declination_of_the_Sun = 0.
	hour_angle = 0.
	sunrise = 0.
	sunset = 0.
	current_Julian_day = 0.
	
	YEAR_current = 0.
	MONTH_current = 0.
	DAY_current = 0.
	J_date = 0.
	
	
	
	
	 !>\b Algorithm \n
    !Default values output variables
    evap = 0.
    cevap = 0.
    evapflows = 0.
	
    !>\b Algorithm \n
    !Start calculate percolation
    firstpercflow = percflow
    cfirstpercflow = cpercflow
	
	
	recession_soilrunoff1 = 0. 
	recession_soilrunoff2 = 0.
	totStorage_removed_1half1  = 0.
	totStorage_removed_2half1 = 0.
	totStorage_removed_1half2  = 0.
	totStorage_removed_2half2 = 0.
	totalfraction1 = 0.
	totalfraction2 = 0.
	Q_integrated = 0.
	
	avail_1_2 = 0.
	avail_2_1 = 0.
	avail_2_2 = 0.
	recession_perc1 = 0.
	recession_perc2 = 0.
	recession_evap1 = 0.
	recession_evap2 = 0.
	
	dayfrac1 = 0.
	dayfrac2 = 0.
	dayfrac3 = 0.
	k_tot_inv = 0.
	
	perc1_day1 = 0.
	perc2_day1 = 0.
	perc1_day2 = 0.
	perc2_day2 = 0.
	perc1_day3 = 0.
	perc2_day3 = 0.	
    runoff1_day1 = 0.
	runoff2_day1 = 0.
	runoff1_day2 = 0.
	runoff2_day2 = 0.
	runoff1_day3 = 0.
	runoff2_day3 = 0.
	runoff3_day1 = 0.
	runoff3_day2 = 0.
	avail1_day = 0.
	avail2_day = 0.
	evap1_day = 0.
	evap1_day0 = 0.
	evap2_day = 0.
	evap2_day0 = 0.
	evap1_lower = 0.
	evap2_lower = 0.
	evap1_default = 0.
    evap2_default = 0.
	extra_perc1 = 0.
	extra_perc2 = 0.
	soil_water_layer1_initial = 0.
	soil_water_layer2_initial = 0.
	soil_water_layer3_initial = 0.
        coef  = 0.
        time_diff_UT = 0.
	eq_time = 0.
        decl = 0.
	year_fractional = 0.
	var1 = 0.
	var2 = 0.
	var3 = 0.
	var4 = 0.
	storage_sunrise = 0.
	storage_sunset = 0.
	storage2_sunrise = 0.
	storage2_sunset = 0.
	evap1 = 0.
	evap2 = 0.
	S1 = 0.
	S2 = 0.
	S3 = 0.
	Arguments_as_real = 0.
	hh = 0
	ktot = 0.
	correctK = 0.
	evapx = 0.
	Strx = 0.
	
	

        if(dayno >=68 .AND. dayno <= 306) then
		time_diff_UT = 5
	else
		time_diff_UT = 4
        end if
		
		
	



	!!https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
       
       
	year_fractional = 2.0*22.0/7.0/365.0*(real(dayno)-1.0)
 
	var1 = 0.001868*cos(year_fractional)
	var2 = 0.032077*sin(year_fractional)
	var3 = 0.014615*cos(year_fractional*2.0)
	var4 = 0.040849*sin(year_fractional*2.0)
	
		eq_time = 229.18*(0.000075+var1- var2-var3-var4)


	decl = (0.006918-0.399912*cos(year_fractional) +					&
	 &		0.070257*sin(year_fractional) - 							&
	 &		0.006758*cos(2.0*year_fractional) + 						&
	 &		0.000907*sin(2.0*year_fractional) - 						&
	 &		0.002697*cos(3.0*year_fractional) + 						&
	 &		0.00148*sin(3.0*year_fractional))



	hour_angle = acos(cos(90.833*22.0/7.0/180.0)/						&
	 &					(cos(latitude_i*22.0/7.0/180.0)*cos(decl))-		&
	 &					tan(latitude_i*22.0/7.0/180.0)*tan(decl))

	sunrise = 720.0 - 4.0*(latitude_i+hour_angle*180/(22.0/7.0))-eq_time

	sunrise = (sunrise/60.0+ time_diff_UT)/24.

	sunset = 720.0 - 4.0*(latitude_i-hour_angle*180/(22.0/7.0))-eq_time

	sunset = (sunset/60.0+ time_diff_UT)/24.												!sunset time as fraction of day

         

	
	!!Record the initial soil water content
	soil_water_layer1_initial  = soilstate%water(1,j,i)
	soil_water_layer2_initial  = soilstate%water(2,j,i)
	soil_water_layer3_initial  = soilstate%water(3,j,i)
	
	!dayfrac = 0.5
	
	dayfrac1 = sunrise
	dayfrac2 = sunset - sunrise
	dayfrac3 = 1. - sunset
	
	
	!>Calculate soil temperature reduction
    soiltemp_reduction=1.
    IF(landpar(m_tredA,classdata(j)%luse).GT.0.)THEN
      DO k=1,2
        IF(soilthick(k,j)>0)THEN
          soiltemp = soilstate%temp(k,j,i)
          IF(soiltemp.GT.landpar(m_ttrig,classdata(j)%luse))THEN
            soiltemp_reduction(k) = 1.-EXP(-landpar(m_tredA,classdata(j)%luse) * (soiltemp-landpar(m_ttrig,classdata(j)%luse))**landpar(m_tredB,classdata(j)%luse))
          ELSE
            soiltemp_reduction(k) = 0.
          ENDIF
        ENDIF
      ENDDO
    ENDIF

      
	
	
	!! Make initial estimate of the evap fluxes assuming that the evap fluxes are the only fluxes
	
	!call calculate_actual_soil_evapotranspiration_no_balance(i,j,temp,epot,wp,fc,  &
    !                                  epotfrac,soilstate,evap1,evap2,barefrac)
	!evap1_default = evap1
    !evap2_default = evap2
	
	!Transformed evap1 formulation
	!evap1 = epot*epotfrac(1)*soiltemp_reduction(1)*Max(1.,((soilstate%water(1,j,i)-wp(1))/(basinlp(i) * fc(1))))
	
	
! !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   !Soil storage1
! 
! !___________________________________________________________________________________
! S1 =soilstate%water(1,j,i)
! S2 = soilstate%water(2,j,i)
! S3 = soilstate%water(3,j,i)
! 
! 
! evap1 = epot*epotfrac(1)*soiltemp_reduction(1)*					&
!   &	Min(1.,(Min(0.,S1-wp(1))/(basinlp(i) * fc(1))))
! 
! 
! !Smoothened HYPE percolation equations
! 
! perc1=Min(soilpar(m_perc1,isoil)*((soilstate%water(1,j,i)-fc(1))/	&
!                                     &	ep(1)),((wp(2)+fc(2)+ep(2)) - soilstate%water(2,j,i)))
! 
! 
! soilrunoff(1) = S1*soilrc(1,j,i)
! 
! !Soil storage2
! !___________________________________________________________________________________
! 
! !HYPE rearranged evapotranspiration equations
! 
! 
! evap2 = epot*epotfrac(2)*soiltemp_reduction(2)*					&
!   &	Min(1.,(Min(0.,S2-wp(2))/(basinlp(i) * fc(2))))
! 
! !Smoothened HYPE percolation equations
! 
! perc2=Min(soilpar(m_perc2,isoil)*((soilstate%water(2,j,i)-fc(2))/	&
!                                     &	ep(2)),((wp(3)+fc(3)+ep(3)) - soilstate%water(3,j,i)))
! 
! 
! soilrunoff(2) = soilstate%water(2,j,i)*soilrc(2,j,i)
! 
! 
! !Soil storage3
! !____________________________________________________________________________________
! 
! !HYPE soil erosion equations
! 
! 
! soilrunoff(3) = S3*soilrc(3,j,i)
! 
! !runoff1_day1 = soilrc(1,j,i)*S1
! !runoff2_day1 = soilrc(2,j,i)*S1
! 
! 



!############################################################################
!#######Begin LNR3_simultaneous_formulation_of_the_state_Space_formulation###
!############################################################################

!********************First day fraction**************************************

					!!!!------------------First Layer------------------------
					


					!!!!------------------First Layer------------------------
					
	!Adjust the concentration of soil layer 1 due to infiltration
	
	CALL add_substance_only(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),infilt,cinfilt)
					
		soilstate%water(1,j,i) = soilstate%water(1,j,i) + infilt
					
					
		perc1 = Min(Max(0.,(soilstate%water(1,j,i)-fc(1)-wp(1))),		&
	 &	(Min(soilpar(m_perc1,isoil)*									&
	 &	Max(0.,Min(1.,(infilt + soilstate%water(1,j,i))/				&
	 &						(fc(1)+ wp(1) + ep(1)))),					&
	 &	Max(0.,infilt + soilstate%water(1,j,i)-fc(1)- wp(1)),			&
     &	Max(0.,(wp(2)+fc(2)+ep(2)) - (soilstate%water(2,j,i))))))*1.0
	 
		
		perc1_day2 = perc1
		
		
		CALL remove_substance_only(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),perc1,cperc,status) 
		CALL add_substance_only(numsubstances,soilstate%water(2,j,i),soilstate%conc(:,2,j,i),perc1,cperc)
		
				
		!Update the first layer storage
		
		
		soilstate%water(1,j,i)= soilstate%water(1,j,i)-	 perc1
	 
	 
	 
	 
		
	 
	 evap1 = Min(Max(0.,soilstate%water(1,j,i)-wp(1)),					&
	 &	(MIN(1.,barefrac)*epot*epotfrac(1)*								&
	 &	soiltemp_reduction(1)*											&
	 &	Min(1.,(Max(0.,soilstate%water(1,j,i)-wp(1))/					&
	 &	(basinlp(i) * fc(1))))))*1.0
	 
	 
	 
	 
		
	 
	 
	 
	 	 IF(numsubstances>0) CALL set_evaporation_concentrations(soilstate%conc(:,1,j,i),cevap1)
      
      !>\li Remove evapotranspiration of soillayer 1
      IF(evap1+realzero<soilstate%water(1,j,i))THEN
        CALL remove_substance_only(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),evap1,cevap1,status)
        IF(status.NE.0) CALL error_remove_water(errstring(7),basin(i)%subid,i,j)
      ELSE
        soilstate%water(1,j,i) = 0.
        IF(numsubstances>0) soilstate%conc(:,1,j,i) = 0.    !remove last traces, safe for wp=0
      ENDIF
	  
	  
		!Update the first layer storage
		
		soilstate%water(1,j,i)= soilstate%water(1,j,i)-	evap1
	 
	 
		soilrunoff(1) = (Max(0.,soilstate%water(1,j,i)-					&
	 &		fc(1)-wp(1))*soilrc(1,j,i))*1.0
	 
	    	!????IF(numsubstances>0.) cpercflow(1,:) = (cperc(:)*perc1+cfirstpercflow(1,:)*firstpercflow(1))/percflow(1)
	IF(numsubstances>0.) cpercflow(1,:) = cperc(:)
	
	IF(numsubstances>0) crunoff = soilstate%conc(:,1,j,i)
    CALL remove_substance_only(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),soilrunoff(1),crunoff,status) 

	 
	 !Update the first layer storage
		
		soilstate%water(1,j,i)= soilstate%water(1,j,i)-	soilrunoff(1)	 
	 
	 
		 
		
		!!!!------------------Second Layer------------------------
					
		soilstate%water(2,j,i) = soilstate%water(2,j,i) + perc1
					
		perc2 = Min(Max(0.,(soilstate%water(2,j,i)-fc(2)-wp(2))),		&
	 &	(Min(soilpar(m_perc2,isoil)*									&
	 &	Max(0.,Min(1.,(perc1+soilstate%water(2,j,i))/					&
	 &						(fc(2)+ wp(2) + ep(2)))),					&
	 &	Max(0.,perc1+soilstate%water(2,j,i)-fc(2)- wp(2)),				&
     &	Max(0.,(wp(3)+fc(3)+ep(3)) - soilstate%water(3,j,i)))))*1.0
	 
	 
	 
	 
	 
	 CALL remove_substance_only(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),perc2,cperc,status) 
    IF(status.NE.0) CALL error_remove_water(errstring(19),subid,i,j)
 	CALL add_substance_only(numsubstances,soilstate%water(3,j,i),soilstate%conc(:,3,j,i),perc2,cperc)
    !???IF(numsubstances>0.) cpercflow(2,:) = (cperc(:)*perc2+cfirstpercflow(2,:)*firstpercflow(2))/percflow(2)
    IF(numsubstances>0.) cpercflow(2,:) = cperc(:)
	
		!Update the second layer storage
		
		soilstate%water(2,j,i)= soilstate%water(2,j,i)- perc2
	
	
	 
	 
	 evap2 = Min(Max(0.,soilstate%water(2,j,i)-wp(2)),					&
	 &	(MIN(1.,barefrac)*epot*epotfrac(2)*								&
	 &	soiltemp_reduction(2)*											&
	 &	Min(1.,(Max(0.,soilstate%water(2,j,i)-wp(2))/					&
	 &	(basinlp(i) * fc(2))))))*1.0

		
	 
	 
	 	 IF(numsubstances>0) CALL set_evaporation_concentrations(soilstate%conc(:,2,j,i),cevap2)
        
        !>\li Remove evapotranspiration of soillayer 2
        IF(evap2+realzero<soilstate%water(2,j,i))THEN
          CALL remove_substance_only(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),evap2,cevap2,status)
          IF(status.NE.0) CALL error_remove_water(errstring(8),basin(i)%subid,i,j)
        ELSE
          soilstate%water(2,j,i) = 0.
          IF(numsubstances>0) soilstate%conc(:,2,j,i) = 0.  !remove last traces, safe for wp=0
        ENDIF
		
		
		!Update the second layer storage
		
		soilstate%water(2,j,i)= soilstate%water(2,j,i)- evap2
	 

		soilrunoff(2) = (Max(0.,soilstate%water(2,j,i)					&
	 &			-fc(2)- wp(2))*soilrc(2,j,i))*1.0
	 
	 
	 
	 
		perc2_day2 = perc2
	
		runoff2_day2 = soilrunoff(2)
	 
	 
			 
	 
	!Update the second layer storage
		
		soilstate%water(2,j,i)= soilstate%water(2,j,i)- soilrunoff(2)
		
		IF(numsubstances>0) crunoff = soilstate%conc(:,2,j,i)
    CALL remove_substance_only(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),soilrunoff(2),crunoff,status) 

	 
	 
	 
		 
	 !!!!-----------------Third layer
		
	soilstate%water(3,j,i)= perc2 + soilstate%water(3,j,i)
		
		soilrunoff(3) = (Max(0.,soilstate%water(3,j,i)- fc(3) -			&
	 &			wp(3))*soilrc(3,j,i))*1.0
	 
	 
	 runoff3_day2 = soilrunoff(3)
	 
	 
		
		
	 
	 
	 
	 
		
		
		

	  
	IF(numsubstances>0) crunoff = soilstate%conc(:,3,j,i)
    CALL remove_substance_only(soilstate%water(3,j,i),numsubstances,soilstate%conc(:,3,j,i),soilrunoff(3),crunoff,status) 


	 !Update the third layer storage
		
		soilstate%water(3,j,i)= soilstate%water(3,j,i)-	soilrunoff(3)
	 
	 
		 
	 
!############################################################################
!#######End LNR3_simultaneous_formulation_of_the_state_Space_formulation#####
!############################################################################
! !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	
	
	
	
	  
	  
	  !!!!!!!@%@%@%@%@%@%@%@%@%@%@
	  
	
	!!Calculate aggregated percolation fluxes
	

		
			
	
      percflow(1) = perc1
	  percflow(2) = perc2
	  
	
	  
	  
	   !>Set output variables
      evap = evap1 + evap2
      evapflows(1) = evap1
      evapflows(2) = evap2
      IF(evap>0 .AND. numsubstances>0)THEN
        cevap(:) = (cevap1(:)*evap1 + cevap2(:)*evap2)/evap
      ENDIF
	  
	  
	  

	  
	  !!!Don't forget to aggregate the evap and pass it to the right variable.
	  
	  
	  crunoff1 = 0.
	  IF(numsubstances>0) crunoff1 = soilstate%conc(:,1,j,i)
	  crunoff2 = 0.
	  IF(numsubstances>0) crunoff2 = soilstate%conc(:,2,j,i)
	  crunoff3 = 0.
	  IF(numsubstances>0) crunoff3 = soilstate%conc(:,3,j,i)
		
		
		csoilrunoff(:,1) = crunoff1
		csoilrunoff(:,2) = crunoff2
		csoilrunoff(:,3) = crunoff3

	
	!write(120041,*)percflow(1),percflow(2),i,j													!!Added by BTW
	!write(120031,*)soilrunoff(1),soilrunoff(2),soilrunoff(3),i,j								!!Added by BTW
	!write(120032,*)avail(1),avail(2),avail(3),i,j								!!Added by BTW
	END SUBROUTINE evaporation_percolation_soilrunoff_linear_reservoir

  !>\brief Drainage level runoff: tile or drainage pipe
  !!
  !> \b Reference ModelDescription Chapter Land routines (Soil water - Runoff through drainage pipes)
  !------------------------------------------------------------------
  SUBROUTINE calculate_tile_drainage(i,j,isoil,subid,wp,fc,ep,&
       sdepth,sthick,rrcscorr,soilstate,runoffd,crunoffd,cweights)

    USE MODVAR, ONLY : soilpar,  &
                       numsubstances,   &
                       maxsoillayers,   &
                       tiledepth
    USE HYPEVARIABLES, ONLY : m_trrcs

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of current subbasin
    INTEGER, INTENT(IN) :: j        !<index of current class
    INTEGER, INTENT(IN) :: isoil    !<soil type index
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers) !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers) !<"field capacity" volume (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers) !<effective porosity volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<Lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<Thickness of soil layers (m)
    REAL, INTENT(IN)    :: rrcscorr    !<correction of recession coefficients
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: runoffd                    !<runoff
    REAL, INTENT(OUT)   :: crunoffd(numsubstances)    !<concentration of runoff 
    REAL, INTENT(OUT)   :: cweights(maxsoillayers)    !<weights for calc. drain.conc from layer.conc (zero or one)
    
    !Local variables 
    INTEGER status    !error status of subroutine call
    REAL    deltah    !groundwater level above tile drainage pipe level (m)
    REAL    trc       !coefficient for runoff recession tile flow (fraction per timestep)

    !>\b Algorithm \n
    !>Set default output values
    runoffd = 0.
    crunoffd(:) = 0.
    cweights(:) = 0.

    IF(tiledepth(j)<=0) RETURN   !no tile drainage

    !>Set local parameters
    trc = soilpar(m_trrcs,isoil)*rrcscorr       !Runoff coefficient for tile runoff
    IF(trc==0)RETURN
    IF(trc>1.) trc = 1.

    !> epending on depth of tile drainage pip calculate:
    IF(tiledepth(j)>0 .AND. trc>0.)THEN
      IF(tiledepth(j)<=sdepth(1))THEN       
        !>\li drainage in uppermost soil layer
        deltah = (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1) - (sdepth(1) - tiledepth(j))  !m
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(1) * ep(1)
          crunoffd(:)=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(15),subid,i,j)
          cweights(1) = 1.
        ENDIF
      ELSEIF(tiledepth(j)<=sdepth(2))THEN   
        !>\li drainage in middle soil layer
        deltah = (soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2) - (sdepth(2) - tiledepth(j))
        IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
          deltah = deltah + (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1)
        ENDIF
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(2) * ep(2)
          IF(runoffd > soilstate%water(2,j,i)-wp(2)-fc(2))  runoffd = soilstate%water(2,j,i)-wp(2)-fc(2)
          crunoffd(:)=soilstate%conc(:,2,j,i)
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(16),subid,i,j)
          cweights(2) = 1.
        ENDIF
      ELSE                                  
        !>\li drainage in deepest soil layer
        deltah = (soilstate%water(3,j,i)-wp(3)-fc(3))/ep(3) * sthick(3) - (sdepth(3) - tiledepth(j))
        IF(soilstate%water(3,j,i)-wp(3)-fc(3)-ep(3)>=0.)THEN
          deltah = deltah + (soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2)
          IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
            deltah = deltah + (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1)
          ENDIF
        ENDIF
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(3) * ep(3)
          IF(runoffd > soilstate%water(3,j,i)-wp(3)-fc(3))  runoffd = soilstate%water(3,j,i)-wp(3)-fc(3)
          crunoffd(:)=soilstate%conc(:,3,j,i)
          CALL remove_water(soilstate%water(3,j,i),numsubstances,soilstate%conc(:,3,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(17),subid,i,j)
          cweights(3) = 1.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_tile_drainage

  !>\brief Soil runoff from all soil layers
  !!
  !> \b Reference ModelDescription Chapter Land routines (Soil water - Groundwater runoff)
  !--------------------------------------------------------------------------------
  SUBROUTINE calculate_soil_runoff(i,j,subid,wp,fc,ep,sdepth,sthick,soilstate,soilrunoff,csoilrunoff,isoil)

    USE MODVAR, ONLY : numsubstances,   &
         maxsoillayers,   &
         streamdepth, &
         modeloption, &
         soilpar, &
         p_frozenrunoff
    USE HYPEVARIABLES, ONLY : soilrc,   &
        m_logsatmp, m_bcosby

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of current subbasin
    INTEGER, INTENT(IN) :: j        !<index of current class
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    INTEGER, INTENT(IN) :: isoil    !<index of soil type
    REAL, INTENT(IN)    :: wp(maxsoillayers) !<volume below wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers) !<"field capacity" volume (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers) !<effective porosity volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<Lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<Thickness of soil layers (m)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(OUT)   :: soilrunoff(maxsoillayers)  !<runoff
    REAL, INTENT(OUT)   :: csoilrunoff(numsubstances,maxsoillayers)    !<concentration of runoff
    
    !Local variables
    INTEGER status    !error status of subroutine call
    INTEGER k        !index of soil layer
    REAL    deltah,deltah1,deltah2    !groundwater level above drainage level (streamdepth) (m)
    REAL    runoff1,runoff2,runoff3
    REAL    crunoff1(numsubstances),crunoff2(numsubstances),crunoff3(numsubstances)
    REAL    liqfrac(maxsoillayers), liqmax(maxsoillayers), porosity, numerator, denominator
    
    SELECT CASE(modeloption(p_frozenrunoff))
    CASE(0) ! Current HYPE soil runoff
      liqfrac(:) = 1.
    CASE(1) ! Reduced runoff due frozen soils
      DO k=1,maxsoillayers
        IF(soilstate%temp(k,j,i) >= 0.) THEN
          liqfrac(k) = 1.
        ELSE
          porosity = ep(k)+wp(k)+fc(k)
          numerator = -334000*soilstate%temp(k,j,i)
          denominator = 9.81*(soilstate%temp(k,j,i)+273.16)*(10**soilpar(m_logsatmp,isoil)/100)
          liqmax(k) = porosity*(numerator/denominator)**(-1/soilpar(m_bcosby,isoil))
          liqfrac(k)= MIN( 1., liqmax(k) / soilstate%water(k,j,i) )
        ENDIF
      END DO
    END SELECT

    !Soil runoff: ditch or local stream
    runoff1 = 0.
    crunoff1(:) = 0.
    runoff2 = 0.
    crunoff2(:) = 0.
    runoff3 = 0.
    crunoff3(:) = 0.
    IF(streamdepth(j)>0)THEN
      IF(streamdepth(j)<=sdepth(1))THEN       !Drainage level in uppermost layer
        deltah = liqfrac(1)*(soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1) - (sdepth(1) - streamdepth(j))  !m
        IF(deltah>0.)THEN
          runoff1 = soilrc(1,j,i) * deltah / sthick(1) * ep(1)
          crunoff1(:) = soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(9),subid,i,j)
        ENDIF
      ELSEIF(streamdepth(j)<=sdepth(2))THEN   !Drainage level in middle layer
        deltah = liqfrac(2)*(soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2) - (sdepth(2) - streamdepth(j))
        deltah1 = liqfrac(1)*(soilstate%water(1,j,i)-wp(1)-fc(1))           
        IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN  
          IF(deltah1>0.) deltah = deltah + deltah1/ep(1) * sthick(1)
        ENDIF
        IF(deltah>0.)THEN
          runoff2 = MIN(liqfrac(2)*(soilstate%water(2,j,i)-wp(2)-fc(2)), soilrc(2,j,i) * deltah / sthick(2) * ep(2))
          crunoff2(:) = soilstate%conc(:,2,j,i)
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoff2,crunoff2,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(10),subid,i,j)
        ENDIF
        IF(deltah1>0.)THEN          
          runoff1 = soilrc(1,j,i) * deltah1   
          crunoff1(:)=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(11),subid,i,j)
        ENDIF
      ELSE                                      !Drainage level in deepest layer
        deltah = liqfrac(3)*(soilstate%water(3,j,i)-wp(3)-fc(3))/ep(3) * sthick(3) - (sdepth(3) - streamdepth(j))
        deltah2 = liqfrac(2)*soilstate%water(2,j,i)-wp(2)-fc(2)        
        deltah1 = liqfrac(1)*soilstate%water(1,j,i)-wp(1)-fc(1)         
        IF(soilstate%water(3,j,i)-wp(3)-fc(3)-ep(3)>=0.)THEN
          IF(deltah2>0.) deltah = deltah + deltah2/ep(2) * sthick(2)
          IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
            IF(deltah1>0.) deltah = deltah + deltah1/ep(1) * sthick(1)
          ENDIF
        ENDIF
        IF(deltah>0.)THEN
          runoff3 = MIN(liqfrac(3)*(soilstate%water(3,j,i)-wp(3)-fc(3)), soilrc(3,j,i) * deltah / sthick(3) * ep(3))
          crunoff3(:) = soilstate%conc(:,3,j,i)
          CALL remove_water(soilstate%water(3,j,i),numsubstances,soilstate%conc(:,3,j,i),runoff3,crunoff3,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(12),subid,i,j)
        ENDIF
        IF(deltah1>0.)THEN           
          runoff1 = soilrc(1,j,i) * deltah1 
          crunoff1(:)=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status)
          IF(status.NE.0) CALL error_remove_water(errstring(13),subid,i,j)
        ENDIF
        IF(deltah2>0.)THEN           
          deltah = deltah+sdepth(2)-streamdepth(j)
          IF(deltah>0.)THEN
            runoff2 = MIN(deltah2,liqfrac(2)*(soilrc(2,j,i) * deltah / sthick(2) * ep(2)))
            crunoff2(:)=soilstate%conc(:,2,j,i)
          ELSE  !(perched watertable)
            runoff2 = soilrc(2,j,i) * deltah2
            crunoff2(:)=soilstate%conc(:,2,j,i)
          ENDIF
        ENDIF
        IF(runoff2>0.)THEN
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoff2,crunoff2,status)
          IF(status.NE.0) CALL error_remove_water(errstring(14),subid,i,j)
        ENDIF
      ENDIF
    ENDIF

    soilrunoff(1) = runoff1
    soilrunoff(2) = runoff2
    soilrunoff(3) = runoff3
    csoilrunoff(:,1) = crunoff1(:)
    csoilrunoff(:,2) = crunoff2(:)
    csoilrunoff(:,3) = crunoff3(:)

  END SUBROUTINE calculate_soil_runoff

  !>\brief Calculate and add infiltration to soil
  !>Includes calculation of surface flow and macropore flow due to
  !>limited infiltration capacity. 
  !>
  !>\b Consequences Module modvar variables soilji and csoilji may change.
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - 
  !> Diversion of surface runoff and macropore flow, Infiltration)
  !-----------------------------------------------------------------------------------
  SUBROUTINE infiltration(i,j,isoil,iluse,wp,fc,ginfilt,cginfilt,infilt,  &
       cinfilt,surfaceflow,csurfaceflow,macroflow,cmacroflow,soilstate)

    USE MODVAR, ONLY : numsubstances,   &
         maxsoillayers,   &
         soilpar, &
         i_in
    USE HYPEVARIABLES, ONLY : m_macrate,m_mactrinf,m_mactrsm, &
         m_srrate

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    INTEGER, INTENT(IN) :: iluse              !<index of landuse
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ginfilt            !<gross infiltration; rain+snowmelt (mm/timestep)
    REAL, INTENT(IN)    :: cginfilt(numsubstances)      !<concentration of gross infiltration
    REAL, INTENT(OUT)   :: infilt             !<infiltration (mm/timestep)
    REAL, INTENT(OUT)   :: cinfilt(numsubstances)       !<concentration of infiltration
    REAL, INTENT(OUT)   :: surfaceflow        !<surface runoff due to limited infiltration capacity (mm/timestep)
    REAL, INTENT(OUT)   :: csurfaceflow(numsubstances)  !<concentration of surface flow
    REAL, INTENT(OUT)   :: macroflow          !<macropore flow (mm/timestep)
    REAL, INTENT(OUT)   :: cmacroflow(numsubstances)    !<concentration of macropore flow
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    REAL help,avail
    REAL macrate,inflowtres,srrate   !infiltration excess runoff and macropore flow parameters
    REAL smtresmm

    !>\b Algorithm \n
    !>Set default output values
    infilt       = ginfilt
    cinfilt      = cginfilt
    surfaceflow  = 0.
    csurfaceflow = 0.
    macroflow    = 0.
    cmacroflow   = 0.

    !>If no incoming water; return
    IF(ginfilt == 0.)  RETURN

    !>Set parameter values
    macrate = soilpar(m_macrate,isoil)         !macropore flow fraction (-)
    inflowtres = soilpar(m_mactrinf,isoil)     !threshold for macropore (and surface runoff) flow (mm/timestep)
    smtresmm = (wp(1)+ fc(1))*soilpar(m_mactrsm,isoil)       !soil moisture threshold for macropore (and surface runoff) flow (mm)
    srrate = soilpar(m_srrate,isoil)           !surface runoff flow fraction (-)
    IF(macrate+srrate>1.)THEN
      help = macrate + srrate
      macrate = macrate/help
      srrate = srrate/help
    ENDIF

    !>Calculate surface flow and macropore flow due to limited infiltration capacity 
    avail = ginfilt - inflowtres 
    IF(avail>0. .AND. soilstate%water(1,j,i) > smtresmm) THEN
      macroflow = macrate * avail
      surfaceflow = srrate * avail
      cmacroflow = cginfilt
      csurfaceflow = cginfilt
    ENDIF

    !>Calculate net infiltration
    infilt = ginfilt - macroflow - surfaceflow
    cinfilt = cginfilt

    !>Add infiltration to the upper soillayer soil, including
    !>transfering of IN in infiltration to solid ON in soil
	
	!#####################################
    !IF(infilt>0)THEN
     !  IF(i_in>0) CALL atmdep_in_loss(iluse,infilt,soilstate%fastN(:,j,i),cinfilt(i_in))   !(Atmospheric) IN moved to fastN 
     !  CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),infilt,cinfilt)
    !ENDIF
	
	!#####################################

  END SUBROUTINE infiltration
       
  !>\brief Calculate and add infiltration to frozen soils
  !>Includes calculation of surface flow and macropore flow due to
  !>limited infiltration capacity, and reduced infiltration under frozen condtions.
  !>
  !>\b Consequences Module modvar variables soilji, csoilji and icelensji may change.
  !>
  !>\b Reference 
  !-----------------------------------------------------------------------------------
  SUBROUTINE infiltrationfrozen(i,j,isoil,iluse,wp,fc,ginfilt,cginfilt,infilt,  &
       cinfilt,surfaceflow,csurfaceflow,macroflow,cmacroflow,soilstate, &
       frozenstate,ep,temp,tmin,tmax)

    USE MODVAR, ONLY : numsubstances,   &
         maxsoillayers,   &
         soilpar, &
         i_in, &
         seconds_per_timestep, &
         missing_value, &
         nclass, nsub
    USE HYPEVARIABLES, ONLY : m_macrate,m_mactrinf,m_mactrsm, &
         m_srrate, m_bfroznsoil

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    INTEGER, INTENT(IN) :: iluse              !<index of landuse
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<effective porosity (mm)
    REAL, INTENT(IN)    :: ginfilt            !<gross infiltration; rain+snowmelt (mm/timestep)
    REAL, INTENT(IN)    :: cginfilt(numsubstances)      !<concentration of gross infiltration
    REAL, INTENT(OUT)   :: infilt             !<infiltration (mm/timestep)
    REAL, INTENT(OUT)   :: cinfilt(numsubstances)       !<concentration of infiltration
    REAL, INTENT(OUT)   :: surfaceflow        !<surface runoff due to limited infiltration capacity (mm/timestep)
    REAL, INTENT(OUT)   :: csurfaceflow(numsubstances)  !<concentration of surface flow
    REAL, INTENT(OUT)   :: macroflow          !<macropore flow (mm/timestep)
    REAL, INTENT(OUT)   :: cmacroflow(numsubstances)    !<concentration of macropore flow
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    TYPE(snowicestatetype),INTENT(IN)  :: frozenstate   !<Frozen states
    REAL, INTENT(IN)  :: temp     !<current class temperature (C)
    REAL, INTENT(IN)  :: tmin     !<current daily min temperature (C)
    REAL, INTENT(IN)  :: tmax     !<current daily max temperature (C)
    
    !Local variables
    REAL help,avail
    REAL macrate,inflowtres,srrate   !infiltration excess runoff and macropore flow parameters
    REAL smtresmm
    REAL tmin_check, tmax_check, t0, potinf ! for frozen soils infiltration

    !>\b Algorithm \n
    !>Set default output values
    infilt       = ginfilt
    cinfilt      = cginfilt
    surfaceflow  = 0.
    csurfaceflow = 0.
    macroflow    = 0.
    cmacroflow   = 0.

    !>If no incoming water; return
    IF(ginfilt == 0.)  RETURN

    !>Set parameter values
    macrate = soilpar(m_macrate,isoil)         !macropore flow fraction (-)
    inflowtres = soilpar(m_mactrinf,isoil)     !threshold for macropore (and surface runoff) flow (mm/timestep)
    smtresmm = (wp(1)+ fc(1))*soilpar(m_mactrsm,isoil)       !soil moisture threshold for macropore (and surface runoff) flow (mm)
    srrate = soilpar(m_srrate,isoil)           !surface runoff flow fraction (-)
    IF(macrate+srrate>1.)THEN
      help = macrate + srrate
      macrate = macrate/help
      srrate = srrate/help
    ENDIF

    !>Calculate surface flow and macropore flow due to limited infiltration capacity 
    avail = ginfilt - inflowtres 
    IF(avail>0. .AND. soilstate%water(1,j,i) > smtresmm) THEN
      macroflow = macrate * avail
      surfaceflow = srrate * avail
      cmacroflow = cginfilt
      csurfaceflow = cginfilt
    ENDIF

    !>Calculate net infiltration
    infilt = ginfilt - macroflow - surfaceflow
    cinfilt = cginfilt
    
    !> Zhao & Gray infiltration into frozen soils
    !> coded by M.K. MacDonald (27 October 2015)
    !> icelens = 1: ice lens present which restrict infiltration; 0 no ice lens
    !> Presence of icelens depends on daily maximum and minimum air temperature. 
    !> If tmin.obs & tmax.obs are not input as model forcing, assume they are +- 5C from daily mean 
    !> (based on winter climate normal data for Winnipeg Richardson Int'l A).
    tmin_check = tmin
    tmax_check = tmax
    IF(tmin_check == missing_value) tmin_check = temp - 5.
    IF(tmax_check == missing_value) tmax_check = temp + 5.
    !> Restricted Case: no infiltration (calculated infilt redirected to macroflow & surfaceflow
    IF((tmin_check < -10. .AND. infilt >= 5.) .OR. (tmax_check < 0. .AND. soilstate%icelens(j,i) == 1)) THEN  
      soilstate%icelens(j,i) = 1 !.TRUE.
      surfaceflow = surfaceflow + infilt*srrate /(macrate+srrate)
      macroflow   = macroflow   + infilt*macrate/(macrate+srrate)
      infilt = 0.
    !> Limited Case: no ice lens, but infiltration limited by presence of frozen soils
    ELSEIF(soilstate%temp(1,j,i) <= 0.) THEN
      soilstate%icelens(j,i) = 0 !.FALSE.
      t0 = MAX(1., 0.65*frozenstate%snow(j,i)-5) ! opportunity time [hours]
      potinf = soilpar(m_bfroznsoil,isoil) * (0.99**2.92) * ((1-soilstate%water(1,j,i)/(ep(1)+wp(1)+fc(1)))**1.64) * ((0.-soilstate%temp(1,j,i))/273.15)**(-0.45) * t0**0.44 ! potential infiltration [mm] 
      potinf = potinf/(t0/(seconds_per_timestep/60/60)) ! potential infiltration [mm/timestep]
      potinf = MAX(potinf, 0.)
      IF(potinf < infilt)THEN ! reduce infilt to restricted potinf and redirected water to macroflow & surfaceflow
       surfaceflow = surfaceflow + (infilt - potinf)*srrate  /(macrate+srrate)
       macroflow   = macroflow   + (infilt - potinf)*macrate /(macrate+srrate)
       infilt = potinf
      ENDIF
    ! Unlimited case: unfrozen soils (regular HYPE infiltration)
    ELSE
      soilstate%icelens(j,i) = 0 !.FALSE.
      infilt = infilt
    ENDIF

    !>Add infiltration to the upper soillayer soil, including
    !>transfering of IN in infiltration to solid ON in soil
    IF(infilt>0)THEN
       IF(i_in>0) CALL atmdep_in_loss(iluse,infilt,soilstate%fastN(:,j,i),cinfilt(i_in))   !(Atmospheric) IN moved to fastN 
       CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),infilt,cinfilt)
    ENDIF

  END SUBROUTINE infiltrationfrozen     

  !>Calculate percolation down through the soil layers 
  !>Includes change in concentration of percolating water. 
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - Percolation)
  !----------------------------------------------------------------------
  SUBROUTINE percolation(i,j,isoil,subid,wp,fc,ep,sthick,percflow,cpercflow,soilstate)

    USE MODVAR, ONLY : numsubstances,   &
                       maxsoillayers,   &
                       genpar,soilpar,landpar, &            
                       classdata
    USE HYPEVARIABLES, ONLY : m_perc1,m_perc2,    &
                              m_crate5,   &
                              m_onpercred, m_pppercred

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    INTEGER, INTENT(IN) :: subid              !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<"effective porosity" volume (mm) (water avaliable for runoff)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<thickness of soil layers (m)
    REAL, INTENT(OUT)   :: percflow(2)        !<percolation (mm/time step)
    REAL, INTENT(OUT)   :: cpercflow(2,numsubstances)        !<concentration of percolation (mm/time step)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    INTEGER status            !error status of subroutine
    REAL perc1,perc2          !percolation soillayer 1 and 2 (mm/timestep)
    REAL maxperc1,maxperc2    !maximum percolation from soillayer 1 and maximum percolation to soillayer 2 (mm/timestep)
    REAL cperc(numsubstances) !concentration of percolation

    !>\b Algorithm \n
    !Start calculate percolation
    percflow = 0.
    IF(numsubstances>0.) cpercflow = 0.
    IF(sthick(2)>0)THEN

      !>Calculate limitations for percolation              
      IF(soilstate%water(1,j,i)-wp(1)>fc(1))THEN
        maxperc1 = MIN((soilstate%water(1,j,i)-wp(1)-fc(1)),soilpar(m_perc1,isoil))      !Maximum percolation from uppermost soil layer
      ELSE
        maxperc1 = 0.
      ENDIF
      maxperc2 = MIN(wp(3)+fc(3)+ep(3)-soilstate%water(3,j,i),soilpar(m_perc2,isoil))   !Maximum percolation deep soil layer can accept

      !>Calculate percolation amount of water
      IF(soilstate%water(2,j,i)+maxperc1<=wp(2)+fc(2))THEN
        perc1=maxperc1
        perc2=0.
      ELSEIF(soilstate%water(2,j,i)+maxperc1<=wp(2)+fc(2)+ep(2))THEN
        perc1=maxperc1
        perc2=MIN(soilstate%water(2,j,i)+perc1-wp(2)-fc(2),maxperc2)
      ELSE
        perc2=MIN(soilstate%water(2,j,i)+maxperc1-wp(2)-fc(2),maxperc2)
        IF(maxperc1+soilstate%water(2,j,i)-perc2<=wp(2)+fc(2)+ep(2))THEN
          perc1=maxperc1
        ELSE  !perc2=maxperc2
          perc1=wp(2)+fc(2)+ep(2)-soilstate%water(2,j,i)+perc2
        ENDIF
      ENDIF

      !>Move percolation water to underlaying soillayer and reduce the concentrations:
      IF(perc1>0)THEN
        cperc(:)=soilstate%conc(:,1,j,i)
        !>\li Reduce OC, ON and PP concentration of water percolating from upper soillayer
        CALL doc_percolation_reduction(numsubstances,cperc,genpar(m_crate5),soilstate%temp(1,j,i),   &
             soilstate%water(1,j,i),wp(1),wp(1)+fc(1)+ep(1),sthick(1))
        CALL onpp_percolation_reduction(numsubstances,cperc,landpar(m_onpercred,classdata(j)%luse),   &
             landpar(m_pppercred,classdata(j)%luse))  
        !>\li Remove water from upper soillayer and add to second soillayer
        CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),perc1,cperc,status) 
        IF(status.NE.0) CALL error_remove_water(errstring(18),subid,i,j)
        CALL add_water(numsubstances,soilstate%water(2,j,i),soilstate%conc(:,2,j,i),perc1,cperc)
        percflow(1) = perc1
        IF(numsubstances>0.) cpercflow(1,:) = cperc(:)
      ENDIF
      IF(perc2>0)THEN
        cperc(:)=soilstate%conc(:,2,j,i)
        !>\li Reduce OC, ON and PP concentration of water percolating from middle soillayer
        CALL doc_percolation_reduction(numsubstances,cperc,genpar(m_crate5),soilstate%temp(2,j,i),   &
             soilstate%water(2,j,i),wp(2),wp(2)+fc(2)+ep(2),sthick(2))
        CALL onpp_percolation_reduction(numsubstances,cperc,landpar(m_onpercred,classdata(j)%luse),   &
             landpar(m_pppercred,classdata(j)%luse))
        !>\li Remove water from middle soillayer and add to third soillayer
        CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),perc2,cperc,status) 
        IF(status.NE.0) CALL error_remove_water(errstring(19),subid,i,j)
        CALL add_water(numsubstances,soilstate%water(3,j,i),soilstate%conc(:,3,j,i),perc2,cperc)
        percflow(2) = perc2
        IF(numsubstances>0.) cpercflow(2,:) = cperc(:)
      ENDIF

    ENDIF

  END SUBROUTINE percolation

  !>\brief Add macropore water flow to soil layer with groundwater level. 
  !!If this soillayer can't take all water the rest is added to the
  !!soillayer(s) above. 
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - Macropore flow)
  !---------------------------------------------------------------------
  SUBROUTINE add_macropore_flow(i,j,macroflow,cmacroflow,wp,fc,ep,pw,sdepth,sthick,slmacroflows,soilstate)

    USE MODVAR, ONLY : numsubstances,   &
         maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    REAL, INTENT(IN)    :: macroflow          !<macropore flow to be added (mm/timestep)
    REAL, INTENT(IN)    :: cmacroflow(numsubstances) !<concentration of macropore flow 
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<"effective porosity" volume (mm) (water avaliable for runoff)
    REAL, INTENT(IN)    :: pw(maxsoillayers)  !<total pore volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<thickness of soil layers (m)
    REAL, INTENT(OUT)   :: slmacroflows(maxsoillayers) !<macropore flow to each soil layer (mm/timestep)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    REAL gwat                   !groundwater table (m) (negative)
    REAL newsoil                
    REAL fill,fill2

    !Start subroutine
    !>\b Algorithm \n
    !>If no macropore flow: return
    slmacroflows = 0
    IF(macroflow<=0) RETURN   

    !>Find soillayer of groundwater table
    CALL calculate_groundwater_table(soilstate%water(:,j,i),wp,fc,ep,pw,sdepth(:),sthick(:),gwat)
    !>If groundwater table in soillayer three:
    IF(-gwat>sdepth(2)) THEN
      !>\li Check if soillayer three has room and add the water there is room for
      newsoil = soilstate%water(3,j,i) + macroflow
      IF(newsoil > pw(3)) THEN
        fill = newsoil - pw(3)
        soilstate%conc(:,3,j,i) = (soilstate%water(3,j,i)*soilstate%conc(:,3,j,i) + (macroflow - fill)*cmacroflow(:))/pw(3)
        soilstate%water(3,j,i) = pw(3)
        slmacroflows(3) = macroflow - fill
      ELSE
        fill = 0.
        soilstate%conc(:,3,j,i) = (soilstate%water(3,j,i)*soilstate%conc(:,3,j,i) + macroflow*cmacroflow(:))/newsoil
        soilstate%water(3,j,i) = newsoil
        slmacroflows(3) = macroflow
      ENDIF
      !>\li If too much water, check if soillayer 2 has room and add the water there is room for
      IF(fill > 0.) THEN
        newsoil = soilstate%water(2,j,i) + fill
        IF(newsoil > pw(2)) THEN
          fill2 = newsoil - pw(2)
          soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + (fill-fill2)*cmacroflow(:))/pw(2)
          soilstate%water(2,j,i) = pw(2)
          slmacroflows(2) = fill - fill2
        ELSE
          fill2 = 0.
          soilstate%conc(:,2,j,i) = ((newsoil-fill)*soilstate%conc(:,2,j,i) + fill*cmacroflow(:))/newsoil
          soilstate%water(2,j,i) = newsoil
          slmacroflows(2) = fill
        ENDIF
        !>\li If still too much water add the rest to soillayer 1
        IF(fill2 > 0.) THEN
          newsoil = soilstate%water(1,j,i) + fill2
          soilstate%conc(:,1,j,i) = (soilstate%water(1,j,i)*soilstate%conc(:,1,j,i) + fill2*cmacroflow(:))/newsoil
          soilstate%water(1,j,i) = newsoil
          slmacroflows(1) = fill2
        ENDIF
      ENDIF

    !>Elseif groundwater table in soillayer two:
    ELSEIF(-gwat>sdepth(1)) THEN
      newsoil = soilstate%water(2,j,i) + macroflow
      !>\li Check if soillayer 2 has room and add the water there is room for
      IF(newsoil > pw(2)) THEN
        fill = newsoil - pw(2)
        soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + (macroflow - fill)*cmacroflow(:))/pw(2)
        soilstate%water(2,j,i) = pw(2)
        slmacroflows(2) = macroflow - fill
      ELSE
        fill = 0.
        soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + macroflow*cmacroflow(:))/newsoil
        soilstate%water(2,j,i) = newsoil
        slmacroflows(2) = macroflow
      ENDIF
      !>\li If too much water add the rest to soillayer 1
      IF(fill > 0.) THEN
        CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),fill,cmacroflow)
        slmacroflows(1) = fill
      ENDIF

    !>Elseif groundwater table in soillayer one:
    ELSE
      !>\li Add macropore flow to soillayer 1
      CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),macroflow,cmacroflow)
      slmacroflows(1) = macroflow
    ENDIF

  END SUBROUTINE add_macropore_flow

  !>Subroutine for calculation of ground water table level (metres
  !>above land surface)
  !!
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_groundwater_table(soil,wpvol,fcvol,&
       epvol,totvol,soildep,thickness,gwat)
    USE MODVAR, ONLY : maxsoillayers

    !Argument declarations
    REAL, INTENT(IN)     :: soil(maxsoillayers)     !<soil moisture (mm)
    REAL, INTENT(IN)     :: wpvol(maxsoillayers)    !<wilting point volume in all layers (mm) 
    REAL, INTENT(IN)     :: fcvol(maxsoillayers)    !<field capacity volume in all layers (mm)
    REAL, INTENT(IN)     :: epvol(maxsoillayers)    !<effective porosity volume in all layers (mm) 
    REAL, INTENT(IN)     :: totvol(maxsoillayers)   !<maximum volume of water in soil layer (pore volume) (mm) 
    REAL, INTENT(IN)     :: soildep(maxsoillayers)  !<depth of soil layers (m)
    REAL, INTENT(IN)     :: thickness(maxsoillayers)  !<thickness of soil layers (m)
    REAL, INTENT(OUT)    :: gwat !<ground water table (m), measured from land surface and up
    
    !Local parameters
    REAL, PARAMETER :: mindiff = 0.000005 !safeguard in choosing groundwater table from soil layer tables
    
    !Local variables
    REAL gwat1,gwat2,gwat3      !groundwater table of each soil layer

    gwat = 0.

    IF(thickness(2)>0)THEN
      IF(soil(1)-wpvol(1)-fcvol(1)>0.0)THEN
        gwat1 = (soil(1) - totvol(1))/epvol(1) * thickness(1)     !negative (m)
        IF(gwat1 > 0) gwat1 = (soil(1) - totvol(1)) * 0.001     !100% porositet above land surface
      ELSE
        gwat1 = -soildep(1)
      ENDIF
      IF(soil(2)-wpvol(2)-fcvol(2)>0.0)THEN
        gwat2 = (soil(2)-totvol(2))/epvol(2) * thickness(2) - soildep(1)     !negative (m)
      ELSE
        gwat2 = -soildep(2)
      ENDIF
      IF(soil(3)-wpvol(3)-fcvol(3)>0.0)THEN
        gwat3 = (soil(3)-totvol(3))/epvol(3) * thickness(3) - soildep(2)     !negative (m)
      ELSE
        gwat3 = -soildep(3)
      ENDIF
      IF(-gwat3>soildep(2)+mindiff)THEN     !Find ground water table as lowest level with filled pores below
        gwat = gwat3
      ELSEIF(-gwat2>soildep(1)+mindiff)THEN
        gwat = gwat2
      ELSE
        gwat = gwat1
      ENDIF
    ELSE
      IF(soil(1)-wpvol(1)-fcvol(1)>0.0)THEN
        gwat1 = (soil(1) - totvol(1))/epvol(1) * thickness(1)     !negative (m)
        IF(gwat1 > 0) gwat1 = (soil(1) - totvol(1)) * 0.001     !100% porositet above land surface
      ELSE
        gwat1 = -soildep(1)
      ENDIF
      gwat = gwat1
    ENDIF

  END SUBROUTINE calculate_groundwater_table

  !>Calculation of snowdepth and age of snow. Snow depth depends on age of snow.
  !>
  !>\b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_snowdepth(snow,oldsnow,snowdens0,snowdensdt,snowdepth,snowage)

    USE MODVAR, ONLY : timesteps_per_day
    
    !Argument declarations
    REAL, INTENT(IN)    :: snow       !<snow water equivalent (mm)
    REAL, INTENT(IN)    :: oldsnow    !<snow water equivalent before addition of snowfall/melt this timestep (mm)
    REAL, INTENT(IN)    :: snowdens0  !<model parameter, snow density at snowfall
    REAL, INTENT(IN)    :: snowdensdt !<model parameter, snow density increase due to ageing (g/cm3.day)
    REAL, INTENT(OUT)   :: snowdepth  !<depth of snow (cm)
    REAL, INTENT(INOUT) :: snowage    !<age of snow pack (timesteps)

    !Local variables
    REAL snowfall    !snowfall today
    REAL snowdens    !Snow density (g/cm3)

    !> \b Algorithm \n
    !>Calculate snowfall
    snowfall = snow - oldsnow     !if>0
    snowdepth  = 0.
    !>If no snow: reset snow variables
    IF(snow==0)THEN
      snowage = 0.
    ELSE
    !>Else: calculate snow age, density and depth
      snowage = snowage + 1.
      IF(snowfall > 0) snowage = snowage * oldsnow / snow
      snowdens = snowdens0 + snowdensdt * snowage / timesteps_per_day
      IF(snowdens>0) snowdepth  = 0.1 * snow / snowdens     !0.1 cm/mm
    ENDIF

  END SUBROUTINE calculate_snowdepth

  !>Calculation of soil temperature in soil layers and deep soil
  !>
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_soiltemp(n,airtemp,snowdepth,soilmemdeep,soilmemlayer,deeptemp,soiltemp)

    USE MODVAR, ONLY : timesteps_per_day
    
    !Argument declarations
    INTEGER, INTENT(IN)  :: n               !<number of soil layers
    REAL, INTENT(IN)     :: airtemp         !<air temperature (degree Celcius) 
    REAL, INTENT(IN)     :: snowdepth       !<snow depth (cm)
    REAL, INTENT(IN)     :: soilmemdeep     !<parameter, temperature memory of deep soil (days)
    REAL, INTENT(IN)     :: soilmemlayer(n) !<parameter, temperature memory of soil layer (timesteps)
    REAL, INTENT(INOUT)  :: deeptemp        !<deep soil temperature (degree Celcius)
    REAL, INTENT(INOUT)  :: soiltemp(n)     !<soil temperature (degree Celcius)
    
    !Local parameters
    REAL, PARAMETER :: spfrost = 10.      !coefficient of equation for weight of air temperature in soil temperature calculation (days/cm)
    REAL, PARAMETER :: weightdeep = 0.001 !weight of deep soil temperature for soil temperature calculation (dimensionless)
    
    !Local variables
    INTEGER k        !layer index
    REAL weightair   !weight of air temperature for soil temperature calculation (dimensionless)

    !> \b Algorithm \n
    !>Calculate deep soil temperature
    weightair = 1./ ((soilmemdeep + spfrost * snowdepth)*timesteps_per_day)
    CALL calculate_weighted_temperature(airtemp,weightair,0.0,0.0,deeptemp)
    !>Calculate soil layer temperature for each soil layer
    DO k = 1,n
      weightair = 1./ (soilmemlayer(k) + spfrost * snowdepth*timesteps_per_day)
      CALL calculate_weighted_temperature(airtemp,weightair,deeptemp,weightdeep,soiltemp(k))
    ENDDO

   END SUBROUTINE calculate_soiltemp

  !>Calculation of soil temperature as an average of
  !>three temperatures: air temperature, deep soil temperature and
  !>previous soil temperature
  !>
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !--------------------------------------------------------------------------------
  SUBROUTINE calculate_weighted_temperature(temp1,weight1,temp2,weight2,soiltemp)

    !Argument declarations
    REAL, INTENT(IN)     :: temp1     !<air temperature (degree Celcius) 
    REAL, INTENT(IN)     :: weight1   !<weight of temp1 (dimensionless)
    REAL, INTENT(IN)     :: temp2     !<temperature of deep soil (ca 1m) (degree Celcius)
    REAL, INTENT(IN)     :: weight2   !<weight of temp2 (dimensionless)
    REAL, INTENT(INOUT)  :: soiltemp  !<soil layer temperature (degree Celcius)

    !> \b Algorithm \n
    !>Calculate weighted temperature
    soiltemp = soiltemp * (1. - weight1 - weight2) + &
         temp1    * weight1  +   &
         temp2    * weight2

  END SUBROUTINE calculate_weighted_temperature

  !>Calculation of soil frost depth depending on temperature of soil
  !>
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !-------------------------------------------------------------------------------------
  SUBROUTINE calculate_frostdepth(fc,cfrost,sfrost,soil,frostdepth,soiltemp,thickness)
  
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    REAL, INTENT(IN)  :: fc           !<water content at field capacity (mm)
    REAL, INTENT(IN)  :: cfrost       !<soil frost coefficient, land use dependent (cm/degree)
    REAL, INTENT(IN)  :: sfrost       !<soil frost coefficient, soil type dependent (cm/degree)
    REAL, INTENT(IN)  :: soil         !<soil water (mm) 
    REAL, INTENT(OUT) :: frostdepth   !<depth of soil frost (cm)
    REAL, INTENT(IN)  :: soiltemp(:)  !<soil temperature (degree Celcius)
    REAL, INTENT(IN)  :: thickness(:) !<soil layer thickness (m)
    
    !Local variables 
    INTEGER dim
    REAL uppertemp  !average temperature upper soil (ca 50 cm)

    !> \b Algorithm \n
    !>If soil frost parameters are set:
    IF(cfrost>0 .AND. sfrost>0)THEN
      dim = SIZE(soiltemp)
      !> \li Calculate average temperature of upper two soil layers
      IF(dim==1)THEN
        uppertemp = soiltemp(1)
      ELSE
        uppertemp = (soiltemp(1)*thickness(1)+soiltemp(2)*thickness(2))/(thickness(1)+thickness(2))
      ENDIF
      !> \li If temperature is negative, calculate soil frost depth
      IF(uppertemp<0)THEN
        frostdepth = cfrost * sfrost * uppertemp * fc / soil
      ELSE
        frostdepth = 0.
      ENDIF
    ELSE
    !Else soil frost is set to missing
      frostdepth = missing_value
    ENDIF

  END SUBROUTINE calculate_frostdepth

  !> \brief Calculation of soil moisture deficit (mm left to field capacity)
  !>in top two layers
  !>
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !---------------------------------------------------------------------
  SUBROUTINE calculate_soil_moisture_deficit(soil,wpvol,fcvol,thickness,smdef)
    USE MODVAR, ONLY : maxsoillayers

    !Argument declaration
    REAL, INTENT(IN)     :: soil(maxsoillayers)       !<soil moisture  (mm)
    REAL, INTENT(IN)     :: wpvol(maxsoillayers)      !<wilting point volume in all layers (mm) 
    REAL, INTENT(IN)     :: fcvol(maxsoillayers)      !<"field capacity" volume in all layers (mm)
    REAL, INTENT(IN)     :: thickness(maxsoillayers)  !<thickness of soil layers (m)
    REAL, INTENT(OUT)    :: smdef                     !<soil moisture deficit (mm)
    
    !Local variables
    REAL smdef1,smdef2      !soil moisture deficit of each soil layer

    !> \b Algorithm \n
    !>Initate soil moisture deficit to zero
    smdef = 0.

    !>Calculate soil moisture deficit in top soil layer, add to total
    smdef1 = fcvol(1)+wpvol(1)-soil(1)
    IF(smdef1>0.) smdef = smdef + smdef1

    !>Calculate soil moisture deficit in second soil layer, add to total
    IF(thickness(2)>0)THEN
      smdef2 = fcvol(2)+wpvol(2)-soil(2)
      IF(smdef2>0.) smdef = smdef + smdef2
    ENDIF

  END SUBROUTINE calculate_soil_moisture_deficit

END MODULE SOIL_PROCESSES
