!> \file sw_proc.f90
!> Contains module surfacewater_processes.

!>Lake and river water related subroutines in HYPE
MODULE SURFACEWATER_PROCESSES

  !Copyright 2012-2016 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !------------------------------------------------------------------------
  USE STATETYPE_MODULE
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water,       &
                                          error_remove_water, &
                                          add_water
  USE SOIL_PROCESSES, ONLY : calculate_snowdepth,  &
                             snowalbedo_function,   &
                             latentheat_tempfunction
       
  !Also uses hypevariables, modvar,general_functions
  IMPLICIT NONE
  PRIVATE
  !--------------------------------------
  !Private procedures 
  !--------------------------------------
  ! set_rivertemp
  ! calc_qbank 
  ! update_qbank
  ! get_current_lake_outflow_parameters
  ! average_flow_rating_curve
  ! riverice_riverwater_interaction
  ! calculate_snow_on_ice
  ! calculate_lakeice_lakewater_interaction
  ! calculate_icedepth
  ! calculate_T2_transfer
  ! calculate_T2_transfer_upper2lower
  !-------------------------------------
  PUBLIC :: calculate_landarea_riverlength,  &
            add_precipitation_to_river, &
            calculate_river_evaporation, &
            calculate_actual_lake_evaporation, &
            sum_upstream_area, &
            set_general_rating_k,  &
            calculate_water_temperature, &
            set_water_temperature,  &
            calculate_river_characteristics, &
            translation_in_river, &
            point_abstraction_from_main_river, &
            point_abstraction_from_outlet_lake, &
            calculate_outflow_from_lake, &
            calculate_flow_from_lake_waterstage, &
            check_outflow_from_lake, &
            remove_outflow_from_lake, &
            calculate_flow_within_lake, &
            calculate_olake_waterstage, &
            calculate_regamp_adjusted_waterstage, &
            calculate_branched_flow,  &
            calculate_lake_volume, &
            T2_processes_in_river, &
            T2_processes_in_lake, &
            ice_processes_in_lake, &
            ice_processes_in_river, &
            add_T2_concentration_in_precipitation_on_water, &
            get_rivertempvol, &
            inline_dam_outflow

  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(8)  !error message for location of remove_water call
  PARAMETER (errstring = (/'evapotranspiration lake, less than lake volume',    &   !1
                           'evapotranspiration lake, more than lake volume',    &   !2
                           'evapotranspiration lake, slowlake part used   ',    &   !3
                           'lake outflow, no NPC simulation               ',    &   !4 
                           'lake outflow, no division in parts (NPC sim)  ',    &   !5
                           'lake outflow, from fastlake part              ',    &   !6
                           'lake outflow, from slowlake part              ',    &   !7
                           'flow between fast- and slowlake parts         ' /))     !8

CONTAINS

  !>\brief Calculate land area of subbasins and determine riverlength 
  !!for local streams and main rivers.
  !>
  !\b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !----------------------------------------------------------------------------
  SUBROUTINE calculate_landarea_riverlength(nsub,landarea,rivlength)

    USE MODVAR, ONLY : basin,classbasin, &
                       slc_ilake,slc_olake, &
                       slc_lriver,slc_mriver
    
    !Argument declarations
    INTEGER, INTENT(IN) :: nsub   !<Number of subbasins
    REAL, INTENT(OUT)   :: landarea(nsub)       !<land area [m2]
    REAL, INTENT(OUT)   :: rivlength(2,nsub)    !<river length [m]
   
    !Local variables
    INTEGER i

    !>\b Algorithm \n
    !>Calculate land area of subbasin
    landarea = basin%area
    IF(slc_ilake>0)  landarea = landarea - basin(:)%area * classbasin(:,slc_ilake)%part
    IF(slc_olake>0)  landarea = landarea - basin%area * classbasin(:,slc_olake)%part
    IF(slc_lriver>0) landarea = landarea - basin%area * classbasin(:,slc_lriver)%part
    IF(slc_mriver>0) landarea = landarea - basin%area * classbasin(:,slc_mriver)%part
    !>Calculate square root of landarea, set local river length
    DO i = 1,nsub
      IF(landarea(i)<0) landarea(i)=0.   !Safe for all lake subbasin (1-1<0)
      rivlength(1,i) = SQRT(landarea(i))
    ENDDO
    !>Set main river length from GeoData, or if zero use same as local river length
    rivlength(2,:) = basin(:)%rivlen  !main river length
    WHERE(rivlength(2,:)==0) rivlength(2,:) = rivlength(1,:)

  END SUBROUTINE calculate_landarea_riverlength

  !>\brief Add precipitation to river, according to volume of watercourse elements.
  !>
  !\b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !----------------------------------------------------------------------------
  SUBROUTINE add_precipitation_to_river(i,pooltype,area,prec,cprec,riverstate)

    USE MODVAR, ONLY : numsubstances
    USE HYPEVARIABLES, ONLY : ttpart,ttstep

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<index of subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<rivertype: 1=lriver, 2=mriver
    REAL, INTENT(IN)    :: area                       !<river area (m2)
    REAL, INTENT(IN)    :: prec                       !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: cprec(numsubstances)       !<concentration of precipitation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state

    !Local variables
    INTEGER l
    REAL precm3
    REAL totvol
    REAL waterfrac

    !>\b Algorithm \n
    precm3 = prec * 1.E-3 * area
    !>Calculate total volume of river to use fractions of river water i different compartment
    totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))

    !>Add precipitation to river watercourse for each compartment in relation to its volume fraction
    IF(totvol>0)THEN
      IF(riverstate%water(pooltype,i)>0)THEN
        waterfrac = riverstate%water(pooltype,i)/totvol
        CALL add_water(numsubstances,riverstate%water(pooltype,i),riverstate%conc(:,pooltype,i),waterfrac*precm3,cprec)
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
          CALL add_water(numsubstances,riverstate%qqueue(l,pooltype,i),riverstate%cqueue(:,l,pooltype,i),waterfrac*precm3,cprec)
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that remaining outflow will be correct
          CALL add_water(numsubstances,riverstate%qqueue(l,pooltype,i),riverstate%cqueue(:,l,pooltype,i),waterfrac*precm3,cprec)
        ENDIF
      ENDIF
    ELSE
      !>If no river volume add all precipitation to river water compartment
      riverstate%water(pooltype,i) = precm3
    ENDIF  

  END SUBROUTINE add_precipitation_to_river

  !>\brief Calculate and remove evaporation from river
  !>
  !\b Reference ModelDescription Chapters Rivers and lakes (Rivers - Common river processes)
  !> and Processes above ground (Evaporation)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_river_evaporation(i,j,pooltype,numsubst,area,temp,epot,evap,cevap,riverstate)

    USE MODVAR, ONLY : basin,classdata, &
                       landpar,   &
                       cwater,    &
                       i_t1,i_t2, &
                       i_sm,i_gm,i_rn,i_li,i_ri,i_si
    USE HYPEVARIABLES, ONLY : m_ttmp,ttpart,ttstep

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: j        !<class index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    INTEGER, INTENT(IN) :: numsubst !<number of substances modelled
    REAL, INTENT(IN)    :: area     !<river area (m2)
    REAL, INTENT(IN)    :: temp     !<air temperature
    REAL, INTENT(IN)    :: epot     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubst) !<concentration in evapotranspiration (eg. mg/L)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<Lake state

    !Local variables
    INTEGER k,l     !loop-variable, substance/queue
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for evaporation (C)
    REAL evapm3     !actual evaporation in m3
    REAL totvol     !total river watercourse volume (m3)
    REAL waterfrac  !fraction of water to be removed
    
    !>\b Algorithm \n
    !>Set default values output variables (zero evaporation)
    evap = 0.
    cevap = 0.

    !Set local parameter
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    IF(temp > tt) THEN
      !>If temperature is above threshold river evaporation is potential
      evapm3 = epot*area*1.E-3
      !>Set concentration of evaporation (T1=O18)
      DO k=1,numsubst
        cevap(k) = 0.
        IF(k==i_t1) cevap(k) = riverstate%conc(k,pooltype,i)  !t1 == O18 in hype (approx. should be different for different waters!)
        !For T2==temperature, set the evaporation concentration = L/c = latent heat of vaporization divided by heat capacity
        IF(k==i_t2) cevap(k) = 1000. * latentheat_tempfunction(riverstate%conc(k,pooltype,i)) / cwater ! 1000 * [MJ/kg] / [KJ/kg/C] = C
       !For water origin model, set evaporation concentration equal riverstate conc
        IF(k==i_sm) cevap(k) = riverstate%conc(k,pooltype,i)  !sm == snowmelt, in water origin trace element model
        IF(k==i_gm) cevap(k) = riverstate%conc(k,pooltype,i)  !gm == glaciermelt, in water origin trace element model
        IF(k==i_rn) cevap(k) = riverstate%conc(k,pooltype,i)  !rn == rainfall, in water origin trace element model
        IF(k==i_li) cevap(k) = riverstate%conc(k,pooltype,i)  !li == initial lake, in water origin trace element model
        IF(k==i_ri) cevap(k) = riverstate%conc(k,pooltype,i)  !ri == initial river, in water origin trace element model
        IF(k==i_si) cevap(k) = riverstate%conc(k,pooltype,i)  !si == initial soil, in water origin trace element model
      ENDDO

      !>Calculate total river volume to use fractions of river water i different compartments
      totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
      IF(totvol<=0.) RETURN

      !>Check if enough water is available for evaporation in each compartment
      IF(evapm3<totvol)THEN
        !>Remove evaporation from river watercourse compartments
        waterfrac = riverstate%water(pooltype,i)/totvol
        CALL remove_water(riverstate%water(pooltype,i),numsubst,riverstate%conc(:,pooltype,i),waterfrac*evapm3,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
        DO l = 1,ttstep(pooltype,i)
          IF(riverstate%qqueue(l,pooltype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
            CALL remove_water(riverstate%qqueue(l,pooltype,i),numsubst,riverstate%cqueue(:,l,pooltype,i),waterfrac*evapm3,cevap,status)
            IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
          ENDIF
        ENDDO
        IF(ttpart(pooltype,i)>0)THEN
          l = ttstep(pooltype,i) + 1
          IF(riverstate%qqueue(l,pooltype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that pool get correct concentration change
            CALL remove_water(riverstate%qqueue(l,pooltype,i),numsubst,riverstate%cqueue(:,l,pooltype,i),waterfrac*evapm3,cevap,status)
            IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
          ENDIF
        ENDIF
        evap = epot
      ELSE
        !>If less water than wanted, remove last traces of substance with the evaporation
        evapm3 = totvol
        riverstate%water(pooltype,i) = 0.
        IF(numsubst>0.) riverstate%conc(:,pooltype,i) = 0.
        DO l = 1,ttstep(pooltype,i)
          riverstate%qqueue(l,pooltype,i) = 0.
          IF(numsubst>0.) riverstate%cqueue(:,l,pooltype,i) = 0.
        ENDDO
        IF(ttpart(pooltype,i)>0)THEN
          l = ttstep(pooltype,i) + 1
          riverstate%qqueue(l,pooltype,i) = 0.
          IF(numsubst>0.) riverstate%cqueue(:,l,pooltype,i) = 0.
        ENDIF
        evap = evapm3/area*1000.
      ENDIF   
    ENDIF

  END SUBROUTINE calculate_river_evaporation

  !>\brief Calculate total volume and mean T2 temperature concentration in river
  !----------------------------------------------------------------------------------------
  SUBROUTINE get_rivertempvol(i,pooltype,riverstate,meanrivertemp,totrivervol)

    USE MODVAR, ONLY : i_t2
    USE HYPEVARIABLES, ONLY : ttpart,ttstep
  
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    TYPE(riverstatetype),INTENT(IN) :: riverstate  !<River state
    REAL, INTENT(OUT)   :: meanrivertemp  !<temperature of river water
    REAL, INTENT(OUT)   :: totrivervol    !<volume of river water
    
    INTEGER l
    
    !Total volume in all river elements (translation boxes and river volume)
    totrivervol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
    
    IF(totrivervol.GT.0.)THEN
      !Weighted average T2 concentration
      meanrivertemp = riverstate%conc(i_t2,pooltype,i) * riverstate%water(pooltype,i)/totrivervol
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          meanrivertemp = meanrivertemp + riverstate%cqueue(i_t2,l,pooltype,i) * riverstate%qqueue(l,pooltype,i)/totrivervol
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          meanrivertemp = meanrivertemp + ttpart(pooltype,i) * riverstate%cqueue(i_t2,l,pooltype,i) * riverstate%qqueue(l,pooltype,i)/totrivervol
        ENDIF
      ENDIF
    ELSE
      meanrivertemp = 0.
    ENDIF

  END SUBROUTINE get_rivertempvol
  
  !>\brief Set a T2 temperature concentration to all river elements
  !----------------------------------------------------------------------------------------
  SUBROUTINE set_rivertemp(i,pooltype,riverstate,meanrivertemp)

    USE MODVAR, ONLY : i_t2
    USE HYPEVARIABLES, ONLY : ttpart,ttstep
  
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state
    REAL, INTENT(IN)   :: meanrivertemp !<temperature of river
    
    INTEGER l
    
    !Riverbox
    IF(riverstate%water(pooltype,i).GT.0.)THEN
      riverstate%conc(i_t2,pooltype,i) = meanrivertemp
    ELSE
      riverstate%conc(i_t2,pooltype,i) = 0.
    ENDIF
      
    !Translation boxes
    DO l = 1,ttstep(pooltype,i)
      IF(riverstate%qqueue(l,pooltype,i)>0)THEN
        riverstate%cqueue(i_t2,l,pooltype,i) = meanrivertemp
      ELSE
        riverstate%cqueue(i_t2,l,pooltype,i) = 0.
      ENDIF
    ENDDO

    IF(ttpart(pooltype,i)>0)THEN
      l = ttstep(pooltype,i) + 1
      IF(riverstate%qqueue(l,pooltype,i)>0)THEN
        riverstate%cqueue(i_t2,l,pooltype,i) = meanrivertemp
      ELSE
        riverstate%cqueue(i_t2,l,pooltype,i) = 0.
      ENDIF
    ENDIF
    
  END SUBROUTINE set_rivertemp
  
  !>\brief Calculate and remove evaporation from lake
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Evaporation)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_actual_lake_evaporation(i,j,itype,numsubst,temp,epot,evap,cevap,lakestate)

    USE MODVAR, ONLY : basin,classdata, &
                       landpar, &
                       cwater,  &
                       i_t1,i_t2, &
                       i_sm,i_gm,i_rn,i_li,i_ri,i_si
    USE HYPEVARIABLES, ONLY : m_ttmp

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: j        !<class index
    INTEGER, INTENT(IN) :: itype    !<lake type (ilake or olake)
    INTEGER, INTENT(IN) :: numsubst !<number of substances modelled
    REAL, INTENT(IN)    :: temp     !<air temperature
    REAL, INTENT(IN)    :: epot     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubst) !<concentration in evapotranspiration (eg. mg/L)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state

    !Local variables
    INTEGER k       !loop-variable, substance
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for evaporation (C)

    !Default values output variables
    evap = 0.
    cevap = 0.

    !Set local parameter
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    IF(temp>tt .AND. epot>0) THEN

      !Calculate actual evaporation, potential evaporation is default for temperature above threshold
      evap = epot
      DO k=1,numsubst
        cevap(k) = 0.
        IF(k==i_t1) cevap(k) = lakestate%conc(k,itype,i)  !t1 == O18 in hype
        !For T2==temperature, set the evaporation concentration = L/c = latent heat of vaporization divided by heat capacity
        IF(k==i_t2) cevap(k) = 1000. * latentheat_tempfunction(lakestate%conc(k,itype,i)) / cwater ! 1000 * [MJ/kg] / [KJ/kg/C] = C
        !For water origin model, set evaporation concentration equal riverstate conc
        IF(k==i_sm) cevap(k) = lakestate%conc(k,itype,i)  !sm == snowmelt, in water origin trace element model
        IF(k==i_gm) cevap(k) = lakestate%conc(k,itype,i)  !gm == glaciermelt, in water origin trace element model
        IF(k==i_rn) cevap(k) = lakestate%conc(k,itype,i)  !rn == rainfall, in water origin trace element model
        IF(k==i_li) cevap(k) = lakestate%conc(k,itype,i)  !li == initial lake, in water origin trace element model
        IF(k==i_ri) cevap(k) = lakestate%conc(k,itype,i)  !ri == initial river, in water origin trace element model
        IF(k==i_si) cevap(k) = lakestate%conc(k,itype,i)  !si == initial soil, in water origin trace element model
      ENDDO

      !Remove evaporation from lake, check if enough water is available              
      IF(evap<lakestate%water(itype,i))THEN
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
      ELSEIF(numsubst==0)THEN
        evap = lakestate%water(itype,i)
        cevap = lakestate%conc(:,itype,i)    !remove last traces of substances when lake dries out
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(2),basin(i)%subid,i,j)
      ELSEIF(lakestate%slowwater(itype,i)>0.)THEN
        !lake divided in two parts for N-, P- and C-simulations
        !mix lake parts and move to slowlake-part, evaporation from mixed lake volume
        CALL add_water(numsubst,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
        lakestate%water(itype,i)=0.
        lakestate%conc(:,itype,i)=0.
        IF(evap>=lakestate%slowwater(itype,i))THEN
          evap = lakestate%slowwater(itype,i)
          cevap = lakestate%concslow(:,itype,i)   !remove last traces of substances when lake dries out
        ENDIF
        CALL remove_water(lakestate%slowwater(itype,i),numsubst,lakestate%concslow(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(3),basin(i)%subid,i,j)
      ELSE
        !no water in slowlake-part, empty lake
        evap = lakestate%water(itype,i)
        cevap = lakestate%conc(:,itype,i)    !remove last traces of substances when lake dries out
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(3),basin(i)%subid,i,j)
      ENDIF
    ENDIF

  END SUBROUTINE calculate_actual_lake_evaporation

  !>Subroutine for summation of the area upstream of the outlet of all
  !>subbasins of the catchment
  !-------------------------------------------------------------------
  SUBROUTINE sum_upstream_area(n,areasum)
  
    USE MODVAR, ONLY : path,        &
                       basin,       &
                       branchdata,  &
                       branchindex

    !Argument declarations
    INTEGER, INTENT(IN)  :: n             !<number of subbasins
    REAL, INTENT(OUT)    :: areasum(n)    !<upstream area (m2)
    
    !Local variables
    INTEGER i,j,k,m         !loop variables
    REAL usarea             !summation variable for upstream area
    INTEGER, DIMENSION(n) :: A
    LOGICAL branchexists    !flag for branchdata available

    branchexists = .FALSE.
    IF(ALLOCATED(branchdata)) branchexists = .TRUE.
    A = 0
    areasum = 0.

    DO i = 1,n
      k = 0
      m = 0 
      USarea = 0. 
      DO j = 1,n
        IF(branchexists)THEN
          IF(branchindex(j)>0)THEN  !branch for this subbasin
            IF(i == path(j)%main)THEN
              m = m + A(j)
              USarea = USarea + areasum(j) * branchdata(branchindex(j))%mainpart
              k = k + 1      
            ENDIF
            IF(i == branchdata(branchindex(j))%branch)THEN
              m = m + a(j)
              usarea = usarea + areasum(j)*(1.-branchdata(branchindex(j))%mainpart)
              k = k + 1      
            ENDIF
          ELSE  !no branch for this subbasin
            IF(i == path(j)%main)THEN
              m = m + A(j)
              usarea = usarea + areasum(j)
              k = k + 1      
            ENDIF
          ENDIF
        ELSE    !no branches at all in the model set-up
           IF(i == path(j)%main)THEN
              m = m + a(j)
              usarea = usarea + areasum(j)
              k = k + 1      
           ENDIF
        ENDIF
      ENDDO
      IF(k==0) THEN                 !no inflows
        A(i) = 1
        areasum(i) = basin(i)%area
      ELSEIF(k==m) THEN             !k inflow, m (all) have their upstream area ready
        A(i) = 1
        areasum(i) = USarea + basin(i)%area
      ELSE                          !not all inflow have their upstream area ready (m<k)
        A(i) = 0                    !this indicates an error in coupling
        areasum(i) = 0.

        WRITE(6,*) 'ERROR in coupling of subbasins, some downstream basin before upstream basin'
        WRITE(6,*) 'i= ',i,' subid= ',basin(i)%subid
        WRITE(6,*) '(k= ',k,' m= ',m,')'
        STOP 1
      ENDIF
    ENDDO

  END SUBROUTINE sum_upstream_area

  !>Subroutine for calculation of general rating curve k-value for each lake
  !!
  !\b Reference ModelDescription Chapter Rivers and lakes (Lakes - Common lake processes)
  !-------------------------------------------------------------------
  SUBROUTINE set_general_rating_k(nl,n,locarea,areasum,rating)
    
    USE HYPEVARIABLES, ONLY : m_grat1,  &
                              m_grat3,   &
                              m_ratcorr, &
                              m_ilrrat1,m_olrrat1
    USE MODVAR, ONLY : genpar,  &
                       regpar,  &
                       basin, &
                       ilregpar,olregpar

    !Argument declarations
    INTEGER, INTENT(IN)  :: nl            !<number of lake types
    INTEGER, INTENT(IN)  :: n             !<number of subbasins
    REAL, INTENT(IN)     :: locarea(n)    !<landarea of subbasin [m2]
    REAL, INTENT(IN)     :: areasum(n)    !<upstream area [m2]
    REAL, INTENT(OUT)    :: rating(nl,n)  !<k-value (?)
    
    !Local variables
    INTEGER i         !loop variables
    REAL ratcorr

    DO i = 1,n
      IF(basin(i)%parregion>0)THEN
        ratcorr = 1. + regpar(m_ratcorr,basin(i)%parregion)   !Correction of general rating curve grat1-parameter
      ELSE
        ratcorr = 1.
      ENDIF
      rating(:,i) = genpar(m_grat1) * ratcorr
      IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar)) rating(1,i) = ilregpar(m_ilrrat1,basin(i)%ilakeregion) * ratcorr !ilake region parameter
      IF(basin(i)%olakeregion.GT.0 .AND. ALLOCATED(olregpar)) rating(2,i) = olregpar(m_olrrat1,basin(i)%olakeregion) * ratcorr !olake region parameter
    ENDDO

    IF(genpar(m_grat3)>0.)THEN
      DO i = 1,n
        IF(locarea(i)*basin(i)%ilakecatch>0.)rating(1,i) = rating(1,i)*(locarea(i)*basin(i)%ilakecatch)**genpar(m_grat3)
        IF(areasum(i)>0.) rating(2,i) = rating(2,i)*(areasum(i))**genpar(m_grat3)
      ENDDO
    ENDIF
      
  END SUBROUTINE set_general_rating_k

  !>\brief Calculates temperature of rivers and lakes and other temperature variables
  !> rivertemp: 20-day moving average of air temperature
  !> laketemp: 5-day moving average of air temperature
  !> Also 10- and 20-day moving average of lake- and river temperature is calculated.
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_water_temperature(i,airtemp,riverstate,lakestate)

    USE MODVAR, ONLY : genpar,        &
                       basin,         &
                       classbasin,    &
                       slc_ilake,     &
                       slc_olake,     &
                       i_in,i_sp,i_oc,  &
                       timesteps_per_day, &
                       ilregpar
    USE HYPEVARIABLES, ONLY : m_laketemp, &
                              m_gldepi, &
                              m_ilrldep

    !Argument declarations
    INTEGER, INTENT(IN) :: i                 !<index of current subbasin
    REAL, INTENT(IN)    :: airtemp           !<air temperature for subbasin
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    
    !Local parameters
    REAL, PARAMETER :: rivertemp_days = 20.     !Number of days for river temperature calculation
    REAL, PARAMETER :: laketemp_days  = 5.      !Number of days for lake temperature calculation
    REAL, PARAMETER :: T10day_parameter = 10.
    REAL, PARAMETER :: T20day_parameter = 20.
    
    !Local variables
    INTEGER watertype               !Internal or main/outlet
    REAL    mtimesteps,mtimesteps2  !Number of timesteps temperature is averaged over
    REAL gldepi
     
    !>\b Algorithm \n
    !>Calculate river temperature, same for local and main river
    IF(timesteps_per_day==1)THEN
      riverstate%temp(:,i) = riverstate%temp(:,i) + ((airtemp - riverstate%temp(:,i)) / rivertemp_days)
    ELSE
      mtimesteps = timesteps_per_day*rivertemp_days
      riverstate%temp(:,i) = riverstate%temp(:,i) + ((airtemp - riverstate%temp(:,i)) / mtimesteps)
    ENDIF

    !>Calculate lake temperature, same for internal and outlet lake (if exist)
    IF(genpar(m_laketemp)==0)THEN
      !>\li If parameter not set: as 5-day moving average
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
          IF(timesteps_per_day==1)THEN
            lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / laketemp_days)
          ELSE
            mtimesteps = timesteps_per_day*laketemp_days
            lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / mtimesteps)
          ENDIF
        ENDIF
      ENDIF
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          IF(timesteps_per_day==1)THEN
            lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / laketemp_days)
          ELSE
            mtimesteps = timesteps_per_day*laketemp_days
            lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / mtimesteps)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      !>\li Elseif parameter set: as a moving average of a period determined by lake depth
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
!          IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))THEN
!            gldepi = ilregpar(m_ilrldep,basin(i)%ilakeregion)
!          ELSE
!            gldepi = genpar(m_gldepi)
!          ENDIF
!          mtimesteps = timesteps_per_day*MIN(MAX(genpar(m_gldepi),5.),5.+genpar(m_laketemp))
!          mtimesteps = timesteps_per_day*MIN(MAX(gldepi,5.),5.+genpar(m_laketemp))
          mtimesteps = timesteps_per_day*MIN(MAX(basin(i)%lakedepth(1),5.),5.+genpar(m_laketemp))
          lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / mtimesteps)
        ENDIF
      ENDIF
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          mtimesteps = timesteps_per_day*MIN(MAX(basin(i)%lakedepth(2),5.),5.+genpar(m_laketemp))
          lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / mtimesteps)
        ENDIF
      ENDIF
    ENDIF

    !>Calculate 10- and 20-day mean of water temperature for N,P or C processes
    IF(i_in>0 .OR. i_sp>0 .OR. i_oc>0)THEN
      mtimesteps = timesteps_per_day*t10day_parameter
      mtimesteps2 = timesteps_per_day*t20day_parameter
      DO watertype = 1,2                   !(1=local/internal, 2=main/outlet)
        lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
        lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
        riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
        riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
      ENDDO
    ENDIF

  END SUBROUTINE calculate_water_temperature

  !>\brief Set temperature of rivers and lakes from T2 and calculate other temperature variables
  !Note: This is for today temperature
  !-----------------------------------------------------------------------
  SUBROUTINE set_water_temperature(waterbody,i,riverstate,lakestate)

    USE MODVAR, ONLY : classbasin,  &
                       slc_ilake,     &
                       slc_olake,     &
                       i_in,i_sp,i_oc,i_t2,  &
                       timesteps_per_day

    !Argument declarations
    INTEGER, INTENT(IN) :: waterbody         !<flag for waterbody
    INTEGER, INTENT(IN) :: i                 !<index of current subbasin
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    
    !Local parameters
    REAL, PARAMETER :: rivertemp_days = 20.     !Number of days for river temperature calculation
    REAL, PARAMETER :: laketemp_days  = 5.      !Number of days for lake temperature calculation
    REAL, PARAMETER :: T10day_parameter = 10.
    REAL, PARAMETER :: T20day_parameter = 20.
    
    !Local variables
    INTEGER watertype               !Internal or main/outlet
    REAL    mtimesteps,mtimesteps2  !Number of timesteps temperature is averaged over

    !>\b Algorithm \n
    !>Set river temperature to T2 temperature
    IF(waterbody==1) riverstate%temp(1,i) = riverstate%conc(i_t2,1,i)
    IF(waterbody==3) riverstate%temp(2,i) = riverstate%conc(i_t2,2,i)

    !>Set lake temperature (if exist)
    IF(waterbody==2)THEN
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
          lakestate%temp(1,i) = lakestate%conc(i_t2,1,i)
        ENDIF
      ENDIF
    ENDIF
    IF(waterbody==4)THEN
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          lakestate%temp(2,i) = lakestate%conc(i_t2,2,i)
        ENDIF
      ENDIF
    ENDIF

    !>Calculate 10- and 20-day mean of water temperature for N,P or C processes
    IF(i_in>0 .OR. i_sp>0 .OR. i_oc>0)THEN
      mtimesteps = timesteps_per_day*t10day_parameter
      mtimesteps2 = timesteps_per_day*t20day_parameter
      IF(waterbody==1.OR.waterbody==2)THEN 
        watertype = 1                   !(1=local/internal, 2=main/outlet)
        IF(waterbody==2) lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
        IF(waterbody==2) lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
        IF(waterbody==1) riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
        IF(waterbody==1) riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
      ENDIF
      IF(waterbody==3.OR.waterbody==4)THEN 
        watertype = 2                   !(1=local/internal, 2=main/outlet)
        IF(waterbody==4) lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
        IF(waterbody==4) lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
        IF(waterbody==3) riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
        IF(waterbody==3) riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
      ENDIF
    ENDIF

  END SUBROUTINE set_water_temperature

  !>\brief Calculate river characteristics
  !>River characteristics include depth, area, bankful flow and 365-day-average-Q
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_river_characteristics(i,itype,flow,calcNP,riverstate,depth,riverarea,qbank)
  
    USE MODVAR, ONLY : basin,   &
                       lregpar, &
                       genpar
    USE HYPEVARIABLES, ONLY : riverlength,  &
                              deadwidth,  &
                              m_velpar1,  &
                              m_velpar2,  &
                              m_velpar3,  &
                              m_widpar1,  &
                              m_widpar2,  &
                              m_widpar3,  &
                              m_maxwidth

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: itype      !<lake type (ilake or olake)
    REAL, INTENT(IN)    :: flow       !<river flow (m3/s) 
    LOGICAL, INTENT(IN) :: calcNP     !<status of NP simulation (to calculate bankful flow)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    REAL, INTENT(OUT)   :: depth      !<river depth (m)
    REAL, INTENT(OUT)   :: riverarea  !<river surface area (m2)
    REAL, INTENT(OUT)   :: qbank      !<flow at bank-ful river channel (m3/s)
    
    !Local variables
    REAL par(6)      !model parameters for velocity and width of river
    REAL rlength     !river length (m)
    REAL velocity    !river velocity (m/s)
    REAL width       !river width (m)

    !>\b Algorithm \n
    !>Set parameter values
    IF(basin(i)%lakeregion>0)THEN
      par = (/lregpar(m_velpar1,basin(i)%lakeregion),lregpar(m_velpar2,basin(i)%lakeregion),lregpar(m_velpar3,basin(i)%lakeregion),lregpar(m_widpar1,basin(i)%lakeregion),lregpar(m_widpar2,basin(i)%lakeregion),lregpar(m_widpar3,basin(i)%lakeregion) /)
    ELSE
      par = (/0,0,0,0,0,0/)   !OK? gives width=1 and depth=flow
    ENDIF

    !>Update state variable 365-day mean river discharge (m3/s)
    riverstate%Qmean(itype,i) = riverstate%Qmean(itype,i) + (flow-riverstate%Qmean(itype,i))/365.

    !>River length,depth and width, depend on velocity
    rlength = riverlength(itype,i)
    depth = 0.020          !low flow default value
    width = depth * 5.     !low flow default value
    IF(riverstate%Qmean(itype,i)>0.01.AND.flow>0.) THEN        
      velocity = (10**par(1)) * (riverstate%Qmean(itype,i)**par(2)) * ((flow/riverstate%Qmean(itype,i))**par(3))
      IF(velocity>0.2) THEN   
        width = (10**par(4)) * (flow/velocity)**(par(5)+par(6)*LOG10(flow/velocity))
        depth = (flow / velocity) / width             
      ENDIF
    ENDIF

    !>River (surface/bottom) area
    IF(genpar(m_maxwidth)>0)THEN
      riverarea = min(max(width,deadwidth(itype,i)),genpar(m_maxwidth)) * rlength
    ELSE  
      riverarea = max(width,deadwidth(itype,i)) * rlength
    ENDIF

    !>Calculate new bankfull flow, stored in Q2max. 
    IF(calcNP) CALL calc_qbank(flow,i,itype,riverstate%Q365(:,itype,i),riverstate%Qdayacc(:,itype,i),Qbank)  !Subroutine also updates Qmax, Qdayacc and riverQ365.
    
  END SUBROUTINE calculate_river_characteristics

  !>Estimates the bank full flow by the second highest q from the
  !>daily values of last year
  !>
  !>\b Consequences Module hypevariables variables qmax, q2mqx, iqmax, and iq2max 
  !> may change.
  !>
  !>\b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !---------------------------------------------------------------------
  SUBROUTINE calc_qbank(flow,i,itype,riverq365,Qdayacc,Qbank)
    
    USE HYPEVARIABLES, ONLY : qmax,q2max,   &  !OUT
                              iqmax,iq2max     !OUT
    USE MODVAR, ONLY : dayno,     &
                       timesteps_per_day, &
                       endofday,  &
                       tsofday

    !Argument declarations
    REAL, INTENT(IN)     :: flow    !<flow current time step (m3/s)
    INTEGER, INTENT(IN)  :: i       !<index of current subbasin
    INTEGER, INTENT(IN)  :: itype   !<river type 1=local, 2=main
    REAL,INTENT(INOUT)   :: riverq365(366)  !<river flow last 365 days (m3/s)
    REAL,INTENT(INOUT)   :: Qdayacc(timesteps_per_day)  !<river flow last day (m3/s)
    REAL, INTENT(OUT)    :: qbank   !<bankfull flow
    
    !local variables
    REAL q        !average flow for day (m3/s)

    !Accumulate flow values for calculation of daily mean
    Qdayacc(tsofday) = flow

    IF(endofday)THEN
      q = SUM(Qdayacc(:))/REAL(timesteps_per_day)
      riverq365(dayno) = q !First year: initial assignment, following years: overwrite oldest value
      !Estimate river bankful flow with second highest flow
      IF(dayno==iqmax(itype,i) .OR. dayno==iq2max(itype,i))THEN !too old values, search whole array for new
        CALL update_qbank(riverq365(:),qmax(itype,i),q2max(itype,i),iqmax(itype,i),iq2max(itype,i))
      ELSEIF(q > qmax(itype,i))THEN
        q2max(itype,i) = qmax(itype,i)     !new highest flow
        iq2max(itype,i) = iqmax(itype,i)
        qmax(itype,i) = q
        iqmax(itype,i) = dayno
      ELSEIF(q > q2max(itype,i))THEN    !new second highest flow
        q2max(itype,i) = q
        iq2max(itype,i) = dayno
      ENDIF
    ENDIF
    qbank = q2max(itype,i)

  END SUBROUTINE calc_qbank

  !>Update highest and second highest flow when one of them reach
  !>retiring age
  !---------------------------------------------------------------------
  SUBROUTINE update_qbank(q_array,qmax,q2,imax,i2)

    !Argument declarations
    REAL, INTENT(IN)     :: q_array(366)  !<flow all days last year
    REAL, INTENT(OUT)    :: qmax          !<highest flow all days last year
    REAL, INTENT(OUT)    :: q2            !<second highest flow all days last year
    INTEGER, INTENT(OUT) :: imax          !<index of highest flow all days last year
    INTEGER, INTENT(OUT) :: i2            !<index of second highest flow all days last year
    
    !Local variables
    INTEGER i

    qmax = 0.
    q2 = 0. 

    DO i = 1, 366
      IF(q_array(i) >= qmax)THEN 
        q2 = qmax
        i2 = imax
        qmax = q_array(i)
        imax = i
      ELSEIF(q_array(i) > q2)THEN
        q2 = q_array(i)
        i2 = i
      ENDIF
    ENDDO

  END SUBROUTINE update_qbank

  !>Translation (delay) in river       
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !-------------------------------------------------------------------
  SUBROUTINE translation_in_river(i,itype,qin,cin,qout,cout,riverstate)
  
    USE MODVAR, ONLY : numsubstances,     &
                       realzero,  &
                       seconds_per_timestep
    USE HYPEVARIABLES, ONLY : transtime,  &
                              ttstep,     &
                              ttpart

    !Argument declaration
    INTEGER, INTENT(IN) :: i                    !<index of current subbasin
    INTEGER, INTENT(IN) :: itype                !<river type (local or main)
    REAL, INTENT(IN)    :: qin                  !<inflow to river train (m3/s)
    REAL, INTENT(IN)    :: cin(numsubstances)   !<concentration of inflow to river train
    REAL, INTENT(OUT)   :: qout                 !<outflow of river train (m3/s)
    REAL, INTENT(OUT)   :: cout(numsubstances)  !<concentration of outflow of river train
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    !Local variables
    INTEGER y     !translation, whole time steps
    REAL    x     !translation, additional part of time step

    !>\b Algoritm \n
    !>Add new inflow to translation variable (river train)
    riverstate%qqueue(0,itype,i) = qin * seconds_per_timestep
    IF(numsubstances>0)riverstate%cqueue(:,0,itype,i) = cin

    IF(transtime(itype,i)>0)THEN
      !>Calculate outflow from river train
      y = ttstep(itype,i)
      x = ttpart(itype,i)
      qout = (1.-x)*riverstate%qqueue(y,itype,i) + x*riverstate%qqueue(y+1,itype,i) !Calculate flow (m3) from river after translation
      IF(qout>realzero)THEN
        cout = ((1.-x)*riverstate%qqueue(y,itype,i)*riverstate%cqueue(:,y,itype,i) + &
                 x*riverstate%qqueue(y+1,itype,i)*riverstate%cqueue(:,y+1,itype,i))/qout
      ELSE
        cout = 0.
      ENDIF
      qout = qout / seconds_per_timestep  !flow (m3/s)

      !>Translate the flows in the river train
      riverstate%qqueue(1:y+1,itype,i) = riverstate%qqueue(0:y,itype,i)
      IF(numsubstances>0)riverstate%cqueue(:,1:y+1,itype,i) = riverstate%cqueue(:,0:y,itype,i)
    ELSE
      !Elseif no delay, outflow = inflow
      qout = qin
      cout = cin
    ENDIF

  END SUBROUTINE translation_in_river


  !>Abstraction of water from main river
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources - Negative point source)
  !-------------------------------------------------------------------
  SUBROUTINE point_abstraction_from_main_river(i,pooltype,riverstate,removedflow)
  
    USE MODVAR, ONLY : basin, load,  &
                       seconds_per_timestep
    USE HYPEVARIABLES, ONLY : ttstep,     &
                              ttpart

    !Argument declaration
    INTEGER, INTENT(IN) :: i                          !<index of current subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<river type (local or main)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<river states
    REAL, INTENT(OUT)   :: removedflow                !<removed flow (m3/timestep)

    !Local variables
    INTEGER l     
    REAL    totvol      !volume in river (m3)
    REAL    absvol      !abstraction volume to be removed (m3)
    REAL    waterfrac   !fraction of water in current pool

    !>\b Algoritm \n
    removedflow = 0.
    IF(load(i)%abstrvol(1)==0) RETURN
    
    !>If abstraction of water: Calculate amount
    totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
    absvol = load(i)%abstrvol(1)*seconds_per_timestep  !m3
    !>Remove abstraction water proportionally from river and queue
    IF(absvol<totvol)THEN
      IF(riverstate%water(pooltype,i)>0.)THEN
        waterfrac = riverstate%water(pooltype,i)/totvol
        riverstate%water(pooltype,i) = riverstate%water(pooltype,i) - waterfrac*absvol
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0.)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that remaining outflow will be correct
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDIF
      removedflow = absvol
    ELSE
      riverstate%water(pooltype,i) = 0.
      riverstate%qqueue(1:ttstep(pooltype,i)+1,pooltype,i) = 0.
      WRITE(6,*) 'Warning: Point source abstraction from river could not be fulfilled, not enough water in river.'
      WRITE(6,*) 'Warning: subbasin ',basin(i)%subid, 'abstracted volume: ',totvol
      removedflow = totvol
    ENDIF

  END SUBROUTINE point_abstraction_from_main_river

  !>Abstraction of water from outlet lake
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources - Negative point source)
  !-------------------------------------------------------------------
  SUBROUTINE point_abstraction_from_outlet_lake(i,pooltype,qunitfactor,lakestate,removedflow)
  
    USE MODVAR, ONLY : load,  &
                       seconds_per_timestep

    !Argument declaration
    INTEGER, INTENT(IN) :: i                    !<index of current subbasin
    INTEGER, INTENT(IN) :: pooltype             !<lake type (local or outlet)
    REAL, INTENT(IN)    :: qunitfactor          !<transformation factor m3/s->mm/timestep
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states
    REAL, INTENT(OUT)   :: removedflow          !<removed flow (m3/timestep)

    !Local variables
    REAL    lakevol     !volume in lake (mm)
    REAL    absvol      !abstraction volume to be removed (mm)
    REAL    waterfrac   !fraction of water in current pool

    !>\b Algoritm \n
    removedflow = 0.
    IF(load(i)%abstrvol(2)==0) RETURN

    !>If abstraction of water: Calculate amount
    lakevol = lakestate%water(pooltype,i)       !mm
    IF(ALLOCATED(lakestate%slowwater)) lakevol = lakevol + lakestate%slowwater(pooltype,i)       !mm
    absvol = load(i)%abstrvol(2) * qunitfactor  !mm/ts
    
    !>Remove abstraction water proportionally from fast and slow lake part
    IF(absvol<lakevol)THEN
      IF(ALLOCATED(lakestate%slowwater))THEN
        waterfrac = lakestate%water(pooltype,i)/lakevol
        lakestate%water(pooltype,i) = lakestate%water(pooltype,i) - absvol*waterfrac
        lakestate%slowwater(pooltype,i) = lakestate%slowwater(pooltype,i) - absvol*(1.-waterfrac)
      ELSE
        lakestate%water(pooltype,i) = lakestate%water(pooltype,i) - absvol
      ENDIF
      removedflow = load(i)%abstrvol(2) * seconds_per_timestep
    ELSE
      lakestate%water(pooltype,i) = 0.
      IF(ALLOCATED(lakestate%slowwater)) lakestate%slowwater(pooltype,i) = 0.
      WRITE(6,*) 'Warning: Abstraction from lake could not be done, not enough water in lake.'
      removedflow = lakevol / qunitfactor * seconds_per_timestep
    ENDIF

  END SUBROUTINE point_abstraction_from_outlet_lake

  !>\brief Subroutine for finding current lake outflow parameters. 
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_lake_outflow_parameters(i,itype,lakeareain,olakewst,   &
                                         ratck,ratcexp,&
                                         upstreamlakebasin,w0Today,wmin,damProd,branchflow,&
                                         qin)
       
    USE HYPEVARIABLES, ONLY : ratingk,      &  
                              m_grat2,      &  
                              m_limprod,    &
                              m_ldlimprod,  &
                              m_gldepi,     &
                              m_krelflood,  &
                              m_kthrflood,  &
                              m_klowflood,  &
                              m_ilrrat2,    &
                              m_olrrat2, &
                              m_ilrldep
    USE MODVAR, ONLY : missing_value,     &
                       dayno,   &
                       pi,      &
                       lake,    &
                       dam,     &
                       lakebasin,  &
                       genpar,  &
                       lakeindex, &
                       damindex,    &
                       lakebasinindex,  &
                       lakedatapar, &
                       lakedataparindex, &
                       ilregpar,olregpar, basin

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2) (from GeoData)
    REAL, INTENT(IN)    :: olakewst      !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: ratck        !<current rating curve parameter rate
    REAL, INTENT(OUT)   :: ratcexp      !<current rating curve parameter exponent
    LOGICAL, INTENT(OUT):: upstreamlakebasin   !<Upstream lake basin?
    REAL, INTENT(OUT)   :: w0Today      !<current water level threshold in w-reference system (m)
    REAL, INTENT(OUT)   :: wmin         !<minimum water level threshold (sänkningsgräns) in w-reference system (m) 
    REAL, INTENT(OUT)   :: damProd      !<current dam production flow (m3/s)
    REAL,OPTIONAL,INTENT(IN) :: qin     !<current inflow to lake (m3/s)
    REAL,INTENT(INOUT):: branchflow !<outflow from dam branch (used subsequently in calculate_branch_flow)
    
    !Local variables
    REAL wlmr                   !water level lake (m)
    REAL deltaw0                !difference in water level threshold period 2 (m)
    REAL qamp,qpha              !parameters for regulation of lake
    REAL qprod1, qprod2, qprodToday   !Production flow (m3/s)
    REAL fracLevel              !Actual reservoir situation, in fraction of (w0 - wmin)
    REAL fillDamThreshold       !Percentage of reservoir capacity bellow which economy regime starts (from "par.txt")
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 
    REAL rating1,rating2        !general rating curve parameters outlet lake
    REAL regrate,regexp         !current parameters for specific rating curve
    INTEGER dayno1, dayno2      !Starting day nr. for production periods 1 and 2
    INTEGER dampurpose          ! Purpose of dam, 1=irrigation, 2=water supply, 3=flood control, 4=hydropower, 5+=special
    REAL qinfmax                ! Max mean monthly inflow
    REAL snowfrac               ! Fraction of precipitaiton falling as snow in catchment upstream of dam
    REAL qprod0                 !Production flow which is equal to mean inflow
    REAL regvol 
    REAL qinftoday              ! Current inflow for 'today'
    REAL qthresh                ! Threshold inflow over which a flood control dam save water
    REAL lthresh                ! Threshold reservoir level over which flood control dam releases extra flow
    REAL qinfmin                ! Min mean monthly inflow
    
    !MH2017:branchflow is not being initialized, candidate for deletion
    !****************************
    branchflow = 0.
    !***************************
    !Initial values
    wlmr=0.
    damProd = 0. 
    qprodToday = 0.
    qprod0 = 0.
    qprod1 = 0.
    qprod2 = 0.
    dayno1 = 0
    dayno2 = 0
    qamp = 0.
    qpha = 0.
    regvol=0.
    qinfmax = 0.
    qinfmin=0.
    snowfrac = 0.
    dampurpose = 0
    regrate = 0.
    regexp = 0.
    ratck = 0.
    ratcexp = 0.
    deltaw0 = 0.
    wmin = missing_value
    lakearea = lakeareain    
    qinftoday=0.  
    upstreamlakebasin = .FALSE.
    lthresh = 0.  
    qthresh = 0.  
    
    

    !Current lake parameter values
    rating1 = ratingk(itype,i)                    !General rating parameters for lake
    rating2 = genpar(m_grat2)                     !-"-
    !ilake and olake region parameters
    IF(itype.EQ.1 .AND. basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))rating2=ilregpar(m_ilrrat2,basin(i)%ilakeregion)
    IF(itype.EQ.2 .AND. basin(i)%olakeregion.GT.0 .AND. ALLOCATED(olregpar))rating2=olregpar(m_olrrat2,basin(i)%olakeregion)

    
    fillDamThreshold = genpar(m_limprod)          !Threshold for fractional height, below which the production rate decreases linearly
    
    IF(PRESENT(qin))THEN
      qinftoday=qin
    ENDIF
     
    IF(ALLOCATED(lakedatapar)) fillDamThreshold = lakedatapar(lakedataparindex(i,itype),m_ldlimprod)

    IF(itype==2)THEN
      IF(ALLOCATED(lakeindex))THEN
        IF(lakeindex(i)>0)THEN
          qprod1 = lake(lakeindex(i))%qprod1          !Production discharge for period 1
          qprod2 = lake(lakeindex(i))%qprod2          !                     for period 2
          dayno1 = lake(lakeindex(i))%datum1          !Starting day nr. for production period 1
          dayno2 = lake(lakeindex(i))%datum2          !                                period 2
          qamp = lake(lakeindex(i))%qamp              !Amplitude of sin-adjustment of qprod
          qpha = lake(lakeindex(i))%qpha              !Phase of sin-adjustment of qprod
          regrate = lake(lakeindex(i))%rate           !Rating curve parameter
          regexp = lake(lakeindex(i))%exp             !Rating curve parameter
          deltaw0 = lake(lakeindex(i))%deltaw0        !difference in lake threshold/"dämningsgräns" period 2
          wmin = lake(lakeindex(i))%wmin              !lake threshold/"sänkningsgräns"
        ENDIF
      ENDIF
      IF(ALLOCATED(lakebasinindex))THEN
        IF(lakebasinindex(i)>0)THEN  
          IF(lakebasin(lakebasinindex(i))%last)THEN          !Recalculate water stage for last lake basin
            qprod1 = lake(lakebasin(lakebasinindex(i))%ilk)%qprod1          !Production discharge for period 1
            qprod2 = lake(lakebasin(lakebasinindex(i))%ilk)%qprod2          !                     for period 2
            dayno1 = lake(lakebasin(lakebasinindex(i))%ilk)%datum1          !Starting day nr. for production period 1
            dayno2 = lake(lakebasin(lakebasinindex(i))%ilk)%datum2          !                                period 2
            qamp = lake(lakebasin(lakebasinindex(i))%ilk)%qamp              !Amplitude of sin-adjustment of qprod
            qpha = lake(lakebasin(lakebasinindex(i))%ilk)%qpha              !Phase of sin-adjustment of qprod
            regrate = lake(lakebasin(lakebasinindex(i))%ilk)%rate           !Rating curve parameter
            regexp = lake(lakebasin(lakebasinindex(i))%ilk)%exp             !Rating curve parameter
            deltaw0 = lake(lakebasin(lakebasinindex(i))%ilk)%deltaw0        !difference in lake threshold/"dämningsgräns" period 2
            wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin              !lake threshold/"sänkningsgräns"
          ELSE
            upstreamlakebasin = .TRUE.
            deltaw0 = lake(lakebasin(lakebasinindex(i))%ilk)%deltaw0        !difference in lake threshold/"dämningsgräns" period 2
            wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin              !lake threshold/"sänkningsgräns"
            dayno1 = lake(lakebasin(lakebasinindex(i))%ilk)%datum1          !Starting day nr. for production period 1
            dayno2 = lake(lakebasin(lakebasinindex(i))%ilk)%datum2          !                                period 2
          ENDIF
        ENDIF
      ENDIF
      IF(ALLOCATED(damindex))THEN
        IF(damindex(i)>0)THEN
          regvol = dam(damindex(i))%regvol          !Regvol
          qprod1 = dam(damindex(i))%qprod1          !Production discharge for period 1 
          qprod2 = dam(damindex(i))%qprod2          !                     for period 2
          dayno1 = dam(damindex(i))%datum1          !Starting day nr. for production period 1
          dayno2 = dam(damindex(i))%datum2          !                                period 2
          snowfrac = dam(damindex(i))%snowfrac      ! Fraction of prec falling as snow upstream of dam
          qamp = dam(damindex(i))%qamp              !Amplitude of sin-adjustment of qprod
          qpha = dam(damindex(i))%qpha              !Phase of sin-adjustment of qprod
          regrate = dam(damindex(i))%rate           !Rating curve parameter
          regexp = dam(damindex(i))%exp             !Rating curve parameter
          wmin = dam(damindex(i))%wmin              !lake threshold/"sänkningsgräns"
          qprod0 = dam(damindex(i))%qinfmed         !natural flows used for determining outflow
          qinfmin = dam(damindex(i))%qinfmin
          qinfmax = dam(damindex(i))%qinfmax
          dampurpose = dam(damindex(i))%purpose              !dam purpose
          lthresh = 0.-genpar(m_klowflood)*regvol*1000000./lakearea       !threshold level for extra flood control releases (typical 1/3 of regvol)
          qthresh = genpar(m_kthrflood)*qinfmax 
          
        ENDIF
      ENDIF
    ENDIF   !itype=2
      
    !Local lake parameter values
    IF(itype==1)THEN
!      w0Today = genpar(m_gldepi)      !ilake depth = threshold (m)
!      IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))w0Today = ilregpar(m_ilrldep,basin(i)%ilakeregion)
      w0Today = basin(i)%lakedepth(1)
      ratck   = rating1
      ratcexp = rating2
      RETURN
    ENDIF

    !Outlet lake/dam parameter values
    IF(itype==2)THEN

      !Water level for outlet lake
      wlmr = olakewst
     
      !Dam in DamData
      IF(damindex(i)>0)THEN    !The below code is repeated first for dams, then for lakes (can be later removed for lakes if reg only occurs on damdata.txt)
        !CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref) !Reads in lakearea for dam and wlmr
              
        !Current production flow for dam with regulation volume (calibrated dams)
        IF(wmin.NE.missing_value)THEN       ! Wmin.NE.missing if RegVol> 0
          IF(qprod1>0.)THEN  !If specific Qprod1 is given in DamData.txt
            qprodToday=qprod1                                !Production rate 1 applies if specified                                        
            IF(dayno1*dayno2 > 0) THEN                       !If both dates for different production regimes are non-zero, ...
              IF(dayno < dayno1 .OR. dayno >= dayno2)THEN
                qprodToday = qprod2   !... and today is not within the periode of production rate 1, then production rate 2 applies
              ENDIF
            ENDIF
          ELSE 
            qprodToday = qprod0            !Otherwise QprodToday = mean inflow
          ENDIF                                 
          IF(fillDamThreshold>0)THEN    !Reduce production flow for low water levels (save water in dam)
            fracLevel = (wlmr - wmin)/(0. - wmin)                         !Current reservoir height, in fraction of height difference between min and top (can be negative if dam is below minimum level)
            IF(fracLevel > 0 .AND. fracLevel < fillDamThreshold) THEN
              qprodToday = fracLevel/fillDamThreshold * qprodToday        !Economy regime, if reservoir fractional filling is lower than threshold
            ENDIF
          ENDIF

          !Calculate dam outflow depending on dam purpose and current production flow
          selectPurpose:IF(dampurpose==4)THEN    !hydroelectric dam
            IF(qamp>0)THEN
              damProd = qprodToday*(1. + qamp * SIN(2.*pi*(dayno+qpha)/365.))  ! if qamp and qpha given, use them
            ELSEIF(snowfrac > 0.35)THEN     !CD2014 This determines if dams seasonally change flow by comparing regulation capacity to inflows (i.e. for dams wtih more than 35 % of precip that is snow
              IF(qpha>0) THEN
                qamp = 0.71                 !CD2014 based on regression from data, used if Qamp not given.                               
                damProd = qprodToday*(1. + qamp * SIN(2.*pi*(dayno+qpha)/365.))   !Seasonal variation according to sin-function
              ENDIF
            ELSE
              damProd=qprodToday              ! If not a snow/seasonal redist dam and qamp not given, damProd=constant (Qinf or Qprod1)
            ENDIF
          ELSEIF(dampurpose==1)THEN selectPurpose  !irrigation dam
            damProd=qprodToday
          ELSEIF(dampurpose==2)THEN selectPurpose  !water supply dam
            damProd=qprodToday
          ELSEIF(dampurpose==3)THEN selectPurpose  !flood control dam (aim is for dam to maintain as empty as possible)
            IF(qinftoday < qthresh)THEN        ! IF inflow today < threshold inflow
              IF(wlmr<lthresh)THEN              ! IF water level today < threshold level
                damProd=Qinftoday                ! Release the inflow
              ELSE      
                IF(qinftoday < qinfmin)THEN
                  damProd=MIN(Qinftoday*genpar(m_krelflood),qthresh)     ! If water level above threshold, release more than inflows (i.e. try empty the dam)
                ENDIF
!                IF(qinftoday > qinfmin)THEN                         !CP201506 testing
!                  damProd=MIN(Qinftoday*genpar(m_krelflood),qthresh)     ! If water level above threshold, release more than inflows (i.e. try empty the dam)
!                ELSE
!                  damProd = Qinftoday
!                ENDIF
              ENDIF
            ELSE                                 ! If inflow today >= threshold inflow
              damProd=qthresh                    ! Release maximum allowable flow
            ENDIF
          ELSEIF(dampurpose==5)THEN selectPurpose !inline dam
            CALL inline_dam_outflow(i,wlmr,damProd)
          ELSEIF(dampurpose==6)THEN selectPurpose !OCBD dam
            CALL OCBD_dam_outflow(i,wlmr,damProd,branchflow)
          ELSE selectPurpose
            damProd=qprodToday                                                         
          ENDIF selectPurpose                                   !End of IF statements for outflow depending on damtype !CD2014
          !damProd=MIN(damProd,regvol*1000000./86400.)          ! Test: Limit damProd to the Regvol for one day
        ENDIF   !Close wmin.ne.missing loop

        !Ordinary outflow threshold
        w0Today = 0.

        !Current threshold for dam with production by rating curves
        IF(deltaw0/=0.)THEN
          IF(dayno1*dayno2 > 0)THEN                             
            IF(dayno < dayno1 .OR. dayno >= dayno2)  w0Today = w0Today + deltaw0   !threshold for period 2 applies
          ENDIF
        ENDIF
        
        !Set rating curve parameters
        IF(wlmr>w0Today)THEN
          IF(regrate>0.)THEN                         !Specific rating curve for lake or dam spill
            ratck = regrate
            ratcexp = regexp
          ELSEIF(wmin.NE.missing_value)THEN         !Dam without rating curve for spill
          ELSE                                      !General rating curve for lake
            ratck = rating1
            ratcexp = rating2
          ENDIF
        ENDIF     

      !Lake/Dam in LakeData or GeoData
      ELSE      ! NOTE ALL? THE FOLLOWING CAN BE REMOVED IF REGULATION ONLY CALCULATED ON DAMDATA file (lakebasinlake regulation??)
!        CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)
        
        !Current production flow for dam with regulation volume
        IF(wmin.NE.missing_value)THEN
          qprodToday = qprod1                                               !Production rate 1 applies in general
          IF(dayno1*dayno2 > 0) THEN                                       !If both dates for different production regimes are non-zero, ...
            IF(dayno < dayno1 .OR. dayno >= dayno2)  qprodToday = qprod2   !... and today is not within the periode of production rate 1, then production rate 2 applies
          ENDIF
          IF(fillDamThreshold>0)THEN
            fracLevel = (wlmr - wmin)/(0. - wmin)                         !Current reservoir height, in fraction of height difference between min and top (can be negative if dam is below minimum level)
            IF(fracLevel > 0. .AND. fracLevel < fillDamThreshold) THEN
              qprodToday = fracLevel/fillDamThreshold * qprodToday        !Economy regime, if reservoir fractional filling is lower than threshold
            ENDIF
          ENDIF
          damProd = qprodToday*(1. + qamp * sin(2.*pi*(dayno+qpha)/365.))   !Seasonal variation according to sin-function
        ENDIF
 
        !Ordinary outflow threshold
        w0Today = 0.                                               !Threshold w0 applies in general

        !Current threshold for dam with production by rating curves
        IF(deltaw0/=0.)THEN
          IF(dayno1*dayno2 > 0) THEN                             
            IF (dayno < dayno1 .OR. dayno >= dayno2)  w0Today = w0Today + deltaw0   !threshold for period 2 applies
          ENDIF
        ENDIF
 
        IF(upstreamlakebasin)THEN                   !Upstream lake basin; all water above threshold
          IF(wmin==missing_value)THEN
            w0Today = w0Today     
          ELSE
            w0Today = wmin      
          ENDIF
        ELSEIF(wlmr>w0Today)THEN
          IF(regrate>0.)THEN                         !Specific rating curve for lake or dam spill
            ratck = regrate
            ratcexp = regexp
          ELSEIF(wmin.NE.missing_value)THEN         !Dam without rating curve for spill
          ELSE                                      !General rating curve for lake
            ratck = rating1
            ratcexp = rating2
          ENDIF
        ENDIF
        
      ENDIF   ! ENDIF for lakes (i.e. damindex =/0)
    ENDIF   !itype==2
    

  END SUBROUTINE get_current_lake_outflow_parameters
 
  !>\brief Subroutine for calculation outflow from lake. 
  !!General rating curve for ilakes. 
  !!For outlet lakes several options exist: 
  !!Specific rating curve, general rating curve, all water above threshold for 
  !!upstream lake basin, regulation with spill by rating curve, 
  !!constant production flow depending on date or two separate rating curves for olake.
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_outflow_from_lake(i,itype,qin,lakeareain,lakewstmm,   &
                                         qunitfactorin,outflowm3s,outflowmm, &
                                         lakestate,branchflow)
       
    USE MODVAR, ONLY : missing_value,    &
                       dam,     &
                       damindex

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    REAL, INTENT(IN)    :: qin           !<inflow of lake (m3/s) 
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2)
    REAL, INTENT(IN)    :: lakewstmm     !<lake water stage (mm)
    REAL, INTENT(IN)    :: qunitfactorin !<factor for transforming flow for lake from m3/s to mm/timestep and back
    REAL, INTENT(OUT)   :: outflowm3s    !<outflow of lake (m3/s)
    REAL, INTENT(OUT)   :: outflowmm     !<outflow of lake (mm)
    REAL,INTENT(INOUT):: branchflow !<outflow from dam branch (used subsequently in calculate_branch_flow)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state

    !Local variables
    LOGICAL upstreamlakebasin   !Upstream lake basin
    REAL wlmr                   !water level lake  (m)
    REAL wmin                   !water levels threshold for production (m)
    REAL w0ref                  !water stage reference level (not used)
    REAL w0Today                !water level threshold  (m)
    REAL ratingc,ratinge        !current rating curve parameters outlet lake
    REAL damProd                !Current dam production flow (m3/s)
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 
    REAL qunitfactor            !factor for transforming flow for lake from m3/s to mm/timestep and back (adjusted for last lakebasin) 

    !Initial values
    outflowm3s = 0.
    wlmr = lakewstmm*0.001              !Water in lake (m) (default=absolute)
    lakearea = lakeareain      
    qunitfactor = qunitfactorin
    
    !Water level for outlet lake
    IF(itype==2) CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)

    !Current parameter values (ilake and olake)
    CALL get_current_lake_outflow_parameters(i,itype,lakeareain,wlmr,   &
                                         ratingc,ratinge,upstreamlakebasin,w0Today,wmin,damProd,branchflow,&
                                         qin)
   
    !Outflow from local lake
    IF(itype==1)THEN
      IF(wlmr>w0Today)  outflowm3s = average_flow_rating_curve(qin,lakearea,wlmr,ratingc,ratinge,w0Today)  !General rating curve (m3/s)
      outflowmm = outflowm3s * qunitfactor        !to mm/timestep
    ENDIF

    IF(itype==2)THEN
!      CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)    !Water level for outlet lake

      !Outflow determination
      IF(upstreamlakebasin)THEN         !Upstream lake basin; all water above threshold
        IF(wlmr-w0Today>0) outflowm3s = (wlmr-w0Today) * 1000. / qunitfactor          
      ELSEIF(wlmr>w0Today)THEN          !Water level above threshold
        IF(ratingc>0)THEN                         !Specific rating curve for lake or dam spill
          outflowm3s = MAX(average_flow_rating_curve(qin,lakearea,wlmr,ratingc,ratinge,w0Today), damProd)
        ELSEIF(wmin.NE.missing_value)THEN         !Dam without rating curve for spill; all water above threshold but at least production
          outflowm3s = MAX((wlmr-w0Today) * 1000. / qunitfactor, damProd)
        ELSE
          WRITE(6,*) 'Error: Ended in else that is not possible. Outflow of outlet lake'
          WRITE(6,*) 'Check input data for this lake.'
          WRITE(6,*) 'i',i,'itype',itype
          WRITE(6,*) 'More info: wlmr',wlmr,'w0Today',w0Today,'wmin',wmin,'ratingc',ratingc
        ENDIF
      ELSEIF(wmin==missing_value)THEN   !No regulated lake flow
        !no damflow
      ELSEIF(wlmr>wmin)THEN             !Production flow to lower threshold
        outflowm3s = MIN((wlmr-wmin) * 1000. / qunitfactor, damProd)
      ENDIF
      ! MH2017: seperate IF statement so that cases with damindex(i)>0, but not special, get recorded
      IF(damindex(i)>0) THEN ! MH2017: damindex(i) can be 0, which isn't allowed for dam()
        IF(dam(damindex(i))%purpose >= 5) outflowm3s = damProd ! MM2016: special dam case  dam(damindex(i))%purpose==5
      ENDIF 
      outflowmm = outflowm3s * qunitfactor        !to mm/d
    
    ENDIF   !itype==2

  END SUBROUTINE calculate_outflow_from_lake

  !>\brief Momentanous flow by rating curve
  !>
  !>Subroutine for calculation momentanous outflow from lake from current lake 
  !>water stage by simple lake rating curve equation. 
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_flow_from_lake_waterstage(i, &
                                                 itype, &
                                                 lakeareain, &
                                                 lakewstmm,   &
                                                 outflowm3s, &
                                                 lakestate, &
                                                 branchflow)
       
    USE GENERAL_FUNCTIONS, ONLY : simple_rating_curve

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2)
    REAL, INTENT(IN)    :: lakewstmm     !<lake water stage (mm)
    REAL, INTENT(OUT)   :: outflowm3s    !<outflow of lake (m3/s)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state
    REAL,INTENT(INOUT)  :: branchflow    !<outflow from dam branch (used subsequently in calculate_branch_flow)
    
    !Local variables
    LOGICAL upstreamlakebasin   !Upstream lake basin; Tappa tillrinning
    REAL wlmr                   !water level lake (m)
    REAL wmin                   !water level threshold in w-reference system (m)
    REAL w0Today                !water level threshold in w-reference system (m)
    REAL w0ref                  !water stage reference level (m) (not used)
    REAL ratingc,ratinge        !general rating curve parameters outlet lake
    REAL damProd                !Dam production flow
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 

    !Initial values
    outflowm3s = 0.
    lakearea = lakeareain     
    wlmr = lakewstmm*0.001              !Water in lake (m) (default=absolute)
    IF(itype==2) CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)
    
    !Current parameter values
    CALL get_current_lake_outflow_parameters(i,itype,lakeareain,wlmr,   &
                                             ratingc,ratinge,upstreamlakebasin,w0Today,wmin,damProd,&
                                             branchflow)

    !Outflow determination
    IF(ratingc>0)THEN
      outflowm3s = simple_rating_curve(wlmr,ratingc,ratinge,w0Today)
    ELSE
      outflowm3s = 0. !Error in indata reaching this else?
    ENDIF

  END SUBROUTINE calculate_flow_from_lake_waterstage

  !>\brief Subroutine for checking lake outflow against lake volume 
  !(bad rating curve parameters or numerical problems)
  !------------------------------------------------------------------------------
  SUBROUTINE check_outflow_from_lake(i,itype,lakeareain,lakewstmm, &
                                     qunitfactor,outflowm3s,outflowmm,  &
                                     lakestate,branchflow,qin)
       
    USE MODVAR, ONLY : missing_value,     &
                       basin,   &
                       lakebasin,  &
                       lakebasinindex, &
                       damindex, dam

    !Argument declarations
    INTEGER, INTENT(IN) :: i            !<index of current subbasin
    INTEGER, INTENT(IN) :: itype        !<lake type (local or main)
    REAL, INTENT(IN)    :: lakeareain   !<lakearea (m2)
    REAL, INTENT(IN)    :: lakewstmm    !<lake water stage (mm)
    REAL, INTENT(IN)    :: qunitfactor  !<factor for transforming flow for lake from m3/s to mm/d and back
    REAL, INTENT(INOUT) :: outflowm3s   !<outflow of lake (m3/s)
    REAL, INTENT(INOUT) :: outflowmm    !<outflow of lake (mm)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state 
    REAL,INTENT(INOUT) ::  branchflow   !<branchflow for special dam (purpose 5)
    REAL,INTENT(INOUT) ::  qin          !<inflow to special dam (purpose 5)
    !REAL branchflow

    !Local variables
    LOGICAL upstreamlakebasin   !Upstream lake basin; Tappa tillrinning
    LOGICAl notspecialdam       !.FALSE. for special dams (dam(i)%purpose>=5), .TRUE. otherwise
    REAL wlm,wlmr               !water level lake (m)
    REAL ldepthm                !depth of lake (to threshold w0) (m)
    REAL wmin                   !water level lower threshold (m)
    REAL w0Today                !water level threshold (m)
    REAL w0ref                  !water stage reference level (m) (not used)
    REAL ratingc,ratinge,damProd
    REAL wthresh                !current water level threshold
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 

    !>\b Algorithm \n
    !Initial values
    wlm = lakewstmm*0.001               !Water in lake (m)
    IF(itype==2) CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)
    
    !>Get current parameter values
    CALL get_current_lake_outflow_parameters(i,itype,lakeareain,wlmr,   &
                                         ratingc,ratinge,upstreamlakebasin,w0Today,wmin,damProd,branchflow,qin)
    
    !>Check outflow from local lake
    IF(itype==1)THEN
      IF(outflowmm*0.001>wlm-w0Today)THEN   !Check for enough water in lake (bad rating curve or numerical problems)
        IF(wlm>w0Today)THEN
          outflowmm = (wlm-w0Today)*1000.
        ELSE
          outflowmm = 0.
        ENDIF
        IF(outflowmm>lakewstmm) outflowmm = lakewstmm   !Safety for rounded wlm used and ldepthm = 0
        outflowm3s = outflowmm/qunitfactor
      ENDIF
    ENDIF

    !special dams can break when using the unnamed IF block, damcheck skips it for special dams
    damcheck:IF(itype==2.AND.damindex(i)>0)THEN
      IF(dam(damindex(i))%purpose>=5)THEN
        notspecialdam=.FALSE.
      ELSE
        notspecialdam=.TRUE.
      ENDIF
    ELSEIF(itype==2)THEN damcheck
      notspecialdam=.TRUE.
    ENDIF damcheck

    !>Check outflow for enough water in olake
    IF(itype==2.AND.notspecialdam)THEN

      !>\i Calculate current threshold
      IF(wmin==missing_value .OR. upstreamlakebasin)THEN
        wthresh = w0Today
      ELSE  
        wthresh = wmin
      ENDIF

      !>\i Calculate current water stage
!      CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)
      IF(ALLOCATED(lakebasinindex))THEN
        IF(lakebasinindex(i)>0)THEN  
          IF(lakebasin(lakebasinindex(i))%last)THEN
            ldepthm = basin(i)%lakedepth(2)          !Olake threshold depth (m)
            wlmr = wlm - ldepthm       !Recalculate water stage for last lake basin to lake basin level, because this is the volume available
          ENDIF
        ENDIF
      ENDIF
      
      !>\i Check against lowest water level allowed
      IF(outflowmm*0.001>wlmr-wthresh)THEN
        IF(wlmr>wthresh)THEN
          outflowmm = (wlmr-wthresh)*1000.
        ELSE
          outflowmm = 0.
          !MM2016: for dampurpose==5
          branchflow = 0.
        ENDIF
        IF(outflowmm>lakewstmm) outflowmm = lakewstmm   !Safety for rounded wlmr used. 
        outflowm3s = outflowmm/qunitfactor
      ENDIF
    
    ENDIF   !itype==2

  END SUBROUTINE check_outflow_from_lake

  !>\brief Removal of outflow from lake and setting of
  !>concentration of outflow.
  !>
  !>\b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !-----------------------------------------------------------------------
  SUBROUTINE remove_outflow_from_lake(i,itype,ns,outflowmm,subid,coutflow,lakestate)

    USE HYPEVARIABLES, ONLY : m_lddeeplake, &
                              m_ldfastlake
    USE MODVAR, ONLY : lakedatapar,     &
                       lakedataparindex,&
                       dam,damindex,basin !MH2017:allows subroutine to skip removing water from special dams

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    INTEGER, INTENT(IN) :: ns            !<number of substances
    REAL, INTENT(IN)    :: outflowmm     !<outflow of lake (mm/timestep)
    INTEGER, INTENT(IN) :: subid         !<subid of current subbasin, for error output
    REAL, INTENT(OUT)   :: coutflow(ns)  !<concentration of outflow of lake
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER status
    REAL fastoutpart  !fraction of outflow from fast lake part (if possible)
    REAL q1,q2        !outflow from lake parts (mm/timestep)

    !>\b Algorithm \n
    !>Preparations: default output
    coutflow = 0.


    
      
    !>Remove outflow and set outflow concentrations:
    !>If lake is not divided:
    IF(ns==0)THEN
      IF(outflowmm>0)THEN
        !>\li Outflow is removed from lake
        coutflow = lakestate%conc(:,itype,i)
        CALL remove_water(lakestate%water(itype,i),ns,lakestate%conc(:,itype,i),outflowmm,coutflow,status)
        IF(status.NE.0) CALL error_remove_water(errstring(4),subid,i,itype)
      ENDIF
    ELSE
      IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN
        !>If divided totally mixed lake (only slowlake):
        !>\li Add all water to slowlake and mix
        IF(lakestate%water(itype,i)>0)THEN   !add lakewater to slowwater and mix
          CALL add_water(ns,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
          lakestate%water(itype,i)=0.
          lakestate%conc(:,itype,i)=0.
        ENDIF
        !>\li Outflow is removed from slowlake
        IF(outflowmm>0)THEN
          coutflow(:) = lakestate%concslow(:,itype,i)
          !MH2017:skip removal of water for special dams
          !intent is to preserve NPC processes
          CALL remove_water(lakestate%slowwater(itype,i),ns,lakestate%concslow(:,itype,i),outflowmm,coutflow,status)
          IF(damindex(i)>0.AND.itype==2)THEN       
            IF(dam(damindex(i))%purpose>=5)THEN !special dam
            !comment out to ENDIF if testing Calibration
              IF(lakestate%water(itype,i)<0)THEN   !if lakewater addition was skipped due to lakestate%water(itype,i) < 0
                CALL add_water(ns,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
                lakestate%water(itype,i)=0.
                lakestate%conc(:,itype,i)=0.
              ENDIF
              lakestate%slowwater(itype,i)=((dam(damindex(i))%stortod - dam(damindex(i))%stord) / dam(damindex(i))%stora) ** (1/dam(damindex(i))%storc) + dam(damindex(i))%storb
              lakestate%slowwater(itype,i)=(lakestate%slowwater(itype,i) - dam(damindex(i))%w0ref + basin(i)%lakedepth(2)) * 1000.
              status = 0
            ENDIF
          ENDIF
          IF(status.NE.0) CALL error_remove_water(errstring(5),subid,i,itype)
        ENDIF
      ELSE
        !>If lake is divided in fast and slow part and outflow occur from both parts:
        !>\li Calculate the outflow fraction from each part
        IF(outflowmm>0)THEN
          fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
          q1 = fastoutpart*outflowmm
          q2 = (1.-fastoutpart)*outflowmm
          IF(q1>lakestate%water(itype,i))THEN
            q2 = q2 + (q1 - lakestate%water(itype,i))
            q1 = lakestate%water(itype,i)
          ENDIF
          IF(q2>lakestate%slowwater(itype,i))THEN
            q1 = q1 + (q2 - lakestate%slowwater(itype,i))
            q2 = lakestate%slowwater(itype,i)
          ENDIF
          !>\li Remove the outflow from both parts
          IF(q1>0.)THEN
            CALL remove_water(lakestate%water(itype,i),ns,lakestate%conc(:,itype,i),q1,lakestate%conc(:,itype,i),status)
            IF(status.NE.0) CALL error_remove_water(errstring(6),subid,i,itype)
          ENDIF
          IF(q2>0)THEN 
            CALL remove_water(lakestate%slowwater(itype,i),ns,lakestate%concslow(:,itype,i),q2,lakestate%concslow(:,itype,i),status)
            IF(status.NE.0) CALL error_remove_water(errstring(7),subid,i,itype)
          ENDIF
          !>\li Calculate the concentration of outflow
          coutflow(:) = (q1*lakestate%conc(:,itype,i)+q2*lakestate%concslow(:,itype,i))/outflowmm
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE remove_outflow_from_lake

  !>\brief Flow between lake parts for divided lake
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !-------------------------------------------------------------------
  SUBROUTINE calculate_flow_within_lake(i,itype,subid,lakestate)

    USE HYPEVARIABLES, ONLY : slowlakeini,  &
                              m_lddeeplake
    USE MODVAR, ONLY : numsubstances,   &
                       lakedatapar, &
                       lakedataparindex
 
    !Argument declarations
    INTEGER, INTENT(IN) :: i     !<index of current subbasin
    INTEGER, INTENT(IN) :: itype !<lake type (local or main)
    INTEGER, INTENT(IN) :: subid !<subid of current subbasin, for error output
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER status
    REAL fill   !flow from lakewi to slowlake

    IF(.NOT.ALLOCATED(lakestate%slowwater)) RETURN !no lake division, could change to if ns=0
    
    !> \b Algorithm \n
    !> If the lake is divided and the slow part has room for more water:
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)>0.)THEN
      IF(lakestate%water(itype,i)>slowlakeini(itype,i)-lakestate%slowwater(itype,i))THEN    !fill up slowlake
        !> If the fast part has enough water; fill up the slow part
        fill=slowlakeini(itype,i)-lakestate%slowwater(itype,i)
        lakestate%concslow(:,itype,i)=(lakestate%conc(:,itype,i)*fill+lakestate%concslow(:,itype,i)*lakestate%slowwater(itype,i))/slowlakeini(itype,i)
        lakestate%slowwater(itype,i)=slowlakeini(itype,i)
        CALL remove_water(lakestate%water(itype,i),numsubstances,lakestate%conc(:,itype,i),fill,lakestate%conc(:,itype,i),status)
        IF(status.NE.0) CALL error_remove_water(errstring(8),subid,i,itype)
      !> Else empty all water from fast part in the slow part
      ELSEIF(lakestate%water(itype,i)>0.AND.slowlakeini(itype,i)-lakestate%slowwater(itype,i)>0)THEN    !move all lakewi-water to slowlake
        CALL add_water(numsubstances,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
        lakestate%water(itype,i)=0.
        lakestate%conc(:,itype,i)=0.
      ENDIF  
    ENDIF

  END SUBROUTINE calculate_flow_within_lake

  !>\brief Flow from rating curve.
  !>Estimates average lake outflow (m3/s) during 1 timestep
  !!Uses rating equation q = k*(w-w0)**p
  !-----------------------------------------------------------------
  REAL FUNCTION average_flow_rating_curve(q_in,l_area,wst,k,p,w0)

    USE MODVAR, ONLY : seconds_per_timestep

    !Argument declarations
    REAL, INTENT(IN) :: q_in   !<inflow (m3/s)
    REAL, INTENT(IN) :: l_area !<lake area (m2)
    REAL, INTENT(IN) :: wst    !<current water level (m)
    REAL, INTENT(IN) :: k      !<rating curve coefficient
    REAL, INTENT(IN) :: p      !<rating curve exponent
    REAL, INTENT(IN) :: w0     !<rating curve threshold (m)
    
    !Local variables
    DOUBLE PRECISION w,h,h0,z,r,qut

    qut = 0.D0
    w = DBLE(wst)-DBLE(q_in)*DBLE(seconds_per_timestep)/DBLE(l_area) !Subtract inflow added in HYPE
    IF(w>DBLE(w0)) THEN
      h0 = w-DBLE(w0)                             !Height over threshold (m)
      r = DBLE(p)*DBLE(k)*(h0**(DBLE(p)-1.))/DBLE(l_area)            !Linearized recession rate (1/sec)
      IF(r>0.D0)THEN
        z = h0+DBLE(q_in)/r/DBLE(l_area)-h0/DBLE(p)   !Auxiliary variable (m)
        h = (h0-z)*EXP(-r*DBLE(seconds_per_timestep))+z  !New height above threshold (m)
        qut = DBLE(q_in)-DBLE(l_area)*(h-h0)/DBLE(seconds_per_timestep)
        IF(qut<0.D0) qut = 0.D0
      ENDIF
    ENDIF
    average_flow_rating_curve = REAL(qut)

  END FUNCTION average_flow_rating_curve

  !>Calculate outlet lake water stage (m) in local reference system and for w-reference system
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Lakes - Outlet lake (olake) as a lake basin)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_olake_waterstage(i,lakewatermm,lakeareain,lakearea,lakewst,lakestate,w0ref)

    USE MODVAR, ONLY : basin,           &
                       classbasin,      &
                       dam,             &
                       damindex,        &
                       lake,            &
                       lakeindex,       &
                       lakebasin,       &
                       lakebasinindex,  &
                       missing_value,   &
                       dayno,           &
                       nsub,            &
                       slc_olake

    !Arguments declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: lakewatermm   !<outlet lake water content (mm)
    REAL, INTENT(IN)    :: lakeareain    !<outlet lake area of subbasin (m2)
    REAL, INTENT(OUT)   :: lakearea      !<outlet lake area (of subbasin or whole lake for last lakebasin) (m2)
    REAL, INTENT(OUT)   :: lakewst       !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: w0ref         !<level to be added for w-ref outlet lake water stage (m)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER isb        !subbasin-loop index
    INTEGER dayno1, dayno2      !Starting day nr. for production periods 1 and 2    
    REAL lakewaterm    !lake water (m)
    REAL wmin          !lake threshold (sänkningsgräns)
    REAL deltaw0       !difference in lake threshold period 2 (m)
    REAL deltaw        !distance between thresholds (m)
    REAL isb_lakewst   !lake water above threshold of lake basins in lake (m)
    REAL lack          !water lacking for lake basin to reach threshold (m3)
    REAL overwmin      !water volume above wmin (to be spread over whole lake area)

    !> \b Algoritm \n
    !>Check for lake existance; return if not found    
    lakearea = lakeareain   !Default output lake area
    w0ref = 0.
    IF(lakeareain==0)THEN
      lakewst = missing_value
      RETURN
    ENDIF

    !Lake water level of current subbasin
    lakewaterm = lakewatermm * 0.001

    !Lake water reference
    IF(ALLOCATED(lakeindex))THEN
      IF(lakeindex(i)>0)THEN
        w0ref = lake(lakeindex(i))%w0ref
      ENDIF
    ENDIF
    IF(ALLOCATED(damindex))THEN
      IF(damindex(i)>0)THEN
        w0ref = dam(damindex(i))%w0ref
      ENDIF
    ENDIF

    !>Calculate lake water stage (for single olake) in w-reference system
    lakewst = lakewaterm - basin(i)%lakedepth(2)

    !Lakebasin lake
    IF(ALLOCATED(lakebasinindex))THEN
      IF(lakebasinindex(i)>0)THEN  
        IF(lakebasin(lakebasinindex(i))%last)THEN
          !>If outlet of lakebasin lake; waterstage is calculated as mean for whole lake
          !>\li Calculate lake water reference at threshold
          w0ref = lake(lakebasin(lakebasinindex(i))%ilk)%w0ref
          deltaw0 = lake(lakebasin(lakebasinindex(i))%ilk)%deltaw0  
          wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin
          dayno1 = lake(lakebasin(lakebasinindex(i))%ilk)%datum1          !Starting day nr. for production period 1
          dayno2 = lake(lakebasin(lakebasinindex(i))%ilk)%datum2          !                                period 2
          IF(wmin==missing_value)THEN
            deltaw = 0
            IF (dayno1*dayno2 > 0) THEN                                       !If both dates for different production regimes are non-zero,
              IF (dayno < dayno1 .OR. dayno >= dayno2)  deltaw = - deltaw0    !and today is not within the period 1, then threshold for period 2 applies
            ENDIF
          ELSE
            deltaw = 0. - wmin
          ENDIF
          !Check for lake basin with water stage below threshold and 
          !>\li Calculate lack of lake water in lakebasins of the lake
          lack = 0
          DO isb=1,nsub
            IF(lakebasinindex(isb)>0)THEN
              IF(lakebasin(lakebasinindex(isb))%ilk == lakebasin(lakebasinindex(i))%ilk)THEN
                isb_lakewst = lakestate%water(2,isb)
                IF(ALLOCATED(lakestate%slowwater)) isb_lakewst = isb_lakewst + lakestate%slowwater(2,isb)
                isb_lakewst = isb_lakewst * 0.001 - (basin(isb)%lakedepth(2) - deltaw)
                IF(isb_lakewst<0)THEN
                  lack = lack - isb_lakewst*classbasin(isb,slc_olake)%part*basin(isb)%area
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          !>\li Calculate average water stage and area for whole lake            
          overwmin = (lakewaterm - (basin(i)%lakedepth(2)-deltaw))*lakearea - lack
          IF(overwmin>0)THEN
            lakewst = overwmin/lake(lakebasin(lakebasinindex(i))%ilk)%area - deltaw   !W>wmin
          ELSE
            overwmin = (lakewaterm - (basin(i)%lakedepth(2)-deltaw))*lakearea
            IF(overwmin>0)THEN
              lakewst = 0. - deltaw   !W=wmin
            ELSE
              lakewst = overwmin/lakearea - deltaw  !W<wmin
            ENDIF
          ENDIF
          lakearea = lake(lakebasin(lakebasinindex(i))%ilk)%area  
        ELSE
          lakewst = lakewaterm - basin(i)%lakedepth(2)
          w0ref = 0.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_olake_waterstage

  !>Calculate outlet lake water stage (m) in local reference system adjusted for "real" regulation amplitude
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_regamp_adjusted_waterstage(i,lakeareain,lakewst,lakewstadj)

    USE MODVAR, ONLY : dam,             &
                       damindex,        &
                       lake,            &
                       lakeindex,       &
                       lakebasin,       &
                       lakebasinindex,  &
                       missing_value

    !Arguments declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: lakeareain    !<outlet lake area of subbasin (m2)
    REAL, INTENT(IN)    :: lakewst       !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: lakewstadj    !<outlet lake water stage adjusted for "real" amplitude of regulation volume (m)
    
    !Local variables
    REAL wfactor       !regulation amplitude scaling factor

    !> \b Algoritm \n
    wfactor = missing_value
    lakewstadj = missing_value
    
    !Check for lake existance; return if not found    
    IF(lakeareain==0) RETURN

    !>Get regulation amplitude adjustment factor
    !Case of single lake:
    IF(ALLOCATED(lakeindex))THEN
      IF(lakeindex(i)>0)THEN
        wfactor = lake(lakeindex(i))%wampcoeff
      ENDIF
    ENDIF
    IF(ALLOCATED(damindex))THEN
      IF(damindex(i)>0)THEN
        wfactor = dam(damindex(i))%wampcoeff
      ENDIF
    ENDIF

    !Case of lakebasin lake (last basin):
    IF(ALLOCATED(lakebasinindex))THEN
      IF(lakebasinindex(i)>0)THEN  
        IF(lakebasin(lakebasinindex(i))%last)THEN
          wfactor = lake(lakebasin(lakebasinindex(i))%ilk)%wampcoeff
        ENDIF
      ENDIF
    ENDIF

    !>Calculate adjusted lake water stage
    lakewstadj = lakewst
    IF(wfactor/=missing_value .AND. lakewst<0.) lakewstadj = lakewst*wfactor

  END SUBROUTINE calculate_regamp_adjusted_waterstage

  !>Calculate subbasin outlet flow division into main channel and branch
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Rivers - Main river)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_branched_flow(i,totflow,mainflow,branchflow)

    USE MODVAR, ONLY : branchdata,    &
         branchindex

    !Argument declaration
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: totflow       !<outflow of subbasin
    REAL, INTENT(OUT)   :: mainflow      !<flow in main channel
    REAL, INTENT(INOUT)   :: branchflow    !<flow in branch

    !Local variables
    REAL part, minQmain, maxQmain, maxQbranch
    
    !> \b Algorithm \n
    
    !>MM2016: for special dam branch; mainflow & branchflow calculated previously in calculate_outflow_from_lake
    IF(branchflow>0. .and. branchindex(i)/=0)THEN
      mainflow = totflow
      branchflow = branchflow
      RETURN
    ENDIF
    
    !>Initialisation, default is all flow in main (single) channel
    mainflow   = totflow
    branchflow = 0.

    !>Check for branch existance and flow>0
    IF(.NOT.ALLOCATED(branchdata)) RETURN
    IF(branchindex(i)==0) RETURN
    IF(totflow == 0) RETURN
   
    !>Set current parameter values
    part = branchdata(branchindex(i))%mainpart
    maxQmain = branchdata(branchindex(i))%maxQ
    minQmain = branchdata(branchindex(i))%minQ
    maxQbranch = branchdata(branchindex(i))%maxQbranch
    
    !>Calculate flow in main channel and in branch
    mainflow = totflow
    IF(totflow>minQmain)THEN
      mainflow = part * (totflow - minQmain) + minQmain
    ENDIF
    IF(maxQmain>0 .AND. mainflow>maxQmain)THEN
      mainflow = maxQmain
    ELSEIF(maxQbranch>0 .AND. (1.-part)*(totflow-minQmain)>maxQbranch)THEN
      mainflow = totflow - maxQbranch
    ENDIF
    branchflow = totflow - mainflow
    
  END SUBROUTINE calculate_branched_flow

  !>Calculate different volumes of lakes for print out 
  !>Volume for ilakes, volume for olakes, and volume for whole lakes (basindivided).
  !------------------------------------------------------------------
  SUBROUTINE calculate_lake_volume(itype,i,a,lakewi,lakebasinvol,lakevol,lakevolsum)
  
    USE MODVAR, ONLY : lakebasin, &
                       lakebasinindex, &
                       nbasinlakes, &
                       missing_value                   

    !Argument declarations
    INTEGER, INTENT(IN)   :: itype            !<lake type; ilake=1, olake=2
    INTEGER, INTENT(IN)   :: i                !<index of current subbasin
    REAL, INTENT(IN)      :: a                !<lake area (m2)
    REAL, INTENT(IN)      :: lakewi           !<lake water stage (mm)
    REAL, INTENT(INOUT)   :: lakebasinvol(2)  !<volume of olake and ilake
    REAL, INTENT(INOUT)   :: lakevol          !<volume of olake/volume for lake with basins in outlet basin
    REAL, INTENT(INOUT)   :: lakevolsum(nbasinlakes)    !<to sum lakebasins to outlet basin
    
    !Local variables
    INTEGER lakeid !lake id, used for lakes with many basins (lakeid i lakedata)

    !> \b Algoritm \n
    !> Calculate lake volume for current lake
    lakebasinvol(itype) = lakewi * 0.001 * a  !volume in lake (m3)

    !> If outlet lake
    IF(itype==2)THEN
      !\li Set volume for olake
      lakevol = lakebasinvol(itype)
      !\li For basin-lakes: calculate volume of outlet, upstream sub-lakebasins volumes set to missing
      IF(ALLOCATED(lakebasinindex))THEN    
        IF(lakebasinindex(i) .NE. 0) THEN   !lakebasin
          lakeid = lakebasin(lakebasinindex(i))%ilk
          lakevol = missing_value
          lakevolsum(lakeid) = lakevolsum(lakeid) + lakebasinvol(itype)
          IF(lakebasin(lakebasinindex(i))%last) THEN  !outlet of lakebasin-lake
            lakevol = lakevolsum(lakeid)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_lake_volume

  
  !>Calculate temperature(T2) processes in rivers
  !----------------------------------------------------------
   SUBROUTINE T2_processes_in_river(i,itype,temp,swrad,riversurft,riverarea,frozenstate,riverstate,freezeupday,freezeuparea)
  
    USE MODVAR, ONLY: genpar, &
                      i_t2, &
                      modeloption, &
                      p_lakeriverice, &
                      cwater
    USE HYPEVARIABLES, ONLY: m_t2trriver, &
                             m_riceTf, &
                             m_tcfriver, &
                             m_scfriver, &
                             m_ccfriver, &
                             m_lcfriver, &
                             m_stbcorr1, &
                             m_stbcorr2, &
                             m_stbcorr3
    
    !Argument variables
    INTEGER, INTENT(IN) :: i               !<index of subbasin
    INTEGER, INTENT(IN) :: itype           !<index of river type (local = 1, main = 2)
    REAL, INTENT(IN)    :: temp            !<air temperature
    REAL, INTENT(IN)    :: swrad           !<solar radiation
    REAL, INTENT(INOUT) :: riversurft(2)   !<water surface temperature
    REAL, INTENT(IN)    :: riverarea       !<river area
    TYPE(snowicestatetype),INTENT(IN)  :: frozenstate   !<Snow and ice states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER,INTENT(OUT) :: freezeupday(2)
    REAL, INTENT(INOUT) :: freezeuparea     !<fraction of riverarea with newice formation
    
    !Local variables    
    REAL t2transfcorr
    REAL watertemp, watervol,icefreefraction

    !Initiate heat deficit and freezeup flag and surface temp variables
    freezeuparea = 0.
    freezeupday(itype) = 0
    riversurft(itype) = 0.  

    !Get total river water volume and mean T2 temperature
    CALL get_rivertempvol(i,itype,riverstate,watertemp,watervol)
    watervol = watervol * 1000. / riverarea    !scale volume [m3] to depth [mm]

    IF(watervol.GT.0.)THEN    !Skip calculations if there is no water in the river

      icefreefraction = 1. - frozenstate%rivericecov(itype,i)      !Fraction of icefree river surface area
      t2transfcorr = 1.      !Seasonal correction of T2 exchange coefficient   

      !River-Atmosphere T2 exchange, only in ice-free conditions and if there is some water in the river
      IF(icefreefraction.GT.0.)THEN    
        !River-atmosphere exchange 
        ! optional models  (will be reduced to one option after some initial testing for EHYPE3.0 and SHYPE2012)
        SELECT CASE(modeloption(p_lakeriverice))
        CASE(2) ! new model based on Piccolroaz et al 2013, with modifications for fractional ice cover, and calculation of fractional freezup area
          CALL calculate_watersurface_heatbalance(temp,swrad,watertemp,watervol*riverarea*0.001,riverarea*icefreefraction, & 
                                                  genpar(m_tcfriver),genpar(m_scfriver),genpar(m_ccfriver),genpar(m_lcfriver), &
                                                  freezeuparea,genpar(m_riceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))
        CASE(1) ! the simple air-water temperature exchange model (Johan/David), with modifications for fractional ice cover, and calculation of fractional freezup area
          CALL calculate_T2_transfer(temp,watertemp,watervol*riverarea*0.001,riverarea*icefreefraction,genpar(m_t2trriver)*t2transfcorr, & 
                                     freezeuparea,genpar(m_riceTf))
        ENDSELECT
        
        !Check the freezeup conditions
        IF(freezeuparea.GT.0.)THEN
          !freezup area is the fraction of previously unfrozen area (riverarea*icefreefraction), where new ice formation is triggered
          !re-scale to a fraction of the entire river area:
          freezeuparea = freezeuparea * icefreefraction
          freezeupday(itype) = 1
        ENDIF
       
        !Assign update values to the riverstate variables
        CALL set_rivertemp(i,itype,riverstate,watertemp)
      
        !Assign river surface temperature if (partly) icefree conditions - it's later rescaled after ice calculations
        riversurft(itype) = riverstate%conc(i_t2,itype,i)
      ENDIF
    ELSE
      !Set T2 temperature to 0. if there is no water
      CALL set_rivertemp(i,itype,riverstate,0.)
    ENDIF
    
  END SUBROUTINE T2_processes_in_river

  !>Calculate ice processes in rivers
  !----------------------------------------------------------
  SUBROUTINE ice_processes_in_river(i,itype,iluse,snowfall,temp,riversurftemp,  &
                  riversnowdepth,riverarea,swrad,frozenstate,riverstate, &
                  freezeupday,breakupday,freezeuparea)
     
    USE MODVAR, ONLY: genpar, &
                      Lfreezing
    USE HYPEVARIABLES, ONLY: m_sndens0,&
                             m_ricesndens, &
                             m_ricetf, &
                             m_ricekika, &
                             m_ricekexp, &
                             m_ricetmelt,  &
                             m_ricewme,  &
                             m_riceTf

    !Argument declaration    
    INTEGER, INTENT(IN) :: i                 !<index of subbasin
    INTEGER, INTENT(IN) :: itype             !<index of lake/river type
    INTEGER, INTENT(IN) :: iluse             !<index of landuse
    REAL,INTENT(IN)     :: snowfall          !<snowfall
    REAL,INTENT(IN)     :: temp              !<air temperature
    REAL,INTENT(INOUT)  :: riversurftemp(2)  !<water surface temperature
    REAL,INTENT(INOUT)  :: riversnowdepth(2) !<snowdepth on lakeice, cm
    REAL,INTENT(IN)     :: riverarea         !<river area
    REAL,INTENT(IN)     :: swrad             !<shortwave radiation
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER, INTENT(IN) :: freezeupday(2)    !<status freeze-up day
    INTEGER, INTENT(INOUT) :: breakupday(2)  !<status ice break-up day
    REAL, INTENT(IN)     :: freezeuparea     !<fraction of river area with newice formation (calculated by surface heat balance function)
    
    !Local variables
    REAL driverwidt, dsnowdt
    REAL oldsnow, melt
    REAL newicesurftemp,newice,newbice,newicesnow,newicesnowdepth,oldsurftemp
    INTEGER newbreakup

    !Local parameters
    REAL, PARAMETER :: L = 3.35E5     !latent heat of freezing, J/kg
    REAL, PARAMETER :: dice = 0.917    !density of ice, fraction of water
    REAL, PARAMETER :: mm2cm = 0.1
    REAL, PARAMETER :: cm2mm = 10.
    
    !RiverIceModel: Initialization of some variables
    breakupday(itype)  = 0
    riversnowdepth(itype) = 0.0
    newicesurftemp = 0.
    newice = 0.
    newbice = 0.
    newicesnow = 0.
    newicesnowdepth = 0.
    newbreakup = 0
    
    !New ice formation on "freezeuparea" (calculated by surface heat balance function)
    IF(freezeuparea>0.)THEN
      CALL calculate_icedepth(newicesurftemp, newice, &
                              newbice,newicesnow,newicesnowdepth, & 
                              temp,driverwidt,dsnowdt,freezeupday(itype),newbreakup, &
                              genpar(m_ricetf),genpar(m_ricekika),genpar(m_ricekexp),genpar(m_ricetmelt))
    ENDIF    

    !Calculate development of the old river ice
    IF(frozenstate%riverice(itype,i)>0.)THEN
       
      !first guess is that the old ice (or snow) is melting at 0 degrees
      oldsurftemp = 0.
      
      !Snow on riverice calculation
      oldsnow = frozenstate%riversnow(itype,i)
      CALL calculate_snow_on_ice(iluse,snowfall,frozenstate%riversnow(itype,i),temp,  &
                                 melt,swrad,frozenstate%riversnowage(itype,i))
         
      !Update snow age and snow depth for snow on ice
      CALL calculate_snowdepth(frozenstate%riversnow(itype,i),oldsnow,genpar(m_sndens0),  &
                               genpar(m_ricesndens),riversnowdepth(itype),frozenstate%riversnowage(itype,i))

      !Ice depth calculation (incl. update of skin temperature)
      CALL calculate_icedepth(oldsurftemp, frozenstate%riverice(itype,i), &
                              frozenstate%riverbice(itype,i),frozenstate%riversnow(itype,i),riversnowdepth(itype), & 
                              temp,driverwidt,dsnowdt,freezeupday(itype),breakupday(itype), &
                              genpar(m_ricetf),genpar(m_ricekika),genpar(m_ricekexp),genpar(m_ricetmelt))

      !If river temperature is above freezing, use the excess heat to melt some river ice from below (see further in corresponding lake routine)
      CALL riverice_riverwater_interaction(i,itype,riverstate,frozenstate,riverarea,breakupday,driverwidt)
    ENDIF

    !Add new ice to the old ice
    IF(newice>0.)THEN
      IF(frozenstate%riverice(itype,i)>0.)THEN
        frozenstate%riversnow(itype,i) = frozenstate%riversnow(itype,i)* frozenstate%rivericecov(itype,i)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        riversnowdepth(itype) = riversnowdepth(itype) * frozenstate%rivericecov(itype,i)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        frozenstate%riverice(itype,i) = (frozenstate%riverice(itype,i)*frozenstate%rivericecov(itype,i) + newice*freezeuparea)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        frozenstate%riverbice(itype,i) = (frozenstate%riverbice(itype,i)*frozenstate%rivericecov(itype,i) + newbice*freezeuparea)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        riversurftemp(itype) = newicesurftemp * freezeuparea + oldsurftemp * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-freezeuparea-frozenstate%rivericecov(itype,i))
        frozenstate%rivericecov(itype,i) = (frozenstate%rivericecov(itype,i)+freezeuparea)      
      ELSE
        frozenstate%riversnow(itype,i) = 0.0
        frozenstate%riversnowage(itype,i) = 0.0
        riversnowdepth(itype) = 0.0
        frozenstate%riverice(itype,i) = newice
        frozenstate%riverbice(itype,i) = newbice
        riversurftemp(itype) = newicesurftemp * freezeuparea + riversurftemp(itype)*(1.-freezeuparea)
        frozenstate%rivericecov(itype,i) = freezeuparea
        !Make sure breakupday is 0 (strange situation with complete meltout of old ice and newice formation at the same time)
        IF(breakupday(itype)==1) breakupday(itype)=0
      ENDIF
    ELSE
      !Only old ice remaining
      IF(frozenstate%riverice(itype,i).GT.0.)THEN
        !weighted surface temperature (oldice and open water surface temperature)
        riversurftemp(itype) = oldsurftemp * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-frozenstate%rivericecov(itype,i))
      ELSE
        !no new snow and no old snow
        !check if there was complete meltout today, in that case make sure all variables are reset
        IF(breakupday(itype).EQ.1)THEN
          frozenstate%riverice(itype,i) = 0.
          frozenstate%riverbice(itype,i) = 0.
          frozenstate%riversnow(itype,i) = 0.
          frozenstate%riversnowage(itype,i) = 0.
          riversnowdepth(itype) = 0.0
          riversurftemp(itype) = genpar(m_riceTf) * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-frozenstate%rivericecov(itype,i))
          frozenstate%rivericecov(itype,i) = 0.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE ice_processes_in_river
  
  !>Calculate interaction between river water and river ice 
  ! - heat from water temperature above freezing is used to melt river ice 
  !   by reducing the fractional area, rather than reducing ice depth
  ! - latent heat correspondning to ice meltwater is also added to the water
  !--------------------------------------------------------------------------
  SUBROUTINE riverice_riverwater_interaction(i, itype, riverstate, frozenstate, riverarea, breakupday, driverwidt)

    USE MODVAR, ONLY : genpar, cwater
    USE HYPEVARIABLES, ONLY : m_riceTf,m_ricewme
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<current subbasin index
    INTEGER, INTENT(IN) :: itype    !<river type
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate      !<River state
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
    REAL, INTENT(IN) :: riverarea   !<river area (m2)
    INTEGER, INTENT(INOUT) :: breakupday(2) !<status of river ice break up
    REAL, INTENT(IN) :: driverwidt
    
    
    !local variables
    REAL watertemp, watervol, icewater, meltheat, waterheat, meltwater, newwatertemp,oldicecover
    
    !parameters
    REAL, PARAMETER :: L = 3.35E5     !latent heat of freezing, J/kg
    REAL, PARAMETER :: dice = .917    !density of ice, fraction of water
    REAL, PARAMETER :: mm2cm = 0.1
    REAL, PARAMETER :: cm2mm = 10.
  
    !Get total river water volume and mean T2 temperature
    CALL get_rivertempvol(i,itype,riverstate,watertemp,watervol)
    
    !scale volume [m3] to depth [mm], volume water per unit area
    watervol = watervol * 1000. / (riverarea)
    
    oldicecover = frozenstate%rivericecov(itype,i)
    
    IF(watervol.GT.0.)THEN
    
      !available heat for melting (C * KG/M2 * 1000 * KJ/KG/C = J/M2), per unit area
      waterheat = (watertemp-genpar(m_riceTf)) * watervol * 1000. * cwater
      
      IF(waterheat.GT.0.)THEN
        !Try to melt some ice, but only if the ice did not already melt completely (breakupday==1)
        IF(breakupday(itype).EQ.0)THEN
          ! !melt the ice, from below, in cm ice
          ! bottommelt = min(frozenstate%riverice(itype,i),waterheat/(L*dice)*mm2cm)
          ! meltheat   = bottommelt * (L*dice) * cm2mm
          ! meltwater = bottommelt * dice *cm2mm
          
          !river ice and snow mass, in mm water, per unit area of ice covered river
          icewater = frozenstate%riverice(itype,i)*dice*cm2mm + frozenstate%riversnow(itype,i)
          
          !ice melt, in mm per unit area of ice-covered river 
          ! - it is thus unly the water below the ice which is interacting with the ice
          ! - the available heat is scaled with the "meltefficiency" parameter 
          meltwater = min(icewater,genpar(m_ricewme)*waterheat/L)
          meltheat = meltwater * L
             
! 3) update the frozen states with bottom melt
          
          !frozenstate%riverice(itype,i)=max(0.,frozenstate%riverice(itype,i)-bottommelt)
          !IF(frozenstate%riverice(itype,i).GT.0.)THEN
          IF((icewater-meltwater).GT.0.)THEN
            !some ice remains, reduce ice content by reducing the fractional area
            frozenstate%rivericecov(itype,i) = min(1.,max(0.,frozenstate%rivericecov(itype,i)*(1.- meltwater/icewater)))           
!            frozenstate%riverbice(itype,i)=max(0.,frozenstate%riverbice(itype,i)-bottommelt)
          ELSE
            !complete melt of the riverice
            frozenstate%riverice(itype,i) =0.
            frozenstate%riverbice(itype,i)=0.
            
            !add heat needed to melt the riversnow to the meltheat
            !meltheat = meltheat + frozenstate%riversnow(itype,i) * L
              
            !add snow to the meltwater
            !meltwater = meltwater + frozenstate%riversnow(itype,i)
              
            !reset the snow states
            frozenstate%riversnow(itype,i)=0.
            frozenstate%riversnowage(itype,i)=0.
            
            !and the ice cover state
            frozenstate%rivericecov(itype,i)=0.
            
            !set breakup flag to 1
            breakupday(itype) = 1
          ENDIF
        ELSE
          !Ice was already melted away by the icedepth function
          meltheat   = 0.
          meltwater  = 0.
        ENDIF
      ELSE
        meltheat = 0.
        meltwater = 0.
      ENDIF
! 4) use any remaining heat and the zero degree melt water to update the river state
        
      !remove melt heat from heat content of the lake water (this is now per unit area previously ice covered river)
      waterheat = waterheat - meltheat
      
      !add any previous surface melt water to the meltwater
      IF(driverwidt.GT.0)THEN
        meltwater = meltwater + driverwidt
      ENDIF
      
      !temperature of water from remaining heat content
      newwatertemp=max(waterheat/(watervol * 1000. * cwater) + genpar(m_riceTf),genpar(m_riceTf))
        
      !dilute with the meltwater, which is at freezing point
      newwatertemp=max(genpar(m_ricetf),watertemp * (watervol - meltwater)/watervol)
       
      !weighted temperature between (previously) ice covered and ice free water
      watertemp = watertemp * (1.-oldicecover) + newwatertemp * oldicecover
       
      !finally, assign update values to the riverstate variables
      CALL set_rivertemp(i,itype,riverstate,watertemp)
   ENDIF
  
  END SUBROUTINE riverice_riverwater_interaction
  
  !>Subroutine for calculation of snow on ice changes; snowfall addition and
  !>snow pack melting
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snow_on_ice(iluse,snowfall,snow,temp,melt,swrad,snowage)
  
    USE MODVAR, ONLY : landpar,   &
                       missing_value,   &
                       modeloption,     &
                       p_snowmelt,      &
                       genpar
    USE HYPEVARIABLES, ONLY : m_ttmp,m_cmlt,m_snalbmin,m_snalbmax,m_snalbkexp,m_cmrad,m_licewcorr,m_cmltcorr

    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: snowfall !<precipitation as snow (mm/timestep) 
    REAL, INTENT(INOUT) :: snow     !<snow pack (mm)
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
    REAL, INTENT(OUT)   :: melt     !<snow melt (mm/timestep)
    REAL, INTENT(IN)    :: swrad    !<shortwave radiation (MJ/m2/day?)
    REAL, INTENT(INOUT) :: snowage  !<snowage (timesteps)

    !Local variables
    REAL tt       !threshold temperature for snow melt (and evaporation) (C)
    REAL cm       !coefficient for snow melt (mm/C/timestep)
    REAL newsnow
    REAL snowalbedo
    REAL snalbmax, snalbmin, snalbkexp
    REAL cmrad     ! radiation index snow melt factor (mm/MJ/timestep)
    REAL snowcover
    !Set parameter values
    tt = landpar(m_ttmp,iluse) !threshold temperature for snow melt
    cm = landpar(m_cmlt,iluse)*genpar(m_cmltcorr) !Coefficient for snow melt
    
    ! just set snowcover = 1., and introduce snowcover calculation on lake ice later...
    snowcover = 1.
    
    !Melting, select model            
    SELECT CASE(modeloption(p_snowmelt))
    CASE(0) ! Original Temperature Index model, WITHOUT snowcover scaling
       IF(snow>0 .AND. temp >= tt) THEN
          melt = cm * (temp - tt)
          melt = MIN(melt, snow)
       ELSE
          melt = 0.
       ENDIF
    CASE(1) ! Original Temperature Index model, WITH snowcover scaling
        IF(snow>0 .AND. temp >= tt) THEN
          melt = cm * (temp - tt)
    
          !Scale melt with fractional snow cover
          melt = melt * snowcover
       
          melt = MIN(melt, snow)
       ELSE
          melt = 0.
       ENDIF
    CASE(2) ! Temperature AND Radiation index model, WITH snowcover scaling
       !Set parameter values
       snalbmin  = landpar(m_snalbmin,iluse)
       snalbmax  = landpar(m_snalbmax,iluse)
       snalbkexp = landpar(m_snalbkexp,iluse)
       cmrad     = landpar(m_cmrad,iluse)
       
       !Get snow albedo           
       snowalbedo = snowalbedo_function(snowage,snalbmin,snalbmax,snalbkexp)
       
       !Radiation component
       melt = cmrad * swrad * (1.-snowalbedo)
       
       !Add Temperature component
       IF(snow>0. .AND. temp >= tt)THEN
          melt = melt + cm * (temp - tt)
       ENDIF
       
       !Subtract refreezing component when temperatures below tt, parameter 0.5 hard coded for the moment
       IF(snow>0. .AND. temp < tt .AND. melt > 0.)THEN
          melt = melt - 0.5 * cm * (tt - temp)
       ENDIF
       
      !Scale melt with fractional snow cover
       melt = melt * snowcover

      !Minimize melt to available snow, and make sure it's positive
       melt = MAX(0.,MIN(melt, snow))
       
    CASE DEFAULT ! Original temperature index model
       IF(snow>0 .AND. temp >= tt) THEN
          melt = cm * (temp - tt)
          melt = MIN(melt, snow)
       ELSE
          melt = 0.
       ENDIF
    END SELECT
    

    !Update the snow with snowfall and melting
    newsnow = max(0.,snow + genpar(m_licewcorr)*snowfall  - melt)

    snow = newsnow

  END SUBROUTINE calculate_snow_on_ice

  !>Calculate lake ice processes
  !----------------------------------------------------------
  SUBROUTINE ice_processes_in_lake(i,itype,iluse,snowfall,temp,lakesurftemp,  &
                                   lakesnowdepth,swrad,frozenstate, &
                                   lakestate,freezeupday,breakupday,lakearea,freezeuparea)
    
    USE MODVAR, ONLY: genpar
    USE HYPEVARIABLES, ONLY: m_sndens0, &
                             m_licesndens, &
                             m_licetf,   &
                             m_licekika, &
                             m_licekexp, &
                             m_licetmelt,  &
                             m_licewme,  &
                             m_liceTf
           
    !Argument declarations
    INTEGER, INTENT(IN) :: i                !<index of subbasin
    INTEGER, INTENT(IN) :: itype            !<index of lake/river type
    INTEGER, INTENT(IN) :: iluse            !<index of landuse
    REAL,INTENT(IN)     :: snowfall         !<snowfall
    REAL,INTENT(IN)     :: temp             !<air temp
    REAL,INTENT(INOUT)  :: lakesurftemp(2)  !<water surface temperature
    REAL,INTENT(OUT)    :: lakesnowdepth(2) !<snowdepth on lakeice, cm
    REAL,INTENT(IN)     :: swrad            !<shortwave radiation
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    INTEGER, INTENT(IN) :: freezeupday(2)   !<status freeze-up day
    INTEGER, INTENT(INOUT) :: breakupday(2)   !<status ice break-up day
    REAL, INTENT(IN)    :: lakearea
    REAL, INTENT(IN)    :: freezeuparea    !<fractional water surface area with newice formation, given by temperature routine
    
    !Local variables
    REAL dlakewidt, dsnowdt
    REAL oldsnow, melt
    REAL newicesurftemp,newice,newbice,newicesnow,newicesnowdepth,oldsurftemp
    INTEGER newbreakupday
   
    !Initialization
    breakupday(itype)  = 0
    lakesnowdepth(itype) = 0.0
    newicesurftemp = 0.
    newice = 0.
    newbice = 0.
    newicesnow = 0.
    newicesnowdepth = 0.
    newbreakupday=0
    
    !Newice formation on "freezeuparea" (calculated by surface heat balance function)
    IF(freezeuparea.GT.0.)THEN
      CALL calculate_icedepth(newicesurftemp, newice, &
                              newbice,newicesnow,newicesnowdepth, & 
                              temp,dlakewidt,dsnowdt,freezeupday(itype),newbreakupday, &
                              genpar(m_licetf),genpar(m_licekika),genpar(m_licekexp),genpar(m_licetmelt))
    ENDIF

    !Calculate development of the old river ice
    IF(frozenstate%lakeice(itype,i).GT.0)THEN
       !FROZEN LAKE

       !first guess is that the ice (or snow) is melting at 0 degrees
       oldsurftemp = 0.0
       
       !snow on lakeice calculation
       oldsnow = frozenstate%lakesnow(itype,i)
       CALL calculate_snow_on_ice(iluse,snowfall,frozenstate%lakesnow(itype,i),temp,melt,swrad,frozenstate%lakesnowage(itype,i))
                  
       !Update snow age and snow depth for snow on ice
       CALL calculate_snowdepth(frozenstate%lakesnow(itype,i),oldsnow,genpar(m_sndens0),genpar(m_licesndens),lakesnowdepth(itype),frozenstate%lakesnowage(itype,i))
      
       !Ice depth calculation (inlc. update of skin temperature)
       CALL calculate_icedepth(oldsurftemp,frozenstate%lakeice(itype,i),frozenstate%lakebice(itype,i),frozenstate%lakesnow(itype,i),lakesnowdepth(itype), & 
                               temp,dlakewidt,dsnowdt,freezeupday(itype),breakupday(itype), &
                               genpar(m_licetf),genpar(m_licekika),genpar(m_licekexp),genpar(m_licetmelt))
       
       !Calculate bottom melt due to heat from lake water temperatures above freezing, as well as influence of surface melt on lake water temperature
       CALL calculate_lakeice_lakewater_interaction(itype,i,lakearea,frozenstate,lakestate,dlakewidt,breakupday)
       
    ENDIF
    
    !Add new ice to the old ice
    IF(newice.GT.0.)THEN
      IF(frozenstate%lakeice(itype,i).GT.0.)THEN
         frozenstate%lakesnow(itype,i) = frozenstate%lakesnow(itype,i)* frozenstate%lakeicecov(itype,i)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         lakesnowdepth(itype) = lakesnowdepth(itype) * frozenstate%lakeicecov(itype,i)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         frozenstate%lakeice(itype,i) = (frozenstate%lakeice(itype,i)*frozenstate%lakeicecov(itype,i) + newice*freezeuparea)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         frozenstate%lakebice(itype,i) = (frozenstate%lakebice(itype,i)*frozenstate%lakeicecov(itype,i) + newbice*freezeuparea)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         lakesurftemp(itype) = newicesurftemp * freezeuparea + oldsurftemp * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1. - freezeuparea - frozenstate%lakeicecov(itype,i))
         frozenstate%lakeicecov(itype,i) = (frozenstate%lakeicecov(itype,i)+freezeuparea)  
      ELSE
         frozenstate%lakesnow(itype,i) = 0.0
         frozenstate%lakesnowage(itype,i) = 0.0
         lakesnowdepth(itype) = 0.0
         frozenstate%lakeice(itype,i) = newice
         frozenstate%lakebice(itype,i) = newbice
         lakesurftemp(itype) = newicesurftemp * freezeuparea + lakesurftemp(itype)*(1. - freezeuparea)
         frozenstate%lakeicecov(itype,i) = freezeuparea
         !Make sure breakupday is 0 (strange situation with complete meltout of old ice and newice formation at the same time)
         IF(breakupday(itype)==1) breakupday(itype)=0
      ENDIF
    ELSE
      !Or just check breakup conditions of old ice, and/or update the lakesurf temperature
      IF(frozenstate%lakeice(itype,i).GT.0.)THEN
        lakesurftemp(itype) = oldsurftemp * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1. - frozenstate%lakeicecov(itype,i))
      ELSE
        !no new snow and no old snow
        !check if there was complete meltout today, in that case make sure all variables are reset
        IF(breakupday(itype).EQ.1)THEN
          frozenstate%lakeice(itype,i) = 0.
          frozenstate%lakebice(itype,i) = 0.
          frozenstate%lakesnow(itype,i) = 0.
          frozenstate%lakesnowage(itype,i) = 0.
          lakesnowdepth(itype) = 0.0
          lakesurftemp(itype) = genpar(m_liceTf) * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1.-frozenstate%lakeicecov(itype,i))
          frozenstate%lakeicecov(itype,i) = 0.
        ENDIF
      ENDIF   
    ENDIF

  END SUBROUTINE ice_processes_in_lake
  
  !>Calculate lake ice melt from heat from lake water, as well as influence of ice surface melt on lake water temperature
  ! - depending on lake type (fast and slow split or not, deep or shallow), a mean water temperature and water volume is
  !   calculated for the interaction with the lake ice. The resulting watertemperature is then assigned to the 
  !   various lake water components
  ! - heat from water temperature above freezing is used to melt lake ice 
  !   by reducing the fractional area, rather than reducing ice depth
  ! - latent heat correspondning to ice meltwater is also added to the water
  SUBROUTINE calculate_lakeice_lakewater_interaction(itype,i,lakearea,frozenstate,lakestate,dlakewidt,breakupday)

    USE MODVAR, ONLY: genpar, lakedataparindex,lakedatapar,i_t2,cwater
    USE HYPEVARIABLES, ONLY: m_lddeeplake, m_ldfastlake, m_licetf, m_licewme

    !Argument declarations
    INTEGER,INTENT(IN) :: i             !<index of subbasin
    INTEGER,INTENT(IN) :: itype         !<index of lake type (ilake = 1, olake = 2)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    REAL, INTENT(IN)   :: lakearea
    REAL, INTENT(IN)   :: dlakewidt
    INTEGER, INTENT(OUT) :: breakupday(2)
    
    !Local variables    
    integer laketype
!    real meantemp,meanwater,waterheat,meltheat,watertemp,watervol,bottommelt,meltwater,newwatertemp
    real meantemp,meanwater,waterheat,meltheat,watertemp,watervol,icewater,meltwater,newwatertemp,oldicecov
    REAL epidepth, fastoutpart

    !parameters
    real, parameter :: L = 3.35E5     ! latent heat of freezing, J/kg
    real, parameter :: dice = .917    ! density of ice, fraction of water
    real, parameter :: mm2cm = 0.1
    real, parameter :: cm2mm = 10.
    
!--------------------------------------------------------------------------------------
! lakewater-lakeice interaction:
!
! 1) find out how much water and at what temperature we have for melting ice from below
! 2) melt corresponding ice (from below: black ice, slush ice, snow)
! 3) update the frozen states
! 4) use remaining heat and heat from melt water to update the lake state
!
! the first and last step is complicated by the various lake water storage configurations
!--------------------------------------------------------------------------------------    
    oldicecov = frozenstate%lakeicecov(itype,i)
    
!1)find out how much water and at what temperature we have for melting ice from below
    
    !Depth to thermocline, function of lake area (REF)
    epidepth = 6.95 * (lakearea / 1000000)**0.185
    
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN
       !--------------------------------------
       !Lake model WITHOUT fast and slow split
       !--------------------------------------
       
       !Calculate total water stage (WATER+SLOWWATER) and average temperature (CONC*WATER+CONCSLOW*SLOWWATER)/(WATER+SLOWWATER)
       meanwater = lakestate%water(itype,i)+lakestate%slowwater(itype,i)
       IF(meanwater.GT.0.)THEN
         meantemp  = (lakestate%water(itype,i)*lakestate%conc(i_t2,itype,i)+lakestate%slowwater(itype,i)*lakestate%concslow(i_t2,itype,i))/meanwater
       
         !Check lake depth, if thermal stratification
         IF(epidepth < meanwater*0.001 .AND. itype==2)THEN !why is this only possible for olakes?
           !Two-layer olake, waterdepth > thermocline, olake
    
           !->derive lake uppertemp(t) from meantemp(t, preliminary) and lowertemp(t-1)
           lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth

           !temperature and water volume interacting with the ice
           watertemp = lakestate%uppertemp(itype,i)
           watervol  = epidepth*1000.
           laketype  = 2
         ELSE
           !one-layer olake or ilakes
           watertemp = meantemp
           watervol  = meanwater
           laketype  = 1
         ENDIF
       ELSE
         !no water in the lake, do nothing
         laketype = 0
       ENDIF
    ELSE
       !--------------------------------------
       !Lake model WITH fast and slow split
       !--------------------------------------
       fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
       
       !areally weighted averaged water stage
       meanwater = lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart
       
       IF(meanwater.GT.0.)THEN
         !calculate weighted average mean temperature
         meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
       
         IF(epidepth < meanwater*0.001 .AND. itype==2)THEN
           !Deep lake with thermal stratification (assume they have a common lowertemp)
           
           !calculate upper temp from meantemp and the lowertemp
           lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
           
           !temperature and water volume interacting with the ice
           watertemp = lakestate%uppertemp(itype,i)
           watervol  = epidepth*1000.
           laketype  = 4
         ELSE
           !One-layer lake
           watertemp = meantemp
           watervol  = meanwater
           laketype  = 3
         ENDIF
       ELSE
         !no water in the lake - do nothing
         laketype = 0
       ENDIF        
    ENDIF
    
! 2) melt corresponding ice (from below: black ice, slush ice, snow), takin fractional ice cover into account
    !watervol = watervol / frozenstate%lakeicecov(itype,i) 
    IF(laketype.GT.0)THEN
      !available heat for melting (C * KG/M2 * 1000 * KJ/KG/C = J/M2)
      waterheat = (watertemp-genpar(m_liceTf)) * watervol * 1000. * cwater 
       
      IF(waterheat.GT.0.)THEN
      
        !Try bottom melt only if there was not already complete meltout (breakupday==1)
        IF(breakupday(itype).EQ.0)THEN
          !!melt the ice from below, in cm ice
          !bottommelt = min(frozenstate%lakeice(itype,i),waterheat/(L*dice)*mm2cm)
          !meltheat   = bottommelt * (L*dice) * cm2mm
          !meltwater  = bottommelt * dice *cm2mm
          
          !lake ice and snow mass, in mm water, per unit area of ice covered lake
          icewater = frozenstate%lakeice(itype,i)*dice*cm2mm + frozenstate%lakesnow(itype,i)
          
          !ice melt, in mm per unit area of ice-covered river
          ! - it is thus unly the water below the ice which is interacting with the ice
          ! - the available heat is scaled with a "Meltefficiency" parameter
          meltwater = min(icewater,genpar(m_licewme)*waterheat/L)
          meltheat = meltwater * L
        
! 3) update the frozen states
!          frozenstate%lakeice(itype,i)=max(0.,frozenstate%lakeice(itype,i)-bottommelt)
        
!          IF(frozenstate%lakeice(itype,i).GT.0.)THEN
          IF((icewater-meltwater).GT.0.)THEN
            !some ice remains, redice icemass by reducing fractional coverage
            frozenstate%lakeicecov(itype,i) = min(1.,max(0.,frozenstate%lakeicecov(itype,i)*(1-meltwater/icewater)))
!            frozenstate%lakebice(itype,i)=max(0.,frozenstate%lakebice(itype,i)-bottommelt)
          ELSE
            !complete melt of the lakeice
            frozenstate%lakeice(itype,i)=0.
            frozenstate%lakebice(itype,i)=0.
        
            !add heat needed to melt the lakesnow to the meltheat
            !meltheat = meltheat + frozenstate%lakesnow(itype,i) * L
            
            !add snow to the meltwater
            !meltwater = meltwater + frozenstate%lakesnow(itype,i)
            
            !reset the snow states
            frozenstate%lakesnow(itype,i)=0.
            frozenstate%lakesnowage(itype,i)=0.
            
            !and ice cover area
            frozenstate%lakeicecov(itype,i) = 0.
            
            !set breakupflag to 1
            breakupday(itype) = 1
          ENDIF
        ELSE
          meltheat = 0.
          meltwater = 0.
        ENDIF
      ELSE
        meltheat = 0.
        meltwater = 0.
      ENDIF
! 4) use any remaining heat and the zero degree melt water to update the lake state
      
      !remove melt heat from heat content of the lake water
      waterheat = waterheat - meltheat
      
      !add any previous surface melt water to the meltwater
      IF(dlakewidt.GT.0)THEN
          meltwater = meltwater + dlakewidt
      ENDIF
      
      !temperature of water from remaining heat content
      newwatertemp=max(waterheat/(watervol * 1000. * cwater) + genpar(m_liceTf),genpar(m_liceTf))
      
      !dilute with the meltwater, which is at freezing point
      newwatertemp = max(genpar(m_licetf),watertemp * (watervol - meltwater)/watervol)
      
      !weighted temperature, between icefree and icecovered water
      watertemp = oldicecov * newwatertemp + (1.-oldicecov)*watertemp
      
      !finally, assign update values to the real state variable
      SELECT CASE(laketype)
      
        CASE(1) !single layer without split
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = watertemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = watertemp
          lakestate%uppertemp(itype,i) = watertemp
          lakestate%lowertemp(itype,i) = watertemp
        
        CASE(2) !two-layer without split
          lakestate%uppertemp(itype,i) = watertemp
          meantemp = (lakestate%uppertemp(itype,i) * epidepth + (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i) )/ (meanwater * 0.001)
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = meantemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = meantemp
          
        CASE(3) !single layer with split
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = watertemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = watertemp
          lakestate%uppertemp(itype,i) = watertemp
          lakestate%lowertemp(itype,i) = watertemp
        
        CASE(4) !two-layer with split
          lakestate%uppertemp(itype,i) = watertemp
          meantemp = (lakestate%uppertemp(itype,i) * epidepth + (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i) )/ (meanwater * 0.001)
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = meantemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = meantemp
      
      END SELECT
    
    ENDIF !if laketype = 0, no water in lake -> do nothing

  END SUBROUTINE calculate_lakeice_lakewater_interaction

  !>Calculate Lake T2 temperature processes
  !
  ! the concept with lakestate%WATER & lakestate%SLOWWATER and the corresponding concentrations lakestate%CONC & lakestate%CONCSLOW
  ! is partly incompatible with the conceptual model for lake temperature and it's vertical distribution (uppertemp and lowertemp)
  !
  ! it's also important to notice that even a lake without partitioning into fastpart and slowpart is still using 
  !    the state variables SLOWWATER and CONCSLOW. The difference is that for a non-splitted lake, everything left in WATER and CONC is transfered
  !    to SLOWWATER and CONCSLOW at the end of the timestep (SUBROUTINE calculate_flow_within_lake). However, this takes
  !    place AFTER the lake temperature (this routine) and ice calculations. Thus, at this point the total water stage and average temperature
  !    must take both WATER and SLOWWATER into account to be correct. If we apply the same concept on the splitted lakes, we can calculate an 
  !    average upper and lower lake temperature, but keeping individual average temperatures in the two parts of the lake.
  !  This comment is also valid for the lakeice subroutine.
  !----------------------------------------------------------
  SUBROUTINE T2_processes_in_lake(i,itype,temp,swrad,lakesurft,lakearea,frozenstate,lakestate,freezeup,freezeuparea)

    USE MODVAR, ONLY: genpar, &
                      lakedataparindex, &
                      lakedatapar, &
                      i_t2,                &
                      modeloption,         &
                      p_lakeriverice
    USE HYPEVARIABLES, ONLY: m_lddeeplake, &
                             m_t2trlake,   &
                             m_ldfastlake, &
                             m_upper2deep, &
                             m_liceTf,     &
                             m_tcflake,    &
                             m_scflake,    &
                             m_ccflake,    &
                             m_lcflake,    &
                             m_stbcorr1,  &
                             m_stbcorr2,  &
                             m_stbcorr3

    !Argument declarations
    INTEGER,INTENT(IN) :: i             !<index of subbasin
    INTEGER,INTENT(IN) :: itype         !<index of lake type (ilake = 1, olake = 2)
    REAL,INTENT(IN)    :: temp          !<air temp
    REAL,INTENT(IN)    :: swrad         !<shortwave radiation, MJ/m2/day
    REAL,INTENT(INOUT) :: lakesurft(2)  !<water surface temperature
    REAL,INTENT(IN)    :: lakearea      !<lake area
    TYPE(snowicestatetype),INTENT(IN)  :: frozenstate   !<Snow and ice states
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    INTEGER, INTENT(INOUT) :: freezeup(2) !<is water cooling below freezing piont (1 yes, 0 no)?
    REAL, INTENT(OUT)  :: freezeuparea  !<fraction of lake area with newice formation

    !Local variables    
    LOGICAL epilimnion
    REAL meantemp, meanwater
    REAL t2transfcorr,epidepth, fastoutpart
    REAL icefreefraction, freezeuparea2

    !0 Some initializations
    freezeup = 0
    epilimnion = .FALSE.
    freezeuparea = 0.
    freezeuparea2 = 0.

    !1 Lake-atmosphere T2 exchange
    ! 
    !1.1 Seasonal correction of T2 exchange coefficient   
    t2transfcorr = 1.  !Modify according to Johans suggestion below?
    
    !1.2 Depth to thermocline, function of lake area (REF)
    epidepth = 6.95 * (lakearea / 1000000)**0.185
    
    !1.3 Lake-Atmosphere T2 exchange, depending on lake distribution type
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN

      !Lake model without fast and slow split
       
      !Calculate total water stage (WATER+SLOWWATER) and average temperature (CONC*WATER+CONCSLOW*SLOWWATER)/(WATER+SLOWWATER)
      meanwater = lakestate%water(itype,i)+lakestate%slowwater(itype,i)
      IF(meanwater.GT.0.)THEN
        meantemp  = (lakestate%water(itype,i)*lakestate%conc(i_t2,itype,i)+lakestate%slowwater(itype,i)*lakestate%concslow(i_t2,itype,i))/meanwater
       
        !Check lake depth, if thermal stratification
        IF(epidepth < meanwater*0.001 .AND. itype==2)THEN !why is this only possible for olakes?
          !!Two-layer olake, waterdepth > thermocline, olake
          epilimnion = .TRUE.
          
          !->derive lake uppertemp(t) from meantemp(t, preliminary) and lowertemp(t-1)
          lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth

          !Introducing fractional ice cover to get smoother transition over the freezing point
          icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
           
          !->exchange with atmosphere - if no ice - updating meantemp(t) and uppertemp(t)
          IF(icefreefraction.GT.0.)THEN
            !temperature flow calculated from (temp-uppertemp), updating the mean temperature (concslow(i_t2,:,:))
            ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
            SELECT CASE(modeloption(p_lakeriverice))
            CASE(2) !new model based on Piccolroaz et al 2013, modified for fractional ice cover and newice formation
              CALL calculate_watersurface_heatbalance(temp-1. * (lakestate%uppertemp(itype,i) - meantemp),swrad,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction, & 
                                                      genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                      freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CASE(1) !original function of Johan/David, modified for fractional ice cover and newice formation
              CALL calculate_T2_transfer(temp-1. * (lakestate%uppertemp(itype,i) - meantemp),meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction, & 
                                         genpar(m_t2trlake)*t2transfcorr,freezeuparea,genpar(m_liceTf)) !JS4
            END SELECT

            !->re-calculate uppertemp
            lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth

            !Check freezeup conditions, indicated by relative freezeuparea
            IF(freezeuparea.GT.0.)THEN
              !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
              !re-scale to a fraction of the entire waterarea:
              freezeuparea = freezeuparea * icefreefraction
              freezeup(itype) = 1
            ENDIF
          ENDIF
        ELSE
          !Otherwise, single-layer, ilake

          !Introducing fractional ice cover to get smoother transition over the freezing point
          icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
          !!->exchange with atmosphere, if no ice, update meantemp(t)
          IF(icefreefraction.GT.0.)THEN
            ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
            SELECT CASE(modeloption(p_lakeriverice))
            CASE(2) ! new model based on Piccolroaz et al 2013
              CALL calculate_watersurface_heatbalance(temp,swrad,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction,genpar(m_tcflake), & 
                                                      genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                      freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CASE(1)
              CALL calculate_T2_transfer(temp,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction,genpar(m_t2trlake)*t2transfcorr, &
                                         freezeuparea,genpar(m_liceTf))
            END SELECT

            !Check freezeup conditions, indicated by relative freezeuparea
            IF(freezeuparea.GT.0.)THEN
              !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
              !re-scale to a fraction of the entire waterarea:
              freezeuparea = freezeuparea * icefreefraction
              freezeup(itype) = 1
            ENDIF
          ENDIF

          lakestate%uppertemp(itype,i) = meantemp
          lakestate%lowertemp(itype,i) = meantemp
        ENDIF
      ELSE
        !no water in the lake, set temperature to 0
        meantemp  = 0.
      ENDIF
 
      !Finally, assign the updated meantemp to the lakestate%conc AND lakestate%concslow 
      lakestate%concslow(i_t2,itype,i) = meantemp
      lakestate%conc(i_t2,itype,i)     = meantemp
             
    ELSE
      !
      !Lakes with split in fast and slow part
      fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
       
      !areally weighted averaged water stage
      meanwater = lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart
       
      IF(epidepth < meanwater*0.001 .AND. itype==2)THEN
        !Deep lake with thermal stratification (assume they have a common lowertemp)
        epilimnion = .TRUE.
          
        !calculate weighted average mean temperature
        IF(meanwater.GT.0.)THEN
          meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
        ELSE
          meantemp=0.
          lakestate%lowertemp(itype,i) = 0.
          lakestate%uppertemp(itype,i) = 0.
        ENDIF        
         
        !calculate upper temp from meantemp and the lowertemp
        lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
         
        !Introducing fractional ice cover to get smoother transition over the freezing point
        icefreefraction = 1. - frozenstate%lakeicecov(itype,i)

        !atmosphere-lake T2 transfer, if there is no ice
        IF(icefreefraction.GT.0.)THEN
          !atmosphere-lake interaction calculated for fastpart and slowpart separately (the same temperature gradient is effectively used here: temp-uppertemp)
          ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
          SELECT CASE(modeloption(p_lakeriverice))
          CASE(2) ! new model based on Piccolroaz et al 2013
            CALL calculate_watersurface_heatbalance(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%conc(i_t2,itype,i)),swrad, & 
                                                    lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                                    lakearea*fastoutpart*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CALL calculate_watersurface_heatbalance(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%concslow(i_t2,itype,i)),swrad, & 
                                                    lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                                    lakearea*(1-fastoutpart)*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    freezeuparea2,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
          CASE(1)
            CALL calculate_T2_transfer(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%conc(i_t2,itype,i)),lakestate%conc(i_t2,itype,i), & 
                                       lakestate%water(itype,i)*lakearea*0.001*fastoutpart,lakearea*fastoutpart*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                                    freezeuparea,genpar(m_liceTf))
            CALL calculate_T2_transfer(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%concslow(i_t2,itype,i)),lakestate%concslow(i_t2,itype,i), & 
                                       lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart),lakearea*(1-fastoutpart)*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                                    freezeuparea2,genpar(m_liceTf))
          END SELECT
            
          !Check freezeup conditions, indicated by relative freezeuparea
          IF(freezeuparea.GT.0. .OR. freezeuparea2.GT.0.)THEN
            !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
            !re-scale to a fraction of the entire waterarea:
            freezeuparea = (freezeuparea*fastoutpart+freezeuparea2*(1-fastoutpart)) * icefreefraction
            freezeup(itype) = 1
          ENDIF

          !->re-calculate averaged meantemp and upper temp, then check for freezeup conditions on the new uppertemp as well
          meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
          lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
          IF(lakestate%uppertemp(itype,i).LT.genpar(m_liceTf))THEN
            freezeup = 1
            lakestate%uppertemp(itype,i) = genpar(m_liceTf)
          ENDIF           
        ENDIF 
      ELSE
        !Shallow lake without thermal stratification

        !Introducing fractional ice cover to get smoother transition over the freezing point
        icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
        IF(icefreefraction.GT.0.)THEN
          ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
          SELECT CASE(modeloption(p_lakeriverice))
          CASE(2) ! new model based on Piccolroaz et al 2013
            CALL calculate_watersurface_heatbalance(temp,swrad, & 
                                                    lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                                    lakearea*fastoutpart*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CALL calculate_watersurface_heatbalance(temp,swrad, & 
                                                    lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                                    lakearea*(1-fastoutpart)*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    freezeuparea2,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
          CASE(1)
            CALL calculate_T2_transfer(temp,lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                       lakearea*fastoutpart*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                       freezeuparea,genpar(m_liceTf))
            CALL calculate_T2_transfer(temp,lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                       lakearea*(1-fastoutpart)*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                       freezeuparea2,genpar(m_liceTf))
          END SELECT 
         
          !Check freezeup conditions, indicated by relative freezeuparea
          IF(freezeuparea.GT.0. .OR. freezeuparea2.GT.0.)THEN
            !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
            !re-scale to a fraction of the entire waterarea:
            freezeuparea = (freezeuparea*fastoutpart+freezeuparea2*(1-fastoutpart)) * icefreefraction
            freezeup(itype) = 1
          ENDIF
       
          !assign new weighted average temperature for the lowertemp variable
          lakestate%lowertemp(itype,i) = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / &
                              (lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart)
           
          !set uppertemp to the same, for the output
          lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i) 
        ENDIF
      ENDIF
    ENDIF
          
    !2: Upper->Lower Lake T2 exchange
                    
    IF(epilimnion)THEN 
      !Autumn circulation
      IF(lakestate%uppertemp(itype,i)< lakestate%lowertemp(itype,i) .AND. lakestate%uppertemp(itype,i) > 3.95)THEN  !autumn circulation
        lakestate%lowertemp(itype,i) = meantemp
        lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i)
      ELSE
      !Spring circulation
        IF(lakestate%uppertemp(itype,i)> lakestate%lowertemp(itype,i) .AND. lakestate%uppertemp(itype,i) < 3.95)THEN  !spring circulation  
          lakestate%lowertemp(itype,i) = meantemp
          lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i)
        ELSE
          !Heat transfer between upper and lower (new function)
          CALL calculate_T2_transfer_upper2lower(lakestate%uppertemp(itype,i),lakestate%lowertemp(itype,i),epidepth*lakearea, & 
                                                (meanwater*0.001-epidepth)*lakearea,lakearea,genpar(m_upper2deep)) 
        ENDIF
      ENDIF
    ENDIF
    
    !Assign lake surface temperature if icefree conditions
    IF((icefreefraction-freezeuparea).GT.0.)lakesurft(itype) = lakestate%uppertemp(itype,i)
    
  END SUBROUTINE T2_processes_in_lake
  
  !>Calculate temperature(T2) "concentration" in lake/river precipitation
  !>due to ice presens
  !Changes the default T2 concentration, set for class
  !Modified to use fractional ice cover
  !----------------------------------------------------------
   SUBROUTINE add_T2_concentration_in_precipitation_on_water(prec,temp,snowfall,rainfall,watertemp,cprec,icecover)
  
    REAL, INTENT(IN)      :: prec       !<precipitation
    REAL, INTENT(IN)      :: temp       !<air temperature
    REAL, INTENT(IN)      :: snowfall   !<snow fall
    REAL, INTENT(IN)      :: rainfall   !<rain fall
    REAL, INTENT(IN)      :: watertemp  !<temperature of water
    REAL, INTENT(INOUT)   :: cprec      !<T2 concentration of precipitation
    REAL, INTENT(IN)      :: icecover   !<ice cover
    
    !This is now much more straight forward, using the fractional ice cover:
    ! Rainfall has always cprec = airtemp (but not lower than freezing point)
    ! Snowfall on the ice-free fraction has cprec = latentheat of freezing + sensible heat content
    ! Snowfall in the ice-covered fraction has cprec = laketemp
    IF(prec.GT.0)THEN
      !Rainfall temperature = max(0,air temp)
      cprec = rainfall * MAX(0.0,temp)
      
      !Snowfalltemp on ice   = watertemp (temporary), negative latent heat is added later when snow is melting in the ice routine
      !Snowfalltemp on water = airtemp + negative latent heat, taking into account diff spec.heat of ice and water
      cprec = cprec + snowfall * (watertemp * icecover + (MIN(temp,0.0)*2.1/4.2 - 335./4.2)*(1.-icecover))

      !Weighting by total precipitation
      cprec = cprec/prec  
    ELSE
      cprec = 0.
    ENDIF
    
  END SUBROUTINE add_T2_concentration_in_precipitation_on_water
  
!>Subroutine to calculate growth of ice on lakes and rivers (only after freezeup has been identified)
! Developers: David Gustafsson(SMHI)
!
! Model is largely based on the review of thermodynamical ice models by Leppäranta (1993)
!
!    Ice growth: "freezing-degree-day"
!    Ice melt:   "positive degree-days" with constant 0.1-0.5 degC day/cm
!
! Snow on ice is considered, however, snowaccumulation and snowmelt is supposed to be calcluated outside of this routine:
!    the routine is calculating freezing of slush ice in case the snowmass is large enough to submerge the ice surface.
!
! Input/output to the model is icedepths, snowmass, snowdepth, and mass rate of lakewater and snowwater transformation from ice
!
! Model was calibrated with Swedish lake ice depth and river ice depth data for the North Hydrology project (Gustafsson et al, 2012)
!---------------------------------------
  SUBROUTINE calculate_icedepth(tsurf,iced,biced,snowm,snowd,tair,dlakewidt,dsnowdt,ifreezeup,ibreakup,tf,kika,kexp,pm)

    !Argument declarations
    REAL, INTENT(inout) :: tsurf       !<lake surface temperature, when the lake is ice and/or snowcovered, Tsurf is back calculated from ice growth, unless its melting, then its set to 0
    REAL, INTENT(inout) :: iced        !<ice depth, cm (black ice + snowice)
    REAL, INTENT(inout) :: biced       !<black ice, cm
    REAL, INTENT(inout) :: snowm       !<snowmass, mm
    REAL, INTENT(inout) :: snowd       !<snowdepth, cm
    REAL, INTENT(in)    :: tair        !<air temperature, C
    REAL, INTENT(out)   :: dlakewidt   !<transformation of lake ice to lake water (positive direction from ice to water)
    REAL, INTENT(out)   :: dsnowdt     !<transformation of snow to lake ice       (positive direction from ice to snow)
    INTEGER, INTENT(in) :: ifreezeup   !<freeze-up day flag (1=yes, 0=no)
    INTEGER, INTENT(out):: ibreakup    !<break-up day flag  (1=yes, 0=no)
    REAL, INTENT(in)    :: tf          !<tf (~0.)  , freezing point temperature of the lake/river water, °C
    REAL, INTENT(in)    :: kika        !<kika(~10.), ratio between thermal conductivity of ice and heat exchange coef in air
    REAL, INTENT(in)    :: kexp        !<kiks(~10.), as above but for snow, actually dependent on snow density, but we use a fixed value
    REAL, INTENT(in)    :: pm          !<pm (~0.5) , degree-day melt factor for ice, cm/°C
    
    !Local variables
    REAL :: slushd    !slush depth,     cm
    REAL :: siced     !snow ice depth, cm
    REAL :: dsnow     !snow density, g/cm3
    REAL :: S         !freezing degree days, Cday
    REAL :: M         !melting degree days, Cday
    REAL :: dHsicedt  !snowice growth, cm/day, potential
    REAL :: dHsicedt2 !snowice growth, cm/day, actual
    REAL :: dHicedt   !blackice growth, cm/day, actual
    REAL :: oldiced   !icedepth at start of calculaton, cm
    !Parameters, calculated in the code from the input parameters
    REAL :: ki       ! ki, thermal conductivity of ice, W/m/C, ki ~= 1.9 (see Leppäranta-93)
    REAL :: ai       ! a, degreeday factor, see Leppäranta(1993)~=3.3 if ki=1.9 W/m/C
    REAL :: ka       ! ka = ki/kika, heat exchange in air
    REAL :: ks       ! ks = ki*(rhosnow/rhoice)^ksexp, thermal conductivity in snow
    REAL :: kiks
               
    !Physical Constants
    REAL, PARAMETER :: L = 335.      !latent heat of freezing, J/g
    REAL, PARAMETER :: dice = 0.917  !density of ice, g/cm3
    REAL, PARAMETER :: mm2cm = 0.1, cm2mm = 10. !parameters for transformation from mm to cm
    
    !Conversion of some parameters
    ki = 2.2                ! standard value 2.2 for fresh water black ice
    ki = ki * 86400. / 100. ! (W/m/oC) -> (J/d/cm/C)
    ai  = (2*ki/dice/L)**0.5      
    ka = ki/kika    ! kika and kexp can be calibration parameters if needed

    !Initialization of some variables
    dHicedt   = 0.
    dHsicedt  = 0.
    dHsicedt2 = 0.
    dlakewidt = 0.
    dsnowdt   = 0.
    oldiced = iced
    siced = iced - biced
    
    !If there is old ice or if freeze-up condition has been met, calculate ice growth and ice melt
    IF(iced.gt.0. .OR. ifreezeup.EQ.1)THEN ! ifreezeup eq. to tsurf < 0...

      !Freezing and melting degree days (actually, HYPE is using a threshold temperature for snow melt)
      S = AMAX1(0.,-Tair)   ! freezing degree days
      M = AMAX1(0.,Tair)     ! melting degree days

      !Accumulation and Melt of snow on ice is treated outside of this function
 
      !Snow density
      IF(snowm.GT.0.0 .AND. snowd.GT.0.0)THEN
        dsnow = snowm * mm2cm / snowd
        ks = ki * (dsnow/dice)**kexp  
      ELSE
        dsnow = 0.0
        ks = ki * (0.1/dice)**kexp  
      ENDIF
      kiks = ki/ks
        
      !Ice growth
      IF(snowm*mm2cm.GT.iced *(1.-dice))THEN
        !Submerged snow on ice, snowmass exceeds floating capacity of the ice
        !slush depth [cm] above ice surface (depends on snow
        !density, snow mass, and ice mass (assuming no capillary rise
        !in snow), limited by snow depth (check density if there is problem):
        slushd = (snowm*mm2cm - iced * (1.-dice))/(dsnow/dice)
        IF(slushd.gt.snowd)THEN
          WRITE(*,*) 'WARNING: slushdepth > snowdepth. slushdepth, snowdepth, dsnow:',slushd,snowd,dsnow
          slushd = snowd
        ENDIF

        ! Snow-ice growth (d(Hsi)/dt), see Leppäranta(1993), eq 21
        IF(Tair.LT.Tf)THEN
          ! height change, of the snow ice, limited by the slush depth
          dHsicedt  = ks * (Tf - Tair)/(snowd+kika)/(dice*L*(1.-dsnow/dice))
          dHsicedt2 = amin1(slushd,dHsicedt)  ! only valid for daily time steps
          
          ! update surface temperature
          tsurf = dHsicedt*(1.-dsnow/dice)*L*dice / ka + Tair
          
          !update ice depths and masses:
          snowd  = amax1(0.,snowd - dHsicedt2)          ! snow depth, cm
          snowm  = amax1(0.,snowm - dHsicedt2 * dsnow * cm2mm)  ! snowmass, mm

          siced  = amax1(0.,siced + dHsicedt2)          ! snow ice depth, cm
          slushd = amax1(0.,slushd - dHsicedt2)         ! slush depth, cm
      
          iced   = amax1(0.,biced + siced)             ! total ice depth, cm

          !how much lake water (mm) and snow mass (mm) is transformed to snow-ice?
          dlakewidt  = dlakewidt - dHsicedt2 * (1.-dsnow/dice) * cm2mm
          dsnowdt    = dsnowdt   - dHsicedt2 * dsnow * cm2mm

          ! if the potential snow-ice growth was larger than the
          ! slushdepth, it means that we have additional heat loss to
          ! freeze also the black ice, which could be used to calculate black ice
          ! growth at this point:
          
        ENDIF
      ELSE
        ! ICE SURFACE ABOVE WATER SURFACE, AND WE MAY ESTIMATE BLACK ICE GROWTH 
        slushd = 0.
        ! (black) ice growth, including insulation of snow on ice (see Leppäranta(1983), dHdt = 0.5*a^2*S/(H + ki/ka + kiks * h)
        dHicedt = 0.5*ai**2 * S /(iced + kika + snowd*kiks)
        iced = iced + dHicedt
        
        ! update surface temperature
        tsurf = dHicedt * L * dice / ka + Tair
        
        ! we do the calculation for the total ice depth (then separate
        ! snow ice from clack ice)
        biced = iced - siced
        
        !how much lake water (mm) and snow mass (mm) is transformed to snow-ice?
        dlakewidt  = dlakewidt - dHicedt * dice * cm2mm
            
      ENDIF

      ! ICE MELT, simple degree day (if there is ice, if there is no snow, and if there is positive degree days)
      IF(M.GT.0. .AND. iced.GT.0. .AND. snowd.LE.0.)THEN
        dHicedt = - AMIN1(iced,M*pm)         ! pm, degree day factor [cm/C/day]
        iced    = AMAX1(0.,iced + dHicedt)   ! total ice depth [cm]
        siced   = AMAX1(0.,siced + dHicedt); ! snow ice [cm], is melted before the black iace
        biced   = AMAX1(0.,iced-siced);      ! black ice [cm]
        
        ! how much lake water is generated?
        dlakewidt  = dlakewidt - dHicedt * dice * cm2mm

        ! set surface temperature to 0
        tsurf = 0.0
      ENDIF
        
      ! BREAK UP DAY
      IF(iced.LE.0. .AND. oldiced .GT. 0.)THEN
        ibreakup = 1
        iced = 0.
        snowd = 0.
        biced = 0.
        snowm=0.
        slushd=0.
        siced = 0.
        tsurf = 0.0
      ENDIF
    ENDIF
      
  END SUBROUTINE calculate_icedepth

!>Subroutine to calculate transfer of heat from air to water
!---------------------------------------
  SUBROUTINE calculate_T2_transfer(airtemp,watertemp,watervol,waterarea,T2transfer,freezeuparea,freezingpoint)

    !Argument declaration
    REAL, INTENT(IN)    :: airtemp       !<air temperature (deg Celsius)
    REAL, INTENT(INOUT) :: watertemp     !<water temperature (deg Celsius)
    REAL, INTENT(IN)    :: watervol      !<surface water volume (m3 or mm)
    REAL, INTENT(IN)    :: waterarea     !<surface water area (m2)
    REAL, INTENT(IN)    :: T2transfer    !<heat transfer parmeter from air to water (deg/m2/day)
    REAL, INTENT(OUT)   :: freezeuparea  !fractional area were ice formation is trigered (fraction, 0-1)
    REAL, INTENT(IN)    :: freezingpoint !freezingpoint temperature, deg C

    !Local variable declarations
    REAL t2_transf              !T2 transfer    
    REAL density
    REAL heatcapacity, thermcond
      
    density = 1000.
    heatcapacity = 4.181 * density * 1000.
    thermcond = T2transfer * 24. * 3600.
    freezeuparea = 0.
      
    IF(airtemp > watertemp)THEN
      t2_transf = MIN((airtemp - watertemp) * watervol * heatcapacity,(airtemp - watertemp)* waterarea * thermcond)
    ELSE
      t2_transf = MAX((airtemp - watertemp) * watervol * heatcapacity,(airtemp - watertemp)* waterarea * thermcond)
    ENDIF
    IF(watervol>0.)THEN
      !evaluate ice formation conditions (new temperature<freezing point)
      IF((watertemp * watervol * heatcapacity + t2_transf) / (watervol * heatcapacity).LT.freezingpoint)THEN
        !estimate a freezup area (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point 
        freezeuparea=  max(0.,min(1.,1. - (freezingpoint * (watervol * heatcapacity) - watertemp * watervol * heatcapacity)/t2_transf))
        watertemp = freezingpoint
      ELSE
        !calculate new temperature, water volume must be in m3!
        watertemp = (watertemp * watervol * heatcapacity + t2_transf) / (watervol * heatcapacity)
      ENDIF
    ENDIF
 
  END SUBROUTINE calculate_T2_transfer
      
!>Subroutine to calculate transfer of heat(temperature) between upper and lower layer in lakes 
!---------------------------------------
  SUBROUTINE calculate_T2_transfer_upper2lower(uppertemp,lowertemp,uppervol,lowervol,waterarea,T2transfer) 

    !Argument declaration
    REAL, INTENT(INOUT) :: uppertemp     !<upper water temperature (deg Celsius)
    REAL, INTENT(INOUT) :: lowertemp     !<lower water temperature (deg Celsius)
    REAL, INTENT(IN)    :: uppervol      !<upepr layer water volume (m3)
    REAL, INTENT(IN)    :: lowervol      !<upepr layer water volume (m3)
    REAL, INTENT(IN)    :: waterarea     !<surface water area (m2)
    REAL, INTENT(IN)    :: T2transfer    !<heat transfer parmeter from air to water (deg/m2/day)

    !Local variable declarations
    REAL t2_transf              !T2 transfer    
    REAL density
    REAL heatcapacity, thermcond
    REAL equiltemp
    
    density = 1000.
    heatcapacity = 4.181 * density * 1000.
    thermcond = T2transfer * 24. * 3600.
    
    !Calculate equilibrium temperature, when heat is evenly distributed
    IF((uppervol+lowervol).GT.0.)THEN
      equiltemp = (uppertemp*uppervol + lowertemp*lowervol)/(uppervol+lowervol)
    
      !calculate heatflow and update temperatures, depending on initial gradient:
      IF(uppertemp > lowertemp)THEN
        !heat flow from upper to lower
        t2_transf = (uppertemp - lowertemp)* waterarea * thermcond
        !Upper and lower temperatures, limited by equilibrium temperature
        uppertemp = MAX(equiltemp,(uppertemp * uppervol * heatcapacity - t2_transf) / (uppervol * heatcapacity))
        lowertemp = MIN(equiltemp,(lowertemp * lowervol * heatcapacity + t2_transf) / (lowervol * heatcapacity))
      ELSE
        !heat flow from lower to upper
        t2_transf = (lowertemp - uppertemp)* waterarea * thermcond
        !Upper and lower temperatures, limited by equilibrium temperature
        uppertemp = MIN(equiltemp,(uppertemp * uppervol * heatcapacity + t2_transf) / (uppervol * heatcapacity))
        lowertemp = MAX(equiltemp,(lowertemp * lowervol * heatcapacity - t2_transf) / (lowervol * heatcapacity))
      ENDIF
    ELSE
      uppertemp=0.
      lowertemp=0.
    ENDIF
 
  END SUBROUTINE calculate_T2_transfer_upper2lower

!>\brief Subroutine to calculate transfer of heat from air to water including a solar radiation term and a residual term.
!>The routine is based on the model sugested by Piccolroaz et al (2013), with modifications to use real (or estimated)
!>shortwave radiation. 
!>Partly ice covered situations can be taken into account by reducing the input waterarea
!>If the heat balance is negative enough to lower temperature below freezing, a reduction in the surface area is
!>estimated, which shows at how large area the ice is forming.
!:TODO make T2 subroutines work for other timestep than day
!---------------------------------------
  SUBROUTINE calculate_watersurface_heatbalance(airtemp,swrad,watertemp,watervol,waterarea,tempcoef,radcoef,constcoef,lincoef,freezeuparea,freezingpoint,stabpar1,stabpar2,stabpar3)

  USE MODVAR, ONLY: cwater,seconds_per_timestep
    
    !Argument declaration
    REAL, INTENT(IN)    :: airtemp       !<air temperature (deg Celsius)
    REAL, INTENT(IN)    :: swrad         !<shortwave radiation (MJ/m2/day)
    REAL, INTENT(INOUT) :: watertemp     !<water temperature (deg Celsius)
    REAL, INTENT(IN)    :: watervol      !<water volume (m3)
    REAL, INTENT(IN)    :: waterarea     !<water surface area (m2)
    REAL, INTENT(IN)    :: tempcoef      !<heat transfer parmeter from air to water (J/m2/s/deg)
    REAL, INTENT(IN)    :: radcoef       !<heat transfer parmeter from radiation to water (fraction, 0-1)
    REAL, INTENT(IN)    :: constcoef     !<heat transfer parmeter, constant residual term (J/m2/s)
    REAL, INTENT(IN)    :: lincoef       !<heat transfer parmeter, linear residualterm (J/m2/s/deg)
    REAL, INTENT(OUT)   :: freezeuparea  !<fractional area were ice formation is trigered (fraction, 0-1)
    REAL, INTENT(IN)    :: freezingpoint !<freezingpoint temperature, deg C
    REAL, INTENT(IN)    :: stabpar1      !<Stability parameter, affects both heating and cooling. No correction if set to zero
    REAL, INTENT(IN)    :: stabpar2      !<Stability parameter, affects cooling. No correction if set to zero
    REAL, INTENT(IN)    :: stabpar3      !<Stability parameter, affects heating. No correction if set to zero
       
    !Local variable declarations
    REAL netheat                  !Net heat flux to the water (J/timestep)    
    REAL density                  !Water density (kg/m3)
    REAL heatcapacity             !heat capacity of water (J/m3/deg)
    REAL timestep                 !timestep, seconds
    REAL tempdiff                 !Temperature difference
    REAL stabfunction             !Stability correction function

    REAL, PARAMETER :: seconds_per_day = 86400.  

    !> \b Algorithm \n
    density      = 1000.                    ! kg/m3, density of water
    heatcapacity = cwater * density * 1000. ! J/m3/deg  [kJ/kg/deg * kg/m3 * 1/k]
    timestep     = seconds_per_timestep     ! 86400 seconds for daily timestep (only working timestep for T2 so far)
           
    freezeuparea = 0.
                          
    !make calculation only if the water has a volume
    IF(watervol>0.)THEN
    
      !initialize the net heat flux, J/timestep
      netheat = 0.
      
      !>Calculate stability correction for heat exchange between air and water
      tempdiff = airtemp - watertemp
      IF(tempdiff>0.) THEN
        stabfunction = 1./(1. + stabpar1 * tempdiff)**stabpar3
      ELSE
        stabfunction = 1./(1. - stabpar1 * tempdiff)**(-stabpar2)
      ENDIF
      
      
      !>Add the air temperature term
      IF(airtemp > watertemp)THEN
        netheat = netheat + MIN((airtemp - watertemp) * watervol * heatcapacity, stabfunction * (airtemp - watertemp)* waterarea * tempcoef * timestep) !J/timestep
      ELSE
        netheat = netheat + MAX((airtemp - watertemp) * watervol * heatcapacity, stabfunction * (airtemp - watertemp)* waterarea * tempcoef * timestep) !J/timestep
      ENDIF
      
      !>Add the radiation term, MJ/m2/day => J/m2/s and then multiplied with timestep in s.
      netheat = netheat + 1.E6 * swrad /seconds_per_day * waterarea * radcoef * timestep
      
      !>Add the residual term, same units as temperature equation
      netheat = netheat + (watertemp*lincoef + constcoef) * waterarea * tempcoef * timestep

      !>Evaluate ice formation conditions (new temperature<freezing point) and calculate new water temperature
      IF((watertemp * watervol * heatcapacity + netheat) / (watervol * heatcapacity).LT.freezingpoint)THEN
        !estimate a freezup area (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point 
        freezeuparea = MAX(0.,MIN(1.,1. - (freezingpoint * (watervol * heatcapacity) - watertemp * watervol * heatcapacity)/netheat))
        watertemp = freezingpoint
      ELSE
        !calculate new temperature, water volume must be in m3!
        watertemp = (watertemp * watervol * heatcapacity + netheat) / (watervol * heatcapacity)
      ENDIF
    ENDIF
 
  END SUBROUTINE calculate_watersurface_heatbalance

!Subroutine to detect lake ice breakup conditions      
!---------------------------------------
!  SUBROUTINE calculate_lakeice_breakup()
!  END SUBROUTINE calculate_lakeice_breakup

!! Subroutine to calculate freezeup of lakes and rivers (two alternative methods)
!! Developers: David Gustafsson, SMHI
!! method 1: (default) dynamic water temperature model (Strömqvist, Lindström, Gustafsson, mfl, SMHI)
!! method 2: moving average air temperature filter as proxy for water surface temperature (e.g. Lindström&Lövenius, 2000)
!!---------------------------------------
!!  SUBROUTINE calculate_lakeice_freezeup(Tair,Twater,iced,TTf,Tf,ifreezeup)
!!
!!    !Argument declarations
!!    REAL, INTENT(IN)     :: Tair          !<air temperature, oC
!!    REAL, INTENT(IN)     :: Twater        !<water surface temperature, oC
!!    REAL, INTENT(IN)     :: iced          !<total ice depth, cm
!!    REAL, INTENT(IN)     :: TTf           !<threshold water temperature for freezeup
!!    REAL, INTENT(IN)     :: Tf            !<freezing point temperature
!!    INTEGER, INTENT(OUT) :: ifreezeup     !<flag for freeze-up and break-up day
!        
!!Parameters
! !       real, intent(IN) :: p1, p2, meandepth ! parameters used for deriving the freezup day filter parameter from lake depth (not used in future)
! !       real             :: pf
!
!! options (for development)
!        INTEGER          :: FreezeupModel    ! option for freezeup model, option 1 or 2        
!        FreezeupModel = 1
!
!        ifreezeup = 0 ! set freezeupday flag initially to 0
!        
!!Check ice conditions: if ice is present(iced>0), no need to calculate freezeupday
!        IF(iced.LE.0.)THEN
!            IF(FreezeupModel.EQ.1)THEN
!!Alt 1: identify freeze-up day with simulated surface water temperature<= freezeup_threshold and air temperature<= freezing_temperature
!                IF(Twater.LE.TTf .AND. Tair.LE.Tf)THEN
!                    ifreezeup = 1
!                ENDIF
!            ELSE
!!Alt 2: identify freeze-up day with air temperature filter
!  !              pf = (p1/amax1(0.5,meandepth))**p2;       ! Delay parameter pf as a function of mean lake depth
!  !              Tfreezeup = Tfreezeup*(1-pf) + pf * Tair ! slowly moving average filter on air temp
!  !              IF(Tfreezeup.lt.0)THEN
!                    ifreezeup = 1
!  !              ENDIF
!            ENDIF
!        ELSE
!            ifreezeup = 0
!        ENDIF
!  END SUBROUTINE calculate_lakeice_freezeup
!      
!!Subroutine to calculate freezeup of rivers
!! Developers: David Gustafsson, SMHI
!! method:    dynamic water temperature model (Strömqvist, Lindström, Gustafsson, mfl, SMHI)
!!---------------------------------------
!  SUBROUTINE calculate_riverice_freezeup(Tair,Twater,iced,TTf,Tf,ifreezeup)
!
!    !Argument declarations
!    REAL, INTENT(IN)     :: iced          !<total ice depth, cm
!    REAL, INTENT(IN)     :: Twater        !<water surface temperature, oC
!    REAL, INTENT(IN)     :: Tair          !<air temperature, oC
!    INTEGER, INTENT(OUT) :: ifreezeup     !<flag for freeze-up and break-up day
!    REAL, INTENT(IN)     :: TTf           !<threshold water temperature for freezeup
!    REAL, INTENT(IN)     :: Tf            !<freezing point temperature
!
!    ifreezeup = 0 ! set freezeupday flag initially to 0
!
!    !Check ice conditions: if ice is present(iced>0), no need to calculate freezeupday
!    IF(iced.LE.0.)THEN
!      !Identify freeze-up day with simulated surface water temperature<= freezeup_threshold and air temperature<= freezing_temperature
!      IF(Twater.LE.TTf .AND. Tair.LE.Tf)THEN
!        ifreezeup = 1
!      ENDIF
!    ELSE
!      ifreezeup = 0
!    ENDIF
!  END SUBROUTINE calculate_riverice_freezeup

  !>\brief Subroutine for finding current lake outflow parameters. 
  !------------------------------------------------------------------------------
  SUBROUTINE inline_dam_outflow(i,wlmr,damProd)
       
    USE MODVAR, ONLY : missing_value,     &
                       dayno,   &
                       pi,      &
                       lake,    &
                       dam,     &
                       lakebasin,  &
                       genpar,  &
                       lakeindex, &
                       damindex,    &
                       lakebasinindex,  &
                       lakedatapar, &
                       lakedataparindex, &
                       ilregpar,olregpar, basin

    !Argument declarations
    INTEGER, INTENT(IN) :: i            !<index of current subbasin
    REAL, INTENT(IN)    :: wlmr         !<water lake level (m)
    REAL, INTENT(OUT)   :: damProd      !<current dam production flow (m3/s)
    
    !Local variable declarations
    INTEGER :: counter  = 0  !<loop counter
    REAL    :: wlmr0    = 0. !<water level in local reference system (m)
    REAL    :: wsl      = 0. !<current water level calculated from storage rule (m)
    REAL    :: opt1m3s  = 0. !<selected outflow
    REAL    :: opt2m3s  = 0. !<outflow after powerhouse restriction
    REAL    :: opt3m3s  = 0. !<outflow after weekly change restriction
    REAL    :: opt4m3s  = 0. !<outflow after daily change restriction
    REAL    :: maxqmult = 0. !<Coefficient to Scale Transition Band Type
    REAL    :: stora    = 0. !<where: S = A * (WSL - B) ^ C + D
    REAL    :: storb    = 0. !<
    REAL    :: storc    = 0. !<
    REAL    :: stord    = 0. !<
    REAL    :: minqa    = 0. !<where: Qo = A * (WSL - B) ^ C + D
    REAL    :: minqb    = 0. !<
    REAL    :: minqc    = 0. !<
    REAL    :: minqd    = 0. !<
    REAL    :: midqa    = 0. !<
    REAL    :: midqb    = 0. !<
    REAL    :: midqc    = 0. !<
    REAL    :: midqd    = 0. !<
    REAL    :: maxqa    = 0. !<
    REAL    :: maxqb    = 0. !<
    REAL    :: maxqc    = 0. !<
    REAL    :: maxqd    = 0. !<
    REAL    :: deltaw   = 0. !<Weekly Flow Change Limit
    REAL    :: deltad   = 0. !<Daily Flow Change Limit
    REAL    :: phqo     = 0. !<Powerhouse Limit Flow Value
    REAL    :: drqo     = 0. !<Minimum Fixed Flow Value
    REAL    :: flqo     = 0. !<Maximum Fixed Flow Value
    INTEGER :: drtype   = 0  !<0 = Dr. Percentile, 1 = Fixed Min.
    INTEGER :: mintype  = 0  !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: midtype  = 0  !<3 = Band Ideal,  0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: maxtype  = 0  !<0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: fltype   = 0  !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    INTEGER :: phtype   = 0  !<0 = No PH Restriction, 1 = Hold Water
     REAL    :: monthly(12,9) !<Used for monthly regression calculations, aka weighted flow, by column
                             !<column 1:minA 2:midA 3:highA 4:natminA 5:natmidA 6: nathighA 7:minB 8:midB 9:highB
    REAL    :: daily(366,18) !<Stores daily extra data by column
    !column 1:drought level 2:lower op.s 3:lower band 4:avg.level 5:upper band 6:upper op.s 7:flood level
    !column 8:(smoothed) drought percentile 9:sm.min.prcntl 10:sm.avg.prcntl 11:sm.max.prcntl 12:sm.flood.prcntl
    !column 13:natural drought percentile 14:nat.min.prcntl 15:nat.avg.prcntl 16:nat.max.prcntl 17:nat.flood.prcntl 18:monthly weighting
    REAL    :: weekly7day(7) !<stores powerhouse restricted flows from past 7 days. weekly7day(1) = today's flow, weekly7day(2) = yesterday's flow, etc.
    REAL    :: dailyyes = 0. !yesterday's daily restricted flow
        
    !time-based parameters
    INTEGER :: idt      = 0  !<used to check if early-model calculations apply
    INTEGER :: idtlast  = 0  !<prevents sabd14day(:), weekly7day(:) updating when get_current_outflow_parameters()
                             !<is called by check_ouflow_from_lake()
    INTEGER :: curmonth = 0  !<current month, extracted from datetype
    INTEGER :: dayofmonth = 0 !<current day of month, extracted from datetype
    INTEGER :: idtlimit = 0  !<equals idt if early iteration, otherwise equals size of qin14day(:), sabd14day(:), weekly7day(:)
    INTEGER :: backmonth = 0 !<if daymonth <= 15, month-1, ow month
    INTEGER :: formonth = 0 !<if daymonth <= 15, month, ow month+1
    REAL    :: mweight  = 0.
    
    !daily water surface level parameters
    REAL    :: drlevel  = 0. !<stores value in daily(dayno,1)
    REAL    :: lowops   = 0. !<daily(dayno,2)
    REAL    :: lowband  = 0. !<daily(dayno,3)
    REAL    :: avglevel = 0. !<daily(dayno,4)
    REAL    :: upband   = 0. !<daily(dayno,5)
    REAL    :: upops    = 0. !<daily(dayno,6)
    REAL    :: fllevel  = 0. !<daily(dayno,7)
    
    !daily smoothed percentile flow (type 0)
    REAL    :: drprcntl = 0. !<daily(dayno,8)
    REAL    :: minprcntl = 0. !<daily(dayno,9)
    REAL    :: midprcntl = 0. !<daily(dayno,10)
    REAL    :: maxprcntl = 0. !<daily(dayno,11)
    REAL    :: flprcntl = 0. !<daily(dayno,12)
    
    !ideal band Q (type 3)
    REAL    :: sabtod   = 0. !<average band storage today
    REAL    :: sabtom   = 0. !<average band storage tomorrow
    REAL    :: sabd     = 0. !<average band storage change
    REAL    :: percentbfm = 0. !<percent band flow multiplier, limited by maxqmult
    REAL    :: idealbandq = 0. !<flow required to achieve average WSL
    REAL    :: qin14day(14)  !<stores inflow from past 14 days. qin14day(1) = today's qin, qin14day(2) = yesterday's qin, etc.
    REAL    :: sabd14day(14) !<stores average band storage change from past 14 days. sabd14day(1) = today's sabd, sabd14day(2) = yesterday's sabd, etc.
       
    !weighted, aka monthly regression flow (type 1)
    REAL    :: weightedlow = 0. !<monthly regressed lower flow
    REAL    :: midbase  = 0. !<weighted middle lower base
    REAL    :: weightedmid = 0. !<monthly regressed middle flow
    REAL    :: highbase = 0. !<weighted high lower base
    REAL    :: weightedhigh = 0. !<monthly regressed high flow
    
    !physical, aka A-B-C-D flow (type 2)
    REAL    :: minphysical = 0. !<A-B-C-D minimum flow
    REAL    :: midphysical = 0. !<A-B-C-D middle flow
    REAL    :: maxphysical = 0. !<A-B-C-D maximum flow
    
    !selected outflows (based on type)
    REAL    :: drflow      = 0. !<0 = Dr. Percentile, 1 = Fixed Min.
    REAL    :: lowflow     = 0. !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: corlowflow  = 0. !<sets lowflow >= droughtflow
    REAL    :: lowtransflow = 0. !<interpolate between lowflow and opflow
    REAL    :: opflow      = 0. !<3 = Band Ideal,  0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: coropflow   = 0. !<sets opflow>=lowflow, <=highflow
    REAL    :: hightransflow = 0. !<interpolate between opflow and highflow
    REAL    :: highflow    = 0. !<0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: corhighflow = 0. !<sets highflow<=flflow
    REAL    :: flflow      = 0. !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    REAL    :: transflow   = 0. !<temporary variable that stores the linear interpolation between flows based on WSL
    
    !calculation of slope of last week's powerhouse flows
    REAL    :: sumx        = 0. !<sum of days
    REAL    :: sumy        = 0. !<sum of ph flows
    REAL    :: sumxx       = 0. !<sum of days squared
    REAL    :: sumxy       = 0. !<sum of days times ph flows
    REAL    :: weeklyslope = 0. !<slope of the powerhouse flows
    
     
    
    !local variable assignment
    maxqmult = dam(damindex(i))%maxqmult
    stora = dam(damindex(i))%stora
    storb = dam(damindex(i))%storb
    storc = dam(damindex(i))%storc
    stord = dam(damindex(i))%stord
    minqa = dam(damindex(i))%minqa
    minqb = dam(damindex(i))%minqb
    minqc = dam(damindex(i))%minqc
    minqd = dam(damindex(i))%minqd
    midqa = dam(damindex(i))%midqa
    midqb = dam(damindex(i))%midqb
    midqc = dam(damindex(i))%midqc
    midqd = dam(damindex(i))%midqd
    maxqa = dam(damindex(i))%maxqa
    maxqb = dam(damindex(i))%maxqb
    maxqc = dam(damindex(i))%maxqc
    maxqd = dam(damindex(i))%maxqd
    deltaw = dam(damindex(i))%deltaw
    deltad = dam(damindex(i))%deltad
    phqo = dam(damindex(i))%phqo
    drqo = dam(damindex(i))%drqo
    flqo = dam(damindex(i))%flqo
    drtype = dam(damindex(i))%drtype
    mintype = dam(damindex(i))%mintype
    midtype = dam(damindex(i))%midtype
    maxtype = dam(damindex(i))%maxtype
    fltype = dam(damindex(i))%fltype
    phtype = dam(damindex(i))%phtype
    monthly = dam(damindex(i))%monthly
    daily = dam(damindex(i))%daily
    idt = dam(damindex(i))%idt
    idtlast = dam(damindex(i))%idtlast
    curmonth = dam(damindex(i))%curmonth
    dayofmonth = dam(damindex(i))%dayofmonth
    mweight = daily(dayno,18)
    qin14day = dam(damindex(i))%qin14day
    sabd14day = dam(damindex(i))%sabd14day
    weekly7day = dam(damindex(i))%weekly7day
    dailyyes = dam(damindex(i))%dailyyes
    
    !determine backwards and forwards month for monthly regression
    IF(curmonth==1.AND.dayofmonth<=15)THEN
      backmonth = 12
    ELSEIF(dayofmonth<=15)THEN
      backmonth = curmonth - 1
    ELSE
      backmonth = curmonth
    ENDIF
    
    IF(curmonth==12.AND.dayofmonth>15)THEN
      formonth = 1
    ELSEIF(dayofmonth>15)THEN
      formonth = curmonth + 1
    ELSE
      formonth = curmonth
    ENDIF
    
    !calculation of WSL
    wlmr0 = wlmr + dam(damindex(i))%w0ref !water level in local reference system
    !IF(wlmr0 < storb .OR. wlmr0 /= wlmr0) wlmr0 = storb + 1.0 !stopgap solution to large negative wlmr values
    dam(damindex(i))%storyes = stora * ((wlmr0 - storb) ** storc) + stord  
    dam(damindex(i))%stortod = dam(damindex(i))%storyes + dam(damindex(i))%qinprev - dam(damindex(i))%qoutprev
    wsl = ((dam(damindex(i))%stortod - stord) / stora) ** (1/storc) + storb
    
    !read in daily WSL bands
    drlevel = daily(dayno,1)
    lowops = daily(dayno,2)
    lowband = daily(dayno,3)
    avglevel = daily(dayno,4)
    upband = daily(dayno,5)
    upops = daily(dayno,6)
    fllevel = daily(dayno,7)
    
    !read in smoothed percentileflows
    drprcntl = daily(dayno,8)
    minprcntl = daily(dayno,9)
    midprcntl = daily(dayno,10)
    maxprcntl = daily(dayno,11)
    flprcntl = daily(dayno,12)
    
    !determine ideal band Q
    updatesabd:IF(idt>idtlast)THEN !only updates sabd14day when called via calculate_outflow_from_lake()
      sabtod = stora * ((avglevel - storb) ** storc) + stord
      IF(dayno+1==367)THEN !goes to next year
        sabtom = stora * ((daily(1,4) - storb) ** storc) + stord
      ELSE
        sabtom = stora * ((daily(dayno+1,4) - storb) ** storc) + stord !daily(dayno+1,4) is tomorrow's avglevel
      ENDIF
      sabd = sabtom - sabtod
      DO counter = 13,1,-1
        sabd14day(counter+1) = sabd14day(counter)
      ENDDO
      sabd14day(1) = sabd
      dam(damindex(i))%sabd14day = sabd14day      
    ENDIF updatesabd
    
    calcpercent:IF(wsl<lowband)THEN
      percentbfm = -1*maxqmult
    ELSEIF(wsl>upband)THEN calcpercent
      percentbfm = maxqmult
    ELSEIF(wsl>avglevel)THEN calcpercent
      percentbfm = maxqmult * (wsl - avglevel) / (upband - avglevel)
    ELSE calcpercent
      percentbfm = -1*maxqmult * (avglevel - wsl) / (avglevel - lowband)
    ENDIF calcpercent
    
    idtlimit = idt
    IF(idtlimit>14) idtlimit=14
    idealbandq = (SUM(qin14day) / idtlimit) + ((1+percentbfm) * (SUM(sabd14day) / idtlimit))
    
    
    !calculate weighted flows (monthly regression)    
    weightedlow = ((monthly(backmonth,1) * (wsl - monthly(backmonth,7))) * mweight) + ((monthly(formonth,1) * (wsl - monthly(formonth,7))) * (1-mweight))
    midbase  = ((monthly(backmonth,1) * (lowops - monthly(backmonth,7))) * mweight) + ((monthly(formonth,1) * (lowops - monthly(formonth,7))) * (1-mweight))
    weightedmid = ((monthly(backmonth,2) * (wsl - monthly(backmonth,8))) * mweight) + ((monthly(formonth,2) * (wsl - monthly(formonth,8))) * (1-mweight)) + midbase
    highbase = ((monthly(backmonth,2) * (upband - monthly(backmonth,8))) * mweight) + ((monthly(formonth,2) * (upband - monthly(formonth,8))) * (1-mweight))
    weightedhigh = ((monthly(backmonth,3) * (wsl - monthly(backmonth,9))) * mweight) + ((monthly(formonth,3) * (wsl - monthly(formonth,9))) * (1-mweight)) + highbase
    
    !calculate physical flows (A-B-C-D)
    minphysical = minqa * ((wsl - minqb) ** minqc) + minqd
    midphysical = midqa * ((wsl - midqb) ** midqc) + midqd
    maxphysical = maxqa * ((wsl - maxqb) ** maxqc) + maxqd
    
    !select outflows (based on type)
    !drought flow
    choosedr:IF(drtype==0)THEN
      drflow = drprcntl
    ELSEIF(drtype==1)THEN choosedr
      drflow = drqo
    ENDIF choosedr
    
    !minimum flow
    chooselow:IF(mintype==0)THEN
      lowflow = minprcntl
    ELSEIF(mintype==1)THEN chooselow
      lowflow = weightedlow
    ELSEIF(mintype==2)THEN chooselow
      lowflow = minphysical
    ENDIF chooselow
    
    !operations flow
    chooseop:IF(midtype==0)THEN
      opflow = midprcntl
    ELSEIF(midtype==1)THEN chooseop
      opflow = weightedmid
    ELSEIF(midtype==2)THEN chooseop
      opflow = midphysical
    ELSEIF(midtype==3)THEN chooseop
      opflow = idealbandq
    ENDIF chooseop
    
    !maximum flow
    choosehigh:IF(maxtype==0)THEN
      highflow = maxprcntl
    ELSEIF(maxtype==1)THEN choosehigh
      highflow = weightedhigh
    ELSEIF(maxtype==2)THEN choosehigh
      highflow = maxphysical
    ENDIF choosehigh
    
    !flood flow
    choosefl:IF(fltype==0)THEN
      flflow = flprcntl
    ELSEIF(fltype==1)THEN choosefl
      flflow = flqo
    ELSEIF(fltype==2)THEN choosefl
      flflow = highflow
    ENDIF choosefl
    
    !correct flows
    !correct minimum flow
    correctlow:IF(lowflow < drflow)THEN
      corlowflow = drflow
    ELSE correctlow
      corlowflow = lowflow
    ENDIF correctlow
    
    !correct maximum flow
    correcthigh:IF(highflow > flflow)THEN
      corhighflow = flflow
    ELSE correcthigh
      corhighflow = highflow
    ENDIF correcthigh
    
    !correct operations flow
    correctop:IF(opflow < corlowflow)THEN
      coropflow = corlowflow
    ELSEIF(opflow > corhighflow)THEN correctop
      coropflow = corhighflow
    ELSE correctop
      coropflow = opflow
    ENDIF correctop
    
    !determine transition flows
    
    !determine low transition flow
    transflow = corlowflow + ((coropflow - corlowflow) * (wsl - lowops) / (lowband - lowops))
    lowtrans:IF(transflow < corlowflow)THEN
      lowtransflow = corlowflow
    ELSEIF(transflow > coropflow)THEN lowtrans
      lowtransflow = coropflow
    ELSE lowtrans
      lowtransflow = transflow
    ENDIF lowtrans
    
    !determine high transition flow
    transflow = coropflow + ((corhighflow - coropflow) * (wsl - upband) / (upops - upband))
    hightrans:IF(transflow < coropflow)THEN
      hightransflow = coropflow
    ELSEIF(transflow > corhighflow)THEN hightrans
      hightransflow = corhighflow
    ELSE hightrans
      hightransflow = transflow
    ENDIF hightrans
    
    !determine outflow before restrictions
    selectflow:IF(wsl > fllevel)THEN
      opt1m3s = flflow
    ELSEIF(wsl > upops)THEN selectflow
      opt1m3s = corhighflow
    ELSEIF(wsl > upband)THEN selectflow
      opt1m3s = hightransflow
    ELSEIF(wsl > lowband)THEN selectflow
      opt1m3s = coropflow
    ELSEIF(wsl > lowops)THEN selectflow
      opt1m3s = lowtransflow
    ELSEIF(wsl > drlevel)THEN selectflow
      opt1m3s = corlowflow
    ELSE selectflow
      opt1m3s = drflow
    ENDIF selectflow
    
    !determine powerhouse restriction
    phrestriction:IF(phtype==1)THEN
      IF(opt1m3s < phqo)THEN
        opt2m3s = opt1m3s
      ELSEIF(wsl > upops)THEN !use spillway?
        opt2m3s = opt1m3s
      ELSE
        opt2m3s = phqo
      ENDIF
    ELSE phrestriction
      opt2m3s = opt1m3s
    ENDIF phrestriction
    
    !determine weekly flow change restriction
    updateweekly:IF(idt>idtlast)THEN !only updates weekly7day when called via calculate_outflow_from_lake()
      DO counter = 6,1,-1
        weekly7day(counter+1) = weekly7day(counter)
      ENDDO
      weekly7day(1) = opt2m3s
      dam(damindex(i))%weekly7day = weekly7day
    ENDIF updateweekly
    
    idtlimit = idt
    IF(idtlimit>7) idtlimit=7
    
    !because weekly7day is in inverse order (recent to older flows)
    !decreasing flows would actually have a positive slope, and vice versa
    DO counter = 1,idtlimit
      sumx  = sumx + REAL(counter)
      sumy  = sumy + weekly7day(counter)
      sumxx = sumxx + (REAL(counter) * REAL(counter))
      sumxy = sumxy + (REAL(counter) * weekly7day(counter))
    ENDDO
    
    weeklyslope = (sumxy - (sumx * sumy / REAL(idtlimit))) / (sumxx - (sumx * sumx / REAL(idtlimit)))
    
    weeklyrestriction:IF(MAXVAL(weekly7day(1:idtlimit)) - MINVAL(weekly7day(1:idtlimit)) < deltaw)THEN
      opt3m3s = opt2m3s
    ELSEIF(weeklyslope < 0)THEN weeklyrestriction !original excel has (slope > 0)
      opt3m3s = MINVAL(weekly7day(1:idtlimit)) + deltaw
    ELSE weeklyrestriction
      opt3m3s = MAXVAL(weekly7day(1:idtlimit)) - deltaw
    ENDIF weeklyrestriction
    
    !determine daily flow change restriction
    dailyrestriction:IF(idt==1)THEN !only one timestep
      opt4m3s = opt3m3s
    ELSEIF(ABS(opt3m3s - dailyyes) < deltad)THEN dailyrestriction
      opt4m3s = opt3m3s
    ELSEIF(opt3m3s > dailyyes)THEN dailyrestriction
      opt4m3s = dailyyes + deltad
    ELSE dailyrestriction
      opt4m3s = dailyyes - deltad
    ENDIF dailyrestriction
    
    IF(idt == idtlast) dam(damindex(i))%dailyyes = opt4m3s !only updates dailyyes when called via check_outflow_from_lake()
    !update yesterday's flow at the end of the timestep, not the beginning
    
    damProd = opt4m3s
    IF(idt > idtlast) dam(damindex(i))%idtlast = idt
    
    
    
  END SUBROUTINE inline_dam_outflow
  
  
  !>\brief Subroutine for finding current lake outflow parameters. 
  !------------------------------------------------------------------------------
  SUBROUTINE OCBD_dam_outflow(i,wlmr,damProd, branchflow)
       
    USE MODVAR, ONLY : missing_value,     &
                       dayno,   &
                       pi,      &
                       lake,    &
                       dam,     &
                       lakebasin,  &
                       genpar,  &
                       lakeindex, &
                       damindex,    &
                       lakebasinindex,  &
                       lakedatapar, &
                       lakedataparindex, &
                       ilregpar,olregpar, basin

    !Argument declarations
    INTEGER, INTENT(IN) :: i            !<index of current subbasin
    REAL, INTENT(IN)    :: wlmr         !<water lake level (m)
    REAL, INTENT(OUT)   :: damProd      !<current dam production flow (m3/s)
    REAL, INTENT(OUT)   :: branchflow   !<current branch flow (m3/s)
    
    !Local variable declarations
    INTEGER :: counter  = 0  !<loop counter
    REAL    :: wlmr0    = 0. !<water level in local reference system (m)
    REAL    :: wsl      = 0. !<current water level calculated from storage rule (m)
    REAL    :: maxqmult = 0. !<Coefficient to Scale Transition Band Type
    REAL    :: stora    = 0. !<where: S = A * (WSL - B) ^ C + D
    REAL    :: storb    = 0. !<
    REAL    :: storc    = 0. !<
    REAL    :: stord    = 0. !<
    INTEGER :: natswitch   = 0 !<0 = No Natural Outflow, 1 = Has Natural Outflow
    INTEGER :: locswitch   = 0 !<0 = No Local Inflow, 1 = Has Local Inflow
    INTEGER :: branchswt   = 0 !0 = reservoir is main, natural is branch, 1 = reservoir is branch, natural is main
    INTEGER :: condswitch  = 0 !<0 = No Flow Condition, 1 = Has Conditioning   
    INTEGER :: condrestr   = 0 !<0 = Low to High, 1 = Drought to Flood
    INTEGER :: condqowsl   = 0 !<0 = Conditioned by Qo, 1 = Conditioned by WSL
    INTEGER :: condsubid   = 0 !<HYPE Basin SubID Used for Conditioning Value
    INTEGER :: condindex   = 0 !<HYPE Basin Index Used for Conditioning Value
    REAL    :: conda       = 0. !<where: Qo = A * (XX - B) ^ C + D
    REAL    :: condb       = 0. !<
    REAL    :: condc       = 0. !<
    REAL    :: condd       = 0. !<
    REAL    :: monthly(12,9) !<Used for monthly regression calculations, aka weighted flow, by column
                             !<column 1:minA 2:midA 3:highA 4:natminA 5:natmidA 6: nathighA 7:minB 8:midB 9:highB
    REAL    :: daily(366,18) !<Stores daily extra data by column
    !column 1:drought level 2:lower op.s 3:lower band 4:avg.level 5:upper band 6:upper op.s 7:flood level
    !column 8:(smoothed) drought percentile 9:sm.min.prcntl 10:sm.avg.prcntl 11:sm.max.prcntl 12:sm.flood.prcntl
    !column 13:natural drought percentile 14:nat.min.prcntl 15:nat.avg.prcntl 16:nat.max.prcntl 17:nat.flood.prcntl 18:monthly weighting
    
    !time-based parameters
    INTEGER :: idt      = 0  !<used to check if early-model calculations apply
    INTEGER :: idtlast  = 0  !<prevents sabd14day(:), weekly7day(:) updating when get_current_outflow_parameters()
                             !<is called by check_ouflow_from_lake()
    INTEGER :: curmonth = 0  !<current month, extracted from datetype
    INTEGER :: dayofmonth = 0 !<current day of month, extracted from datetype
    INTEGER :: idtlimit = 0  !<equals idt if early iteration, otherwise equals size of qin14day(:), sabd14day(:), weekly7day(:)
    INTEGER :: backmonth = 0 !<if daymonth <= 15, month-1, ow month
    INTEGER :: formonth = 0 !<if daymonth <= 15, month, ow month+1
    REAL    :: mweight  = 0.
    
    !reservoir variable declarations
    REAL    :: opt1m3s  = 0. !<selected outflow
    REAL    :: opt2m3s  = 0. !<outflow after powerhouse restriction
    REAL    :: opt3m3s  = 0. !<outflow after weekly change restriction
    REAL    :: opt4m3s  = 0. !<outflow after daily change restriction   
    REAL    :: minqa    = 0. !<where: Qo = A * (WSL - B) ^ C + D
    REAL    :: minqb    = 0. !<
    REAL    :: minqc    = 0. !<
    REAL    :: minqd    = 0. !<
    REAL    :: midqa    = 0. !<
    REAL    :: midqb    = 0. !<
    REAL    :: midqc    = 0. !<
    REAL    :: midqd    = 0. !<
    REAL    :: maxqa    = 0. !<
    REAL    :: maxqb    = 0. !<
    REAL    :: maxqc    = 0. !<
    REAL    :: maxqd    = 0. !<
    REAL    :: deltaw   = 0. !<Weekly Flow Change Limit
    REAL    :: deltad   = 0. !<Daily Flow Change Limit
    REAL    :: phqo     = 0. !<Powerhouse Limit Flow Value
    REAL    :: drqo     = 0. !<Minimum Fixed Flow Value
    REAL    :: flqo     = 0. !<Maximum Fixed Flow Value
    INTEGER :: drtype   = 0  !<0 = Dr. Percentile, 1 = Fixed Min.
    INTEGER :: mintype  = 0  !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: midtype  = 0  !<3 = Band Ideal,  0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: maxtype  = 0  !<0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: fltype   = 0  !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    INTEGER :: phtype   = 0  !<0 = No PH Restriction, 1 = Hold Water
    
    REAL    :: weekly7day(7) !<stores powerhouse restricted flows from past 7 days. weekly7day(1) = today's flow, weekly7day(2) = yesterday's flow, etc.
    REAL    :: dailyyes = 0. !yesterday's daily restricted flow
    
    
    !daily water surface level parameters
    REAL    :: drlevel  = 0. !<stores value in daily(dayno,1)
    REAL    :: lowops   = 0. !<daily(dayno,2)
    REAL    :: lowband  = 0. !<daily(dayno,3)
    REAL    :: avglevel = 0. !<daily(dayno,4)
    REAL    :: upband   = 0. !<daily(dayno,5)
    REAL    :: upops    = 0. !<daily(dayno,6)
    REAL    :: fllevel  = 0. !<daily(dayno,7)
    
    !daily smoothed percentile flow (type 0)
    REAL    :: drprcntl = 0. !<daily(dayno,8)
    REAL    :: minprcntl = 0. !<daily(dayno,9)
    REAL    :: midprcntl = 0. !<daily(dayno,10)
    REAL    :: maxprcntl = 0. !<daily(dayno,11)
    REAL    :: flprcntl = 0. !<daily(dayno,12)
    
    !ideal band Q (type 3)
    REAL    :: sabtod   = 0. !<average band storage today
    REAL    :: sabtom   = 0. !<average band storage tomorrow
    REAL    :: sabd     = 0. !<average band storage change
    REAL    :: percentbfm = 0. !<percent band flow multiplier, limited by maxqmult
    REAL    :: idealbandq = 0. !<flow required to achieve average WSL
    REAL    :: qin14day(14)  !<stores inflow from past 14 days. qin14day(1) = today's qin, qin14day(2) = yesterday's qin, etc.
    REAL    :: sabd14day(14) !<stores average band storage change from past 14 days. sabd14day(1) = today's sabd, sabd14day(2) = yesterday's sabd, etc.
           
    !weighted, aka monthly regression flow (type 1)
    REAL    :: weightedlow = 0. !<monthly regressed lower flow
    REAL    :: midbase  = 0. !<weighted middle lower base
    REAL    :: weightedmid = 0. !<monthly regressed middle flow
    REAL    :: highbase = 0. !<weighted high lower base
    REAL    :: weightedhigh = 0. !<monthly regressed high flow
    
    !physical, aka A-B-C-D flow (type 2)
    REAL    :: minphysical = 0. !<A-B-C-D minimum flow
    REAL    :: midphysical = 0. !<A-B-C-D middle flow
    REAL    :: maxphysical = 0. !<A-B-C-D maximum flow
    
    !conditioned flow (reservoir)
    REAL    :: wsl14day(15)  !<holds past 14 day's WSL for use by conditioned dams, size 15 is intentional
    REAL    :: qo14day(15)   !<holds past 14 day's Qo for use by conditioned dams,
                             !size 15 is intentional, taken from conditioning dam
    INTEGER :: orderswitch = 0 !<0 = conditioning dam calculated first, 1 = conditioned dam calculated first
                             !<0 = use ...14day(2:15), 1 = use ...14day(1:14)
    REAL    :: xx14dayavg  = 0. !<average of wsl14day or qo14day
    REAL    :: condflow    = 0. !<= A * (xx14dayavg - B) ^ C + D
    REAL    :: condop      = 0. !<condswitch = 0: opflow, condswitch = 1: condflow
    
    !selected outflows (based on type)
    REAL    :: drflow      = 0. !<0 = Dr. Percentile, 1 = Fixed Min.
    REAL    :: lowflow     = 0. !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: corlowflow  = 0. !<sets lowflow >= droughtflow
    REAL    :: lowtransflow = 0. !<interpolate between lowflow and opflow
    REAL    :: opflow      = 0. !<3 = Band Ideal,  0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: coropflow   = 0. !<sets opflow>=lowflow, <=highflow
    REAL    :: hightransflow = 0. !<interpolate between opflow and highflow
    REAL    :: highflow    = 0. !<0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: corhighflow = 0. !<sets highflow<=flflow
    REAL    :: flflow      = 0. !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    REAL    :: transflow   = 0. !<temporary variable that stores the linear interpolation between flows based on WSL
    
    !calculation of slope of last week's powerhouse flows
    REAL    :: sumx        = 0. !<sum of days
    REAL    :: sumy        = 0. !<sum of ph flows
    REAL    :: sumxx       = 0. !<sum of days squared
    REAL    :: sumxy       = 0. !<sum of days times ph flows
    REAL    :: weeklyslope = 0. !<slope of the powerhouse flows
    
    !natural variable declarations
    REAL    :: natopt1m3s  = 0. !<selected outflow
    REAL    :: natopt2m3s  = 0. !<outflow after powerhouse restriction
    REAL    :: natopt3m3s  = 0. !<outflow after weekly change restriction
    REAL    :: natopt4m3s  = 0. !<outflow after daily change restriction
    REAL    :: natminqa    = 0. !<where: Qo = A * (WSL - B) ^ C + D
    REAL    :: natminqb    = 0. !<
    REAL    :: natminqc    = 0. !<
    REAL    :: natminqd    = 0. !<
    REAL    :: natmidqa    = 0. !<
    REAL    :: natmidqb    = 0. !<
    REAL    :: natmidqc    = 0. !<
    REAL    :: natmidqd    = 0. !<
    REAL    :: natmaxqa    = 0. !<
    REAL    :: natmaxqb    = 0. !<
    REAL    :: natmaxqc    = 0. !<
    REAL    :: natmaxqd    = 0. !<
    REAL    :: natdeltaw   = 0. !<Weekly Flow Change Limit
    REAL    :: natdeltad   = 0. !<Daily Flow Change Limit
    REAL    :: natphqo     = 0. !<Powerhouse Limit Flow Value
    REAL    :: natdrqo     = 0. !<Minimum Fixed Flow Value
    REAL    :: natflqo     = 0. !<Maximum Fixed Flow Value
    INTEGER :: natdrtype   = 0  !<0 = Dr. Percentile, 1 = Fixed Min.
    INTEGER :: natmintype  = 0  !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natmidtype  = 0  !<0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natmaxtype  = 0  !<3 = Storage Change, 0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natfltype   = 0  !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    INTEGER :: natphtype   = 0  !<0 = No PH Restriction, 1 = Hold Water
    REAL    :: natweekly7day(7) !<stores powerhouse restricted flows from past 7 days. weekly7day(1) = today's flow, weekly7day(2) = yesterday's flow, etc.
    REAL    :: natdailyyes = 0. !yesterday's daily restricted flow       
    
    !daily smoothed percentile flow (type 0)
    REAL    :: natdrprc    = 0. !<daily(dayno,8)
    REAL    :: natminprc   = 0. !<daily(dayno,9)
    REAL    :: natmidprc   = 0. !<daily(dayno,10)
    REAL    :: natmaxprc   = 0. !<daily(dayno,11)
    REAL    :: natflprc    = 0. !<daily(dayno,12)
       
    !weighted, aka monthly regression flow (type 1)
    REAL    :: natweightedlow = 0. !<monthly regressed lower flow
    REAL    :: natmidbase  = 0. !<weighted middle lower base
    REAL    :: natweightedmid = 0. !<monthly regressed middle flow
    REAL    :: nathighbase = 0. !<weighted high lower base
    REAL    :: natweightedhigh = 0. !<monthly regressed high flow
    
    !physical, aka A-B-C-D flow (type 2)
    REAL    :: natminphysical = 0. !<A-B-C-D minimum flow
    REAL    :: natmidphysical = 0. !<A-B-C-D middle flow
    REAL    :: natmaxphysical = 0. !<A-B-C-D maximum flow
    
    !storage delta high flow type (type 3)   
    REAL    :: resqo7day(7)     !<holds past 7 days' reservoir Qo for use by Storage Delta High Flow Type
    REAL    :: qin7day(7)       !<holds past 7 days' Qin for use by Storage Delta High Flow Type
                                !<note: timing is different from qin14day(), DO NOT interchange them
    REAL    :: stordeltaavg = 0. !<average inflow 7 days - average reservoir outflow 7 days
    REAL    :: stordeltahft = 0. !<storage delta high flow type
    
    !selected outflows (based on type)
    REAL    :: natdrflow      = 0. !<0 = Dr. Percentile, 1 = Fixed Min.
    REAL    :: natlowflow     = 0. !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: natcorlowflow  = 0. !<sets lowflow >= droughtflow
    REAL    :: natlowtransflow = 0. !<interpolate between lowflow and opflow
    REAL    :: natopflow      = 0. !<0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: natcoropflow   = 0. !<sets opflow>=lowflow, <=highflow
    REAL    :: nathightransflow = 0. !<interpolate between opflow and highflow
    REAL    :: nathighflow    = 0. !<3 = Storage Change, 0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    REAL    :: natcorhighflow = 0. !<sets highflow<=flflow
    REAL    :: natflflow      = 0. !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    REAL    :: nattransflow   = 0. !<temporary variable that stores the linear interpolation between flows based on WSL
    
    !calculation of slope of last week's powerhouse flows
    REAL    :: natsumx        = 0. !<sum of days
    REAL    :: natsumy        = 0. !<sum of ph flows
    REAL    :: natsumxx       = 0. !<sum of days squared
    REAL    :: natsumxy       = 0. !<sum of days times ph flows
    REAL    :: natweeklyslope = 0. !<slope of the powerhouse flows
     
    
    !local variable assignment
    maxqmult = dam(damindex(i))%maxqmult
    stora = dam(damindex(i))%stora
    storb = dam(damindex(i))%storb
    storc = dam(damindex(i))%storc
    stord = dam(damindex(i))%stord
    natswitch = dam(damindex(i))%natswitch
    locswitch = dam(damindex(i))%locswitch
    branchswt = dam(damindex(i))%branchswt
    condswitch = dam(damindex(i))%condswitch
    condrestr = dam(damindex(i))%condrestr
    condqowsl = dam(damindex(i))%condqowsl
    condsubid = dam(damindex(i))%condsubid
    condindex = dam(damindex(i))%condindex
    conda = dam(damindex(i))%conda
    condb = dam(damindex(i))%condb
    condc = dam(damindex(i))%condc
    condd = dam(damindex(i))%condd
    minqa = dam(damindex(i))%minqa
    minqb = dam(damindex(i))%minqb
    minqc = dam(damindex(i))%minqc
    minqd = dam(damindex(i))%minqd
    midqa = dam(damindex(i))%midqa
    midqb = dam(damindex(i))%midqb
    midqc = dam(damindex(i))%midqc
    midqd = dam(damindex(i))%midqd
    maxqa = dam(damindex(i))%maxqa
    maxqb = dam(damindex(i))%maxqb
    maxqc = dam(damindex(i))%maxqc
    maxqd = dam(damindex(i))%maxqd
    natminqa = dam(damindex(i))%natminqa
    natminqb = dam(damindex(i))%natminqb
    natminqc = dam(damindex(i))%natminqc
    natminqd = dam(damindex(i))%natminqd
    natmidqa = dam(damindex(i))%natmidqa
    natmidqb = dam(damindex(i))%natmidqb
    natmidqc = dam(damindex(i))%natmidqc
    natmidqd = dam(damindex(i))%natmidqd
    natmaxqa = dam(damindex(i))%natmaxqa
    natmaxqb = dam(damindex(i))%natmaxqb
    natmaxqc = dam(damindex(i))%natmaxqc
    natmaxqd = dam(damindex(i))%natmaxqd
    deltaw = dam(damindex(i))%deltaw
    deltad = dam(damindex(i))%deltad
    natdeltaw = dam(damindex(i))%natdeltaw
    natdeltad = dam(damindex(i))%natdeltad
    phqo = dam(damindex(i))%phqo
    drqo = dam(damindex(i))%drqo
    flqo = dam(damindex(i))%flqo
    natphqo = dam(damindex(i))%natphqo
    natdrqo = dam(damindex(i))%natdrqo
    natflqo = dam(damindex(i))%natflqo
    drtype = dam(damindex(i))%drtype
    mintype = dam(damindex(i))%mintype
    midtype = dam(damindex(i))%midtype
    maxtype = dam(damindex(i))%maxtype
    fltype = dam(damindex(i))%fltype
    phtype = dam(damindex(i))%phtype
    natdrtype = dam(damindex(i))%natdrtype
    natmintype = dam(damindex(i))%natmintype
    natmidtype = dam(damindex(i))%natmidtype
    natmaxtype = dam(damindex(i))%natmaxtype
    natfltype = dam(damindex(i))%natfltype
    natphtype = dam(damindex(i))%natphtype
    monthly = dam(damindex(i))%monthly
    daily = dam(damindex(i))%daily
    idt = dam(damindex(i))%idt
    idtlast = dam(damindex(i))%idtlast
    curmonth = dam(damindex(i))%curmonth
    dayofmonth = dam(damindex(i))%dayofmonth
    mweight = daily(dayno,18)
    qin14day = dam(damindex(i))%qin14day
    sabd14day = dam(damindex(i))%sabd14day
    weekly7day = dam(damindex(i))%weekly7day
    natweekly7day = dam(damindex(i))%natweekly7day
    dailyyes = dam(damindex(i))%dailyyes
    natdailyyes = dam(damindex(i))%natdailyyes
    wsl14day = dam(damindex(i))%wsl14day
    qo14day = dam(damindex(i))%qo14day
    orderswitch = dam(damindex(i))%orderswitch
    resqo7day = dam(damindex(i))%resqo7day
    qin7day = dam(damindex(i))%qin7day
    
    !determine backwards and forwards month for monthly regression
    IF(curmonth==1.AND.dayofmonth<=15)THEN
      backmonth = 12
    ELSEIF(dayofmonth<=15)THEN
      backmonth = curmonth - 1
    ELSE
      backmonth = curmonth
    ENDIF
    
    IF(curmonth==12.AND.dayofmonth>15)THEN
      formonth = 1
    ELSEIF(dayofmonth>15)THEN
      formonth = curmonth + 1
    ELSE
      formonth = curmonth
    ENDIF
    
    !calculation of WSL
    wlmr0 = wlmr + dam(damindex(i))%w0ref !water level in local reference system
    !IF(wlmr0 < storb .OR. wlmr0 /= wlmr0) wlmr0 = storb + 1.0 !stopgap solution to large negative wlmr values
    dam(damindex(i))%storyes = stora * ((wlmr0 - storb) ** storc) + stord
    dam(damindex(i))%stortod = dam(damindex(i))%storyes + dam(damindex(i))%qinprev - dam(damindex(i))%qoutprev
    wsl = ((dam(damindex(i))%stortod - stord) / stora) ** (1/storc) + storb
    
    !read in daily WSL bands
    drlevel = daily(dayno,1)
    lowops = daily(dayno,2)
    lowband = daily(dayno,3)
    avglevel = daily(dayno,4)
    upband = daily(dayno,5)
    upops = daily(dayno,6)
    fllevel = daily(dayno,7)
    
    !read in smoothed percentile flows
    drprcntl = daily(dayno,8)
    minprcntl = daily(dayno,9)
    midprcntl = daily(dayno,10)
    maxprcntl = daily(dayno,11)
    flprcntl = daily(dayno,12)
    natdrprc = daily(dayno,13)
    natminprc = daily(dayno,14)
    natmidprc = daily(dayno,15)
    natmaxprc = daily(dayno,16)
    natflprc = daily(dayno,17)
    
    !determine ideal band Q
    updatesabd:IF(idt>idtlast)THEN !only updates sabd14day when called via calculate_outflow_from_lake()
      sabtod = stora * ((avglevel - storb) ** storc) + stord
      IF(dayno+1==367)THEN !goes to next year
        sabtom = stora * ((daily(1,4) - storb) ** storc) + stord
      ELSE
        sabtom = stora * ((daily(dayno+1,4) - storb) ** storc) + stord !daily(dayno+1,4) is tomorrow's avglevel
      ENDIF
      sabd = sabtom - sabtod
      DO counter = 13,1,-1
        sabd14day(counter+1) = sabd14day(counter)
      ENDDO
      sabd14day(1) = sabd
      dam(damindex(i))%sabd14day = sabd14day      
    ENDIF updatesabd
    
    calcpercent:IF(wsl<lowband)THEN
      percentbfm = -1*maxqmult
    ELSEIF(wsl>upband)THEN calcpercent
      percentbfm = maxqmult
    ELSEIF(wsl>avglevel)THEN calcpercent
      percentbfm = maxqmult * (wsl - avglevel) / (upband - avglevel)
    ELSE calcpercent
      percentbfm = -1*maxqmult * (avglevel - wsl) / (avglevel - lowband)
    ENDIF calcpercent
    
    idtlimit = idt
    IF(idtlimit>14) idtlimit=14
    idealbandq = (SUM(qin14day) / idtlimit) + ((1+percentbfm) * (SUM(sabd14day) / idtlimit))
    
    
    !calculate weighted flows (monthly regression)
    weightedlow = ((monthly(backmonth,1) * (wsl - monthly(backmonth,7))) * mweight) + ((monthly(formonth,1) * (wsl - monthly(formonth,7))) * (1-mweight))
    midbase  = ((monthly(backmonth,1) * (lowops - monthly(backmonth,7))) * mweight) + ((monthly(formonth,1) * (lowops - monthly(formonth,7))) * (1-mweight))
    weightedmid = ((monthly(backmonth,2) * (wsl - monthly(backmonth,8))) * mweight) + ((monthly(formonth,2) * (wsl - monthly(formonth,8))) * (1-mweight)) + midbase
    highbase = ((monthly(backmonth,2) * (upband - monthly(backmonth,8))) * mweight) + ((monthly(formonth,2) * (upband - monthly(formonth,8))) * (1-mweight))
    weightedhigh = ((monthly(backmonth,3) * (wsl - monthly(backmonth,9))) * mweight) + ((monthly(formonth,3) * (wsl - monthly(formonth,9))) * (1-mweight)) + highbase
    natweightedlow = ((monthly(backmonth,4) * (wsl - monthly(backmonth,7))) * mweight) + ((monthly(formonth,4) * (wsl - monthly(formonth,7))) * (1-mweight))
    natmidbase  = ((monthly(backmonth,4) * (lowops - monthly(backmonth,7))) * mweight) + ((monthly(formonth,4) * (lowops - monthly(formonth,7))) * (1-mweight))
    natweightedmid = ((monthly(backmonth,5) * (wsl - monthly(backmonth,8))) * mweight) + ((monthly(formonth,5) * (wsl - monthly(formonth,8))) * (1-mweight)) + natmidbase
    nathighbase = ((monthly(backmonth,5) * (upband - monthly(backmonth,8))) * mweight) + ((monthly(formonth,5) * (upband - monthly(formonth,8))) * (1-mweight))
    natweightedhigh = ((monthly(backmonth,6) * (wsl - monthly(backmonth,9))) * mweight) + ((monthly(formonth,6) * (wsl - monthly(formonth,9))) * (1-mweight)) + nathighbase
    
    !calculate physical flows (A-B-C-D)
    minphysical = minqa * ((wsl - minqb) ** minqc) + minqd
    midphysical = midqa * ((wsl - midqb) ** midqc) + midqd
    maxphysical = maxqa * ((wsl - maxqb) ** maxqc) + maxqd
    natminphysical = natminqa * ((wsl - natminqb) ** natminqc) + natminqd
    natmidphysical = natmidqa * ((wsl - natmidqb) ** natmidqc) + natmidqd
    natmaxphysical = natmaxqa * ((wsl - natmaxqb) ** natmaxqc) + natmaxqd
    
    !calculate conditioned flows
    idtlimit = idt - 1
    IF(idtlimit>=14) idtlimit=14
    IF(idtlimit==0)THEN !skip conditioning on first timestep
      condswitch = 0
    ELSE
      IF(condqowsl==0)THEN !use qo
        xx14dayavg = SUM(qo14day(1+orderswitch:idtlimit+orderswitch)) / idtlimit
      ELSEIF(condqowsl==1)THEN !use wsl
        xx14dayavg = SUM(wsl14day(1+orderswitch:idtlimit+orderswitch)) / idtlimit
      ENDIF !Qo = A * (XX - B) ^ C + D
      condflow = conda * ((xx14dayavg - condb) ** condc) + condd
    ENDIF
    
    !select outflows (based on type)
    !drought flow
    choosedr:IF(drtype==0)THEN
      drflow = drprcntl
    ELSEIF(drtype==1)THEN choosedr
      drflow = drqo
    ENDIF choosedr
    
    !minimum flow
    chooselow:IF(mintype==0)THEN
      lowflow = minprcntl
    ELSEIF(mintype==1)THEN chooselow
      lowflow = weightedlow
    ELSEIF(mintype==2)THEN chooselow
      lowflow = minphysical
    ENDIF chooselow
    
    !operations flow
    chooseop:IF(midtype==0)THEN
      opflow = midprcntl
    ELSEIF(midtype==1)THEN chooseop
      opflow = weightedmid
    ELSEIF(midtype==2)THEN chooseop
      opflow = midphysical
    ELSEIF(midtype==3)THEN chooseop
      opflow = idealbandq
    ENDIF chooseop
    
    !conditioned flow
    choosecond:IF(condswitch==0)THEN
      condop = opflow
    ELSEIF(condswitch==1)THEN choosecond
      condop = condflow
    ENDIF choosecond
    
    !maximum flow
    choosehigh:IF(maxtype==0)THEN
      highflow = maxprcntl
    ELSEIF(maxtype==1)THEN choosehigh
      highflow = weightedhigh
    ELSEIF(maxtype==2)THEN choosehigh
      highflow = maxphysical
    ENDIF choosehigh
    
    !flood flow
    choosefl:IF(fltype==0)THEN
      flflow = flprcntl
    ELSEIF(fltype==1)THEN choosefl
      flflow = flqo
    ELSEIF(fltype==2)THEN choosefl
      flflow = highflow
    ENDIF choosefl
    
    !correct flows
    
    !correct minimum flow
    correctlow:IF(lowflow < drflow)THEN
      corlowflow = drflow
    ELSE correctlow
      corlowflow = lowflow
    ENDIF correctlow
    
    !correct maximum flow
    correcthigh:IF(highflow > flflow)THEN
      corhighflow = flflow
    ELSE correcthigh
      corhighflow = highflow
    ENDIF correcthigh
    
    !correct operations flow (conditioned version)
    correctop:IF(condrestr==0)THEN
      IF(condop < corlowflow)THEN
        coropflow = corlowflow
      ELSEIF(condop > corhighflow)THEN
        coropflow = corhighflow
      ELSE
        coropflow = condop
      ENDIF
    ELSEIF(condrestr==1)THEN correctop
      IF(condop < drflow)THEN
        coropflow = drflow
      ELSEIF(condop > flflow)THEN
        coropflow = flflow
      ELSE
        coropflow = condop
      ENDIF
    ENDIF correctop
     
    !determine transition flows
    
    !determine low transition flow
    transflow = corlowflow + ((coropflow - corlowflow) * (wsl - lowops) / (lowband - lowops))
    lowtrans:IF(transflow < corlowflow)THEN
      lowtransflow = corlowflow
    ELSEIF(transflow > coropflow)THEN lowtrans
      lowtransflow = coropflow
    ELSE lowtrans
      lowtransflow = transflow
    ENDIF lowtrans
    
    !determine high transition flow
    transflow = coropflow + ((corhighflow - coropflow) * (wsl - upband) / (upops - upband))
    hightrans:IF(transflow < coropflow)THEN
      hightransflow = coropflow
    ELSEIF(transflow > corhighflow)THEN hightrans
      hightransflow = corhighflow
    ELSE hightrans
      hightransflow = transflow
    ENDIF hightrans
    
    !determine outflow before restrictions
    selectflow:IF(wsl > fllevel)THEN
      opt1m3s = flflow
    ELSEIF(wsl > upops)THEN selectflow
      opt1m3s = corhighflow
    ELSEIF(wsl > upband)THEN selectflow
      opt1m3s = hightransflow
    ELSEIF(wsl > lowband)THEN selectflow
      opt1m3s = coropflow
    ELSEIF(wsl > lowops)THEN selectflow
      opt1m3s = lowtransflow
    ELSEIF(wsl > drlevel)THEN selectflow
      opt1m3s = corlowflow
    ELSE selectflow
      opt1m3s = drflow
    ENDIF selectflow
    
    !determine powerhouse restriction
    phrestriction:IF(phtype==1)THEN
      IF(opt1m3s < phqo)THEN
        opt2m3s = opt1m3s
      ELSEIF(wsl > upops)THEN !use spillway?
        opt2m3s = opt1m3s
      ELSE
        opt2m3s = phqo
      ENDIF
    ELSE phrestriction
      opt2m3s = opt1m3s
    ENDIF phrestriction
    
    !determine weekly flow change restriction
    updateweekly:IF(idt>idtlast)THEN !only updates weekly7day when called via calculate_outflow_from_lake()
      DO counter = 6,1,-1
        weekly7day(counter+1) = weekly7day(counter)
      ENDDO
      weekly7day(1) = opt2m3s
      dam(damindex(i))%weekly7day = weekly7day
    ENDIF updateweekly
    
    idtlimit = idt
    IF(idtlimit>7) idtlimit=7
    
    !because weekly7day is in inverse order (recent to older flows)
    !decreasing flows would actually have a positive slope, and vice versa
    DO counter = 1,idtlimit
      sumx  = sumx + REAL(counter)
      sumy  = sumy + weekly7day(counter)
      sumxx = sumxx + (REAL(counter) * REAL(counter))
      sumxy = sumxy + (REAL(counter) * weekly7day(counter))
    ENDDO
    
    weeklyslope = (sumxy - (sumx * sumy / REAL(idtlimit))) / (sumxx - (sumx * sumx / REAL(idtlimit)))
    
    weeklyrestriction:IF(MAXVAL(weekly7day(1:idtlimit)) - MINVAL(weekly7day(1:idtlimit)) < deltaw)THEN
      opt3m3s = opt2m3s
    ELSEIF(weeklyslope < 0)THEN weeklyrestriction !original excel has (slope > 0)
      opt3m3s = MINVAL(weekly7day(1:idtlimit)) + deltaw
    ELSE weeklyrestriction
      opt3m3s = MAXVAL(weekly7day(1:idtlimit)) - deltaw
    ENDIF weeklyrestriction
    
    !determine daily flow change restriction
    dailyrestriction:IF(idt==1)THEN !only one timestep
      opt4m3s = opt3m3s
    ELSEIF(ABS(opt3m3s - dailyyes) < deltad)THEN dailyrestriction
      opt4m3s = opt3m3s
    ELSEIF(opt3m3s > dailyyes)THEN dailyrestriction
      opt4m3s = dailyyes + deltad
    ELSE dailyrestriction
      opt4m3s = dailyyes - deltad
    ENDIF dailyrestriction
    
    IF(idt == idtlast) dam(damindex(i))%dailyyes = opt4m3s !only updates dailyyes when called via check_outflow_from_lake()
    !update yesterday's flow at the end of the timestep, not the beginning
    
    
    
    
    
    
    !***********************************************
    !calculate natcorlowflow for use by Storage delta High Flow Type
    !drought flow
    natchoosedr:IF(natdrtype==0)THEN
      natdrflow = natdrprc
    ELSEIF(drtype==1)THEN natchoosedr
      natdrflow = natdrqo
    ENDIF natchoosedr
    
    !minimum flow
    natchooselow:IF(natmintype==0)THEN
      natlowflow = natminprc
    ELSEIF(natmintype==1)THEN natchooselow
      natlowflow = natweightedlow
    ELSEIF(natmintype==2)THEN natchooselow
      natlowflow = natminphysical
    ENDIF natchooselow
    
    !correct minimum flow
    natcorrectlow:IF(natlowflow < natdrflow)THEN
      natcorlowflow = natdrflow
    ELSE natcorrectlow
      natcorlowflow = natlowflow
    ENDIF natcorrectlow
    
    !calculate Storage delta High Flow Type   
    idtlimit = idt - 1
    sdhft:IF(idtlimit == 0)THEN !use defined value on first timestep
      stordeltahft = natdeltad
    ELSE sdhft
      IF(idtlimit > 7) idtlimit = 7
      stordeltaavg = (SUM(qin7day) / idtlimit) - (SUM(resqo7day) / idtlimit)
      IF(stordeltaavg < natcorlowflow)THEN
        stordeltahft = natcorlowflow
      ELSE
        stordeltahft = stordeltaavg
      ENDIF
    ENDIF sdhft
    
    
    
    
    !select outflows (based on type)
    
    
    !operations flow
    natchooseop:IF(natmidtype==0)THEN
      natopflow = natmidprc
    ELSEIF(natmidtype==1)THEN natchooseop
      natopflow = natweightedmid
    ELSEIF(midtype==2)THEN natchooseop
      natopflow = natmidphysical
    ENDIF natchooseop
    
    !maximum flow
    natchoosehigh:IF(natmaxtype==0)THEN
      nathighflow = natmaxprc
    ELSEIF(natmaxtype==1)THEN natchoosehigh
      nathighflow = natweightedhigh
    ELSEIF(natmaxtype==2)THEN natchoosehigh
      nathighflow = natmaxphysical
    ELSEIF(natmaxtype==3)THEN natchoosehigh
      nathighflow = stordeltahft
    ENDIF natchoosehigh
    
    !flood flow
    natchoosefl:IF(natfltype==0)THEN
      natflflow = natflprc
    ELSEIF(natfltype==1)THEN natchoosefl
      natflflow = natflqo
    ELSEIF(natfltype==2)THEN natchoosefl
      natflflow = nathighflow
    ENDIF natchoosefl
    
    !correct flows
    
    
    !correct maximum flow
    natcorrecthigh:IF(nathighflow > natflflow)THEN
      natcorhighflow = natflflow
    ELSE natcorrecthigh
      natcorhighflow = nathighflow
    ENDIF natcorrecthigh
    
    !correct operations flow
    natcorrectop:IF(natopflow < natcorlowflow)THEN
      natcoropflow = natcorlowflow
    ELSEIF(natopflow > natcorhighflow)THEN natcorrectop
      natcoropflow = natcorhighflow
    ELSE natcorrectop
      natcoropflow = natopflow
    ENDIF natcorrectop
    
    !determine transition flows
    
    !determine low transition flow
    transflow = natcorlowflow + ((natcoropflow - natcorlowflow) * (wsl - lowops) / (lowband - lowops))
    natlowtrans:IF(transflow < natcorlowflow)THEN
      natlowtransflow = natcorlowflow
    ELSEIF(transflow > natcoropflow)THEN natlowtrans
      natlowtransflow = natcoropflow
    ELSE natlowtrans
      natlowtransflow = transflow
    ENDIF natlowtrans
    
    !determine high transition flow
    transflow = natcoropflow + ((natcorhighflow - natcoropflow) * (wsl - upband) / (upops - upband))
    nathightrans:IF(transflow < natcoropflow)THEN
      nathightransflow = natcoropflow
    ELSEIF(transflow > natcorhighflow)THEN nathightrans
      nathightransflow = natcorhighflow
    ELSE nathightrans
      nathightransflow = transflow
    ENDIF nathightrans
    
    !determine outflow before restrictions
    natselectflow:IF(wsl > fllevel)THEN
      natopt1m3s = natflflow
    ELSEIF(wsl > upops)THEN natselectflow
      natopt1m3s = natcorhighflow
    ELSEIF(wsl > upband)THEN natselectflow
      natopt1m3s = nathightransflow
    ELSEIF(wsl > lowband)THEN natselectflow
      natopt1m3s = natcoropflow
    ELSEIF(wsl > lowops)THEN natselectflow
      natopt1m3s = natlowtransflow
    ELSEIF(wsl > drlevel)THEN natselectflow
      natopt1m3s = natcorlowflow
    ELSE natselectflow
      natopt1m3s = natdrflow
    ENDIF natselectflow
    
    !determine powerhouse restriction
    natphrestriction:IF(natphtype==1)THEN
      IF(natopt1m3s < natphqo)THEN
        natopt2m3s = natopt1m3s
      ELSEIF(wsl > upops)THEN !use spillway?
        natopt2m3s = natopt1m3s
      ELSE
        natopt2m3s = natphqo
      ENDIF
    ELSE natphrestriction
      natopt2m3s = natopt1m3s
    ENDIF natphrestriction
    
    !determine weekly flow change restriction
    natupdateweekly:IF(idt>idtlast)THEN !only updates weekly7day when called via calculate_outflow_from_lake()
      DO counter = 6,1,-1
        natweekly7day(counter+1) = natweekly7day(counter)
      ENDDO
      natweekly7day(1) = natopt2m3s
      dam(damindex(i))%natweekly7day = natweekly7day
    ENDIF natupdateweekly
    
    idtlimit = idt
    IF(idtlimit>7) idtlimit=7
    
    !because weekly7day is in inverse order (recent to older flows)
    !decreasing flows would actually have a positive slope, and vice versa
    DO counter = 1,idtlimit
      natsumx  = natsumx + REAL(counter)
      natsumy  = natsumy + natweekly7day(counter)
      natsumxx = natsumxx + (REAL(counter) * REAL(counter))
      natsumxy = natsumxy + (REAL(counter) * natweekly7day(counter))
    ENDDO
    
    natweeklyslope = (sumxy - (sumx * sumy / REAL(idtlimit))) / (sumxx - (sumx * sumx / REAL(idtlimit)))
    
    natweeklyrestriction:IF(MAXVAL(natweekly7day(1:idtlimit)) - MINVAL(natweekly7day(1:idtlimit)) < natdeltaw)THEN
      natopt3m3s = natopt2m3s
    ELSEIF(natweeklyslope < 0)THEN natweeklyrestriction !original excel has (slope > 0)
      natopt3m3s = MINVAL(natweekly7day(1:idtlimit)) + natdeltaw
    ELSE natweeklyrestriction
      natopt3m3s = MAXVAL(natweekly7day(1:idtlimit)) - natdeltaw
    ENDIF natweeklyrestriction
    
    !determine daily flow change restriction
    natdailyrestriction:IF(idt==1)THEN !only one timestep
      natopt4m3s = natopt3m3s
    ELSEIF(ABS(natopt3m3s - natdailyyes) < natdeltad)THEN natdailyrestriction
      natopt4m3s = natopt3m3s
    ELSEIF(natopt3m3s > natdailyyes)THEN natdailyrestriction
      natopt4m3s = natdailyyes + natdeltad
    ELSE natdailyrestriction
      natopt4m3s = natdailyyes - natdeltad
    ENDIF natdailyrestriction
    
    
    
    IF(idt == idtlast) dam(damindex(i))%natdailyyes = natopt4m3s !only updates dailyyes when called via check_outflow_from_lake()
    !update yesterday's flow at the end of the timestep, not the beginning
    
    IF(natswitch==0) natopt4m3s = 0 !natural flow is turned off
    
    
    !***************************************
    IF(idt > idtlast) dam(damindex(i))%idtlast = idt
    
    choosebranchmain:IF(branchswt==0)THEN !reservoir = mainflow, natural = branchflow
      damProd = opt4m3s
      branchflow = natopt4m3s
    ELSEIF(branchswt==1)THEN choosebranchmain !natural = mainflow, reservoir = branchflow
      damProd = natopt4m3s
      branchflow = opt4m3s
    ENDIF choosebranchmain
    
    
  END SUBROUTINE OCBD_dam_outflow

      

END MODULE SURFACEWATER_PROCESSES
