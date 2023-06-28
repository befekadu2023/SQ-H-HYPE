!> \file hypetypes.f90
!> Contains module statetype_module.

!>Types for state variables for the HYPE model (HYdrological Predictions for the Environment)
!>
!>Five types are distinguished; four for different hydrological features in the 
!>landscape (snow and ice, soil, river and lake), and one for miscellaneous states.
!>Not all state variables are states of the hydrological model, but hold memory of 
!>older values or averages that is needed if the model is restarted.
!>Each statetype has seven subroutines that handle the state variables as a group; 
!>allocation, deallocation, initiation (to zero), initiation for submodel simulation, 
!>get size of states/array, write to array and read from array.
!!
MODULE STATETYPE_MODULE

!Copyright 2013-2015 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
!Used by modules
!------------------------------------------------------------------------------
!"main"
!"optim"
!datamodule 
!modelmodule
!soilmodel_default
!glacier_soilmodel
!irrigation
!modelmodule
!soil_processes
!surfacewater_processes
!npc_soil_processes
!npc_surfacewater_processes
!regional_groundwater
!------------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

!(No) Parameter declarations

!Type declarations for model state variables
!> \brief Type for state variables related to snow and ice
  !
  ! For consistency, we should sort the content of the STATETYPES in the same order as they are read/written/allocated in the rest of the code
  !
  TYPE SNOWICESTATETYPE
    REAL,ALLOCATABLE :: snow(:,:)         !<snow on soil, snow water equivalent per soil-landuse combination  (class,subbasin)        !snowji
    REAL,ALLOCATABLE :: csnow(:,:,:)      !<concentration of snow on soil, snow water equivalent per soil-landuse combination  (substance,class,subbasin)        !csnowji
    REAL,ALLOCATABLE :: snowage(:,:)      !<age of snow pack on soil per land-soil combination (class,subbasin) !snowageji
    REAL,ALLOCATABLE :: snowcov(:,:)      !<degree of snow cover on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: snowmax(:,:)      !<max snow pack on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: glacvol(:)        !<glacier ice volume, one glacier class per subbasin (m3)      !glac_vol
    REAL,ALLOCATABLE :: riversnow(:,:)    !<snow on river, snow water equivalent per river type (type,subbasin)        !riversnowi
    REAL,ALLOCATABLE :: riversnowage(:,:) !<age of snow on riverper river type (type,subbasin)        !riversnowagei
    REAL,ALLOCATABLE :: riverice(:,:)     !<ice on river, total (cm) per river type (type,subbasin)        !rivericei
    REAL,ALLOCATABLE :: riverbice(:,:)    !<black ice on river (cm) per river type (type,subbasin)        !riverbicei
    REAL,ALLOCATABLE :: lakesnow(:,:)     !<snow on lake, snow water equivalent per lake type (type,subbasin) !!OR snow(slc_ilake,i)?        !lakesnowi
    REAL,ALLOCATABLE :: lakesnowage(:,:)  !<age of snow on lakeper lake type (type,subbasin) !!OR snow(slc_ilake,i)?        !lakesnowi
    REAL,ALLOCATABLE :: lakeice(:,:)      !<ice on lake, total (cm) per lake type (type,subbasin)        !lakeicei
    REAL,ALLOCATABLE :: lakebice(:,:)     !<black ice on lake (cm) per lake type (type,subbasin)        !lakebicei
    REAL,ALLOCATABLE :: rivericecov(:,:)  !<degree of ice cover on river (fraction) per river type (type,subbasin)
    REAL,ALLOCATABLE :: lakeicecov(:,:)   !<degree of ice cover on lake (fraction) per lake type (type,subbasin)
    !liquid water content in snow - for later applications with Remote Sensing data
    !REAL,ALLOCATABLE :: wcsnow(:,:)       !<snow liquid water content, mm (class, subbasin)
  END TYPE SNOWICESTATETYPE
!> \brief Type for state variables related to soil
  TYPE SOILSTATETYPE
    REAL,ALLOCATABLE :: water(:,:,:)  !<soil moisture per soil-landuse combination  (soillayer,class,subbasin)        !soilji,soil2ji,soil3ji
    REAL,ALLOCATABLE :: temp(:,:,:)   !<temperature of soil (degree Celcius) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: deeptemp(:,:) !<temperature of deep soil (degree Celcius) (class,subbasin)
    REAL,ALLOCATABLE :: conc(:,:,:,:) !<concentration of soil moisture per land-soil combination (substance,soillayer,rivertype,subbasin)  !csoilji,csoil2ji,csoil3ji
    REAL,ALLOCATABLE :: humusN(:,:,:) !<humusN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastN(:,:,:)  !<fastN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: partP(:,:,:)  !<partP pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastP(:,:,:)  !<fastP pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: humusP(:,:,:) !<humusN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastC(:,:,:)  !<fastC pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: humusC(:,:,:) !<humusC pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: PPrelpool(:,:)!<pool of delayed PP in runoff (mg/m2) (class,subbasin)  
    REAL,ALLOCATABLE :: oldgrw(:,:)   !<class ground water table yesterday (m) (class,subbasin)
    !LOGICAL,ALLOCATABLE:: icelens(:,:)!<boolean indicating presence of ice lens for infiltration into frozen soils (nclass,subbasin)
    INTEGER,ALLOCATABLE:: icelens(:,:)!<integer indicating presence of ice lens for infiltration into frozen soils (nclass,subbasin) (0 no ice lens, 1 icelens)
  END TYPE SOILSTATETYPE
!> \brief Type for state variables related to groundwater/aquifers
  TYPE AQUIFERSTATETYPE
    REAL,ALLOCATABLE :: water(:)  !<water volume of aquifer (m3 or m) (aquifer)
    REAL,ALLOCATABLE :: conc(:,:) !<concentration of aquifer water (substance,aquifer)
    REAL,ALLOCATABLE :: lastrecharge(:)    !<recharge flow last timestep (m3 or m) (aquifer)
    REAL,ALLOCATABLE :: clastrecharge(:,:) !<concentration of recharge flow last time step (substance,aquifer)
    REAL,ALLOCATABLE :: nextoutflow(:)     !<outflow of aquifer to be applied next time step (m3) (aquifer)
    REAL,ALLOCATABLE :: cnextoutflow(:,:)  !<concentration of outflow of aquifer (substance,aquifer)
  END TYPE AQUIFERSTATETYPE
!> \brief Type for state variables related to rivers
  TYPE RIVERSTATETYPE
    REAL,ALLOCATABLE :: water(:,:)     !<water volume of river (rivertype,subbasin) (m3)    !riverboxi
    REAL,ALLOCATABLE :: temp(:,:)   !<river temperature (rivertype,subbasin)                !rivertempi
    REAL,ALLOCATABLE :: conc(:,:,:)   !<concentration in river box (substance,rivertype,subbasin)   !criverboxi
    REAL,ALLOCATABLE :: TPmean(:,:)  !<365-day mean Total Phosphorus concentration of river (rivertype,subbasin)  !riverTPmean
    REAL,ALLOCATABLE :: temp10(:,:)    !<10-day mean river temperature (rivertype,subbasin)   !riverT10
    REAL,ALLOCATABLE :: temp20(:,:)    !<20-day mean river temperature (rivertype,subbasin)   !riverT20
    REAL,ALLOCATABLE :: Psed(:,:)     !<phosphorus sediment in river (kg) (rivertype,subbasin)  !riversed
    REAL,ALLOCATABLE :: qqueue(:,:,:)       !<q of water translated in river (m3) !riverq
    REAL,ALLOCATABLE :: cqueue(:,:,:,:)     !<c of water translated in river (mg/l) !criver
    REAL,ALLOCATABLE :: cwetland(:,:,:)  !<concentration of river wetland volume (mg/l)  !cwetlandvol
    REAL,ALLOCATABLE :: Qdayacc(:,:,:)      !<Q on time step<1 day, used in calculation of daily mean (m3/s) (timestep,rivertype,subbasin)    !Qdayacc
    REAL,ALLOCATABLE :: Q365(:,:,:)    !<Q river last 365 days (m3/s) (dayno,rivertype,subbasin)      !riverQ365
    REAL,ALLOCATABLE :: Qmean(:,:)   !<average discharge in river (365-days-MA) (m3/s) (rivertype,subbasin)   !Qmeani
  END TYPE RIVERSTATETYPE
!> \brief Type for state variables related to lake
  TYPE LAKESTATETYPE
    REAL,ALLOCATABLE :: water(:,:)      !<water stage of lake (mm) (laketype,subbasin)                             !lakewi
    REAL,ALLOCATABLE :: temp(:,:)       !<lake temperature (laketype,subbasin)                                     !laketempi
    REAL,ALLOCATABLE :: conc(:,:,:)     !<concentration in lakes (mg/L etc) (substances,laketype,subbasin)         !clakei
    REAL,ALLOCATABLE :: TPmean(:,:)     !<365-day mean Total Phosphorus concentration of lake (laketype,subbasin)  !lakeTPmean
    REAL,ALLOCATABLE :: slowwater(:,:)  !<volume or level of long turnover lake part (mm) (laketype,subbasin)             !slowlakei
    REAL,ALLOCATABLE :: concslow(:,:,:) !<concentrations of long turnover lake part (mg/l) (substance,laketype,subbasin)  !cslowlakei
    REAL,ALLOCATABLE :: temp10(:,:)     !<10-day mean lake temperature (laketype,subbasin)   !lakeT10
    REAL,ALLOCATABLE :: temp20(:,:)     !<20-day mean lake temperature (laketype,subbasin)   !lakeT20
    REAL,ALLOCATABLE :: uppertemp(:,:)  !<upper layer (epilimnion) temperature (laketype,subbasin)   !lakeuppertemp
    REAL,ALLOCATABLE :: lowertemp(:,:)  !<lower layer (hypolimnion) temperature (laketype,subbasin)  !deeplaketemp
    !MH2017: conditioning variables for use by OCBD dams
    INTEGER,ALLOCATABLE :: useswitch(:) !<0 = not used for conditioning, 1 = used for conditioning (subbasin)
    REAL,ALLOCATABLE :: wsl14day(:,:)  !<holds past 14 days' WSL for use by conditioned dams, size 15 is intentional (subbasin,day)
    REAL,ALLOCATABLE :: qo14day(:,:)   !<holds past 14 days' Qo for use by conditioned dams, size 15 is intentional (subbasin,day)
  END TYPE LAKESTATETYPE
!> \brief Type for miscellaneous state variables
  TYPE MISCSTATETYPE
    REAL,ALLOCATABLE :: temp5(:)      !<5-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp10(:)     !<10-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp20(:)     !<20-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp30(:)     !<30-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: gdd(:,:,:)    !<degree days for growing season start
    INTEGER, ALLOCATABLE :: gsbegin(:,:,:)      !<growing seasons start dayno
    REAL,ALLOCATABLE :: nextirrigation(:,:)     !<irrigation water to be applied next timestep (mm) (class,subbasin)
    REAL,ALLOCATABLE :: cnextirrigation(:,:,:)  !<concentration of the irrigation water 
    REAL,ALLOCATABLE :: updatestationsarcorr(:) !<AR-error of each subbasin (update function qar)      
  END TYPE MISCSTATETYPE

  PUBLIC SNOWICESTATETYPE
  PUBLIC SOILSTATETYPE
  PUBLIC AQUIFERSTATETYPE
  PUBLIC RIVERSTATETYPE
  PUBLIC LAKESTATETYPE
  PUBLIC MISCSTATETYPE
  PUBLIC allocate_model_states                !Eventuellt ska dessa flyttas ut i HYSS, i en data-module med starttillstånd-hantering
  PUBLIC deallocate_model_states
  PUBLIC initiate_state_zero
  PUBLIC initiate_frozenstate_submodel
  PUBLIC initiate_soilstate_submodel
  PUBLIC initiate_aquiferstate_submodel
  PUBLIC initiate_riverstate_submodel
  PUBLIC initiate_lakestate_submodel
  PUBLIC initiate_miscstate_submodel
  PUBLIC get_frozenstate_variables_arraysize
  PUBLIC set_frozenstate_variables_to_array
  PUBLIC set_frozenstate_variables_from_array
  PUBLIC get_soilstate_variables_arraysize
  PUBLIC set_soilstate_variables_to_array
  PUBLIC set_soilstate_variables_from_array
  PUBLIC get_aquiferstate_variables_arraysize
  PUBLIC set_aquiferstate_variables_to_array
  PUBLIC set_aquiferstate_variables_from_array
  PUBLIC get_riverstate_variables_arraysize
  PUBLIC set_riverstate_variables_to_array
  PUBLIC set_riverstate_variables_from_array
  PUBLIC get_lakestate_variables_arraysize
  PUBLIC set_lakestate_variables_to_array
  PUBLIC set_lakestate_variables_from_array
  PUBLIC get_miscstate_variables_arraysize
  PUBLIC set_miscstate_variables_to_array
  PUBLIC set_miscstate_variables_from_array

CONTAINS

  !>Allocate state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_model_states(n,ns,nc,na,nsl,nr,nl,ml,mt,  &
                                   isN,isP,isC,ist2,iswetl,isglac,islrice,isirr,isqar,isgsm,  &
                                   frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    INTEGER,INTENT(IN) :: ml      !<Maximum river lag steps
    INTEGER,INTENT(IN) :: mt      !<timestep per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Soil states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL allocate_frozenstate_variables(n,ns,nc,nr,nl,isglac,islrice,frozenstate)
    CALL allocate_soilstate_variables(n,ns,nc,nsl,isN,isP,isC,soilstate)
    CALL allocate_aquiferstate_variables(na,ns,aquiferstate)
    CALL allocate_riverstate_variables(n,ns,nr,ml,mt,isN,isP,iswetl,riverstate)
    CALL allocate_lakestate_variables(n,ns,nl,ist2,lakestate)
    CALL allocate_miscstate_variables(n,ns,nc,isirr,isqar,isgsm,miscstate)

  END SUBROUTINE allocate_model_states

  !---------------------------------------------------------------
  !>Deallocate model states 
  !---------------------------------------------------------------
  SUBROUTINE deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL deallocate_frozenstate_variables(frozenstate)
    CALL deallocate_soilstate_variables(soilstate)
    CALL deallocate_aquiferstate_variables(aquiferstate)
    CALL deallocate_riverstate_variables(riverstate)
    CALL deallocate_lakestate_variables(lakestate)
    CALL deallocate_miscstate_variables(miscstate)

  END SUBROUTINE deallocate_model_states

  !---------------------------------------------------------------
  !>Initiate states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_state_zero(ns,na,ist2,iswetl,isglac,islrice,isirr,isqar,isgsm,   &
                                     frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL initiate_frozenstate_zero(ns,isglac,islrice,frozenstate)
    CALL initiate_soilstate_zero(ns,soilstate)
    CALL initiate_aquiferstate_zero(ns,na,aquiferstate)
    CALL initiate_riverstate_zero(ns,iswetl,riverstate)
    CALL initiate_lakestate_zero(ns,ist2,lakestate)
    CALL initiate_miscstate_zero(ns,isirr,isqar,isgsm,miscstate)
  
  END SUBROUTINE initiate_state_zero

  !>Allocate snow and ice state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_frozenstate_variables(n,ns,nc,nr,nl,isglac,islrice,frozenstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states

    IF(.NOT.ALLOCATED(frozenstate%snow))      ALLOCATE(frozenstate%snow(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowage))   ALLOCATE(frozenstate%snowage(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowcov))   ALLOCATE(frozenstate%snowcov(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowmax))   ALLOCATE(frozenstate%snowmax(nc,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(frozenstate%csnow))   ALLOCATE(frozenstate%csnow(ns,nc,n))
    ENDIF  
    IF(isglac)THEN
      IF(.NOT.ALLOCATED(frozenstate%glacvol))   ALLOCATE(frozenstate%glacvol(n))
    ENDIF
    IF(islrice)THEN
      IF(.NOT.ALLOCATED(frozenstate%lakesnow))     ALLOCATE(frozenstate%lakesnow(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakesnowage))  ALLOCATE(frozenstate%lakesnowage(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakeice))      ALLOCATE(frozenstate%lakeice(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakebice))     ALLOCATE(frozenstate%lakebice(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakeicecov))   ALLOCATE(frozenstate%lakeicecov(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%riversnow))    ALLOCATE(frozenstate%riversnow(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riversnowage)) ALLOCATE(frozenstate%riversnowage(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riverice))     ALLOCATE(frozenstate%riverice(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riverbice))    ALLOCATE(frozenstate%riverbice(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%rivericecov))  ALLOCATE(frozenstate%rivericecov(nr,n))
    ENDIF  
    !IF(.NOT.ALLOCATED(frozenstate%wcsnow))      ALLOCATE(frozenstate%wcsnow(nc,n))
    
  END SUBROUTINE allocate_frozenstate_variables

  !>Deallocate snow and ice states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_frozenstate_variables(frozenstate)

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<River states

    IF(ALLOCATED(frozenstate%snow))    DEALLOCATE(frozenstate%snow)
    IF(ALLOCATED(frozenstate%csnow))   DEALLOCATE(frozenstate%csnow)
    IF(ALLOCATED(frozenstate%snowage)) DEALLOCATE(frozenstate%snowage)
    IF(ALLOCATED(frozenstate%snowcov)) DEALLOCATE(frozenstate%snowcov)
    IF(ALLOCATED(frozenstate%snowmax)) DEALLOCATE(frozenstate%snowmax)
    IF(ALLOCATED(frozenstate%glacvol)) DEALLOCATE(frozenstate%glacvol)
    IF(ALLOCATED(frozenstate%lakesnow))  DEALLOCATE(frozenstate%lakesnow)
    IF(ALLOCATED(frozenstate%lakesnowage)) DEALLOCATE(frozenstate%lakesnowage)
    IF(ALLOCATED(frozenstate%lakeice))   DEALLOCATE(frozenstate%lakeice)
    IF(ALLOCATED(frozenstate%lakebice))  DEALLOCATE(frozenstate%lakebice)
    IF(ALLOCATED(frozenstate%lakeicecov))  DEALLOCATE(frozenstate%lakeicecov)
    IF(ALLOCATED(frozenstate%riversnow)) DEALLOCATE(frozenstate%riversnow)
    IF(ALLOCATED(frozenstate%riversnowage)) DEALLOCATE(frozenstate%riversnowage)
    IF(ALLOCATED(frozenstate%riverice))  DEALLOCATE(frozenstate%riverice)
    IF(ALLOCATED(frozenstate%riverbice)) DEALLOCATE(frozenstate%riverbice)
    IF(ALLOCATED(frozenstate%rivericecov))  DEALLOCATE(frozenstate%rivericecov)
   ! IF(ALLOCATED(frozenstate%wcsnow))    DEALLOCATE(frozenstate%wcsnow)
    
  END SUBROUTINE deallocate_frozenstate_variables

  !>Initiate submodel snow and ice states by temporary storage of initial state
  !---------------------------------------------------------------
  SUBROUTINE initiate_frozenstate_submodel(ns,n,indexarray,isglac,islrice,frozenstate,frozenstate2)
  
    INTEGER, INTENT(IN) :: ns             !<number of substances
    INTEGER, INTENT(IN) :: n              !<number of subbasins in states to be set
    INTEGER, INTENT(IN) :: indexarray(n)  !<index for basemodel
    LOGICAL,INTENT(IN)  :: isglac         !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice        !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice state, submodel
    TYPE(snowicestatetype),INTENT(IN) :: frozenstate2    !<Snow and ice state, whole model setup

    frozenstate%snow(:,:) = frozenstate2%snow(:,indexarray(:))
    frozenstate%snowage(:,:) = frozenstate2%snowage(:,indexarray(:))
    frozenstate%snowcov(:,:) = frozenstate2%snowcov(:,indexarray(:))
    frozenstate%snowmax(:,:) = frozenstate2%snowmax(:,indexarray(:))
    IF(ns>0)THEN
      frozenstate%csnow(:,:,:) = frozenstate2%csnow(:,:,indexarray(:))
    ENDIF
    IF(isglac) frozenstate%glacvol(:) = frozenstate2%glacvol(indexarray(:))
    IF(islrice)THEN
      frozenstate%lakesnow(:,:)     = frozenstate2%lakesnow(:,indexarray(:))
      frozenstate%lakesnowage(:,:)  = frozenstate2%lakesnowage(:,indexarray(:))
      frozenstate%lakeice(:,:)      = frozenstate2%lakeice(:,indexarray(:))
      frozenstate%lakebice(:,:)     = frozenstate2%lakebice(:,indexarray(:))
      frozenstate%lakeicecov(:,:)   = frozenstate2%lakeicecov(:,indexarray(:))
      frozenstate%riversnow(:,:)    = frozenstate2%riversnow(:,indexarray(:))
      frozenstate%riversnowage(:,:) = frozenstate2%riversnowage(:,indexarray(:))
      frozenstate%riverice(:,:)     = frozenstate2%riverice(:,indexarray(:))
      frozenstate%riverbice(:,:)    = frozenstate2%riverbice(:,indexarray(:))
      frozenstate%rivericecov(:,:)  = frozenstate2%rivericecov(:,indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_frozenstate_submodel

  !>Initiate snow and ice states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_frozenstate_zero(ns,isglac,islrice,frozenstate)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    LOGICAL,INTENT(IN)  :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice   !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states

    frozenstate%snow = 0.
    frozenstate%snowage = 0.
    frozenstate%snowcov = 0.
    frozenstate%snowmax = 0.
    IF(ns>0)THEN
      frozenstate%csnow = 0.
    ENDIF 
    IF(isglac) frozenstate%glacvol = 0.
    IF(islrice)THEN
      frozenstate%lakesnow(:,:)     = 0.
      frozenstate%lakesnowage(:,:)  = 0.
      frozenstate%lakeice(:,:)      = 0.
      frozenstate%lakebice(:,:)     = 0.
      frozenstate%lakeicecov(:,:)   = 0.
      frozenstate%riversnow(:,:)    = 0.
      frozenstate%riversnowage(:,:) = 0.
      frozenstate%riverice(:,:)     = 0.
      frozenstate%riverbice(:,:)    = 0.
      frozenstate%rivericecov(:,:)  = 0.
    ENDIF
  
  END SUBROUTINE initiate_frozenstate_zero

  !>Allocate soil state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_soilstate_variables(n,ns,nc,nsl,isN,isP,isC,soilstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    IF(.NOT.ALLOCATED(soilstate%water))   ALLOCATE(soilstate%water(nsl,nc,n))
    IF(.NOT.ALLOCATED(soilstate%temp))    ALLOCATE(soilstate%temp(nsl,nc,n))
    IF(.NOT.ALLOCATED(soilstate%deeptemp))ALLOCATE(soilstate%deeptemp(nc,n))
    IF(.NOT.ALLOCATED(soilstate%icelens))ALLOCATE(soilstate%icelens(nc,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(soilstate%conc))      ALLOCATE(soilstate%conc(ns,nsl,nc,n))
      IF(.NOT.ALLOCATED(soilstate%PPrelpool)) ALLOCATE(soilstate%PPrelpool(nc,n))
      IF(isN)THEN
        IF(.NOT.ALLOCATED(soilstate%humusN))  ALLOCATE(soilstate%humusN(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastN))   ALLOCATE(soilstate%fastN(nsl,nc,n))
      ENDIF
      IF(isP)THEN
        IF(.NOT.ALLOCATED(soilstate%humusP))  ALLOCATE(soilstate%humusP(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastP))   ALLOCATE(soilstate%fastP(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%partP))   ALLOCATE(soilstate%partP(nsl,nc,n))
      ENDIF
      IF(isC)THEN
        IF(.NOT.ALLOCATED(soilstate%humusC))  ALLOCATE(soilstate%humusC(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastC))   ALLOCATE(soilstate%fastC(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%oldgrw))  ALLOCATE(soilstate%oldgrw(nc,n))
      ENDIF
    ENDIF

  END SUBROUTINE allocate_soilstate_variables

  !>Deallocate soil states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_soilstate_variables(soilstate)

    !Argument declarations
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    IF(ALLOCATED(soilstate%water))    DEALLOCATE(soilstate%water)
    IF(ALLOCATED(soilstate%temp))     DEALLOCATE(soilstate%temp)
    IF(ALLOCATED(soilstate%deeptemp)) DEALLOCATE(soilstate%deeptemp)
    IF(ALLOCATED(soilstate%icelens))  DEALLOCATE(soilstate%icelens)
    IF(ALLOCATED(soilstate%conc))     DEALLOCATE(soilstate%conc)
    IF(ALLOCATED(soilstate%PPrelpool))DEALLOCATE(soilstate%PPrelpool)
    IF(ALLOCATED(soilstate%humusN))   DEALLOCATE(soilstate%humusN)
    IF(ALLOCATED(soilstate%fastN))    DEALLOCATE(soilstate%fastN)
    IF(ALLOCATED(soilstate%humusP))   DEALLOCATE(soilstate%humusP)
    IF(ALLOCATED(soilstate%fastP))    DEALLOCATE(soilstate%fastP)
    IF(ALLOCATED(soilstate%partP))    DEALLOCATE(soilstate%partP)
    IF(ALLOCATED(soilstate%humusC))   DEALLOCATE(soilstate%humusC)
    IF(ALLOCATED(soilstate%fastC))    DEALLOCATE(soilstate%fastC)
    IF(ALLOCATED(soilstate%oldgrw))   DEALLOCATE(soilstate%oldgrw)

  END SUBROUTINE deallocate_soilstate_variables

  !>Initiate submodel soil states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_soilstate_submodel(ns,n,indexarray,soilstate,soilstate2)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: n         !<number of subbasins in states to be set
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil state, submodel
    TYPE(soilstatetype),INTENT(IN) :: soilstate2    !<Soil state, whole model setup

    soilstate%water(:,:,:)  = soilstate2%water(:,:,indexarray(:))
    soilstate%temp(:,:,:)   = soilstate2%temp(:,:,indexarray(:))
    soilstate%deeptemp(:,:) = soilstate2%deeptemp(:,indexarray(:))
    soilstate%icelens(:,:)  = soilstate2%icelens(:,indexarray(:))
    IF(ns>0)THEN
      soilstate%conc(:,:,:,:) = soilstate2%conc(:,:,:,indexarray(:))
      soilstate%PPrelpool(:,:) = soilstate2%PPrelpool(:,indexarray(:))
      IF(ALLOCATED(soilstate%humusN))THEN
        soilstate%humusN(:,:,:) = soilstate2%humusN(:,:,indexarray(:))
        soilstate%fastN(:,:,:)  = soilstate2%fastN(:,:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%humusP))THEN
        soilstate%humusP(:,:,:) = soilstate2%humusP(:,:,indexarray(:))
        soilstate%fastP(:,:,:)  = soilstate2%fastP(:,:,indexarray(:))
        soilstate%partP(:,:,:)  = soilstate2%partP(:,:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%humusC))THEN
        soilstate%humusC(:,:,:) = soilstate2%humusC(:,:,indexarray(:))
        soilstate%fastC(:,:,:)  = soilstate2%fastC(:,:,indexarray(:))
        soilstate%oldgrw(:,:)   = soilstate2%oldgrw(:,indexarray(:))
      ENDIF  
    ENDIF
  
  END SUBROUTINE initiate_soilstate_submodel

  !>Initiate soil states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_soilstate_zero(ns,soilstate)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    soilstate%water = 0.
    soilstate%temp = 0.
    soilstate%deeptemp = 0.
    soilstate%icelens = 0 !.FALSE.
    IF(ns>0)THEN
      soilstate%conc = 0.
      soilstate%PPrelpool = 0.
      IF(ALLOCATED(soilstate%humusN))THEN
        soilstate%humusN(:,:,:) = 0.
        soilstate%fastN(:,:,:)  = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%humusP))THEN
        soilstate%humusP(:,:,:) = 0.
        soilstate%fastP(:,:,:)  = 0.
        soilstate%partP(:,:,:)  = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%humusC))THEN
        soilstate%humusC(:,:,:) = 0.
        soilstate%fastC(:,:,:)  = 0.
        soilstate%oldgrw(:,:)   = 0.
      ENDIF  
    ENDIF
  
  END SUBROUTINE initiate_soilstate_zero

  !>Allocate aquifer state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_aquiferstate_variables(na,ns,aquiferstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(na>0)THEN
      IF(.NOT.ALLOCATED(aquiferstate%water))   ALLOCATE(aquiferstate%water(na))
      IF(.NOT.ALLOCATED(aquiferstate%lastrecharge))   ALLOCATE(aquiferstate%lastrecharge(na))
      IF(.NOT.ALLOCATED(aquiferstate%nextoutflow))   ALLOCATE(aquiferstate%nextoutflow(na))
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(aquiferstate%conc))  ALLOCATE(aquiferstate%conc(ns,na))
        IF(.NOT.ALLOCATED(aquiferstate%clastrecharge))  ALLOCATE(aquiferstate%clastrecharge(ns,na))
        IF(.NOT.ALLOCATED(aquiferstate%cnextoutflow))  ALLOCATE(aquiferstate%cnextoutflow(ns,na))
      ENDIF
    ENDIF

  END SUBROUTINE allocate_aquiferstate_variables

  !>Deallocate aquifer states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_aquiferstate_variables(aquiferstate)

    !Argument declarations
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(ALLOCATED(aquiferstate%water))         DEALLOCATE(aquiferstate%water)
    IF(ALLOCATED(aquiferstate%lastrecharge))  DEALLOCATE(aquiferstate%lastrecharge)
    IF(ALLOCATED(aquiferstate%nextoutflow))   DEALLOCATE(aquiferstate%nextoutflow)
    IF(ALLOCATED(aquiferstate%conc))          DEALLOCATE(aquiferstate%conc)
    IF(ALLOCATED(aquiferstate%clastrecharge)) DEALLOCATE(aquiferstate%clastrecharge)
    IF(ALLOCATED(aquiferstate%cnextoutflow))  DEALLOCATE(aquiferstate%cnextoutflow)

  END SUBROUTINE deallocate_aquiferstate_variables

  !>Initiate submodel aquifer states by temporary storage of initial
  !>state (copy)
  !---------------------------------------------------------------
  SUBROUTINE initiate_aquiferstate_submodel(ns,na,n,indexarray,aquiferstate,aquiferstate2)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: na        !<Number of aquifers
    INTEGER, INTENT(IN) :: n         !<number of subbasins in states to be set
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer state, submodel
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate2    !<aquifer state, whole model setup

    IF(na>0.)THEN
      aquiferstate%water(:)        = aquiferstate2%water(:)
      aquiferstate%lastrecharge(:) = aquiferstate2%lastrecharge(:)
      aquiferstate%nextoutflow(:)  = aquiferstate2%nextoutflow(:)
      IF(ns>0)THEN
        aquiferstate%conc(:,:)          = aquiferstate2%conc(:,:)
        aquiferstate%clastrecharge(:,:) = aquiferstate2%clastrecharge(:,:)
        aquiferstate%cnextoutflow(:,:)  = aquiferstate2%cnextoutflow(:,:)
      ENDIF
    ENDIF
  
  END SUBROUTINE initiate_aquiferstate_submodel

  !>Initiate aquifer states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_aquiferstate_zero(ns,na,aquiferstate)
  
    INTEGER, INTENT(IN) :: ns      !<number of substances
    INTEGER, INTENT(IN) :: na      !<number of aquifers
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(na>0)THEN
      aquiferstate%water = 0.
      aquiferstate%lastrecharge = 0.
      aquiferstate%nextoutflow = 0.
      IF(ns>0)THEN
        aquiferstate%conc = 0.
        aquiferstate%clastrecharge = 0.
        aquiferstate%cnextoutflow = 0.
      ENDIF
    ENDIF
  
  END SUBROUTINE initiate_aquiferstate_zero

  !>Allocate river state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_riverstate_variables(n,ns,nr,ml,mt,isN,isP,iswetl,riverstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nr      !<number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum river lag steps
    INTEGER,INTENT(IN) :: mt      !<timestep per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    IF(.NOT.ALLOCATED(riverstate%water))  ALLOCATE(riverstate%water(nr,n))
    IF(.NOT.ALLOCATED(riverstate%temp))   ALLOCATE(riverstate%temp(nr,n))
    IF(.NOT.ALLOCATED(riverstate%qqueue)) ALLOCATE(riverstate%qqueue(0:ml,nr,n))
    IF(.NOT.ALLOCATED(riverstate%Qmean))  ALLOCATE(riverstate%Qmean(nr,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(riverstate%TPmean)) ALLOCATE(riverstate%TPmean(nr,n))
      IF(.NOT.ALLOCATED(riverstate%Psed))   ALLOCATE(riverstate%Psed(nr,n))
      IF(.NOT.ALLOCATED(riverstate%conc))   ALLOCATE(riverstate%conc(ns,nr,n))
      IF(.NOT.ALLOCATED(riverstate%temp10))   ALLOCATE(riverstate%temp10(nr,n))
      IF(.NOT.ALLOCATED(riverstate%temp20))   ALLOCATE(riverstate%temp20(nr,n))
      IF(.NOT.ALLOCATED(riverstate%cqueue)) ALLOCATE(riverstate%cqueue(ns,0:ml,nr,n))
    ENDIF  
    IF(isN.OR.isP)THEN
      IF(.NOT.ALLOCATED(riverstate%Qdayacc)) ALLOCATE(riverstate%Qdayacc(mt,nr,n))
      IF(.NOT.ALLOCATED(riverstate%Q365)) ALLOCATE(riverstate%Q365(366,nr,n))
    ENDIF
    IF(iswetl)THEN
      IF(.NOT.ALLOCATED(riverstate%cwetland)) ALLOCATE(riverstate%cwetland(ns,nr,n))
    ENDIF

  END SUBROUTINE allocate_riverstate_variables

  !>Deallocate river states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_riverstate_variables(riverstate)

    !Argument declarations
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    IF(ALLOCATED(riverstate%water))    DEALLOCATE(riverstate%water)
    IF(ALLOCATED(riverstate%temp))     DEALLOCATE(riverstate%temp)
    IF(ALLOCATED(riverstate%qqueue))   DEALLOCATE(riverstate%qqueue)
    IF(ALLOCATED(riverstate%Qmean))    DEALLOCATE(riverstate%Qmean)
    IF(ALLOCATED(riverstate%TPmean))   DEALLOCATE(riverstate%TPmean)
    IF(ALLOCATED(riverstate%Psed))     DEALLOCATE(riverstate%Psed)
    IF(ALLOCATED(riverstate%conc))     DEALLOCATE(riverstate%conc)
    IF(ALLOCATED(riverstate%temp10))   DEALLOCATE(riverstate%temp10)
    IF(ALLOCATED(riverstate%temp20))   DEALLOCATE(riverstate%temp20)
    IF(ALLOCATED(riverstate%cqueue))   DEALLOCATE(riverstate%cqueue)
    IF(ALLOCATED(riverstate%Qdayacc))  DEALLOCATE(riverstate%Qdayacc)
    IF(ALLOCATED(riverstate%Q365))     DEALLOCATE(riverstate%Q365)
    IF(ALLOCATED(riverstate%cwetland)) DEALLOCATE(riverstate%cwetland)

  END SUBROUTINE deallocate_riverstate_variables

  !>Initiate submodel river states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_riverstate_submodel(ns,n,isN,isP,iswetl,indexarray,  &
                                          riverstate,riverstate2)
  
    INTEGER,INTENT(IN) :: ns        !<number of substances
    INTEGER,INTENT(IN) :: n         !<number of subbasins in states to be set
    LOGICAL,INTENT(IN) :: isN       !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP       !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: iswetl    !<Status of wetland simulation
    INTEGER,INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state, submodel
    TYPE(riverstatetype),INTENT(IN) :: riverstate2    !<River state, whole model setup

    riverstate%water(:,:) = riverstate2%water(:,indexarray(:))
    riverstate%temp(:,:)  = riverstate2%temp(:,indexarray(:))
    riverstate%qqueue(:,:,:) = riverstate2%qqueue(:,:,indexarray(:))
    riverstate%Qmean(:,:) = riverstate2%Qmean(:,indexarray(:))
    IF(ns>0)THEN
      riverstate%TPmean(:,:) = riverstate2%TPmean(:,indexarray(:))
      riverstate%Psed(:,:)   = riverstate2%Psed(:,indexarray(:))
      riverstate%conc(:,:,:) = riverstate2%conc(:,:,indexarray(:))
      riverstate%temp10(:,:) = riverstate2%temp10(:,indexarray(:))
      riverstate%temp20(:,:) = riverstate2%temp20(:,indexarray(:))
      riverstate%cqueue(:,:,:,:) = riverstate2%cqueue(:,:,:,indexarray(:))
    ENDIF
    IF(isN.OR.isP)THEN
      riverstate%Qdayacc(:,:,:) = riverstate2%Qdayacc(:,:,indexarray(:))
      riverstate%Q365(:,:,:) = riverstate2%Q365(:,:,indexarray(:))
    ENDIF
    IF(iswetl)THEN
      riverstate%cwetland(:,:,:) = riverstate2%cwetland(:,:,indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_riverstate_submodel

  !>Initiate river states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_riverstate_zero(ns,iswetl,riverstate)
  
    INTEGER,INTENT(IN) :: ns        !<number of substances
    LOGICAL,INTENT(IN) :: iswetl    !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    riverstate%water = 0.
    riverstate%temp = 0.
    riverstate%qqueue = 0.
    riverstate%Qmean = 0.
    IF(ns>0)THEN
      riverstate%TPmean = 0.
      riverstate%Psed = 0.
      riverstate%conc = 0.
      riverstate%temp10 = 0.
      riverstate%temp20 = 0.
      riverstate%cqueue = 0.
    ENDIF
!    IF(isN.OR.isP)THEN
!      riverstate%Qdayacc = 0.   !initialized to 0.0001 later in HYPE
!      riverstate%Q365 = 0.   !initialized to 0.0001 later in HYPE
!    ENDIF
    IF(iswetl) riverstate%cwetland = 0.
  
  END SUBROUTINE initiate_riverstate_zero

  !>Allocate lake state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_lakestate_variables(n,ns,nl,ist2,lakestate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nl      !<number of lake types
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states

    IF(.NOT.ALLOCATED(lakestate%water)) ALLOCATE(lakestate%water(nl,n))
    IF(.NOT.ALLOCATED(lakestate%temp))  ALLOCATE(lakestate%temp(nl,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(lakestate%TPmean))    ALLOCATE(lakestate%TPmean(nl,n))
      IF(.NOT.ALLOCATED(lakestate%conc))      ALLOCATE(lakestate%conc(ns,nl,n))
      IF(.NOT.ALLOCATED(lakestate%slowwater)) ALLOCATE(lakestate%slowwater(nl,n))
      IF(.NOT.ALLOCATED(lakestate%concslow))  ALLOCATE(lakestate%concslow(ns,nl,n))
      IF(.NOT.ALLOCATED(lakestate%temp10))    ALLOCATE(lakestate%temp10(nl,n))
      IF(.NOT.ALLOCATED(lakestate%temp20))    ALLOCATE(lakestate%temp20(nl,n))
    ENDIF  
    IF(ist2)THEN
      IF(.NOT.ALLOCATED(lakestate%uppertemp)) ALLOCATE(lakestate%uppertemp(nl,n))
      IF(.NOT.ALLOCATED(lakestate%lowertemp)) ALLOCATE(lakestate%lowertemp(nl,n))
    ENDIF  

  END SUBROUTINE allocate_lakestate_variables

  !>Deallocate lake states used for temporary storage of initial
  !>state for submodel simulation
  !---------------------------------------------------------------
  SUBROUTINE deallocate_lakestate_variables(lakestate)

    !Argument declarations
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states

    IF(ALLOCATED(lakestate%water))     DEALLOCATE(lakestate%water)
    IF(ALLOCATED(lakestate%temp))      DEALLOCATE(lakestate%temp)
    IF(ALLOCATED(lakestate%TPmean))    DEALLOCATE(lakestate%TPmean)
    IF(ALLOCATED(lakestate%conc))      DEALLOCATE(lakestate%conc)
    IF(ALLOCATED(lakestate%slowwater)) DEALLOCATE(lakestate%slowwater)
    IF(ALLOCATED(lakestate%concslow))  DEALLOCATE(lakestate%concslow)
    IF(ALLOCATED(lakestate%temp10))    DEALLOCATE(lakestate%temp10)
    IF(ALLOCATED(lakestate%temp20))    DEALLOCATE(lakestate%temp20)
    IF(ALLOCATED(lakestate%uppertemp)) DEALLOCATE(lakestate%uppertemp)
    IF(ALLOCATED(lakestate%lowertemp)) DEALLOCATE(lakestate%lowertemp)

  END SUBROUTINE deallocate_lakestate_variables

  !>Initiate submodel lake states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_lakestate_submodel(ns,n,ist2,indexarray,lakestate,lakestate2)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns       !<number of substances
    INTEGER, INTENT(IN) :: n        !<number of subbasins in states to be set
    LOGICAL, INTENT(IN) :: ist2     !<Status of T2 simulation
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state, submodel
    TYPE(lakestatetype),INTENT(IN) :: lakestate2    !<Lake state, whole model setup

    lakestate%water(:,:) = lakestate2%water(:,indexarray(:))
    lakestate%temp(:,:) = lakestate2%temp(:,indexarray(:))
    IF(ns>0)THEN
      lakestate%TPmean(:,:) = lakestate2%TPmean(:,indexarray(:))
      lakestate%conc(:,:,:) = lakestate2%conc(:,:,indexarray(:))
      lakestate%slowwater(:,:) = lakestate2%slowwater(:,indexarray(:))
      lakestate%concslow(:,:,:) = lakestate2%concslow(:,:,indexarray(:))
      lakestate%temp10(:,:) = lakestate2%temp10(:,indexarray(:))
      lakestate%temp20(:,:) = lakestate2%temp20(:,indexarray(:))
    ENDIF
    IF(ist2)THEN
      lakestate%uppertemp(:,:) = lakestate2%uppertemp(:,indexarray(:))
      lakestate%lowertemp(:,:) = lakestate2%lowertemp(:,indexarray(:))
    ENDIF  
  
  END SUBROUTINE initiate_lakestate_submodel

  !>Initiate lake states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_lakestate_zero(ns,ist2,lakestate)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns                       !<number of substances
    LOGICAL, INTENT(IN) :: ist2                     !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state, submodel

    lakestate%water = 0.
    lakestate%temp = 0.
    IF(ns>0)THEN
      lakestate%TPmean = 0.
      lakestate%conc = 0.
      lakestate%slowwater = 0.
      lakestate%concslow = 0.
      lakestate%temp10 = 0.
      lakestate%temp20 = 0.
    ENDIF
    IF(ist2)THEN
      lakestate%uppertemp(:,:) = 0.
      lakestate%lowertemp(:,:) = 0.
    ENDIF  

  END SUBROUTINE initiate_lakestate_zero

  !>Allocate miscellaneous state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_miscstate_variables(n,ns,nc,isirr,isqar,isgsm,miscstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nc      !<number of classes
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state

    IF(.NOT.ALLOCATED(miscstate%temp5))    ALLOCATE(miscstate%temp5(n))
    IF(.NOT.ALLOCATED(miscstate%temp10))   ALLOCATE(miscstate%temp10(n))
    IF(.NOT.ALLOCATED(miscstate%temp20))   ALLOCATE(miscstate%temp20(n))
    IF(.NOT.ALLOCATED(miscstate%temp30))   ALLOCATE(miscstate%temp30(n))
    IF(isgsm)THEN
      IF(.NOT.ALLOCATED(miscstate%gdd))      ALLOCATE(miscstate%gdd(2,nc,n))
      IF(.NOT.ALLOCATED(miscstate%gsbegin))  ALLOCATE(miscstate%gsbegin(2,nc,n))
    ENDIF
    IF(isirr)THEN
      IF(.NOT.ALLOCATED(miscstate%nextirrigation)) ALLOCATE(miscstate%nextirrigation(nc,n))
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(miscstate%cnextirrigation)) ALLOCATE(miscstate%cnextirrigation(ns,nc,n))
      ENDIF
    ENDIF
    IF(isqar)THEN
      IF(.NOT.ALLOCATED(miscstate%updatestationsarcorr)) ALLOCATE(miscstate%updatestationsarcorr(n))
    ENDIF

  END SUBROUTINE allocate_miscstate_variables

  !>Deallocate miscellaneous states used for temporary storage of initial
  !>state for submodel simulation
  !---------------------------------------------------------------
  SUBROUTINE deallocate_miscstate_variables(miscstate)

    !Argument declarations
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state

    IF(ALLOCATED(miscstate%temp5))  DEALLOCATE(miscstate%temp5)
    IF(ALLOCATED(miscstate%temp10)) DEALLOCATE(miscstate%temp10)
    IF(ALLOCATED(miscstate%temp20)) DEALLOCATE(miscstate%temp20)
    IF(ALLOCATED(miscstate%temp30)) DEALLOCATE(miscstate%temp30)
    IF(ALLOCATED(miscstate%gdd))    DEALLOCATE(miscstate%gdd)
    IF(ALLOCATED(miscstate%gsbegin))         DEALLOCATE(miscstate%gsbegin)
    IF(ALLOCATED(miscstate%nextirrigation))  DEALLOCATE(miscstate%nextirrigation)
    IF(ALLOCATED(miscstate%cnextirrigation)) DEALLOCATE(miscstate%cnextirrigation)
    IF(ALLOCATED(miscstate%updatestationsarcorr)) DEALLOCATE(miscstate%updatestationsarcorr)

  END SUBROUTINE deallocate_miscstate_variables

  !>Initiate submodel miscellaneous states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_miscstate_submodel(ns,n,isirr,isqar,isgsm,indexarray,miscstate,miscstate2)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: n         !<number of subbasins in states to be set
    LOGICAL,INTENT(IN)  :: isirr     !<Status of irrigation simulation
    LOGICAL,INTENT(IN)  :: isqar     !<Status of updating q with AR
    LOGICAL,INTENT(IN)  :: isgsm     !<Status of growth season model
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state, submodel
    TYPE(miscstatetype),INTENT(IN) :: miscstate2    !<misc state, whole model setup

    miscstate%temp5(:)  = miscstate2%temp5(indexarray(:))
    miscstate%temp10(:) = miscstate2%temp10(indexarray(:))
    miscstate%temp20(:) = miscstate2%temp20(indexarray(:))
    miscstate%temp30(:) = miscstate2%temp30(indexarray(:))
    IF(isgsm)THEN
      miscstate%gdd(:,:,:) = miscstate2%gdd(:,:,indexarray(:))
      miscstate%gsbegin(:,:,:) = miscstate2%gsbegin(:,:,indexarray(:))
    ENDIF
    IF(isirr)THEN
      miscstate%nextirrigation(:,:) = miscstate2%nextirrigation(:,indexarray(:))
      IF(ns>0)THEN
        miscstate%cnextirrigation(:,:,:) = miscstate2%cnextirrigation(:,:,indexarray(:))
      ENDIF
    ENDIF
    IF(isqar)THEN
      miscstate%updatestationsarcorr(:) = miscstate2%updatestationsarcorr(indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_miscstate_submodel

  !>Initiate miscellaneous states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_miscstate_zero(ns,isirr,isqar,isgsm,miscstate)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns        !<number of substances
    LOGICAL,INTENT(IN)  :: isirr     !<Status of irrigation simulation
    LOGICAL,INTENT(IN)  :: isqar     !<Status of updating q with AR
    LOGICAL,INTENT(IN)  :: isgsm     !<Status of growth season model
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state, submodel

    miscstate%temp5  = 0.
    miscstate%temp10 = 0.
    miscstate%temp20 = 0.
    miscstate%temp30 = 0.
    IF(isgsm)THEN
      miscstate%gdd = 0.
      miscstate%gsbegin = 0
    ENDIF
    IF(isirr)THEN
      miscstate%nextirrigation = 0.
      IF(ns>0)THEN
        miscstate%cnextirrigation = 0.
      ENDIF
    ENDIF
    IF(isqar)THEN
      miscstate%updatestationsarcorr(:) = 0.
    ENDIF
    
  END SUBROUTINE initiate_miscstate_zero

!------------------------------------
!Subroutines for state file handling
!------------------------------------

  !>Get number of frozen states
  !-------------------------------------------------
  SUBROUTINE get_frozenstate_variables_arraysize(n,ns,nc,nr,nl,isglac,islrice,dim)

    !Argument declarations
    INTEGER,INTENT(IN)  :: n       !<Number of subbasins
    INTEGER,INTENT(IN)  :: ns      !<Number of substances
    INTEGER,INTENT(IN)  :: nc      !<Number of classes
    INTEGER,INTENT(IN)  :: nr      !<Number of river types
    INTEGER,INTENT(IN)  :: nl      !<Number of lake types
    LOGICAL,INTENT(IN)  :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice !<Status of lakeriver ice model
    INTEGER,INTENT(OUT) :: dim     !<Number of array elements

    dim = 0
    dim = dim + n*nc   !snow
    dim = dim + n*nc   !snowage
    dim = dim + n*nc   !snowcov
    dim = dim + n*nc   !snowmax
    IF(ns>0) dim = dim + n*nc*ns    !csnow
    IF(isglac) dim = dim + n    !glacvol
    IF(islrice)THEN
      dim = dim + nl*n   !lakesnow
      dim = dim + nl*n   !lakesnowage
      dim = dim + nl*n   !lakeice
      dim = dim + nl*n   !lakebice
      dim = dim + nl*n   !lakeicecov
      dim = dim + nr*n   !riversnow
      dim = dim + nr*n   !riversnowage
      dim = dim + nr*n   !riverice
      dim = dim + nr*n   !riverbice
      dim = dim + nr*n   !rivericecov
    ENDIF  

  END SUBROUTINE get_frozenstate_variables_arraysize

  !>Write snow and ice state variables to array
  !-------------------------------------------------
  SUBROUTINE set_frozenstate_variables_to_array(n,ns,nc,nr,nl,isglac,islrice,frozenstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(IN) :: frozenstate  !<Snow and ice states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snow(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowage(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowcov(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowmax(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%csnow(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isglac)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%glacvol(i)
      ENDDO
    ENDIF
    IF(islrice)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakesnow(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakesnowage(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakeice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakebice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakeicecov(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riversnow(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riversnowage(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riverice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riverbice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%rivericecov(j,i)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_frozenstate_variables_to_array

  !>Read snow and ice state variables from array
  !-------------------------------------------------
  SUBROUTINE set_frozenstate_variables_from_array(n,ns,nsfil,nc,nr,nl,isglac,islrice,frozenstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snow(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowage(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowcov(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowmax(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%csnow(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n * nc * nsfil  !skip saved concentrations
    ENDIF  
    IF(isglac)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%glacvol(i) = array(iarr-arrshift)
      ENDDO
    ENDIF
    IF(islrice)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakesnow(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakesnowage(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakeice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakebice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakeicecov(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riversnow(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riversnowage(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riverice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riverbice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%rivericecov(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_frozenstate_variables_from_array

  !>Get number of soil states
  !-------------------------------------------------
  SUBROUTINE get_soilstate_variables_arraysize(n,ns,nc,nsl,isN,isP,isC,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = nsl*nc*n        !water
    dim = dim + nsl*nc*n  !temp
    dim = dim + nc*n      !deeptemp
    dim = dim + nc*n      !icelens
    IF(ns>0)THEN
      dim = dim + ns*nsl*nc*n  !conc
      dim = dim + nc*n         !PPrelpool
      IF(isN)THEN
        dim = dim + nsl*nc*n   !humusN
        dim = dim + nsl*nc*n   !fastN
      ENDIF
      IF(isP)THEN
        dim = dim + nsl*nc*n   !humusP
        dim = dim + nsl*nc*n   !fastP
        dim = dim + nsl*nc*n   !partP
      ENDIF
      IF(isC)THEN
        dim = dim + nsl*nc*n   !humusC
        dim = dim + nsl*nc*n   !fastC
        dim = dim + nc*n       !oldgrw
      ENDIF
    ENDIF

  END SUBROUTINE get_soilstate_variables_arraysize

  !>Write soil state variables to array
  !-------------------------------------------------
  SUBROUTINE set_soilstate_variables_to_array(n,ns,nc,nsl,isN,isP,isC,soilstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(soilstatetype),INTENT(IN) :: soilstate  !<Soil states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%water(k,j,i)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%temp(k,j,i)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%deeptemp(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = REAL(soilstate%icelens(j,i))
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,nsl
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%conc(l,k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(iarr>iarrlast) RETURN    !Escape
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%PPrelpool(j,i)
        ENDDO
      ENDDO
      IF(isN)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusN(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastN(k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(iarr>iarrlast) RETURN    !Escape
      IF(isP)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%partP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(iarr>iarrlast) RETURN    !Escape
      IF(isC)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusC(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastC(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%oldgrw(j,i)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE set_soilstate_variables_to_array

  !>Read soil state variables from array
  !-------------------------------------------------
  SUBROUTINE set_soilstate_variables_from_array(n,ns,nsfil,nc,nsl,isN,isP,isC,soilstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    iarr = 0
    arrshift = iarrfirst - 1 
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%water(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%temp(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%deeptemp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%icelens(j,i) = INT(array(iarr-arrshift))
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,nsl
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%conc(l,k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%PPrelpool(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nc*nsl*nsfil + n*nc
    ENDIF
    IF(ns>0)THEN
      IF(isN)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusN(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastN(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(isP)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%partP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(isC)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusC(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastC(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%oldgrw(j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF(nsfil>0)THEN
    ENDIF

  END SUBROUTINE set_soilstate_variables_from_array

  !>Get number of aquifer states
  !-------------------------------------------------
  SUBROUTINE get_aquiferstate_variables_arraysize(na,ns,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 3*na        !water,lastrecharge,nextoutflow
    IF(ns>0)THEN
      dim = dim + 3*ns*na  !conc,clastrecharge,cnextoutflow
    ENDIF

  END SUBROUTINE get_aquiferstate_variables_arraysize

  !>Write aquifer state variables to array
  !-------------------------------------------------
  SUBROUTINE set_aquiferstate_variables_to_array(na,ns,aquiferstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate  !<aquifer states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%water(i)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%lastrecharge(i)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%nextoutflow(i)
    ENDDO
    IF(ns>0)THEN
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%conc(l,i)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%clastrecharge(l,i)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%cnextoutflow(l,i)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_aquiferstate_variables_to_array

  !>Read aquifer state variables from array
  !-------------------------------------------------
  SUBROUTINE set_aquiferstate_variables_from_array(na,ns,nsfil,aquiferstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,l
    INTEGER iarr,arrshift
    
    iarr = 0
    arrshift = iarrfirst - 1 
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%water(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%lastrecharge(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%nextoutflow(i) = array(iarr-arrshift)
    ENDDO
    IF(ns>0)THEN
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%conc(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%clastrecharge(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%cnextoutflow(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + 2*na*nsfil
    ENDIF

  END SUBROUTINE set_aquiferstate_variables_from_array

  !>Get number of river states
  !-------------------------------------------------
  SUBROUTINE get_riverstate_variables_arraysize(n,ns,nr,ml,mt,isN,isP,iswetl,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nr      !<number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    dim = dim + n*nr        !water
    dim = dim + n*nr        !temp
    dim = dim + (ml+1)*n*nr !qqueue
    dim = dim + n*nr        !Qmean
    IF(ns>0)THEN
      dim = dim + n*nr            !TPmean
      dim = dim + n*nr            !Psed
      dim = dim + ns*n*nr         !conc
      dim = dim + n*nr            !temp10
      dim = dim + n*nr            !temp20
      dim = dim + ns*(ml+1)*nr*n  !cqueue
    ENDIF
    IF(isN.OR.isP)THEN
      dim = dim + mt*n*nr   !Qdayacc
      dim = dim + 366*n*nr  !Q365
    ENDIF
    IF(iswetl)THEN
      dim = dim + ns*n*nr   !cwetland
    ENDIF

  END SUBROUTINE get_riverstate_variables_arraysize

  !>Write river state variables to array
  !-------------------------------------------------
  SUBROUTINE set_riverstate_variables_to_array(n,ns,nr,ml,mt,isN,isP,iswetl,riverstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(IN) :: riverstate  !<River states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%water(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        DO k = 0,ml
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%qqueue(k,j,i)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Qmean(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%TPmean(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Psed(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%conc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp10(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp20(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 0,ml
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%cqueue(l,k,j,i)
            ENDDO
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isN.OR.isP)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,mt
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Qdayacc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,366
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Q365(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(iswetl)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%cwetland(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  

  END SUBROUTINE set_riverstate_variables_to_array

  !>Read river state variables from array
  !-------------------------------------------------
  SUBROUTINE set_riverstate_variables_from_array(n,ns,nsfil,nr,ml,mt,isN,isP,iswetl, &
                                                 Nfil,Pfil,riverstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: ml      !<Maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: Nfil    !<Status of nitrogen simulation of file
    LOGICAL,INTENT(IN) :: Pfil    !<Status of phosphorus simulation of file
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%water(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        DO k = 0,ml
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%qqueue(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Qmean(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%TPmean(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Psed(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%conc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp10(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp20(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 0,ml
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%cqueue(l,k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nr*2+n*nr*nsfil+n*nr*2+n*nr*(ml+1)*nsfil
    ENDIF  
    IF(isN.OR.isP)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,mt
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Qdayacc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,366
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Q365(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(Nfil.OR.Pfil)THEN
      iarr = iarr + n*nr*mt
      iarr = iarr + n*nr*366
    ENDIF  
    IF(iswetl)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%cwetland(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  

  END SUBROUTINE set_riverstate_variables_from_array

  !>Get number of lake states
  !-------------------------------------------------
  SUBROUTINE get_lakestate_variables_arraysize(n,ns,nl,ist2,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    dim = dim + n*nl      !water
    dim = dim + n*nl      !temp
    IF(ns>0)THEN
      dim = dim + n*nl    !TPmean
      dim = dim + ns*n*nl !conc
      dim = dim + n*nl    !slowwater
      dim = dim + ns*n*nl !concslow
      dim = dim + n*nl    !temp10
      dim = dim + n*nl    !temp20
    ENDIF
    IF(ist2)THEN
      dim = dim + n*nl    !uppertemp
      dim = dim + n*nl    !lowertemp
    ENDIF  

  END SUBROUTINE get_lakestate_variables_arraysize

  !>Write lake state variables to array
  !-------------------------------------------------
  SUBROUTINE set_lakestate_variables_to_array(n,ns,nl,ist2,lakestate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n             !<Number of subbasins
    INTEGER,INTENT(IN) :: ns            !<Number of substances
    INTEGER,INTENT(IN) :: nl            !<Number of lake types
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<lake states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim           !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%water(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%TPmean(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%conc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%slowwater(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%concslow(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp10(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp20(j,i)
        ENDDO
      ENDDO
    ENDIF
    IF(ist2)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%uppertemp(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%lowertemp(j,i)
        ENDDO
      ENDDO
    ENDIF  

  END SUBROUTINE set_lakestate_variables_to_array

  !>Read lake state variables from array
  !-------------------------------------------------
  SUBROUTINE set_lakestate_variables_from_array(n,ns,nsfil,nl,ist2,lakestate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation (of file)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<lake states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%water(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%TPmean(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%conc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%slowwater(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%concslow(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp10(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp20(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nl + n*nl*nsfil
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%water(j,i) = lakestate%water(j,i) + array(iarr-arrshift)
        ENDDO
      ENDDO
      iarr = iarr + n*nl*nsfil + n*nl*2
    ENDIF  
    IF(ns==0.AND.ist2)THEN
      iarr = iarr + n*nl*2
    ELSEIF(ist2)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%uppertemp(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%lowertemp(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDIF  

  END SUBROUTINE set_lakestate_variables_from_array

  !>Get number of miscellaneous states
  !-------------------------------------------------
  SUBROUTINE get_miscstate_variables_arraysize(n,ns,nc,isirr,isqar,isgsm,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    dim = dim + n   !temp5
    dim = dim + n   !temp10
    dim = dim + n   !temp20
    dim = dim + n   !temp30
    IF(isgsm)THEN
      dim = dim + 2*n*nc   !gdd
      dim = dim + 2*n*nc   !gsbegin
    ENDIF
    IF(isirr)THEN
      dim = dim + n*nc   !nextirrigation
      IF(ns>0)   dim = dim + n*nc*ns  !cnextirrigation
    ENDIF
    IF(isqar)THEN
      dim = dim + n   !updatestationsarcorr
    ENDIF  
    
  END SUBROUTINE get_miscstate_variables_arraysize

  !>Write miscellaneous state variables to array
  !-------------------------------------------------
  SUBROUTINE set_miscstate_variables_to_array(n,ns,nc,isirr,isqar,isgsm,miscstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    TYPE(miscstatetype),INTENT(IN) :: miscstate  !<misc states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp5(i)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp10(i)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp20(i)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp30(i)
    ENDDO
    IF(isgsm)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%gdd(k,j,i)
          ENDDO
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = REAL(miscstate%gsbegin(k,j,i))
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isirr)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%nextirrigation(j,i)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%cnextirrigation(k,j,i)
            ENDDO
          ENDDO  
        ENDDO
      ENDIF
    ENDIF  
    IF(isqar)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%updatestationsarcorr(i)
      ENDDO
    ENDIF

  END SUBROUTINE set_miscstate_variables_to_array

  !>Read misc state variables from array
  !-------------------------------------------------
  SUBROUTINE set_miscstate_variables_from_array(n,ns,nsfil,nc,isirr,isqar,isqar2,isgsm,miscstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isqar2  !<Status of updating q with AR for file
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp5(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp10(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp20(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,n
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp30(i) = array(iarr-arrshift)
    ENDDO
    IF(isgsm)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%gdd(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%gsbegin(k,j,i) = INT(array(iarr-arrshift))
          ENDDO
        ENDDO  
      ENDDO
    ENDIF
    IF(isirr)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%nextirrigation(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%cnextirrigation(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO  
        ENDDO
      ELSEIF(nsfil>0)THEN
        iarr = iarr + n*nc*nsfil  
      ENDIF
    ENDIF
    IF(isqar.AND.isqar2)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%updatestationsarcorr(i) = array(iarr-arrshift)
      ENDDO
    ELSEIF(isqar2)THEN
      iarr = iarr + n
    ELSEIF(isqar)THEN
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%updatestationsarcorr = 0.    !no error from previous simulation
    ENDIF    

  END SUBROUTINE set_miscstate_variables_from_array

END MODULE STATETYPE_MODULE
