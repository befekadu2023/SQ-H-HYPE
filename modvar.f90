!> \file modvar.f90
!> Contains module modvar.

!>\brief Declaration of model data types and variables
!>
!>This is the "interface" between the model and the simulation system
!!"HYSS". NOTE: everything is public
!>
!>The module modvar contains data types, constants, variables, and procedures 
!>that are HYSS provide for the model. It also holds procedures 
!>for allocating and setting its variables and to locate data in them.
!>
!>The constants define specific numbers, dimensions of arrays, and codes 
!>for accumulation period, model parameter dependence and updating methods.
!>
!>The variables are for holding input data on model set-up, substance modelled,
!>forcing data of current time step, state variables that are common for hydrological 
!>models, variables of (current) time, size of model set-up, settings for 
!>the model run (e.g. output, model parameter values, updating information), 
!>and variables for temporary state transformation.

MODULE MODVAR

!Copyright 2011-2015 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------
!Procedures in this module (all public)
!-----------------------------------------------------------------------------------------
! allocate_state_variables
! allocate_other_modvar
! allocate_outvar
! allocate_substances
! deallocate_modvar
! allocate_modelparameters
! allocate_basinvariables
! initiate_xobsinformation
! save_xobsinformation
! set_soillayerdepth
! set_soildepth
! set_coded_classes
! set_update_nutrientcorr
! allocate_temp_state_variables
! allocate_temp_other_modvar
! allocate_temp_substances
! deallocate_temp_modvar_state
! deallocate_modvar_state
! deallocate_modelmodule_states
! allocate_temp2_modelmodule_state
! deallocate_temp2_modelmodule_states
! initiate_cropdata
! find_cropdata
! find_croppart
!-----------------------------------------------------------------------------------------
  USE LIBDATE, ONLY : DateType
  
  IMPLICIT NONE
  SAVE

!> \name Miscellaneous constant parameters
!> \{
  INTEGER, PARAMETER :: max_classoutvar = 8     !<Number of class load output variables
  INTEGER, PARAMETER :: max_basinoutvar = 20    !<Number of basin load output variables
  INTEGER, PARAMETER :: maxsoillayers = 3       !<Maximum number of soil layers
  INTEGER, PARAMETER :: max_pstype = 3          !<Maximum number of point sources types
  REAL, PARAMETER    :: missing_value = -9999.  !<Value used for missing value
  REAL, PARAMETER    :: realzero = 1.E-37       !<Value of zero in real
  REAL, PARAMETER    :: pi = 3.1415927          !<Mathematical constant pi
  REAL, PARAMETER    :: solar = 0.0820          !<Solar constant, MJ/m2/min
  REAL, PARAMETER    :: cwater = 4.2            !<Heat capacity of water, kJ/kg/C
  REAL, PARAMETER    :: cice   = 2.1            !<Heat capacity of ice, kJ/kg/C
  REAL, PARAMETER    :: Lfreezing = 335.        !<Latent heat of freezing, kJ/kg
!> \}
!> \name Model parameter dependence type parameters
!> These codes are used for identifying which model parameter-array to find a specific parameter in.
!> \{
  INTEGER, PARAMETER :: m_gpar   = 1        !<Code for general parameter type
  INTEGER, PARAMETER :: m_spar   = 2        !<Code for soil type dependent parameter type
  INTEGER, PARAMETER :: m_lpar   = 3        !<Code for land use dependent parameter type
  INTEGER, PARAMETER :: m_bpar   = 4        !<Code for basin dependent parameter type
  INTEGER, PARAMETER :: m_rpar   = 5        !<Code for parameter region dependent parameter type
  INTEGER, PARAMETER :: m_lrpar  = 6        !<Code for lake region dependent parameter type
  INTEGER, PARAMETER :: m_ldpar  = 7        !<Code for lake specific parameter imported from LakeData.txt
  INTEGER, PARAMETER :: m_mpar   = 8        !<Code för monthly dependent parameters
  INTEGER, PARAMETER :: m_wqrpar = 9        !<Code for water quality parameter region dependent parameter type
  INTEGER, PARAMETER :: m_ilrpar =10        !<Code for ilake parameter region dependent parameter type
  INTEGER, PARAMETER :: m_olrpar =11        !<Code for olake parameter region dependent parameter type
!> \}
!> \name Output variable type parameters, code for accumulation of data
!> These codes are used to define how a period value of a variable is calculated.
!> \{
  INTEGER, PARAMETER :: i_sum    = 0    !<Type of output variable, sum over period
  INTEGER, PARAMETER :: i_mean   = 1    !<Type of output variable, mean over period
  INTEGER, PARAMETER :: i_wmean  = 2    !<Type of output variable, weigted mean over period
!> \}
!> \name Updating parameters, code and name of updating method
!> Seven updating methods exist. These can be used separately or together. 
!> \{
  INTEGER, PARAMETER :: i_quseobs   = 1          !<Code for updating method, subbasin outflow Q-station updating
  INTEGER, PARAMETER :: i_tpcorr    = 2          !<Code for updating method, updating of total phosphorus out of subbasin
  INTEGER, PARAMETER :: i_wendupd   = 3          !<Code for updating method, updating of outlet lake water stage at end of time step
  INTEGER, PARAMETER :: i_qar       = 4          !<Code for updating method, subbasin outflow Q-station AR-updating
  INTEGER, PARAMETER :: i_war       = 5          !<Code for updating method, subbasin outflow Q updating based on AR on w-error
  INTEGER, PARAMETER :: i_tncorr    = 6          !<Code for updating method, updating of total nitrogen out of subbasin
  INTEGER, PARAMETER :: i_tnloccorr = 7          !<Code for updating method, updating of total nitrogen locally out of subbasin
  INTEGER, PARAMETER :: i_tploccorr = 8          !<Code for updating method, updating of total phosphorus locally out of subbasin
  INTEGER, PARAMETER :: dim_update  = 8          !<Number of updating methods (for array dimension)
  CHARACTER(LEN=9) :: updatename(dim_update)     !<Name of update functions
  PARAMETER(updatename = (/'quseobs  ','tpcorr   ','wendupd  ','q-AR     ','w-AR     ','tncorr   ','tnloccorr','tploccorr'/))
!> \}

  !Type declarations
!> \brief Type for holding data about subbasin characteristics
  TYPE BASINTYPE
    INTEGER :: subid           !<identification number for subbasin
    REAL    :: area            !<area (m2)
    REAL    :: elev            !<elevation above sea level (m)
    REAL    :: selev           !<spread of elevation (m)
    REAL    :: slope           !<slope (?)
    REAL    :: sslope          !<spread of slope (?)
    REAL    :: xcoord          !<x-coordinate
    REAL    :: ycoord          !<y-coordinate
    REAL    :: longitude       !<longitude, degrees E [-180,180]
    REAL    :: latitude        !<latitude, degrees N [-90,90]
    INTEGER :: region          !<vegetation region used by CropData
    INTEGER :: lakeregion = 1  !<lake region for lake parameters
    INTEGER :: parregion = 1   !<parameter region
    INTEGER :: wqparregion = 1 !<water quality parameter region
    REAL    :: rivlen = 0.     !<length of main river (m)
    REAL    :: ilakecatch = 0. !<part of catchment/runoff that runs to ilake
    REAL    :: lakedepth(1:2)       !<olake water depth at threshold (m)
    REAL    :: closewater = 0. !<part of area (possibly agricultural area) close to water courses
    REAL    :: buffer = 0.     !<part of river length with buffer zones (of river in the area above (e.g. agricultural area))
    ! regions for olakes and ilakes
    INTEGER :: olakeregion = 0 !<ilake parameter region (region 0 uses general parameters, regions>0 uses ilakeparreg parameters)
    INTEGER :: ilakeregion = 0 !<olake parameter region (region 0 uses general parameters, regions>0 uses olakeparreg parameters, and (as usual) all cases replaced by lakedata or damdata parameters) 
    ! MM2016:non-contributing fraction
    REAL    :: noncontr = 0.
  END TYPE BASINTYPE
    
!> \brief Type for holding information about dams
  TYPE DAMTYPE
    REAL    :: qprod1 = 0.   !<regulated dam, discharge between wmin and w0 (m3/s) for period 1
    REAL    :: qprod2 = 0.   !                                                      for period 2
    REAL    :: qamp = 0.     !<regulated lake, amplitude of variation of average discharge between wmin and w0
    REAL    :: qpha = 102.   !<regulated lake, phase of variation of average discharge
    REAL    :: w0ref = 0.    !<threshold for dam, "dämningsgräns" (m in w-reference-system)
    REAL    :: wmin = 0      !<threshold for dam, "sänkningsgräns" (m) in local system ("dämningsgräns"=0)
    REAL    :: wampcoeff     !<adjustment factor for waterstage to "real" regulation amplitude
    REAL    :: rate = 0.     !<rate parameter of outflow equation (m2/s)
    REAL    :: exp = 0.      !<exponent parameter of outflow equation
    REAL    :: regvol = 0.   !<Volume of lake used for regulation
    REAL    :: qinfmed = 0   !<Observed or Simulated inflow to dam (from natural run), mean of Year
    REAL    :: qinfmin = 0   !<Observed or Simulated inflow to dam (from natural run), minimum monthly mean
    REAL    :: qinfmax = 0   !<Observed or Simulated inflow to dam (from natural run), maximum monthly mean
    REAL    :: snowfrac = 0  !<Observed of Simulated fraction of precip that is snow
    INTEGER :: purpose = 0   !<Purpose of dam, 1=irrigation, 2=water supply, 3=flood control, 4=hydroelectricity, 0 = natural
    INTEGER :: datum1  = 0   !<date for two date based prod rate
    INTEGER :: datum2  = 0   !<date for two date based prod rate
    ! MH2017: extradata v2 parameters
    INTEGER :: identify = 0  !<Unique 4-digit code used to open DailyDataXXXX.txt files
    INTEGER :: purpose2 = 0  !<purpose of dam prior to being made a special dam;
                             !<allows turning special dams 'off' so they function as normal dams
    REAL    :: maxqmult = 0. !<Coefficient to Scale Transition Band Type
    REAL    :: wslinit  = 0. !<Lake Water Surface Level used to Initialize Model
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
    REAL, ALLOCATABLE :: monthly(:,:) !<Used for monthly regression calculations, aka weighted flow, by column
                                      !<column 1:minA 2:midA 3:highA 4:natminA 5:natmidA 6: nathighA 7:minB 8:midB 9:highB
    REAL, ALLOCATABLE :: daily(:,:) !<Stores daily extra data by column
    !column 1:drought level 2:lower op.s 3:lower band 4:avg.level 5:upper band 6:upper op.s 7:flood level 
    !column 8:(smoothed) drought percentile 9:sm.min.prcntl 10:sm.avg.prcntl 11:sm.max.prcntl 12:sm.flood.prcntl
    !column 13:natural drought percentile 14:nat.min.prcntl 15:nat.avg.prcntl 16:nat.max.prcntl 17:nat.flood.prcntl 18:monthly weighting
    INTEGER :: idt      = 0  !<used to check if early-model calculations apply
    INTEGER :: idtlast  = 0  !<prevents sabd14day(:), weekly7day(:) updating when get_current_outflow_parameters()
    INTEGER :: curmonth = 0  !<current month, extracted from datetype
    INTEGER :: dayofmonth = 0 !<current day of month, extracted from datetype
                             !<is called by check_ouflow_from_lake()
    REAL    :: storyes  = 0. !<storage (m3) yesterday, = A * (lakewst - B) ^ C + D
    REAL    :: stortod  = 0. !<storage (m3) today, = storyes + qinprev - qoutprev
    REAL    :: qinprev  = 0. !<holds today's inflow (m3/s) for next timestep
    REAL    :: qoutprev = 0. !<...outflow...
    REAL, ALLOCATABLE :: qin14day(:) !<stores inflow from past 14 days. qin14day(1) = today's qin, qin14day(2) = yesterday's qin, etc.
    REAL, ALLOCATABLE :: sabd14day(:) !<stores average band storage change from past 14 days. sabd14day(1) = today's sabd, sabd14day(2) = yesterday's sabd, etc.
    REAL, ALLOCATABLE :: weekly7day(:) !<stores powerhouse restricted flows from past 7 days. weekly7day(1) = today's flow, weekly7day(2) = yesterday's flow, etc.
    REAL    :: dailyyes = 0. !<yesterday's daily restricted flow
    !additional OCBD dam parameters
    INTEGER :: natswitch = 0 !<0 = No Natural Outflow, 1 = Has Natural Outflow
    INTEGER :: locswitch = 0 !<0 = No Local Inflow, 1 = Has Local Inflow
    INTEGER :: condswitch = 0 !<0 = No Flow Condition, 1 = Has Conditioning
    INTEGER :: condrestr = 0 !<0 = Low to High, 1 = Drought to Flood
    INTEGER :: condqowsl = 0 !<0 = Conditioned by Qo, 1 = Conditioned by WSL
    INTEGER :: condsubid = 0 !<HYPE Basin SubID Used for Conditioning Value
    INTEGER :: condindex = 0 !<HYPE Basin Index Used for Conditioning Value
    REAL    :: conda    = 0. !<where: Qo = A * (XX - B) ^ C + D
    REAL    :: condb    = 0. !<
    REAL    :: condc    = 0. !<
    REAL    :: condd    = 0. !<
    REAL    :: natminqa = 0. !<where: Qo = A * (WSL - B) ^ C + D
    REAL    :: natminqb = 0. !<
    REAL    :: natminqc = 0. !<
    REAL    :: natminqd = 0. !<
    REAL    :: natmidqa = 0. !<
    REAL    :: natmidqb = 0. !<
    REAL    :: natmidqc = 0. !<
    REAL    :: natmidqd = 0. !<
    REAL    :: natmaxqa = 0. !<
    REAL    :: natmaxqb = 0. !<
    REAL    :: natmaxqc = 0. !<
    REAL    :: natmaxqd = 0. !<
    REAL    :: natdeltaw = 0. !<Weekly Flow Change Limit
    REAL    :: natdeltad = 0. !<Daily Flow Change Limit
    REAL    :: natphqo  = 0. !<Powerhouse Limit Flow Value
    REAL    :: natdrqo  = 0. !<Minimum Fixed Flow Value
    REAL    :: natflqo  = 0. !<Maximum Fixed Flow Value
    INTEGER :: natdrtype = 0  !<0 = Dr. Percentile, 1 = Fixed Min.
    INTEGER :: natmintype = 0  !<0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natmidtype = 0  !<0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natmaxtype = 0  !<3 = Storage Change, 0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
    INTEGER :: natfltype = 0  !<0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
    INTEGER :: natphtype = 0  !<0 = No PH Restriction, 1 = Hold Water
    REAL    :: natdailyyes = 0. !<yesterday's daily restricted flow
    REAL,ALLOCATABLE :: wsl14day(:)  !<holds past 14 days' WSL for use by conditioned dams, size 15 is intentional
    REAL,ALLOCATABLE :: qo14day(:)   !<holds past 14 days' Qo for use by conditioned dams, size 15 is intentional 
    INTEGER :: orderswitch = 0 !<0 = conditioned dam calculated first, 1 = conditioning dam calculated first
                             !<0 = use ...14day(1:14), 1 = use ...14day(2:15)
    REAL,ALLOCATABLE :: resqo7day(:) !<holds past 7 days' reservoir Qo for use by Storage Delta High Flow Type
    REAL,ALLOCATABLE :: qin7day(:) !<holds past 7 days' Qin for use by Storage Delta High Flow Type
                                   !<note: timing is different from qin14day(), DO NOT interchange them
    REAL, ALLOCATABLE :: natweekly7day(:) !<stores powerhouse restricted flows from past 7 days. weekly7day(1) = today's flow, weekly7day(2) = yesterday's flow, etc.                               
    INTEGER :: branchswt = 0 !<0 = reservoir is main, natural is branch, 1 = reservoir is branch, natural is main
  END TYPE DAMTYPE
  
!> \brief Type for holding information about lakes
  TYPE LAKETYPE
    REAL    :: area          !<(whole) lake area (m2)
    REAL    :: qprod1 = 0.   !<regulated lake, discharge between wmin and w0 (m3/s) for period 1
    REAL    :: qprod2 = 0.   !                                                      for period 2
    INTEGER :: datum1 = 0    !<Start production regime period 1
    INTEGER :: datum2 = 0    !<Start production regime period 2
    REAL    :: qamp = 0.     !<regulated lake, amplitude of variation of average discharge between wmin and w0
    REAL    :: qpha = 102.   !<regulated lake, phase of variation of average discharge
    REAL    :: rate = 0.     !<rate parameter of outflow equation (m2/s)
    REAL    :: exp = 0.      !<exponent parameter of outflow equation
    REAL    :: w0ref = 0.    !<threshold for lake or "dämningsgräns" for regulated lake (m in w-reference-system)
    REAL    :: deltaw0 = 0.  !<difference in threshold for lake rating curve production period 2, waterstage parameter of outflow equation (m)
    REAL    :: wmin          !<threshold for lake or "sänkningsgräns" for regulated lake, waterstage parameter of produktion flow (m in local system ("dämningsgräns"=0)
    REAL    :: wampcoeff     !<adjustment factor for waterstage to "real" regulation amplitude
  END TYPE LAKETYPE
!> \brief Type for holding information about lakebasins
  TYPE LAKEBASINTYPE
    INTEGER :: ilk = 0       !<index for lake in lake-variable in which this lake basin is a part
    LOGICAL :: last = .FALSE.!<flag for lakebasin with outlet of (whole) lake
  END TYPE LAKEBASINTYPE
!> \brief Type for holding information about classes in a subbasin
  TYPE CLASSBASINTYPE
    REAL    :: part         !<part of area for landuse-soil combinations
    REAL    :: part2cr = 0. !<part of class-area that has also a secondary crop
    REAL    :: deltah  = 0. !<deviation of class from mean elevation (m)
  END TYPE CLASSBASINTYPE
!> \brief Type for holding information about glacier class
  TYPE GLACIERTYPE
    INTEGER :: gtype         !<glacier type (0=glacier, 1=icecap)
    REAL    :: volcorr       !<correction of volume-area relationship (glacvcoef) due to combining of several glacier into one class
    INTEGER :: glacdate      !<date of slc area information
    REAL    :: glacinimb     !<annual mass balance for correction of initial volume (mm/year)
  END TYPE GLACIERTYPE
!> \brief Type for holding information about aquifers
  TYPE AQUIFERTYPE
    INTEGER :: parregion     !<parameter region
    REAL    :: area          !<aquifer horizontal area (m2)
    REAL    :: porosity      !<(average) porosity of aquifer (-)
    REAL    :: maxvol        !<maximum volume of aquifer (m3) 
    REAL    :: basedepth     !<base depth (m) of aquifer
    REAL    :: inivol        !<initial volume of aquifer (m3)
    REAL    :: temperature = 0.  !<temperature of aquifer water (degree Celsius)
    REAL    :: conc_IN = 0.  !<initial IN concentration of aquifer (mg/L)
    REAL    :: conc_SP = 0.  !<initial SP concentration of aquifer (mg/L)
    REAL    :: percdelay     !<delay of deep percolation (days)
    REAL    :: retrate       !<rate parameter of return flow 
  END TYPE AQUIFERTYPE
!> \brief Type for holding information about flow path between subbasins
  TYPE PATHTYPE
    INTEGER :: main = 0      !<receiving subbasin main flow !!kan lägga in en bakåtkoppling här om behövs
    INTEGER :: grw1 = 0      !<receiving subbasin regional groundwater flow TODO this can be combined with aquid (recieveid)
    REAL    :: grwtolake = 0 !<part of groundwater flow that goes to the subbasins olake TODO:THIS can be combined with recievefraction
    INTEGER :: aquid            !<id of aquifer this subbasin belongs to, 1..naquifers
    LOGICAL :: rechargebasin    !<delivering subbasin
    REAL    :: recievefraction  !<fraction of aquifer outflow to this subbasin
  END TYPE PATHTYPE
!> \brief Type for holding information about branching flow
  TYPE BRANCHTYPE
    INTEGER :: source        !<subbasin with branch flow
    INTEGER :: branch        !<subbasin receiving branch flow
    REAL    :: mainpart      !<main part of outflow
    REAL    :: maxQ          !<maximum flow in main part
    REAL    :: minQ          !<minimum flow in main part
    REAL    :: maxQbranch    !<maximum flow in branch
  END TYPE BRANCHTYPE
!> \brief Type for holding information about the forcing data
  TYPE FORCINGTYPE
    REAL :: tobselev      !<elevation above sea level for Tobs-grid (m)
  END TYPE FORCINGTYPE
!> \brief Type for holding applied nutrient load data
  TYPE NPCLOADTYPE  
    REAL    :: inwetdep  = 0.  !<IN concentration of wet deposition (mg/l)
    REAL    :: indrydep(3) = 0. !<IN load of dry deposition on land use group 1-3(vegtype in GeoClass) (kg/km2/timestep)
    REAL    :: tpconcloc = 0.  !<P concentration of local diffuse source (mg/l)
    REAL    :: tnconcloc = 0.  !<N concentration of local diffuse source (mg/l)
    REAL    :: volloc    = 0.  !<Water discharge from local diffuse source (m3/day)
    REAL    :: inpartloc = 0.  !<Fraction of inorganic nitrogen of local diffuse source
    REAL    :: sppartloc = 0.  !<Fraction of SRP of local diffuse source
    REAL    :: abstrvol(2) = 0.  !<Abstraction volume from 1=main river,2=olake (m3/s)
    REAL,ALLOCATABLE :: psvol(:)     !<Water discharge from point sources (max_pstype) (m3/s)
    REAL,ALLOCATABLE :: psload(:,:)  !<Nutrient load from point sources (max_pstype,numsubstances) (kg/timestep)
  END TYPE NPCLOADTYPE
!> \brief Type for holding general substance load to surface water
  TYPE T1LOADTYPE  
    INTEGER :: subindex   !<Subbasin index (i, not subid)
    INTEGER :: sw_code    !<Recieving surface water code; 1=local river,2=ilake,3=main river,4=olake
    REAL    :: psvol      !<Water discharge from point source (m3/s)
    REAL    :: psconc     !<T1 load from point source (conc)
    REAL    :: pstemp     !<T2-water temperature of point source (degree C)
  END TYPE T1LOADTYPE
!> \brief Type for holding information about irrigation
  TYPE IRRIGATIONTYPE
    INTEGER :: subid                  !<identification number for subbasin
    INTEGER :: regsourceid = 0        !<identification number for regional source for irrigation of subbasin
    INTEGER :: demandtype             !<type of irrigation water demand equation
    LOGICAL :: dam         = .FALSE.  !<status of irrigation dam?
    REAL    :: local_eff   = 1.       !<efficiency of irrigation in the subbasin (-)
    REAL    :: reg_eff     = 1.       !<efficiency of irrigation in the regional irrigation net (-)
    REAL    :: gw_part                !<part of irrigation from (local) groundwater (-)
    REAL    :: sw_part                !<part of irrigation from (local) surface water (-)
  END TYPE IRRIGATIONTYPE
!> \brief Type for holding information about classes
  TYPE CLASSTYPE
    INTEGER :: luse          !<landuse code
    INTEGER :: soil          !<soil type code
    INTEGER :: crop          !<main crop id number
    INTEGER :: crop2         !<secondary crop id number
  END TYPE CLASSTYPE
!> \brief Type for holding information about vegetation; nutrient uptake, fertilization etc.
  TYPE CROPDATATYPE
    REAL    :: fertnamount1  !<fertiliser amount N (kg/km2)
    REAL    :: fertpamount1  !<fertiliser amount P (kg/km2)
    INTEGER :: fertday1      !<day for fertilizing 1 (dayno)
    REAL    :: fertdown1     !<part of fertilizer amount ploughed down to soillayer 2
    REAL    :: mannamount1   !<manure amount N (kg/km2)
    REAL    :: manpamount1   !<manure amount P (kg/km2)
    INTEGER :: manday1       !<day for manureing 1 (dayno)
    REAL    :: mandown1      !<part of manure amount ploughed down to soillayer 2
    REAL    :: fertnamount2  !<fertiliser amount N (kg/km2)
    REAL    :: fertpamount2  !<fertiliser amount P (kg/km2)
    INTEGER :: fertday2      !<day for fertilizing 2 (dayno)
    REAL    :: fertdown2     !<part of fertilizer amount ploughed down to soillayer 2
    REAL    :: mannamount2   !<manure amount N (kg/km2)
    REAL    :: manpamount2   !<manure amount P (kg/km2)
    INTEGER :: manday2       !<day for manureing 2 (dayno)
    REAL    :: mandown2      !<part of manure amount ploughed down to soillayer 2
    REAL    :: resnamount    !<plant residual amount of N (kg/km2/yr)
    REAL    :: respamount    !<plant residual amount of P (kg/km2/yr)
    REAL    :: rescamount    !<plant resudual (litterfall) of C (kg/km2/yr)
    INTEGER :: resdayno      !<day for residual (dayno)
    REAL    :: resdown       !<part of residual amount ploughed down to soillayer 2
    REAL    :: resfast       !<part of residual amount to fast pool (rest to humus)
    REAL    :: uptake1       !<parameter 1 for plant uptake of nutrient
    REAL    :: uptake2       !<parameter 2 for plant uptake of nutrient
    REAL    :: uptake3       !<parameter 3 for plant uptake of nutrient
    REAL    :: uptakeupper   !<fraction of plant uptake in upper soil layer
    REAL    :: PNuptakeRatio !<phosphorus plant uptake as a factor of nitrogen uptake
    INTEGER :: baredayno1    !<day for beginning of first period with bare soil (dayno), typically 1
    INTEGER :: baredayno2    !<day for end of first period with bare soil (dayno), typically sawing date + a few days
    INTEGER :: baredayno3    !<day for beginning of second period with bare soil (dayno), typically ploughing date
    INTEGER :: baredayno4    !<day for end of second period with bare soil (dayno), typically 365
    INTEGER :: baredayno5    !<day for end of second period with bare soil (dayno), typically saw date of autumn crop
    REAL    :: ccmax1        !<crop cover during summer
    REAL    :: ccmax2        !<crop cover during winter and all year for year-round-crops (e.g. forest)
    REAL    :: gcmax1        !<ground cover during summer
    REAL    :: gcmax2        !<ground cover during winter and all year for year-round-crops (e.g. forest)
    REAL    :: gddsow        !<accumulated growing degree days needed for sowing
    REAL    :: daylength     !<minimum daylength required for start to accumulate gdd (hours)
    REAL    :: basetemp      !<temperature which is deducted from air temperature when calculating gdd (degrees)
    INTEGER :: firstday      !<first possible day for gdd accumulation (dayno)
  END TYPE CROPDATATYPE
!> \brief Type for holding information irrigation for vegetation
  TYPE CROPIRRDATATYPE
    INTEGER :: season_end         !<end of irrigation season (may be larger than 365)
    INTEGER :: plantingdayno      !<day of planting (dayno)/irrigation season start
    INTEGER :: lengthini          !<length of initial crop growth state (days)
    INTEGER :: lengthdev          !<length of development crop growth state (days)
    INTEGER :: lengthmid          !<length of middle crop growth state (days)
    INTEGER :: lengthlate         !<length of late crop growth state (days)
    INTEGER :: imm_start = 0      !<day for start of immersion season (dayno)
    INTEGER :: imm_end            !<day for end of immersion season (dayno)
    REAL    :: kcbini             !<parameter for plant transpiration at initial growth state
    REAL    :: kcbmid             !<parameter for plant transpiration at middle growth state
    REAL    :: kcbend             !<parameter for plant transpiration at end of late growth state
    REAL    :: dlref              !<reference depletion level (-)
  END TYPE CROPIRRDATATYPE
!> \brief Type for holding information about wetlands
  TYPE WETLANDTYPE
    REAL    :: area              !<river wetland surface area (m2)
    REAL    :: depth             !<river wetland depth (m)
    REAL    :: part              !<river wetland part of riverflow as inflow (-)
  END TYPE WETLANDTYPE
!> \brief Type for holding information about output variables
  TYPE IDTYPE
    CHARACTER(LEN=4)  :: shortname      !<short name of output variable
    CHARACTER(LEN=20) :: longname       !<long name of output variable
    INTEGER           :: vartype        !<type of output variable; state,flow or conc
    INTEGER           :: water          !<Corresponding water to concentration (outvar number)
    CHARACTER(LEN=6)  :: shortunit      !<short unit of output variable
    CHARACTER(LEN=20) :: longunit       !<unit of output varible
    CHARACTER(LEN=12) :: mapfilename    !<file name of map-file of output variable
    CHARACTER(LEN=12) :: timefilename   !<file name of time-file of output variable
  END TYPE IDTYPE
!> \brief Type for holding definition of model parameters
  TYPE MODPARIDTYPE
    CHARACTER(LEN=10) :: shortname       !<Name of parameter
    INTEGER           :: deptype         !<Type of dependence for the parameter
    INTEGER           :: parno           !<Parameter number
  END TYPE MODPARIDTYPE

!> \name Variables for configuration of simulation
!> \brief HYSS make available variables to hold information on the configuration of the simulation.
!> \{
  !Substances available for modelling, index of state variables
  INTEGER :: i_in     !<index of state variable, inorganic nitrogen
  INTEGER :: i_on     !<index of state variable, organic nitrogen
  INTEGER :: i_sp     !<index of state variable, soluble (reactive) phosphorus, i.e. phosphate
  INTEGER :: i_pp     !<index of state variable, particulate phosphorus
  INTEGER :: i_t1     !<index of state variable, conservative substance 1, evaporates!
  INTEGER :: i_t2     !<index of state variable, tracer 2, water temperature
  INTEGER :: i_t3     !<index of state variable, conservative substance 3
  INTEGER :: i_oc     !<index of state variable, (dissolved) organic carbon
  !Substances for keeping track of water sources
  INTEGER :: i_sm     !<index of state variable, water from snowmelt
  INTEGER :: i_gm     !<index of state variable, water from glacier melt
  INTEGER :: i_rn     !<index of state variable, water from rain
  INTEGER :: i_li     !<index of state variable, water from lake initial state
  INTEGER :: i_ri     !<index of state variable, water from river initial state
  INTEGER :: i_si     !<index of state variable, water from soil initial state
  INTEGER :: i_sr     !<index of state variable, water from surface runoff
  INTEGER :: i_dr     !<index of state variable, water from drainage
  

  !Model simulation status variables (make a simulationoption like modeloption?)
  LOGICAL :: conductN             !<status of nitrogen simulation
  LOGICAL :: conductP             !<status of phosphorus simulation
  LOGICAL :: conductC             !<status of organic carbon simulation
  LOGICAL :: conductT             !<status of tracer simulation; T1, T2 or T3
  LOGICAL :: conductT2            !<status of water temperature simulation
  LOGICAL :: doirrigation         !<Status if irrigation included in model set-up, till hypevar använd i model if(doirrigation) CALL subrooutine...
  LOGICAL :: irrunlimited         !<If true, simulate unlimited irrigation, else withdraw from actual sources
  LOGICAL :: soiliniwet           !<If true, initiate soil water to porosity, else filed capacity
  LOGICAL :: conductxoms          !<Status if XobsXOMN or XobsXOSN-files should be read
  LOGICAL :: conductwb            !<Status if water balance output should be calculated and written
  LOGICAL :: conductregest        !<Flag for activation of regional parameter estimation
  logical :: conductws            !<Flag for activation of water origin trace model

  !Lotta suggest to do this in a new way in info.txt (model submodeloption value, for instance "model snowfall [0-1]")
  ! however, it's still not ready, so I use the existing structure for now
  !
  !Optional Model Structure Settings (maybe these should be HYPEVAR, together with modeloptionname for info)
  !Processes with model options in HYPE. 0 is always the default.
  INTEGER, PARAMETER :: p_snowfall       = 1  !Snowfall process options, (0) threshold temperatures, (1) snowfall fraction input (SFobs.txt) 
  INTEGER, PARAMETER :: p_snowmelt       = 2  !Snowmelt process options, (0) temperature index, (1) temperature index with snowcover scaling, (2) temperature+radiation index
  INTEGER, PARAMETER :: p_lakeriverice   = 3  !Lake river ice model options, (0) off, (1) version1, (2) version2
  INTEGER, PARAMETER :: p_petmodel       = 4  !PET process options, (0) current HYPE, (1) current HYPE without replacement for observations, (2) (modified) Jensen, (3) (modified) Hargreaves-Samani, (4) Priestly-Taylor, (5) FAO PenmanMonteith
  INTEGER, PARAMETER :: p_deepgroundwater= 5  !Regional groundwater options, (0) none, (1) regional groundwater flow delivered instantly, (2)aquifer modelled explicitly
  INTEGER, PARAMETER :: p_swtemperature  = 6  !Surface water temperature model, (0) old (MA-temp), (1) use T2 for lake and river processes
  INTEGER, PARAMETER :: p_snowevap       = 7  !Snowevaporation model, (0) off, (1) epotsnow = epot * fepotsnow (landuse dependent)
  INTEGER, PARAMETER :: p_growthstart    = 8  !Growth season starting model, (0) use bd2 from CropData, (1) use degree day calculation
  INTEGER, PARAMETER :: p_infiltration   = 9  !Infiltration scheme (0) current HYPE, (1) restricted infiltration into frozen soils
  INTEGER, PARAMETER :: p_frozenrunoff   = 10 !Groundwater runoff option (0) current HYPE, (1) reduced due to ice content in subsurface
  INTEGER, PARAMETER :: num_modelprocess = 10 !Number of processes with options in HYPE
  INTEGER,ALLOCATABLE :: modeloption(:)
!> \}

!> \name Variables for input data
!> \brief HYSS make available variables to hold input data for simulation, 
!> both static geographical data and time series data for forcing the model 
!> and for evaluation. Variables ending with suffix i are dependent on subbasin.
!> \{
  !Observed variable declarations
  REAL,ALLOCATABLE :: tempi(:)             !<Temperature data for current time step
  REAL,ALLOCATABLE :: preci(:)             !<Precipitaion data for current time step
  REAL,ALLOCATABLE :: qobsi(:)             !<Runoff data for current time step
  REAL,ALLOCATABLE :: xobsi(:)             !<Other observations for current time step
  INTEGER, ALLOCATABLE :: xobsindex(:,:)   !<Index for finding variables in xobsi
  REAL,ALLOCATABLE :: snowfraci(:)         !<Precipitation snowfall fraction [-] for current time step
  REAL,ALLOCATABLE :: shortwavei(:)        !<Shortwave (downwards) radiation [MJ/m2/day] for current time step
  REAL,ALLOCATABLE :: windi(:)             !<Winds speed [m/s] for current time step
  REAL,ALLOCATABLE :: humidi(:)            !<Relative humidity [0-1] for current time step
  REAL,ALLOCATABLE :: tmini(:)             !<Daily minimum air temperature for current time step
  REAL,ALLOCATABLE :: tmaxi(:)             !<Daily maximum air temperature for current time step

  !Other (static) indata
  TYPE(BASINTYPE),ALLOCATABLE :: basin(:)         !<Basin characteristics
  TYPE(DAMTYPE),ALLOCATABLE :: dam(:)             !<Dam characteristics 
  TYPE(LAKETYPE),ALLOCATABLE :: lake(:)           !<Lake characteristics
  TYPE(LAKEBASINTYPE),ALLOCATABLE :: lakebasin(:) !<Lakebasin characteristics
  TYPE(PATHTYPE),ALLOCATABLE  :: path(:)          !<Subbasin coupling with index
  TYPE(PATHTYPE),ALLOCATABLE  :: pathsubid(:)     !<Subbasin coupling with subid (read from file)
  TYPE(BRANCHTYPE),ALLOCATABLE  :: branchdata(:)  !<Subbasin branch coupling with index
  TYPE(BRANCHTYPE),ALLOCATABLE  :: branchsubid(:) !<Subbasin branch coupling with subid (read from file)
  TYPE(GLACIERTYPE),ALLOCATABLE  :: glacier(:)    !<Glacier characteristics
  TYPE(AQUIFERTYPE),ALLOCATABLE  :: aquifer(:)    !<Aquifer characteristics
  TYPE(FORCINGTYPE),ALLOCATABLE :: forcing(:)     !<Forcing data informations
  TYPE(NPCLOADTYPE),ALLOCATABLE  :: load(:)       !<Loads of nitrogen and phosphorus from point sources and local diffuse load
  TYPE(T1LOADTYPE),ALLOCATABLE  :: tload(:)       !<Load of general substance T1 from point source
  TYPE(CLASSBASINTYPE),ALLOCATABLE :: classbasin(:,:) !<Basin characteristics that is class dependent
  TYPE(CLASSTYPE),ALLOCATABLE :: classdata(:)     !<Class characteristics
  TYPE(CROPDATATYPE),ALLOCATABLE :: cropdata(:)   !<Crop characteristics
  TYPE(CROPIRRDATATYPE),ALLOCATABLE :: cropirrdata(:)  !<Crop characteristics regarding irrigation
  TYPE(IRRIGATIONTYPE),ALLOCATABLE :: irrigationsystem(:) !<Irrigation characteristics
  TYPE(WETLANDTYPE), ALLOCATABLE :: wetland(:,:)  !<River wetland characteristics (local and main river)
  INTEGER, ALLOCATABLE :: cropindex(:,:)          !<Index find cropdata for class/region
  INTEGER, ALLOCATABLE :: damindex(:)             !<Index find dam for subbasin
  INTEGER, ALLOCATABLE :: lakeindex(:)            !<Index find lake for subbasin
  INTEGER, ALLOCATABLE :: glacierindex(:)         !<Index find glacier for subbasin
  INTEGER, ALLOCATABLE :: lakebasinindex(:)       !<Index find lakebasin for subbasin
  INTEGER, ALLOCATABLE :: branchindex(:)          !<Index find branched subbasins
  LOGICAL, ALLOCATABLE :: tloadexist(:)           !<Index find tload for subbasin
!> \}

!> \name Variables for configuration of model set-up
!> \brief HYSS make available variables to hold information on the 
!> configuration of the model set-up.
!> \{
  !Information on model classes
  REAL    soillayerdepth(maxsoillayers)           !<Lower border of soil layers (used by NHMG2)
  REAL    soillayerthick(maxsoillayers)           !<Thickness of soil layers
  REAL, ALLOCATABLE ::    soildepth(:,:)          !<Lower border of soil layers (used by HYPE) (m) (layer,class)
  REAL, ALLOCATABLE ::    soilthick(:,:)          !<Thickness of soil layers (m) (layer,class)
  REAL, ALLOCATABLE ::    tiledepth(:)            !<Depth of drainage pipe or tile
  REAL, ALLOCATABLE ::    streamdepth(:)          !<Depth to ditch or local stream
  INTEGER, ALLOCATABLE :: vegtype(:)              !<Vegetation type (groups land uses into 1=open, 2=forest, 3=water)
  INTEGER, ALLOCATABLE :: rotation(:)             !<Crop rotation group of classes
  LOGICAL :: wetlandexist = .TRUE.                !<Status of wetlands in model set-up (default is wetland present)
  LOGICAL :: glacierexist = .FALSE.               !<Status of glacier class in model set-up (default is no glacier present)
  !Other model variables - to config type? /KN
  INTEGER nsub                    !<Number of subbasins in current run
  INTEGER nsub_basemodel          !<Number of subbasins in original model set up
  INTEGER numsubstances           !<Number of substances currently modelled
  INTEGER ncrop                   !<Number of crops in CropData.txt
  INTEGER nluse                   !<Number of land uses in GeoClass. (max number)
  INTEGER nsoil                   !<Number of soil types in GeoClass. (max number)
  INTEGER nclass                  !<Number of slc-classes GeoClass
  INTEGER nregions                !<Number of parameter regions in parregion
  INTEGER nwqregions              !<Number of water quality parameter regions
  INTEGER nlakeregions            !<Number of lake parameter regions in lakeregion
  INTEGER nbasinlakes             !<Number of lakes divided into basins
  INTEGER naquifers               !<Number of aquifers
  INTEGER nglaciers               !<Number of glaciers
  INTEGER nilakeregions           !<Number of ilake parameter regions
  INTEGER nolakeregions           !<Number of olake parameter regions
  !Coded soil-landuse classes
  INTEGER,ALLOCATABLE :: classmodel(:)   !<Soil model for class
  INTEGER :: slc_olake     !<Class index for outlet lake
  INTEGER :: slc_ilake     !<Class index for local lake
  INTEGER :: slc_lriver    !<Class index for local river
  INTEGER :: slc_mriver    !<Class index for main river
  INTEGER dimriverlag      !<Dimension, river translation lag
!>\}

!> \name Time variables
!> \brief HYSS make available variables to hold information on time  
!> during simulation run and for model set-up.
!> \{
  TYPE(DateType) currentdate      !<Current date in DateType format Year yyyy, Month 01-12, Day 01-31, Hour 00-23, Minute 00-59, not used by HYPE
  INTEGER month                   !<Current month, not used by HYPE
  INTEGER dayno                   !<Current day number of year
  INTEGER prevdoy                 !<Number of days the previous year
  INTEGER tsofday                 !<Current time step of day (e.g. second ts of 24 this day)
  INTEGER seconds_per_timestep    !<Number of seconds for one time step
  INTEGER timesteps_per_day       !<Number of time steps per day
  LOGICAL firstoutstep            !<Flag for first timestep with print out
  LOGICAL endofday                !<Flag for timestep that is the last of the day
!> \}

!> \name Output variables
!> \brief HYSS make available variables to hold output values and 
!> information on output variables.
!> \{
  LOGICAL          :: conductload      !<flag for calculation of NP load
  INTEGER          :: max_outvar       !<number of defined output variables, set in initiate_output_variables
  INTEGER          :: max_outvar_trunk !<max_outvar in the trunk branch (outvarid of first additional in the CryoProcess is max_outvar_trunk+1)
  TYPE(IDTYPE), ALLOCATABLE :: outvarid(:)    !<information of output variables for this model
  REAL,ALLOCATABLE :: outvar(:,:)             !<variable to save output for print out in
  REAL,ALLOCATABLE :: outvar_classload(:,:,:,:)  !<variable to save class and substance dependent load outputs for yearly accumulation
  REAL,ALLOCATABLE :: outvar_basinload(:,:,:)    !<variable to save substance dependent load outputs for yearly accumulation
  INTEGER,ALLOCATABLE :: changecritvar(:,:)   !<variable to tell which variables to use for criteria calculation instead of those given as input
  CHARACTER (LEN=6),DIMENSION(max_classoutvar+max_basinoutvar+1) :: loadheadings  !<heading for load-files
  NAMELIST /HEADNML/ loadheadings
!> \}

!> \name Model parameter variables
!> \brief HYSS make available variables to hold information on and values for model parameters.
!> \{
  INTEGER          :: max_par           !<maximum number of defined model parameters
  TYPE(MODPARIDTYPE), ALLOCATABLE :: modparid(:)  !<information of model parameters for this model
  REAL,ALLOCATABLE :: genpar(:)         !<parameter values, general for whole model set up
  REAL,ALLOCATABLE :: basinpar(:,:)     !<parameter values for each subbasin
  REAL,ALLOCATABLE :: landpar(:,:)      !<parameter values for each land use
  REAL,ALLOCATABLE :: soilpar(:,:)      !<parameter values for each soil type
  REAL,ALLOCATABLE :: regpar(:,:)       !<parameter values for each parregion
  REAL,ALLOCATABLE :: wqregpar(:,:)     !<parameter values for each wqparregion
  REAL,ALLOCATABLE :: ilregpar(:,:)     !<parameter values for each ilakeregion
  REAL,ALLOCATABLE :: olregpar(:,:)     !<parameter values for each olakeregion
  REAL,ALLOCATABLE :: lregpar(:,:)      !<parameter values for each lake region
  REAL,ALLOCATABLE :: monthpar(:,:)     !<parameter values for each month
  REAL,ALLOCATABLE :: lakedatapar(:,:)  !<parameter values from LakeData.txt [FC, 14.04.11]
  INTEGER,ALLOCATABLE :: lakedataparindex(:,:)  !<table giving for each subbasin ilake and olake the row where the lake-related parameters are to be found in lakedatapar
!> \}

!TODO: All these updating variables can be collected in a type (except doupdate).
!> \name Variables for updating
!> \brief HYSS make available variables to hold information on updating. Several different updating functions exist.
!> \{
  LOGICAL :: doupdate(dim_update)               !<status of update functions
  INTEGER :: wobsvar                            !<variable (outvarid-index) to be used for updating wend and war
  LOGICAL, ALLOCATABLE :: updatestations(:)     !<update status of each subbasin (update function quseobs)
  LOGICAL, ALLOCATABLE :: updatestationsqar(:)  !<update status of each subbasin (update function qar)
  LOGICAL, ALLOCATABLE :: updatestationswar(:)  !<update status of each subbasin (update function war)
  REAL,    ALLOCATABLE :: updatestationsarfact(:) !<AR-factor of each subbasin (update function qar and war)
  REAL,    ALLOCATABLE :: updatetpcorr(:)       !<update parameter value for each subbasin (update function tpcorr)
  REAL,    ALLOCATABLE :: updatetncorr(:)       !<update parameter value for each subbasin (update function tncorr)
  REAL,    ALLOCATABLE :: updatetploccorr(:)    !<update parameter value for each subbasin (update function tploccorr)
  REAL,    ALLOCATABLE :: updatetnloccorr(:)    !<update parameter value for each subbasin (update function tnloccorr)
  LOGICAL, ALLOCATABLE :: wendupdstations(:)    !<update status of each subbasin (update function wendupd)
!> \}

CONTAINS

  !>Initiate the variable for output
  !>
  !> \b Consequences Module variables outvar,outvar_classload and 
  !>outvar_basinload may be allocated.
  !----------------------------------------------------------------
  SUBROUTINE allocate_outvar(n,nj,ns)

    INTEGER, INTENT(IN) :: n      !<number of subbasins
    INTEGER, INTENT(IN) :: nj     !<number of classes
    INTEGER, INTENT(IN) :: ns     !<number of substances

    IF(.NOT.ALLOCATED(outvar)) ALLOCATE(outvar(n,max_outvar))
    IF(conductload) THEN
       IF(.NOT.ALLOCATED(outvar_classload)) ALLOCATE(outvar_classload(nj,max_classoutvar,ns,n))
       IF(.NOT.ALLOCATED(outvar_basinload)) ALLOCATE(outvar_basinload(ns,max_basinoutvar,n))
    ENDIF

  END SUBROUTINE allocate_outvar

  !>Deallocate modvar arrays
  !>
  !> \b Consequences A lot of module variables are deallocated.
  !---------------------------------------
  SUBROUTINE deallocate_modvar(n)
  
  !Argument declarations
  INTEGER, INTENT(IN) :: n    !<Number of subbasins

  !Local variables
  INTEGER i
  
    IF(ALLOCATED(xobsi))      DEALLOCATE(xobsi)
    IF(ALLOCATED(xobsindex))  DEALLOCATE(xobsindex)
    IF(ALLOCATED(basin))      DEALLOCATE(basin)
    IF(ALLOCATED(lake))       DEALLOCATE(lake)  
    IF(ALLOCATED(lakebasin))  DEALLOCATE(lakebasin)
    IF(ALLOCATED(path))       DEALLOCATE(path)
    IF(ALLOCATED(aquifer))    DEALLOCATE(aquifer)
    IF(ALLOCATED(glacier))    DEALLOCATE(glacier)
    IF(ALLOCATED(load))THEN
      DO i = 1,n
        DEALLOCATE(load(i)%psvol)
        DEALLOCATE(load(i)%psload)
      ENDDO
      DEALLOCATE(load)
    ENDIF
    IF(ALLOCATED(outvar))     DEALLOCATE(outvar)
    IF(ALLOCATED(basinpar))   DEALLOCATE(basinpar)
    IF(ALLOCATED(landpar))    DEALLOCATE(landpar)
    IF(ALLOCATED(soilpar))    DEALLOCATE(soilpar)
    IF(ALLOCATED(regpar))     DEALLOCATE(regpar)
    IF(ALLOCATED(wqregpar))   DEALLOCATE(wqregpar)
    IF(ALLOCATED(ilregpar))   DEALLOCATE(ilregpar)
    IF(ALLOCATED(olregpar))   DEALLOCATE(olregpar)
    IF(ALLOCATED(lregpar))    DEALLOCATE(lregpar)
    IF(ALLOCATED(updatestations))       DEALLOCATE(updatestations)
    IF(ALLOCATED(updatestationsqar))    DEALLOCATE(updatestationsqar)
    IF(ALLOCATED(updatestationswar))    DEALLOCATE(updatestationswar)
    IF(ALLOCATED(updatestationsarfact)) DEALLOCATE(updatestationsarfact)
    IF(ALLOCATED(lakeindex))            DEALLOCATE(lakeindex)     
    IF(ALLOCATED(lakebasinindex))       DEALLOCATE(lakebasinindex)

  END SUBROUTINE deallocate_modvar

  !>Deallocate glacier arrays
  !>
  !> \b Consequences Module variables glacier is deallocated.
  !---------------------------------------
  SUBROUTINE deallocate_glacier()
  
    IF(ALLOCATED(glacier)) DEALLOCATE(glacier)

  END SUBROUTINE deallocate_glacier

  !>Allocate the model parameter variables
  !>
  !> \b Consequences Six module variables for model parameters are allocated.
  !------------------------------------------------------------------
  SUBROUTINE allocate_modelparameters(ns)

    INTEGER, INTENT(IN) :: ns     !<number of subbasins
    
    !Local variables
    INTEGER ngenpar,nsoilpar,nlandpar,nbasinpar,nregpar,nlregpar,nmonthpar,nwqregpar,nolregpar,nilregpar
    INTEGER i

    !> \b Algorithm \n
    !> Calculate number of model parameters of each dependence type
    ngenpar = 0; nsoilpar = 0; nlandpar = 0
    nbasinpar = 0; nregpar = 0; nlregpar = 0
    nmonthpar = 0; nwqregpar = 0 ;
    nilregpar=0; nolregpar=0
    DO i = 1,max_par
      IF(modparid(i)%deptype==m_gpar)   ngenpar   = MAX(ngenpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_spar)   nsoilpar  = MAX(nsoilpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_lpar)   nlandpar  = MAX(nlandpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_bpar)   nbasinpar = MAX(nbasinpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_rpar)   nregpar   = MAX(nregpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_wqrpar) nwqregpar = MAX(nwqregpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_ilrpar) nilregpar = MAX(nilregpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_olrpar) nolregpar = MAX(nolregpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_lrpar)  nlregpar  = MAX(nlregpar,modparid(i)%parno)
      IF(modparid(i)%deptype==m_mpar)   nmonthpar = MAX(nmonthpar,modparid(i)%parno)
    ENDDO

    !>Allocate variables for holding model parameter values
    IF(.NOT.ALLOCATED(genpar))   ALLOCATE(genpar(ngenpar))
    IF(.NOT.ALLOCATED(soilpar))  ALLOCATE(soilpar(nsoilpar,nsoil))
    IF(.NOT.ALLOCATED(landpar))  ALLOCATE(landpar(nlandpar,nluse))
    IF(.NOT.ALLOCATED(basinpar)) ALLOCATE(basinpar(nbasinpar,ns))
    IF(.NOT.ALLOCATED(lregpar))  ALLOCATE(lregpar(nlregpar,nlakeregions))
    IF(.NOT.ALLOCATED(regpar))   ALLOCATE(regpar(nregpar,nregions))
    IF(.NOT.ALLOCATED(wqregpar)) ALLOCATE(wqregpar(nwqregpar,nwqregions))
    IF(.NOT.ALLOCATED(ilregpar).AND.nilakeregions.GT.0) ALLOCATE(ilregpar(nilregpar,nilakeregions))
    IF(.NOT.ALLOCATED(olregpar).AND.nolakeregions.GT.0) ALLOCATE(olregpar(nolregpar,nolakeregions))
    IF(.NOT.ALLOCATED(monthpar)) ALLOCATE(monthpar(nmonthpar,12))

    !>Initialize variables to zero
    IF(ALLOCATED(genpar))   genpar   = 0.
    IF(ALLOCATED(soilpar))  soilpar  = 0.
    IF(ALLOCATED(landpar))  landpar  = 0.
    IF(ALLOCATED(basinpar)) basinpar = 0.
    IF(ALLOCATED(lregpar))  lregpar  = 0.
    IF(ALLOCATED(regpar))   regpar   = 0.
    IF(ALLOCATED(wqregpar)) wqregpar = 0.
    IF(ALLOCATED(ilregpar)) ilregpar = 0.
    IF(ALLOCATED(olregpar)) olregpar = 0.
    IF(ALLOCATED(monthpar)) monthpar = 0.

  END SUBROUTINE allocate_modelparameters

  !>Allocate variables for subbasin information
  !>
  !> \b Consequences Module variables for indata are allocated.
  !--------------------------------------------------------------
  SUBROUTINE allocate_basinvariables(ns,nsubst)

    INTEGER, INTENT(IN) :: ns     !<number of subbasins
    INTEGER, INTENT(IN) :: nsubst !<number of substances

    INTEGER i 
    
    IF(.NOT.ALLOCATED(basin))      ALLOCATE(basin(ns))
    IF(.NOT.ALLOCATED(pathsubid))  ALLOCATE(pathsubid(ns))
    ALLOCATE(load(ns))
    DO i = 1,ns
      ALLOCATE(load(i)%psvol(max_pstype))
      ALLOCATE(load(i)%psload(max_pstype,nsubst))
    ENDDO
    IF(.NOT.ALLOCATED(wetland))    ALLOCATE(wetland(ns,2))
    IF(.NOT.ALLOCATED(classbasin)) ALLOCATE(classbasin(ns,nclass))

  END SUBROUTINE allocate_basinvariables

  !>Initiate information about other observation series than PTQ
  !>
  !> \b Consequences Module variable xobsindex is allocated and initialized.
  !---------------------------------------------------------------
  SUBROUTINE initiate_xobsinformation(ns)

    INTEGER, INTENT(IN) :: ns         !<Number of subbasins

    IF(.NOT.ALLOCATED(xobsindex)) ALLOCATE(xobsindex(max_outvar,ns))
    xobsindex = 0

  END SUBROUTINE initiate_xobsinformation

  !>Collect information about finding other observation series
  !>
  !> \b Consequences Module variable xobsindex is set.
  !-------------------------------------------------------
  SUBROUTINE save_xobsinformation(n,var2,ns)

    INTEGER, INTENT(IN) :: n          !<Number of columns in Xobs
    INTEGER, INTENT(IN) :: var2(n,2)  !<Variables and subid for Xobs columns
    INTEGER, INTENT(IN) :: ns         !<Number of subbasins
    
    !Local variables
    INTEGER ivar,isub,icol

    DO ivar = 1, max_outvar
       DO isub = 1, ns
          DO icol = 1,n
             IF(var2(icol,1)==ivar .AND. var2(icol,2)==basin(isub)%subid) THEN
                xobsindex(ivar,isub) = icol
                EXIT
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE save_xobsinformation

  !>Calculate the soil layer depth and thicknesses.
  !>
  !>\b Consequences Module variable soillayerdepth and soillayerthick are set.
  !------------------------------------------------------------------------------
  INTEGER FUNCTION set_soillayerdepth(dim,depth)

    INTEGER, INTENT(IN) :: dim        !<Number of soil layers
    REAL, INTENT(IN)    :: depth(dim) !<Soil layer depth (m)
    !< \retval set_soillayerdepth error status of function
    
    !Local variables
    INTEGER i

    set_soillayerdepth = 0
    soillayerdepth(1:dim) = depth
    IF(dim<maxsoillayers)THEN
       soillayerdepth(dim+1:maxsoillayers) = soillayerdepth(dim)
    ENDIF
    soillayerthick = 0
    soillayerthick(1) = depth(1)
    DO i = 2,dim
       soillayerthick(i) = soillayerdepth(i) - soillayerdepth(i-1)
    ENDDO

  END FUNCTION set_soillayerdepth

  !>Calculate the soil layer depth and thicknesses.
  !>
  !>\b Consequences Module variable soildepth and soilthick are allocated and set.
  !----------------------------------------------------------------------------------
  INTEGER FUNCTION set_soildepth(nlayer,slc,depth)

    INTEGER, INTENT(IN) :: nlayer         !<number of soil layers
    INTEGER, INTENT(IN) :: slc            !<class index
    REAL, INTENT(IN)    :: depth(nlayer)  !<soil layer depth (m)
    !< \retval set_soildepth error status of function

    !Local variables
    INTEGER layer

    set_soildepth = 0
    IF(.NOT.ALLOCATED(soildepth)) ALLOCATE(soildepth(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(soilthick)) ALLOCATE(soilthick(maxsoillayers,nclass))

    soildepth(1:nlayer,slc) = depth
    DO layer = 1,nlayer
       IF(layer==1)THEN
          soilthick(layer,slc) = depth(layer)
       ELSE
          soilthick(layer,slc) = depth(layer)-soildepth(layer-1,slc)
       ENDIF
    ENDDO
    IF(nlayer>0.AND.nlayer<maxsoillayers)THEN
       soildepth(nlayer+1:maxsoillayers,slc) = soildepth(nlayer,slc)
       soilthick(nlayer+1:maxsoillayers,slc) = 0.
    ENDIF

  END FUNCTION set_soildepth

  !---------------------------------------------------------------
  !>\brief Find special coded classes
  !>
  !>\b Consequences Module variable slc_ilake, slc_olake, slc_lriver
  !> and slc_mriver are set, as is classmodel.
  !---------------------------------------------------------------
  SUBROUTINE set_coded_classes(dim,array,ilakearg,olakearg,lriverarg,mriverarg)

    INTEGER, INTENT(IN)  :: dim                 !<nclass
    REAL   , INTENT(IN)  :: array(dim)          !<special class column from GeoClass
    INTEGER, INTENT(OUT),OPTIONAL :: ilakearg   !<class index for ilake
    INTEGER, INTENT(OUT),OPTIONAL :: olakearg   !<clss index for olake
    INTEGER, INTENT(OUT),OPTIONAL :: lriverarg  !<class index for local river
    INTEGER, INTENT(OUT),OPTIONAL :: mriverarg  !<clss index for main river

    !Variable declaration
    INTEGER j                       !loop variable, class

    !>\b Algorithm \n
    !>Set special model code for all classes
    IF(.NOT.ALLOCATED(classmodel)) ALLOCATE(classmodel(dim))
    classmodel = NINT(array)

    slc_ilake = 0
    slc_olake = 0
    slc_lriver = 0
    slc_mriver = 0
    DO j = 1,dim
      IF(array(j)==1) slc_olake = j
      IF(array(j)==2) slc_ilake = j
      IF(array(j)==11) slc_lriver = j
      IF(array(j)==12) slc_mriver = j
    ENDDO

    !>Set coded lake classes for SourceApp
    IF(PRESENT(ilakearg)) ilakearg = slc_ilake
    IF(PRESENT(olakearg)) olakearg = slc_olake
    IF(PRESENT(lriverarg)) lriverarg = slc_lriver
    IF(PRESENT(mriverarg)) mriverarg = slc_mriver

  END SUBROUTINE set_coded_classes

  !---------------------------------------------------------------
  !>\brief Set variable for phosphorus and nitrogen correction
  !>
  !>\b Consequences Module variable doupdate is set.
  !---------------------------------------------------------------
  SUBROUTINE set_update_nutrientcorr(dim,ns,nrow,icorr,dataarray,updarray)

    INTEGER, INTENT(IN)  :: dim         !<number of subbasins of base model
    INTEGER, INTENT(IN)  :: ns          !<number of subbasins
    INTEGER, INTENT(IN)  :: nrow        !<number of data
    INTEGER, INTENT(IN)  :: icorr       !<index for type of correction
    REAL   , INTENT(IN)  :: dataarray(dim,2)        !<read data
    REAL,ALLOCATABLE,INTENT(INOUT)  :: updarray(:)  !<update correction value

    !Variable declaration
    INTEGER i,j

    !>\b Algorithm \n
    !>For each subbasin and data point: If equal set correction value for this subbasin
    IF(.NOT.ALLOCATED(updarray)) ALLOCATE(updarray(ns))
    updarray = 0.
    DO i = 1,ns
      DO j = 1, nrow
        IF(basin(i)%subid == NINT(dataarray(j,1))) THEN
          updarray(i)=dataarray(j,2)
          EXIT
        ENDIF
      ENDDO
    ENDDO
    !>If no correction value are found: turn off updating and deallocate array
    IF(.NOT.ANY(updarray/=0))THEN
      WRITE(6,*) 'WARNING: missing data in update.txt file, no ',TRIM(updatename(icorr)),' updating simulated'
      doupdate(icorr) = .FALSE.
      DEALLOCATE(updarray)
    ENDIF

  END SUBROUTINE set_update_nutrientcorr

  !>Initiate CropData variables to zero
  !>
  !>\b Consequences Module variable cropdata are set.
  !------------------------------------------------------------------------
  SUBROUTINE initiate_cropdata()

    cropdata(:)%fertnamount1  = 0      !fertiliser amount N (kg/km2)
    cropdata(:)%fertpamount1  = 0      !fertiliser amount P (kg/km2)
    cropdata(:)%fertday1      = 0      !day for fertilizing 1 (dayno)
    cropdata(:)%fertdown1     = 0      !part of fertilizer amount ploughed down to soillayer 2
    cropdata(:)%mannamount1   = 0      !manure amount N (kg/km2)
    cropdata(:)%manpamount1   = 0      !manure amount P (kg/km2)
    cropdata(:)%manday1       = 0      !day for manureing 1 (dayno)
    cropdata(:)%mandown1      = 0      !part of manure amount ploughed down to soillayer 2
    cropdata(:)%fertnamount2  = 0      !fertiliser amount N (kg/km2)
    cropdata(:)%fertpamount2  = 0      !fertiliser amount P (kg/km2)
    cropdata(:)%fertday2      = 0      !day for fertilizing 2 (dayno)
    cropdata(:)%fertdown2     = 0      !part of fertilizer amount ploughed down to soillayer 2
    cropdata(:)%mannamount2   = 0      !manure amount N (kg/km2)
    cropdata(:)%manpamount2   = 0      !manure amount P (kg/km2)
    cropdata(:)%manday2       = 0      !day for manureing 2 (dayno)
    cropdata(:)%mandown2      = 0      !part of manure amount ploughed down to soillayer 2
    cropdata(:)%resnamount    = 0      !residual amount of N (kg/km2)
    cropdata(:)%respamount    = 0      !residual amount of P (kg/km2)
    cropdata(:)%rescamount    = 0      !litter fall of C (kg/km2)
    cropdata(:)%resdayno      = 0      !day for residual (dayno), code 0 give residuals every day
    cropdata(:)%resdown       = 0      !part of residual amount ploughed down to soillayer 2
    cropdata(:)%resfast       = 0      !part of residual amount to fastN/P pool (rest to humusN/P)
    cropdata(:)%uptake1       = 0      !parameter for plant uptake of nutrient 1
    cropdata(:)%uptake2       = 0      !parameter for plant uptake of nutrient 2
    cropdata(:)%uptake3       = 0      !parameter for plant uptake of nutrient 3
    cropdata(:)%baredayno1    = 0      !day for beginning of first period with bare soil (dayno), typically 1
    cropdata(:)%baredayno2    = 0      !day for end of first period with bare soil (dayno), typically sawing date + a few days
    cropdata(:)%baredayno3    = 0      !day for beginning of second period with bare soil (dayno), typically ploughing date
    cropdata(:)%baredayno4    = 0      !day for end of second period with bare soil (dayno), typically 365
    cropdata(:)%baredayno5    = 0      !day for end of second period with bare soil (dayno), typically saw date of autumn crop
    cropdata(:)%ccmax1        = 0      !crop cover during summer
    cropdata(:)%ccmax2        = 0      !crop cover during winter and all year for year-round-crops (e.g. forest)
    cropdata(:)%gcmax1        = 0      !ground cover during summer
    cropdata(:)%gcmax2        = 0      !ground cover during winter and all year for year-round-crops (e.g. forest)
    cropdata(:)%uptakeupper   = 1.     !plant uptake from upper soillayer
    cropdata(:)%PNuptakeRatio = 0      !plant uptake P:N ratio
    cropdata(:)%daylength     = 0      !minimum daylength required for start to accumulate gdd
    cropdata(:)%gddsow        = 0      !accumulated growing degree days needed for sowing
    cropdata(:)%basetemp      = 0      !temperature which is deducted from air temperature when calculating gdd
    cropdata(:)%firstday      = 0      !first possible day for gdd accumulation
    cropirrdata(:)%plantingdayno   = 0    !dayno for irrigation season start

  END SUBROUTINE initiate_cropdata

  !>Set the index variable for CropData
  !>
  !>\b Consequences Module variable cropindex may be allocated.
  !---------------------------------------------------------------------
  SUBROUTINE set_cropdataindex(n,id,reg,cindex)

    INTEGER, INTENT(IN)  :: n         !<number of crops/rows in CropData
    INTEGER, INTENT(IN)  :: id(n)     !<cropid
    INTEGER, INTENT(IN)  :: reg(n)    !<region
    INTEGER, ALLOCATABLE, INTENT(INOUT)  :: cindex(:,:) !<cropindex
    
    !Local variables
    INTEGER irow
    INTEGER maxreg,maxid

    !Allocate and initialise cropindex      
    maxreg = MAXVAL(reg)
    maxid  = MAXVAL(id)
    IF(.NOT.ALLOCATED(cindex)) ALLOCATE(cindex(maxid,maxreg))
    cindex = 0

    DO irow = 1,n
      cindex(id(irow),reg(irow)) = irow
    ENDDO

  END SUBROUTINE set_cropdataindex

  !>Finds the corresponding row in CropData for the given crop and the
  !>fraction of the classarea its occupying.
  !---------------------------------------------------------------------
  SUBROUTINE find_croppart(i,j,kcrop,row,part)

    INTEGER, INTENT(IN)  :: i       !<index of current subbasin
    INTEGER, INTENT(IN)  :: j       !<index of current class
    INTEGER, INTENT(IN)  :: kcrop   !<main (1) or secondary (2) crop
    INTEGER, INTENT(OUT) :: row     !<row in CropData for this crop
    REAL,    INTENT(OUT) :: part    !<fraction of classarea for this crop
    
    !Local variables
    INTEGER c,r
    row = 0

    !> \b Algorithm \n
    !>For main crop: get the row from cropindex, areafraction is 1
    IF(kcrop==1)THEN
       c = classdata(j)%crop
       r = basin(i)%region
       IF(c>0.AND.r>0) row = cropindex(c,r)
       part = 1.
    !>For secondary crop: get the row from cropindex and areafraction from classbasin
    ELSE
       c = classdata(j)%crop2
       r = basin(i)%region
       IF(c>0.AND.r>0) row = cropindex(c,r)
       part = classbasin(i,j)%part2cr
    ENDIF

  END SUBROUTINE find_croppart

END MODULE MODVAR
