!> \file compout.f90
!> Contains module compout.

!>\brief The module compout contains procedures relating to calculating output and criteria.
!>
!>Procedures for for accumulation output data and for preparing output for file writing. 
!>Output are accumulated to the time period wanted for print out and for criteria calculations. 
!>Several different criteria is calculated.
MODULE COMPOUT

!Copyright 2011-2015 SMHI
!
!This file is part of HYPE.
!
!HYPE is free software: you can redistribute it and/or modify it under
!the terms of the Lesser GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or (at
!your option) any later version. HYPE is distributed in the hope that
!it will be useful, but WITHOUT ANY WARRANTY; without even the implied
!warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
!the Lesser GNU General Public License for more details. You should
!have received a copy of the Lesser GNU General Public License along
!with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !Uses modvar, worldvar, libdate, timeroutines, and convert

  IMPLICIT NONE
  PRIVATE 
  !----------------------------------------------------------------------
  ! Private procedures
  !----------------------------------------------------------------------
  ! compute_period 
  ! acc_one_data 
  ! stop_acc_data 
  ! check_and_acc_one_data 
  ! check_and_outmean_one_data 
  ! compute_yearfact 
  ! compute_sr2limit 
  ! find_acrit_corresponding_to_crit
  ! find_crit_corresponding_to_acrit 
  ! save_variables_for_later_critcalc 
  ! extend_ktvariables
  ! accumulate_periodvalues
  ! reset_periodvalues 
  ! kendallstau 
  ! calculate_fit 
  ! calculate_ra 
  ! calculate_median
  ! calculate_mean_criterion
  ! calculate_meanabs_criterion
  ! calculate_subbasin_mean_var
  ! calculate_nash_sutcliffe
  ! calculate_regional_nash_sutcliffe
  ! calculate_spatial_nash_sutcliffe
  ! calculate_all_subbasins_tau 
  ! calculate_all_subbasins_ra 
  ! calculate_regional_ra
  ! calculate_relative_error
  ! calculate_regional_relative_error
  ! calculate_errors
  ! calculate_regional_errors
  ! calculate_variance_based_criteria
  ! calculate_kling_gupta
  !-----------------------------------------------------------------------
  PUBLIC :: compute_outvar,&
            compute_mapvar,&
            compute_timevar,&
            compute_outloads,&
            prepare_to_compute_crit,&
            calculate_criteria,&
            find_variable_index, &
            calculate_median

  CONTAINS
  
  !>Compute number of timesteps in period so far, and check if period end.
  !--------------------------------------------------------------------
  SUBROUTINE compute_period(date,pcode,nper,pend,idt,ndt)

    USE WORLDVAR, ONLY : i_t,i_h,i_d,i_w,i_m,i_y,i_s,dtskip
    USE MODVAR, ONLY : seconds_per_timestep
    USE TIMEROUTINES
    USE LIBDATE

    !Argument declaration
    TYPE(DateType), INTENT(IN)  :: date  !<Current date
    INTEGER, INTENT(IN)  :: pcode        !<Code for period
    INTEGER, INTENT(OUT) :: nper         !<Number of values in current period
    LOGICAL, INTENT(OUT) :: pend         !<Flag for end of summation period
    INTEGER, INTENT(IN)  :: idt          !<Current time step number
    INTEGER, INTENT(IN)  :: ndt          !<Number of time steps in run
    
    !Variable declarations
    INTEGER h,d,m,y
    INTEGER dayofem,dayofew,dayofey           !Last day of month, week, and year
    INTEGER ntsperday                         !Number of timesteps per day
    LOGICAL endday,endweek,endmonth,endyear   !Indicators for end of the available accumulation periods
    
    INTEGER, PARAMETER :: seconds_per_hour = 3600
    INTEGER, PARAMETER :: seconds_per_day  = 86400 
    
    !> \b Algorithm \n
    !>Calculate day of year, day of month, last day of year etc.
    !>Calculate end of year, month, week and day.
    ntsperday=seconds_per_day/seconds_per_timestep
    h=date%Hour
    d=date%Day
    dayofew=MOD(idt,7*ntsperday)
    m=date%Month
    y=date%Year
    dayofem = day_of_month(y,m)
    endday = (h==(24-seconds_per_timestep/seconds_per_hour))
    endweek = (dayofew==0)
    endmonth = (dayofem==d .AND. endday)
    dayofey = day_of_month(y,12)
    endyear = (dayofey==d .AND. m==12 .AND. endday)
    
    !>Depending on current period: Calculate number of time steps in period and check if end of period
    pend=.FALSE.
    nper = 0
    SELECT CASE(pcode)
    CASE(i_h) !hourly output
      nper=1 !hour is currently the shortest timestep supported here
      pend=.TRUE.
    CASE(i_t) !timesteply output
      nper=1 
      pend=.TRUE.
    CASE(i_d) !daily output
      IF(endday)THEN
        nper=ntsperday  !Number of time steps per day
        pend=.TRUE.
      ENDIF
    CASE(i_w) !weekly output
      nper=7*ntsperday
      IF(endweek)THEN
        pend=.TRUE.
      ENDIF
    CASE(i_m) !monthly output
      nper=d*ntsperday
      IF(endmonth)THEN
        pend=.TRUE.
      ENDIF
    CASE(i_y) !yearly output
      IF(endyear)THEN
        pend=.TRUE.
        nper=NINT(numdays_of_year(date))*ntsperday
      ENDIF
    CASE(i_s) !output once for simulation period
      IF(idt==ndt)THEN
        nper=ndt-dtskip
        pend=.TRUE.
      ENDIF
    END SELECT
    
  END SUBROUTINE compute_period
  
  !>\brief Accumulate one data point
  !> Depending on type of output variable accumulate data
  !--------------------------------------------------------------------
  SUBROUTINE acc_one_data(sumtype,y,yhelp,acc,acch)
    
    USE MODVAR, ONLY : i_sum,   &
                       i_wmean, &
                       i_mean

    !Argument declaration
    INTEGER, INTENT(IN)  :: sumtype            !<Type of output variable
    REAL, INTENT(IN)     :: y                  !<Data
    REAL, INTENT(IN)     :: yhelp              !<Help data (water for weighting)
    REAL, INTENT(INOUT)  :: acc                !<Variable for accumulated data
    REAL, INTENT(INOUT)  :: acch               !<Variable for accumulated help data
    
    SELECT CASE(sumtype)
    CASE(i_wmean)
      acc = acc + y*yhelp
      acch = acch + yhelp
    CASE(i_mean)
      acc = acc + y
    CASE(i_sum)
      acc = acc + y
    END SELECT
    
  END SUBROUTINE acc_one_data
      
  !>Missing data, stop further accumulation for this period.
  !--------------------------------------------------------------------
  SUBROUTINE stop_acc_data(miss,acc,accok)
    
    !Argument declaration
    REAL, INTENT(IN)     :: miss                 !<Missing value
    REAL, INTENT(OUT)    :: acc                  !<Variable for accumulation of data
    LOGICAL, INTENT(OUT) :: accok                !<Variable for status of accumulation
    
    accok =.FALSE.
    acc = miss
    
  END SUBROUTINE stop_acc_data
  
  !>Check if data ok and then accumulate data.
  !TODO: This is not working correctly (what is that?) for observations and simulated values.
  !They should maybe be treated differently. 
  !TODO: No need to acc_one_data when simulation timestep is equal to print out meanperiod. 
  !Intoduce special case for that?
  !--------------------------------------------------------------------
  SUBROUTINE check_and_acc_one_data(simple,sumtype,y,yhelp,miss,acc,acchelp,accok)
    
    USE MODVAR, ONLY : i_sum,   &
                       i_wmean, &
                       i_mean

    !Argument declaration
    LOGICAL, INTENT(IN)    :: simple             !<Flag for timesteply output
    INTEGER, INTENT(IN)    :: sumtype            !<Accumulation type of data
    REAL, INTENT(IN)       :: y                  !<Data
    REAL, INTENT(IN)       :: yhelp              !<Help data (water for weighting)
    REAL, INTENT(IN)       :: miss               !<Missing value
    REAL, INTENT(INOUT)    :: acc                !<Variable for accumulated data
    REAL, INTENT(INOUT)    :: acchelp            !<Variable for accumulated help data
    LOGICAL, INTENT(INOUT) :: accok              !<Variable for status of accumulation
    
    !> \b Algorithm \n
    IF(simple)THEN
      !>If timsteply output use the data as is (weight equal to 1)
      CALL acc_one_data(sumtype,y,1.,acc,acchelp)
    ELSE  
      !>Else: Depending on data accumulation type: accumulate or end accumulation
      SELECT CASE(sumtype)
      CASE(i_mean,i_sum)
        IF(accok.AND.y/=miss)THEN
          CALL acc_one_data(sumtype,y,yhelp,acc,acchelp)
        ELSE
          CALL stop_acc_data(miss,acc,accok)
        ENDIF
      CASE(i_wmean)
        IF(accok)THEN
          IF(yhelp==0.)THEN
           !do nothing
          ELSEIF(y==miss)THEN
            CALL stop_acc_data(miss,acc,accok)
          ELSEIF(yhelp==miss.OR.yhelp==0.)THEN
            CALL acc_one_data(sumtype,y,1.,acc,acchelp)     !Use weight 1 if no observation of discharge or zero discharge
          ELSE
            CALL acc_one_data(sumtype,y,yhelp,acc,acchelp)
          ENDIF
        ENDIF
      END SELECT
    ENDIF
    
  END SUBROUTINE check_and_acc_one_data
      
  !>Compute the correct output average at the end of mean period for accumulated data
  !------------------------------------------------------------------------------
  SUBROUTINE check_and_outmean_one_data(sumtype,pcode,nper,miss,acc,acch,x)

    USE WORLDVAR, ONLY : i_s,           &
                         meandaysofyear
    USE MODVAR, ONLY : i_sum,   &
                       i_mean,  &
                       i_wmean, &
                       seconds_per_timestep
    
    !Argument declaration
    INTEGER, INTENT(IN) :: sumtype                 !<Code for variable average type
    INTEGER, INTENT(IN) :: pcode                   !<Code for period
    INTEGER, INTENT(IN) :: nper                    !<Period length
    REAL, INTENT(IN)    :: miss                    !<Missing value
    REAL, INTENT(IN)    :: acc                     !<Variable for accumulated data
    REAL, INTENT(IN)    :: acch                    !<Variable for accumulated help data
    REAL, INTENT(OUT)   :: x                       !<Average or sum of accumulated data
    
    !Variable declarations
    REAL div        !Denominator for average
    LOGICAL nodiv   !Flag for no division necessary
    
    !>\b Algorithm \n
    !>Calculate denominator for average depending on accumulation type
    nodiv = .FALSE.
    SELECT CASE(sumtype)
    CASE(i_mean)
      div = REAL(nper)
    CASE(i_wmean)
      div = acch
    CASE(i_sum)
      IF(pcode==i_s)THEN
        div = REAL(nper)/REAL(86400./seconds_per_timestep)/meandaysofyear
      ELSE
        div = 1.
        nodiv = .TRUE.
      ENDIF
    END SELECT
    IF(nper==1) nodiv=.TRUE.    !Timesteply output, no need for weighting
    
    !>Calculate average          
    x = miss
    IF(nodiv)THEN
      x = acc
    ELSE
      IF(acc.NE.miss .AND. div/=0.) THEN
        x = acc/div
      ENDIF
    ENDIF
    
  END SUBROUTINE check_and_outmean_one_data

  !>Accumulate and calculate the mean of period for timeserie output for specific subbasin
  !>
  !>\b Consequences Module variables accdata, accdatahelp, and accdata ok may change.
  !----------------------------------------------------------------------------------------------
  SUBROUTINE compute_outvar(d,isb,dim,x,pend,idt,ndt)

    USE WORLDVAR, ONLY : accdata,       &    !OUT
                         accdatahelp,   &    !OUT
                         dtskip,        &
                         o_nout,        &
                         o_wperiod,     &
                         outvarinfo,    &
                         outvarindex,   &
                         outvarbasinindex,    &
                         accdataok           !OUT
    USE MODVAR, ONLY : outvar,        &
                       outvarid,      &
                       missing_value, &
                       i_wmean
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    TYPE(DateType), INTENT(IN)  :: d    !<Current date
    INTEGER, INTENT(IN)  :: isb         !<Current subbasin
    INTEGER, INTENT(IN)  :: dim         !<Data dimension, number of variables for print out
    REAL, INTENT(OUT)    :: x(dim)      !<Data
    LOGICAL, INTENT(OUT) :: pend        !<Flag for end of summation period
    INTEGER, INTENT(IN)  :: idt         !<Current time step number
    INTEGER, INTENT(IN)  :: ndt         !<Number of time steps in run
    
    !Variable declarations
    INTEGER iad               !Subbasin index in accdata
    INTEGER ivar              !Loop index, variable
    REAL yhelp                !Corresponding water
    INTEGER nper              !Number of values in current period
    LOGICAL simple            !Flag for timesteply output
    
    !>\b Algorithm \n
    !>Skip warmup period
    IF(idt<=dtskip)THEN
      pend=.FALSE.
      x=0.0
      RETURN
    ENDIF
    
    !>Calculate period end and number of time steps included
    CALL compute_period(d,outvarinfo(1,o_wperiod),nper,pend,idt,ndt)
    simple = ((nper==1).AND. pend)
    yhelp = 0.
    iad = outvarbasinindex(isb)
    
    !>For every variable for print out:
    DO ivar = 1, outvarinfo(1,o_nout)
       
      !>\li Accumulate timestep data      
      IF(outvarindex(1,ivar,2)==i_wmean) yhelp = outvar(isb,outvarid(outvarindex(1,ivar,1))%water)
      CALL check_and_acc_one_data(simple,outvarindex(1,ivar,2),outvar(isb,outvarindex(1,ivar,1)),yhelp,missing_value,accdata(iad,ivar),accdatahelp(iad,ivar),accdataok(1,isb,ivar))
       
      !> \li If end of period: calculate output average
      IF(pend)THEN
        CALL check_and_outmean_one_data(outvarindex(1,ivar,2),outvarinfo(1,o_wperiod),nper,missing_value,accdata(iad,ivar),accdatahelp(iad,ivar),x(ivar))
        accdataok(1,isb,ivar) =.TRUE.
        accdata(iad,ivar)=0.0
        accdatahelp(iad,ivar)=0.0
      ELSE
        x(ivar)=0.0
      ENDIF
       
    ENDDO
    
  END SUBROUTINE compute_outvar
      
  !>Accumulate and calculate the mean of period for map output
  !>
  !>\b Consequences Module variables mapdata,maptime,tmap,mapdatahelp 
  !>and accdataok may change.
  !--------------------------------------------------------------------
  SUBROUTINE compute_mapvar(d,idt,ndt,dim,intformat,nmapperiod)

    USE WORLDVAR, ONLY : mapdata,      &    !OUT, Accumulated data for mapping
                         maptime,      &    !OUT, Time for data, stored as character
                         tmap,         &    !OUT
                         mapdatahelp,  &    !OUT
                         maxmap,       &
                         i_t,i_d,i_w,i_m,i_y,i_s,   &
                         dtskip,       &
                         o_wperiod,    &
                         outvarinfo,   &
                         outstartdate, &
                         outvarindex,  &
                         accdataok          !OUT
    USE MODVAR, ONLY : outvar,       &
                       outvarid,     &
                       nsub,         &
                       i_wmean,      &
                       missing_value
    USE LIBDATE, ONLY : DateType, format_date, OPERATOR(-)

    !Argument declaration
    TYPE(DateType), INTENT(IN) :: d       !<Current date
    INTEGER, INTENT(IN)  :: idt           !<Current time step number
    INTEGER, INTENT(IN)  :: ndt           !<Number of time steps in run
    INTEGER, INTENT(IN)  :: dim           !<Dimension data, number of variables for print out
    LOGICAL, INTENT(IN)  :: intformat     !<Write matlab format
    INTEGER, INTENT(OUT) :: nmapperiod    !<Number of periods for map print out
    
    !Variable declarations
    INTEGER ivar         !Loop index, variable
    INTEGER isb          !Loop index, subbasin
    REAL x               !Calculated average
    REAL yhelp           !Corresponding water
    INTEGER nper         !Number of values in current period
    LOGICAL pend         !Flag for end of summation period
    LOGICAL simple       !Flag for timesteply output
    TYPE(DateType) :: aweek
    
    !Initial check if calculation needed
    IF(dim==0) RETURN  !No mapfile print out

    !>\b Algorithm \n
    !>Skip warmup period
    IF(idt<=dtskip)THEN
      RETURN
    ENDIF
    
    !>Calculate period end and number of time steps included
    CALL compute_period(d,outvarinfo(2,o_wperiod),nper,pend,idt,ndt)
    simple = ((nper==1).AND. pend)
    yhelp = 0.0
    
    !>For every variable and subbasin for print out:
    DO ivar = 1,dim
      DO isb = 1,nsub
       !>\li Accumulate timestep data      
        IF(outvarindex(2,ivar,2)==i_wmean) yhelp = outvar(isb,outvarid(outvarindex(2,ivar,1))%water)
        CALL check_and_acc_one_data(simple,outvarindex(2,ivar,2),outvar(isb,outvarindex(2,ivar,1)),yhelp,missing_value,mapdata(isb,ivar,tmap),mapdatahelp(isb,ivar,tmap),accdataok(2,isb,ivar))
       !> \li If end of period: calculate output average
        IF(pend)THEN
          accdataok(2,isb,ivar) =.TRUE.
          CALL check_and_outmean_one_data(outvarindex(2,ivar,2),outvarinfo(2,o_wperiod),nper,missing_value,mapdata(isb,ivar,tmap),mapdatahelp(isb,ivar,tmap),x)
          mapdata(isb,ivar,tmap)=x
        ENDIF
      ENDDO
    ENDDO
    
    !>If end of period: Set the map time and prepare for next accumulation period
    IF(pend)THEN
      IF(intformat)THEN
        SELECT CASE(outvarinfo(2,o_wperiod))
        CASE(i_t)
          CALL format_date(d,'yyyymmddHHMM',maptime(tmap)) 
        CASE(i_d)
          CALL format_date(d,'yyyymmdd',maptime(tmap)) 
        CASE(i_w)
          aweek = DateType(0,0,6,0,0)
          CALL format_date(d-aweek,'yyyymmdd',maptime(tmap))  !beginning of week
        CASE(i_m)
          WRITE(maptime(tmap),'(I4,I2.2)') d%Year,d%Month
        CASE(i_y)
          WRITE(maptime(tmap),'(I4)') d%Year
        CASE(i_s)
          WRITE(maptime(tmap),'(I4,I4)') outstartdate%Year,d%Year  !eg. 19601999
        END SELECT
      ELSE
        SELECT CASE(outvarinfo(2,o_wperiod))
        CASE(i_t)
          CALL format_date(d,'yyyy-mm-dd HH:MM',maptime(tmap)) 
        CASE(i_d)
          CALL format_date(d,'yyyy-mm-dd',maptime(tmap)) 
        CASE(i_w)
          aweek = DateType(0,0,6,0,0)
          CALL format_date(d-aweek,'yyyy-mm-dd',maptime(tmap))  !beginning of week
        CASE(i_m)
          WRITE(maptime(tmap),'(I4,A,I2.2)') d%Year,'-',d%Month
        CASE(i_y)
          WRITE(maptime(tmap),'(I4)') d%Year
        CASE(i_s)
          WRITE(maptime(tmap),'(I4,A,I4)') outstartdate%Year,'-',d%Year  !eg. 1960-1999
        END SELECT
      ENDIF
      tmap=tmap+1
      IF(tmap<=maxmap)THEN
        mapdata(:,:,tmap)=0.0
        mapdatahelp(:,:,tmap)=0.0
      ENDIF
    ENDIF
    !>Calculate number of periods for map print out (so far)
    nmapperiod=tmap-1
    
  END SUBROUTINE compute_mapvar

  !>Accumulate and calculate the mean of period for timeserie output for one variable 
  !>
  !>\b Consequences Worldvar module variables taccdata,taccdatahelp 
  !>and accdataok may change.
  !--------------------------------------------------------------------
  SUBROUTINE compute_timevar(d,ivar,dim,x,pend,idt,ndt)

    USE WORLDVAR, ONLY : taccdata,       &    !OUT, Accumulated data
                         taccdatahelp,   &    !OUT
                         dtskip,         &
                         o_wperiod,      &
                         outvarinfo,     &
                         outvarindex,    &
                         accdataok            !OUT
    USE MODVAR, ONLY : outvar,         &
                       outvarid,       &
                       missing_value,  &
                       i_wmean
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    TYPE(DateType), INTENT(IN)  :: d      !<Current date
    INTEGER, INTENT(IN)  :: ivar          !<Current output variable
    INTEGER, INTENT(IN)  :: dim           !<Data dimension, number of subbasins
    REAL, INTENT(OUT)    :: x(dim)        !<Data
    LOGICAL, INTENT(OUT) :: pend          !<Flag for end of summation period
    INTEGER, INTENT(IN)  :: idt           !<Current time step number
    INTEGER, INTENT(IN)  :: ndt           !<Number of time steps in run
    
    !Variable declarations
    INTEGER isb          !Loop index, subbasins
    REAL yhelp           !Corresponding water for concentrations
    INTEGER nper         !Number of values in current period
    LOGICAL simple       !Flag for simple accumulation; timesteply output
    
    !>\b Algorithm \n
    !>Skip warmup period
    IF(idt<=dtskip)THEN
      pend=.FALSE.
      x=0.0
      RETURN
    ENDIF
    
    !>Calculate period end and number of time steps included
    CALL compute_period(d,outvarinfo(3,o_wperiod),nper,pend,idt,ndt)
    simple = ((nper==1).AND. pend)
    yhelp = 0.
    !>For every subbasin:
    DO isb = 1, dim
      !>\li Accumulate timestep data      
      IF(outvarindex(3,ivar,2)==i_wmean) yhelp = outvar(isb,outvarid(outvarindex(3,ivar,1))%water) !Set flow weight if flow weighted
      CALL check_and_acc_one_data(simple,outvarindex(3,ivar,2),outvar(isb,outvarindex(3,ivar,1)),yhelp,missing_value,taccdata(isb,ivar),taccdatahelp(isb,ivar),accdataok(3,isb,ivar))
      !> \li If end of period: calculate output average
      IF(pend)THEN
        accdataok(3,isb,ivar) =.TRUE. 
        CALL check_and_outmean_one_data(outvarindex(3,ivar,2),outvarinfo(3,o_wperiod),nper,missing_value,taccdata(isb,ivar),taccdatahelp(isb,ivar),x(isb))
        taccdata(isb,ivar)=0.0
        taccdatahelp(isb,ivar)=0.0
      ELSE
        x(isb)=0.0
      ENDIF
    ENDDO
    
  END SUBROUTINE compute_timevar
  
  !>Accumulate yearly totals of pollutant loads from daily values
  !>
  !>\b Consequences Module variables accdata_classload and accdata_basinload 
  !>may change.
  !--------------------------------------------------------------------
  SUBROUTINE compute_outloads(d,pend,idt,ndt)   

    USE WORLDVAR, ONLY : accdata_classload, &     !OUT, Accumulated data (class dependent)
                         accdata_basinload, &     !OUT, Accumulated data (non-class dependent)
                         dtskip,        &
                         i_y        
    USE MODVAR,   ONLY : outvar_classload,  &
                         outvar_basinload,  &  
                         max_classoutvar,   &
                         max_basinoutvar
    USE LIBDATE,  ONLY : DateType
    
    !Argument declaration
    TYPE(DateType), INTENT(IN)  :: d    !<Current date
    LOGICAL, INTENT(OUT) :: pend        !<Flag for end of summation period
    INTEGER, INTENT(IN)  :: idt         !<Current time step number
    INTEGER, INTENT(IN)  :: ndt         !<Number of time steps in run
    
    !Variable declarations
    INTEGER nper                    !Number of values in current period
    
    !>\b Algorithm \n
    !>Skip warmup period
    IF(idt<=dtskip)THEN
      pend=.FALSE.
      RETURN
    ENDIF

    !>Calculate period end and number of time steps included
    CALL compute_period(d,i_y,nper,pend,idt,ndt) 

    !>Accumulate data
    accdata_classload(:,:,:,:) = accdata_classload(:,:,:,:) + outvar_classload(:,:,:,:)
    accdata_basinload(:,:,:)   = accdata_basinload(:,:,:) + outvar_basinload(:,:,:)  
    
  END SUBROUTINE compute_outloads
  
  !>Accumulate data for criteria calculation 
  !>
  !>\b Consequences Worldvar module variables critvec,rs,sc and ts may 
  !>be allocated and/or changed
  !--------------------------------------------------------------------
  SUBROUTINE prepare_to_compute_crit(d,idt,ndt)

    USE WORLDVAR, ONLY : nacrit,           &
                         acccalvar,        &
                         calvarper,        &
                         dtskip,           &
                         critvec,          &  !OUT
                         rs,cs,ts             !OUT
    USE MODVAR, ONLY : nsub,               &
                       missing_value,    &
                       outvar
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    TYPE(DateType), INTENT(IN) :: d  !<Current date 
    INTEGER, INTENT(IN) :: idt       !<Current time step number
    INTEGER, INTENT(IN) :: ndt       !<Number of time steps in run
    
    !Variable declarations
    INTEGER i,isb
    INTEGER p
    INTEGER nper                    !Number of values in current period
    LOGICAL pend                    !Flag for end of summation period
    REAL    r(nsub),c(nsub)         !Recorded and computed value of variable to be used for criteria
    
    !>\b Algorithm \n
    !>Skip warmup period
    IF(idt<=dtskip)THEN
      RETURN
    ENDIF
    !>Allocate criteria accumulation variables if not already done
    IF(.NOT.ALLOCATED(critvec)) ALLOCATE(critvec(9,nsub,nacrit))  !DG additional optimization criterias
    IF(.NOT.ALLOCATED(rs)) ALLOCATE(rs(nsub,nacrit))
    IF(.NOT.ALLOCATED(cs)) ALLOCATE(cs(nsub,nacrit))
    IF(.NOT.ALLOCATED(ts)) ALLOCATE(ts(nsub,nacrit))
    
    !>If first time step of accumulation period: initialise accumulation variables
    IF(idt==1.OR.idt==dtskip+1)THEN
      critvec = 0.D0
      CALL reset_periodvalues(rs,cs,ts,nsub,nacrit)
    ENDIF
    
    !>Calculate period end and number of time steps included
    p=calvarper
    CALL compute_period(d,p,nper,pend,idt,ndt)
    !>For every criteria with unique variables:
    DO i = 1,nacrit
      r(1:nsub) = outvar(1:nsub,acccalvar(i)%rec)
      c(1:nsub) = outvar(1:nsub,acccalvar(i)%comp)
       
      !>Accumulate data for period
      CALL accumulate_periodvalues(missing_value,r(:),c(:),rs(1:nsub,i),&
            cs(1:nsub,i),ts(1:nsub,i),nsub)
      !>If end of period: 
      IF(pend)THEN
        !>\li Upscale data for periods with not all observations
        WHERE(ts(1:nsub,i)>0 .AND. ts(1:nsub,i)<nper)
          rs(1:nsub,i) = rs(1:nsub,i) * nper / ts(1:nsub,i)
          cs(1:nsub,i) = cs(1:nsub,i) * nper / ts(1:nsub,i)
          ts(1:nsub,i) = nper
        ENDWHERE
        CALL reform_accumulated_periodvalues(nsub,acccalvar(i)%vartype,rs(1:nsub,i), &  
            cs(1:nsub,i),ts(1:nsub,i))     !Reform some variables to be average instead of sum
        !>\li If needed: save data for later (for Kendalls Tau and RA the variables need to be saved)
        IF(acccalvar(i)%saveend)THEN
          CALL save_variables_for_later_critcalc(nsub,i,cs(1:nsub,i),rs(1:nsub,i),ts(1:nsub,i))
        ENDIF
        !>\li For every subbasin: accumulate data period data for later criteria calculations
        DO isb = 1,nsub
          IF(ts(isb,i)>0)THEN                 !Accumulations for computation of critera
            critvec(1,isb,i) = critvec(1,isb,i) + 1.D0
            critvec(2,isb,i) = critvec(2,isb,i) + DBLE(rs(isb,i))
            critvec(3,isb,i) = critvec(3,isb,i) + DBLE(rs(isb,i) * rs(isb,i))
            critvec(4,isb,i) = critvec(4,isb,i) + DBLE(cs(isb,i) - rs(isb,i))
            critvec(5,isb,i) = critvec(5,isb,i) + DBLE(rs(isb,i) * cs(isb,i))
            critvec(6,isb,i) = critvec(6,isb,i) + DBLE(ABS(rs(isb,i)-cs(isb,i)))
            critvec(7,isb,i) = critvec(7,isb,i) + DBLE(cs(isb,i) * cs(isb,i))
            critvec(8,isb,i) = dmax1(critvec(8,isb,i),ABS(DBLE(cs(isb,i)))) ! Max Computed for NRMSE [David, 20130221]
            critvec(9,isb,i) = dmax1(critvec(9,isb,i),ABS(DBLE(rs(isb,i)))) ! Max Recorded for NRMSE
          ENDIF
        ENDDO
        !>\li (Re-)initialize accumulation variables
        CALL reset_periodvalues(rs(1:nsub,i),cs(1:nsub,i),ts(1:nsub,i),nsub,1)
      ELSE
      ENDIF
    ENDDO
    
  END SUBROUTINE prepare_to_compute_crit

  !>Calculates objective function (chosen multi-criteria) for the last time step.
  !>Calculates optional other criteria 
  !>
  !>\b Consequences Worldvar module variables ktcomp,ktrec and ktnum etc. may be 
  !>allocated and changed (sort of?).
  !--------------------------------------------------------------------
  SUBROUTINE calculate_criteria(crit0,basincrit,simcrit,crit1,thres1)

    USE WORLDVAR, ONLY : calvar,         &
                         ncrit,          &
                         nacrit,         &
                         maxperf,        &
                         maxsubass,      &
                         i_rnse,i_snse,  &
                         i_mnse,         &
                         i_rmae,i_sre,   &  
                         i_rre,i_mre,    &
                         i_rra,i_sra,    &
                         i_mra,i_tau,    &
                         i_mdnse,i_mdra, &
                         i_mstdre,i_mcc, &
                         i_mdkg,i_mabsre, &
                         i_mnrmse, &
                         ktcomp, ktrec, ktnum,     &    !OUT, sort of
                         ktcomp2,ktrec2,ktnum2,    &    !OUT, sort of
                         ktcomp3,ktrec3,ktnum3,    &    !OUT, sort of
                         ktcomp4,ktrec4,ktnum4          !OUT, sort of
    USE MODVAR, ONLY : nsub,               &
                       i_sum,              &
                       missing_value
    
    !Argument declaration
    REAL, INTENT(OUT)   :: crit0                  !<Value of chosen objective function
    REAL, OPTIONAL, INTENT(OUT)   :: basincrit(nsub,maxsubass,nacrit)    !<Performance criteria per subbasin and variable
    REAL, OPTIONAL, INTENT(OUT)   :: simcrit(maxperf,nacrit)             !<Performance criteria per variable
    REAL, OPTIONAL, INTENT(OUT)   :: crit1        !<Conditional criteria
    REAL, OPTIONAL, INTENT(OUT)   :: thres1       !<Threshold of conditional criteria

    !Parameter declaration
    INTEGER, PARAMETER :: sr2limityear = 5   !Minimum years with data for calculation of QR,QC and SR2
    INTEGER, PARAMETER :: sr2limitbasin = 5  !Minimum number of subbasins for calculation of SR2 etc.
    
    !Variable declarations
    LOGICAL calcsubass,calcsimass
    CHARACTER(LEN=3) str
    INTEGER i             !index criteria variables (ncrit)
    INTEGER iac           !index accumulation criteria variables (nacrit)
    REAL weight         !weigth of criterion
    REAL spatialcrit(3)                    !Spatial criteria: R2, RMSE, bias
    REAL   :: critadd            !criterion value to be added
    REAL   :: stdr(nsub),stdc(nsub),stderr(nsub)
    REAL   :: locrs(nsub)        !Relative error of standard deviation
    REAL   :: loccc(nsub)        !Correlation coefficient
    REAL   :: locmrs             !Relative error of standard deviation, mean of subbasins
    REAL   :: locrve             !Relative error of mean value, regional
    REAL   :: locwve             !Relative error of mean value, mean of subbasins
    REAL   :: locmare            !Absolute error in percent of recorded series, mean of subbasins DG130130
    REAL   :: locmkg             !Kling-Gupta Efficiecy, median of subbasins, DG130130
    REAL   :: locr2(nsub)           !R2 per subbasin
    REAL   :: locve(nsub)           !Error in percent of recorded series
    REAL   :: locqr(nsub)           !Recorded total (discharge)
    REAL   :: locqc(nsub)           !Computed total (discharge)
    REAL   :: locmae(nsub)          !Mean absolute error
    REAL   :: locrmse(nsub)         !Root mean square error
    REAL   :: loctau(nsub)          !Kendalls Tau, non-parametric criteria
    REAL   :: locra(nsub)           !RA per subbasin
    REAL   :: lockge(nsub)          !Kling-Gupta Efficiency per subbasin
    REAL   :: lockge1(nsub)         !Kling-Gupta Efficiency part 1 (cc) per subbasin
    REAL   :: lockge2(nsub)         !Kling-Gupta Efficiency part 2 (std-quotient) per subbasin
    REAL   :: lockge3(nsub)         !Kling-Gupta Efficiency part 3 (mean-quotient) per subbasin
    REAL   :: locnrmse(nsub)        !Normalized RMSE per subbasin
    REAL   :: thres                 !Conditional criteria threshold    
    
    !>\b Algorithm \n
    !>Initializations
    !Find which calculations is asked for
    IF(PRESENT(basincrit))THEN
      calcsubass = .TRUE.
    ELSE
      calcsubass = .FALSE.
    ENDIF
    IF(PRESENT(simcrit))THEN
      calcsimass = .TRUE.
    ELSE
      calcsimass = .FALSE.
    ENDIF
    !Set optional output variables to missing_value
    IF(calcsimass) simcrit = missing_value
    
    !Calculate subbasin and simulation assessment
    !--------------------------------------------
    !>For every criteria with unique variables:
    DO iac = 1, nacrit
       
      !>If assessment print out: Calculate diagnostic criteria (R2, RE etc) for each subbasin and subbasin average value of variables
      IF(calcsubass.OR.calcsimass)THEN
        CALL find_crit_corresponding_to_acrit(iac,i)
        CALL calculate_nash_sutcliffe(iac,nsub,locr2)       !R2
        CALL calculate_relative_error(iac,nsub,locve)       !Relative error
        CALL calculate_variance_based_criteria(iac,nsub,stdr,stdc,stderr,locrs,loccc)
        CALL calculate_kling_gupta(iac,nsub,lockge,lockge1,lockge2,lockge3)         !KGE
        CALL calculate_mean_and_errors(iac,nsub,locqr,locqc,locmae,locrmse,locnrmse) !MAE, RMSE, and NRMSE
      ENDIF
       
      IF(calcsubass)THEN
        !>If subbasin assessment: Set optional subbasin assessment criteria dummy argument
        basincrit(:,1,iac) = locr2(:)
        basincrit(:,2,iac) = loccc(:)
        basincrit(:,3,iac) = locve(:) 
        WHERE(basincrit(:,3,iac)/=missing_value) basincrit(:,3,iac) = basincrit(:,3,iac) * 100.    !%
        basincrit(:,4,iac) = locrs(:)
        WHERE(basincrit(:,4,iac)/=missing_value) basincrit(:,4,iac) = basincrit(:,4,iac) * 100.    !%
        basincrit(:,5,iac) = locqc(:)
        basincrit(:,6,iac) = locqr(:)
        basincrit(:,7,iac) = stdc(:)
        basincrit(:,8,iac) = stdr(:)
        basincrit(:,9,iac) = locmae(:)
        basincrit(:,10,iac) = locrmse(:)
        basincrit(:,11,iac) = locqc(:)-locqr(:)
        basincrit(:,12,iac) = stderr(:)
        basincrit(:,13,iac) = lockge(:)
        basincrit(:,14,iac) = lockge2(:)
        basincrit(:,15,iac) = lockge3(:)
        basincrit(:,16,iac) = locnrmse(:)
      ENDIF
       
      IF(calcsimass)THEN
        !>If simulation assessment: Calculate and set simulation assessment criteria; average and median of subbasin criteria
        CALL calculate_mean_criterion(nsub,locr2,simcrit(i_mnse,iac))            !Aritmetic mean of NSE
        CALL calculate_mean_criterion(nsub,locve,simcrit(i_mre,iac))             !Aritmetic mean of relative error
        CALL calculate_mean_criterion(nsub,locrs,simcrit(i_mstdre,iac))          !Aritmetic mean of relative standarddeviation error
        CALL calculate_mean_criterion(nsub,loccc,simcrit(i_mcc,iac))             !Aritmetic mean of correlation coefficients
        CALL calculate_meanabs_criterion(nsub,locve,simcrit(i_mabsre,iac))       !Aritmetic mean of Absolute Relative Error
        CALL calculate_median(nsub,locr2,missing_value,simcrit(i_mdnse,iac))     !Median NSE
        CALL calculate_median(nsub,lockge,missing_value,simcrit(i_mdkg,iac))     !Median KGE
        CALL calculate_median(nsub,locnrmse,missing_value,simcrit(i_mnrmse,iac)) !Median of Normalized RMSE
        !<If simulation assessment: Calculate and set simulation assessment criteria, spatial and combined ("regional")
        CALL calculate_regional_errors(iac,simcrit(i_rmae,iac))
        CALL calculate_regional_relative_error(iac,simcrit(i_rre,iac))
        CALL calculate_regional_nash_sutcliffe(iac,simcrit(i_rnse,iac))
        CALL calculate_subbasin_mean_var(i,iac,sr2limityear,nsub,locqr,locqc)
        CALL calculate_fit(nsub,locqc,locqr,missing_value,sr2limitbasin,spatialcrit)          !Spatial R2,RMSE and bias
        !>If simulation assessment: Set optional simulation assessment criteria dummy argument
        simcrit(i_snse,iac) = spatialcrit(1)
        simcrit(i_sre,iac)  = spatialcrit(3)
      ENDIF
    ENDDO

    !Calculate the multi-criteria optimization criterion, the objective function
    crit0 = 0.
    !>For every criterion included in the objective function:
    DO i = 1, ncrit
      IF(.NOT.calvar(i)%cond)THEN
        str = calvar(i)%crit
        weight = calvar(i)%weight
        CALL find_acrit_corresponding_to_crit(i,iac)
        !>\li If TAU: Calculate TAU for each subbasin, calculate mean TAU, add to objective function
        IF(str=='TAU')THEN                      !Kendalls Tau, non-parametric criteria
          CALL calculate_all_subbasins_tau(iac,nsub,loctau)
          CALL calculate_mean_criterion(nsub,loctau,critadd)     !Aritmetic mean of tau
          IF(calcsimass) simcrit(i_tau,iac) = critadd
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If RRA: Calculate regional RA, add to objective function
        IF(str=='RRA')THEN                      !Regional RA
          CALL calculate_regional_ra(i,iac,critadd)
          IF(calcsimass) simcrit(i_rra,iac) = critadd
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If SRA: Calculate variable average for each subbasin, calculate RA on these, add to objective function
        IF(str=='SRA')THEN                      !Spatial RA
          CALL calculate_subbasin_mean_var(i,iac,sr2limityear,nsub,locqr,locqc)
          CALL calculate_ra(nsub,locqc(:),locqr(:),missing_value,calvar(i)%coeff,sr2limitbasin,critadd)
          IF(calcsimass) simcrit(i_sra,iac) = critadd
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If SR2: Calculate spatial NSE if not already done, add to objective function
        IF(str=='SR2')THEN                      !Spatial R2 (NSE)
          IF(calcsimass)THEN
            critadd = simcrit(i_snse,iac)
          ELSE
            CALL calculate_spatial_nash_sutcliffe(i,iac,sr2limityear,sr2limitbasin,critadd)
          ENDIF
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If RR2: Calculate regional NSE if not already done, add to objective function
        IF(str=='RR2')THEN                      !Regional R2 (NSE)
          IF(calcsimass)THEN
            critadd = simcrit(i_rnse,iac)
          ELSE
            CALL calculate_regional_nash_sutcliffe(iac,critadd)
          ENDIF
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If MR2: Calculate NSE for every subbasin and their average if not already done, add to objective function
        IF(str=='MR2')THEN                      !R2 (NSE), mean of all subbasins
          IF(calcsimass)THEN
            critadd = simcrit(i_mnse,iac)
          ELSE
            CALL calculate_nash_sutcliffe(iac,nsub,locr2)       
            CALL calculate_mean_criterion(nsub,locr2,critadd)   
          ENDIF
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If MD2: Calculate NSE for every subbasin and their median if not already done, add to objective function
        IF(str=='MD2')THEN                      !R2 (NSE), median of all subbasins
          IF(calcsimass)THEN
            critadd = simcrit(i_mdnse,iac)
          ELSE
            CALL calculate_nash_sutcliffe(iac,nsub,locr2)    
            CALL calculate_median(nsub,locr2,missing_value,critadd)   
          ENDIF
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If MRA: Calculate RA for every subbasin, calculate their average, add to objective function
        IF(str=='MRA')THEN                      !RA, mean of all subbasins
          CALL calculate_all_subbasins_ra(i,iac,nsub,locra)
          CALL calculate_mean_criterion(nsub,locra,critadd)   
          IF(calcsimass) simcrit(i_mra,iac) = critadd
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If MDA: Calculate RA for every subbasin, calculate their median, add to objective function
        IF(str=='MDA')THEN                      !RA, median of all subbasins
          CALL calculate_all_subbasins_ra(i,iac,nsub,locra)
          CALL calculate_median(nsub,locra,missing_value,critadd)
          IF(calcsimass) simcrit(i_mdra,iac) = critadd
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If RRE: Calculate regional RE if not already done, add to objective function
        IF(str=='RRE')THEN                      !Regional Relative error of mean
          IF(calcsimass)THEN
            locrve = simcrit(i_rre,iac)
          ELSE
            CALL calculate_regional_relative_error(iac,locrve)
          ENDIF
          crit0 = crit0 + ABS(locrve) * weight
        ENDIF
        !>\li If MRE: Calculate RE for every subbasin and their average if not already done, add to objective function
        IF(str=='MRE')THEN                      !Relative error of mean, mean over all subbasins
          IF(calcsimass)THEN
            locwve = simcrit(i_mre,iac)
          ELSE
            CALL calculate_relative_error(iac,nsub,locve)        !Relative error
            CALL calculate_mean_criterion(nsub,locve,locwve)     !Aritmetic mean of RE
          ENDIF
          crit0 = crit0 + ABS(locwve) * weight
        ENDIF
        !>\li If MAR: Calculate RE for every subbasin and their average absolute value if not already done, add to objective function
        IF(str=='MAR')THEN                      !Mean Absolute Relative error, DG130130/CP130212
          IF(calcsimass)THEN
            locmare = simcrit(i_mabsre,iac)
          ELSE
            CALL calculate_relative_error(iac,nsub,locve)           !Relative error
            CALL calculate_meanabs_criterion(nsub,locve,locmare)    !Aritmetic mean of Absolute RE
          ENDIF
          crit0 = crit0 + locmare * weight
        ENDIF
        !>\li If MCC: Calculate correlation coefficient for every subbasin and their average if not already done, add to objective function
        IF(str=='MCC')THEN                      !Correlation coefficient, mean over all subbasins
          IF(calcsimass)THEN
            critadd = simcrit(i_mcc,iac)
          ELSE
            CALL calculate_variance_based_criteria(iac,nsub,stdr,stdc,stderr,locrs,loccc)
            CALL calculate_mean_criterion(nsub,loccc,critadd)  
          ENDIF
          crit0 = crit0 - critadd * weight
        ENDIF
        !>\li If MKG: Calculate Kling-Gupta Efficiency for every subbasin and their median if not already done, add to objective function
        IF(str=='MKG')THEN                      !Median Kling-Gupta Efficiency, DG130130/CP130212
          IF(calcsimass)THEN
            locmkg = simcrit(i_mdkg,iac)
          ELSE
            CALL calculate_kling_gupta(iac,nsub,lockge,lockge1,lockge2,lockge3)             !KGE
            CALL calculate_median(nsub,lockge,missing_value,locmkg)    !Median of KGE
          ENDIF
          crit0 = crit0 - locmkg * weight
        ENDIF
        !>\li If MRS: Calculate error in standard deviations for every subbasin and their average if not already done, add to objective function
        IF(str=='MRS')THEN                      !Relative error of standard deviation, mean of all subbasins
          IF(calcsimass)THEN
            locmrs = simcrit(i_mstdre,iac)
          ELSE
            CALL calculate_variance_based_criteria(iac,nsub,stdr,stdc,stderr,locrs,loccc)
            CALL calculate_mean_criterion(nsub,locrs,locmrs)   
          ENDIF
          crit0 = crit0 + ABS(locmrs) * weight
        ENDIF
        !>\li If MNR: Calculate normalised RMSE for every subbasin and their median if not already done, add to objective function
        IF(str=='MNR')THEN                      !Median normalized RMSE 
          IF(calcsimass)THEN
            critadd = simcrit(i_mnrmse,iac) 
          ELSE
            CALL calculate_mean_and_errors(iac,nsub,locqr,locqc,locmae,locrmse,locnrmse)   !MAE and RMSE
            CALL calculate_median(nsub,locnrmse,missing_value,critadd)  !Median of Normalized RMSE 
            !CALL calculate_mean_criterion(nsub,locnrmse(:),critadd)        !Aritmetic mean of Normalized RMSE
          ENDIF
          crit0 = crit0 + critadd * weight
        ENDIF
      ENDIF
    ENDDO
    
    !>Calculation of conditional criterion
    IF(PRESENT(thres1))THEN
      crit1 = 0.
      thres1 = 0.
       
      DO i = 1, ncrit
        IF(calvar(i)%cond)THEN

          str = calvar(i)%crit
          weight = calvar(i)%weight
          thres = calvar(i)%thres
          CALL find_acrit_corresponding_to_crit(i,iac)
               
          IF(str=='TAU')THEN                      !Kendalls Tau, non-parametric criteria
            CALL calculate_all_subbasins_tau(iac,nsub,loctau)
            CALL calculate_mean_criterion(nsub,loctau,critadd)     !Aritmetic mean of tau
            IF(calcsimass) simcrit(i_tau,iac) = critadd
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF

          IF(str=='RRA')THEN                      !Regional RA
            CALL calculate_regional_ra(i,iac,critadd)
            IF(calcsimass) simcrit(i_rra,iac) = critadd
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF

          IF(str=='SRA')THEN                      !Spatial RA
            CALL calculate_subbasin_mean_var(i,iac,sr2limityear,nsub,locqr,locqc)
            CALL calculate_ra(nsub,locqc,locqr,missing_value,calvar(i)%coeff,sr2limitbasin,critadd)
            IF(calcsimass) simcrit(i_sra,iac) = critadd
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
                
          IF(str=='SR2')THEN                      !Spatial R2 (NSE)
            IF(calcsimass)THEN
              critadd = simcrit(i_snse,iac)
            ELSE
              CALL calculate_spatial_nash_sutcliffe(i,iac,sr2limityear,sr2limitbasin,critadd)
            ENDIF
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='RR2')THEN                      !Regional R2 (NSE)
            IF(calcsimass)THEN
              critadd = simcrit(i_rnse,iac)
            ELSE
              CALL calculate_regional_nash_sutcliffe(iac,critadd)
            ENDIF
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MR2')THEN                      !R2 (NSE), mean of all subbasins
            IF(calcsimass)THEN
              critadd = simcrit(i_mnse,iac)
            ELSE
              CALL calculate_nash_sutcliffe(iac,nsub,locr2)       
              CALL calculate_mean_criterion(nsub,locr2,critadd)   
            ENDIF
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MD2')THEN                      !R2 (NSE), median of all subbasins
            IF(calcsimass)THEN
              critadd = simcrit(i_mdnse,iac)
            ELSE
              CALL calculate_nash_sutcliffe(iac,nsub,locr2)    
              CALL calculate_median(nsub,locr2,missing_value,critadd)   
            ENDIF
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MRA')THEN                      !RA, mean of all subbasins
            CALL calculate_all_subbasins_ra(i,iac,nsub,locra)
            CALL calculate_mean_criterion(nsub,locra,critadd)   
            IF(calcsimass) simcrit(i_mra,iac) = critadd
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MDA')THEN                      !RA, median of all subbasins
            CALL calculate_all_subbasins_ra(i,iac,nsub,locra)
            CALL calculate_median(nsub,locra,missing_value,critadd)
            IF(calcsimass) simcrit(i_mdra,iac) = critadd
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='RRE')THEN                      !Regional Relative error of mean
            IF(calcsimass)THEN
              locrve = simcrit(i_rre,iac)
            ELSE
              CALL calculate_regional_relative_error(iac,locrve)
            ENDIF
            crit1 = crit1 + ABS(locrve) * weight
            thres1 = thres1 + thres * weight
          ENDIF
               
          IF(str=='MRE')THEN                      !Relative error of mean, mean over all subbasins
            IF(calcsimass)THEN
              locwve = simcrit(i_mre,iac)
            ELSE
              CALL calculate_relative_error(iac,nsub,locve)        !Relative error
              CALL calculate_mean_criterion(nsub,locve,locwve)     !Aritmetic mean of RE
            ENDIF
            crit1 = crit1 + ABS(locwve) * weight
            thres1 = thres1 + thres * weight
          ENDIF
               
          IF(str=='MAR')THEN                      !Mean Absolute Relative error, DG130130/CP130212
            IF(calcsimass)THEN
              locmare = simcrit(i_mabsre,iac)
            ELSE
              CALL calculate_relative_error(iac,nsub,locve)           !Relative error
              CALL calculate_meanabs_criterion(nsub,locve,locmare)    !Aritmetic mean of Absolute RE
            ENDIF
            crit1 = crit1 + locmare * weight
            thres1 = thres1 + thres * weight
          ENDIF
               
          IF(str=='MCC')THEN                      !Correlation coefficient, mean over all subbasins
            IF(calcsimass)THEN
              critadd = simcrit(i_mcc,iac)
            ELSE
              CALL calculate_variance_based_criteria(iac,nsub,stdr,stdc,stderr,locrs,loccc)
              CALL calculate_mean_criterion(nsub,loccc,critadd)  
            ENDIF
            crit1 = crit1 - critadd * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MKG')THEN                      !Median Kling-Gupta Efficiency
            IF(calcsimass)THEN
              locmkg = simcrit(i_mdkg,iac)
            ELSE
              CALL calculate_kling_gupta(iac,nsub,lockge,lockge1,lockge2,lockge3)             !KGE
              CALL calculate_median(nsub,lockge,missing_value,locmkg)    !Median of KGE
            ENDIF
            crit1 = crit1 - locmkg * weight
            thres1 = thres1 - thres * weight
          ENDIF
               
          IF(str=='MRS')THEN                      !Relative error of standard deviation, mean of all subbasins
            IF(calcsimass)THEN
              locmrs = simcrit(i_mstdre,iac)
            ELSE
              CALL calculate_variance_based_criteria(iac,nsub,stdr,stdc,stderr,locrs,loccc)
              CALL calculate_mean_criterion(nsub,locrs,locmrs)   
            ENDIF
            crit1 = crit1 + ABS(locmrs) * weight
            thres1 = thres1 + thres * weight
          ENDIF
               
          IF(str=='MNR')THEN                      !Median normalised RMSE 
            IF(calcsimass)THEN
              critadd = simcrit(i_mnrmse,iac) 
            ELSE
              CALL calculate_mean_and_errors(iac,nsub,locqr,locqc,locmae,locrmse,locnrmse)   !MAE and RMSE
              CALL calculate_median(nsub,locnrmse,missing_value,critadd)  !Median of Normalized RMSE 
            ENDIF
            crit1 = crit1 + critadd * weight
            thres1 = thres1 + thres * weight
          ENDIF         
        ENDIF
      ENDDO
    ENDIF
    
    !>Deallocate help variables for calculating Kendalls Tau or RA
    IF(ALLOCATED(ktcomp))  DEALLOCATE(ktcomp, ktrec, ktnum)
    IF(ALLOCATED(ktcomp2)) DEALLOCATE(ktcomp2,ktrec2,ktnum2)
    IF(ALLOCATED(ktcomp3)) DEALLOCATE(ktcomp3,ktrec3,ktnum3)
    IF(ALLOCATED(ktcomp4)) DEALLOCATE(ktcomp4,ktrec4,ktnum4)
    
  END SUBROUTINE calculate_criteria

  !>Finds the index for an output variable short name string
  !--------------------------------------------------------
  INTEGER FUNCTION find_variable_index(str,iout,flow)
    
    USE MODVAR, ONLY : max_outvar,outvarid
    USE CONVERT, ONLY : lower_case, scalar_lower_case
    
    !Argument declarations
    CHARACTER(LEN=4), INTENT(IN) :: str !<variable short name
    INTEGER, INTENT(OUT) :: iout        !<index of variable
    INTEGER, INTENT(OUT) :: flow        !<type of variable (0=state,1=flow,2=conc)
    !> \retval find_variable_index error status of function

    !Local variables
    INTEGER i
    CHARACTER(LEN=4) :: lowname(max_outvar)
    CHARACTER(LEN=4) :: str2
    
    !>\b Algorithm \n
    !>Transform all variable name to lower case.
    find_variable_index = 0
    lowname = outvarid(:)%shortname
    CALL lower_case(lowname,max_outvar)
    str2 = str
    CALL scalar_lower_case(str2)
    
    !>Find index of current output variable
    DO i = 1,max_outvar
      IF(str2==lowname(i))EXIT
    ENDDO
    !>If found: Set output dummy argument, else: return missing.
    IF(i<=max_outvar)THEN
      iout=i
      flow=outvarid(i)%vartype
    ELSE
      find_variable_index = 1   !variable missing
      iout  = 0
      flow = 0
      RETURN
    ENDIF
    
    RETURN
  END FUNCTION find_variable_index
  
  !>Calculate the scaling factor to get the right unit for yearly mean
  !>values based on current accumulation period
  !--------------------------------------------------------------------
  REAL FUNCTION compute_yearfact(vartype,per)
    
    USE WORLDVAR, ONLY : i_t,i_d,i_w,i_m,i_y, &
                         meandaysofyear
    USE MODVAR,   ONLY : i_sum,  &
                         seconds_per_timestep
    
    !Argument declarations
    INTEGER, INTENT(IN) :: vartype    !<type of accumulation for variable
    INTEGER, INTENT(IN) :: per        !<accumulation period code
    !> \retval compute_yearfact scaling factor
    
    !Local parameter
    INTEGER, PARAMETER  :: seconds_per_day = 86400
    
    !Local variables
    REAL yearfact, num_ts_per_day
    
    !>\b Algorithm \n
    !>If accumulation type is sum over period: Calculate scaling factor as number of time steps in period
    IF(vartype==i_sum)THEN    !yearly sum (yearfact=num_ts_per_period)
      IF(per == i_t)THEN
        yearfact = meandaysofyear * REAL(seconds_per_day/seconds_per_timestep)
      ELSEIF(per == i_d)THEN
        yearfact = meandaysofyear  
      ELSEIF(per == i_w)THEN
        yearfact = meandaysofyear/7.
      ELSEIF(per == i_m)THEN
        yearfact = 12.
      ELSEIF(per == i_y)THEN
        yearfact = 1.
      ENDIF
    !>If accumulation type is average over period: Calculate scaling factor as reciprocal of number of time steps in period
    ELSE                      !yearly mean (yearfact=1/num_ts_per_period)
      num_ts_per_day = REAL(seconds_per_day/seconds_per_timestep)
      IF(per == i_t)THEN  !For i_mean and i_wmean to give yearly mean locqr and locqc
        yearfact = 1.
      ELSEIF(per == i_d)THEN  !For i_mean and i_wmean to give yearly mean locqr and locqc
        yearfact = 1./num_ts_per_day
      ELSEIF(per == i_w)THEN
        yearfact = 1./(7.*num_ts_per_day)
      ELSEIF(per == i_m)THEN
        yearfact = 12./(meandaysofyear*num_ts_per_day)
      ELSEIF(per == i_y)THEN
        yearfact = 1./(meandaysofyear*num_ts_per_day)
      ENDIF
    ENDIF
    
    compute_yearfact = yearfact
    
  END FUNCTION compute_yearfact
  
  !>Calculate minimum number of data points necessary to have minimum
  !>number of years with data for SR2-criterium calculation
  !--------------------------------------------------------------------
  REAL FUNCTION compute_sr2limit(yearlim,per)

    USE WORLDVAR, ONLY : i_t,i_d,i_w,i_m,meandaysofyear
    USE MODVAR, ONLY : seconds_per_timestep
    
    !Argument declarations
    INTEGER, INTENT(IN) :: yearlim      !<Mimimum number of years needed
    INTEGER, INTENT(IN) :: per          !<Accumulation period code
    !> \retval compute:sr2limit number of data points needed

    !Local parameter
    INTEGER, PARAMETER :: seconds_per_day = 86400
    
    !Local variables
    REAL sr2limit
                  
    sr2limit = yearlim
    IF(per == i_t)THEN
      sr2limit = NINT(meandaysofyear) * yearlim * (seconds_per_day/seconds_per_timestep)
    ELSEIF(per == i_d)THEN
      sr2limit = NINT(meandaysofyear) * yearlim
    ELSEIF(per == i_w)THEN
      sr2limit = NINT(meandaysofyear * REAL(yearlim)/7.)
    ELSEIF(per == i_m)THEN
      sr2limit = 12 * yearlim
    ENDIF
    
    compute_sr2limit = sr2limit
    
  END FUNCTION compute_sr2limit

  !>Find what accumulated variables corresponds to current criteria variable
  !--------------------------------------------------------------------
  SUBROUTINE find_acrit_corresponding_to_crit(i,out)

    USE WORLDVAR, ONLY : acccalvar,  &
                         calvar,     &
                         nacrit

    !Argument declarations
    INTEGER, INTENT(IN)  :: i     !<Current criterion index
    INTEGER, INTENT(OUT) :: out   !<Accumulated variable index
    
    ! Local variables
    INTEGER j
          
    DO j = 1,nacrit
      IF(calvar(i)%comp==acccalvar(j)%comp .AND. calvar(i)%rec==acccalvar(j)%rec)EXIT 
    ENDDO
    out = j
    
  END SUBROUTINE find_acrit_corresponding_to_crit
    
  !>Find what criteria variable corresponds to current accumulated variable
  !--------------------------------------------------------------------
  SUBROUTINE find_crit_corresponding_to_acrit(j,out)
    
    USE WORLDVAR, ONLY : acccalvar,  &
         calvar,     &
         ncrit
    
    !Argument declarations
    INTEGER, INTENT(IN)  :: j     !<Accumulated variable index
    INTEGER, INTENT(OUT) :: out   !<Current criterion index

    !Local variables
    INTEGER i
    
    DO i = 1,ncrit
      IF(calvar(i)%comp==acccalvar(j)%comp .AND. calvar(i)%rec==acccalvar(j)%rec)EXIT 
    ENDDO
    out = i
    
  END SUBROUTINE find_crit_corresponding_to_acrit
    
  !>Save the variables needed to calculate Kendalls tau or RA after the simulation
  !>
  !>\b Consequences Worldvar module variables ktcomp,ktrec and ktnum etc
  !>may be allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE save_variables_for_later_critcalc(n,icrit,cs,rs,ts)

    USE WORLDVAR, ONLY : ktcomp, ktrec, ktnum,      & !OUT
                         ktcomp2,ktrec2,ktnum2,     & !OUT
                         ktcomp3,ktrec3,ktnum3,     & !OUT
                         ktcomp4,ktrec4,ktnum4        !OUT
    
    !Argument declaration
    INTEGER, INTENT(IN) :: n               !<Number of values (subbasins)
    INTEGER, INTENT(IN) :: icrit           !<Current accumulated variable index
    REAL, INTENT(IN)    :: cs(n)           !<Computed value
    REAL, INTENT(IN)    :: rs(n)           !<Recorded value
    INTEGER, INTENT(IN) :: ts(n)           !<Number of values used for cs,rs

    ! Local variables
    INTEGER i
    
    !>\b Algorithm \n
    !>Depending on index of accumulated variable:
    SELECT CASE(icrit)
    CASE(1)
      !>\li Allocate variables for holding values until end of simulation if necessary
      IF(.NOT.(ALLOCATED(ktcomp)))THEN 
        ALLOCATE(ktcomp(n,100))
        ALLOCATE(ktrec(n,100))
        ALLOCATE(ktnum(n))
        ktnum = 1
      ENDIF
      !>\li For every subbasin: If values exist set variable. If variable full increase its size.
      DO i = 1,n
        IF(ts(i)>0)THEN
          ktcomp(i,ktnum(i)) = cs(i)
          ktrec(i,ktnum(i))  = rs(i)
          ktnum(i) = ktnum(i) + 1
          !Increase the size of variables if necessary
          IF(ktnum(i).GT.SIZE(ktcomp,2))  &
            CALL extend_ktvariables(n,SIZE(ktcomp,2),icrit)
        ENDIF
      ENDDO
    CASE(2)
      !Allocate if necessary
      IF(.NOT.(ALLOCATED(ktcomp2))) THEN 
        ALLOCATE(ktcomp2(n,100))
        ALLOCATE(ktrec2(n,100))
        ALLOCATE(ktnum2(n))
        ktnum2 = 1
      ENDIF
      !Save if values exist
      DO i = 1,n
        IF(ts(i)>0)THEN
          ktcomp2(i,ktnum2(i)) = cs(i)
          ktrec2(i,ktnum2(i))  = rs(i)
          ktnum2(i) = ktnum2(i) + 1
          !Increase the size of variables if necessary
          IF(ktnum2(i).GT.SIZE(ktcomp2,2)) &
            CALL extend_ktvariables(n,SIZE(ktcomp2,2),icrit)
        ENDIF
      ENDDO
    CASE(3)
      !Allocate if necessary
      IF(.NOT.(ALLOCATED(ktcomp3))) THEN 
        ALLOCATE(ktcomp3(n,100))
        ALLOCATE(ktrec3(n,100))
        ALLOCATE(ktnum3(n))
        ktnum3 = 1
      ENDIF
      !Save if values exist
      DO i = 1,n
        IF(ts(i)>0)THEN
          ktcomp3(i,ktnum3(i)) = cs(i)
          ktrec3(i,ktnum3(i))  = rs(i)
          ktnum3(i) = ktnum3(i) + 1
          !Increase the size of variables if necessary
          IF(ktnum3(i).GT.SIZE(ktcomp3,2)) &
            CALL extend_ktvariables(n,SIZE(ktcomp3,2),icrit)
        ENDIF
      ENDDO
    CASE(4)
      !Allocate if necessary
      IF(.NOT.(ALLOCATED(ktcomp4))) THEN 
        ALLOCATE(ktcomp4(n,100))
        ALLOCATE(ktrec4(n,100))
        ALLOCATE(ktnum4(n))
        ktnum4 = 1
      ENDIF
      !Save if values exist
      DO i = 1,n
        IF(ts(i)>0)THEN
          ktcomp4(i,ktnum4(i)) = cs(i)
          ktrec4(i,ktnum4(i))  = rs(i)
          ktnum4(i) = ktnum4(i) + 1
          !Increase the size of variables if necessary
          IF(ktnum4(i).GT.SIZE(ktcomp4,2)) &
            CALL extend_ktvariables(n,SIZE(ktcomp4,2),icrit)
        ENDIF
      ENDDO
    END SELECT
    
  END SUBROUTINE save_variables_for_later_critcalc
  
  !>Subroutine that extends two allocatable arrays by 100
  !>
  !>\b Consequences Module worldvar variables are increased in size
  !--------------------------------------------------------------------
  SUBROUTINE extend_ktvariables(n,m,i)

    USE WORLDVAR, ONLY : ktcomp, ktrec,     &  !OUT
                         ktcomp2,ktrec2,    &  !OUT
                         ktcomp3,ktrec3,    &  !OUT
                         ktcomp4,ktrec4        !OUT

    !Argument declaration
    INTEGER, INTENT(IN) :: n              !<Number of values (subbasins)
    INTEGER, INTENT(IN) :: m              !<Second dimension, to be increased
    INTEGER, INTENT(IN) :: i              !<Current accumulated variable index
    
    !Variable declaration
    REAL, ALLOCATABLE :: loccomp(:,:),locrec(:,:)
    
    SELECT CASE(i)
    CASE(1)
      ALLOCATE(loccomp(n,m+100))
      loccomp(1:n,1:m) = ktcomp
      CALL move_alloc(loccomp,ktcomp)
      ALLOCATE(locrec(n,m+100))
      locrec(1:n,1:m)  = ktrec
      CALL move_alloc(locrec,ktrec)
    CASE(2)
      ALLOCATE(loccomp(n,m+100))
      loccomp(1:n,1:m) = ktcomp2
      CALL move_alloc(loccomp,ktcomp2)
      ALLOCATE(locrec(n,m+100))
      locrec(1:n,1:m)  = ktrec2
      CALL move_alloc(locrec,ktrec2)
    CASE(3)
      ALLOCATE(loccomp(n,m+100))
      loccomp(1:n,1:m) = ktcomp3
      CALL move_alloc(loccomp,ktcomp3)
      ALLOCATE(locrec(n,m+100))
      locrec(1:n,1:m)  = ktrec3
      CALL move_alloc(locrec,ktrec3)
    CASE(4)
      ALLOCATE(loccomp(n,m+100))
      loccomp(1:n,1:m) = ktcomp4
      CALL move_alloc(loccomp,ktcomp4)
      ALLOCATE(locrec(n,m+100))
      locrec(1:n,1:m)  = ktrec4
      CALL move_alloc(locrec,ktrec4)
    END SELECT
    
  END SUBROUTINE extend_ktvariables
  
  !>Accumulate variables to period used by criterion
  !--------------------------------------------------------------------
  SUBROUTINE accumulate_periodvalues(m,r,c,rs,cs,ts,n)
    
    !Argument declaration
    REAL, INTENT(IN)       :: m               !<Mask value, missing value
    INTEGER, INTENT(IN)    :: n               !<Size of variables, number of subbasins
    REAL, INTENT(IN)       :: r(n)            !<Recorded value of variable to be used for criteria
    REAL, INTENT(IN)       :: c(n)            !<Computed value of variable to be used for criteria
    REAL, INTENT(INOUT)    :: rs(n)           !<Accumulation of recorded value
    REAL, INTENT(INOUT)    :: cs(n)           !<Accumulation of computed value
    INTEGER, INTENT(INOUT) :: ts(n)           !<Accumulation of number of values
    
    WHERE(r(:)/=m .AND. c(:)/=m)
      rs(:) = rs(:) + r(:)
      cs(:) = cs(:) + c(:)
      ts(:) = ts(:) + 1
    ENDWHERE
    
  END SUBROUTINE accumulate_periodvalues

  !>Reform period accumulate variables to average value for some variables
  !--------------------------------------------------------------------
  SUBROUTINE reform_accumulated_periodvalues(n,acctype,rs,cs,ts)
    
    USE MODVAR, ONLY : i_mean,i_wmean
    
    !Argument declaration
    INTEGER, INTENT(IN)    :: n               !<Dimension of variables, number of subbasins
    INTEGER, INTENT(IN)    :: acctype         !<Type of accumulation for variable; sum, mean, weighted mean
    REAL, INTENT(INOUT)    :: rs(n)           !<Accumulation of recorded value
    REAL, INTENT(INOUT)    :: cs(n)           !<Accumulation of computed value
    INTEGER, INTENT(INOUT) :: ts(n)           !<Accumulation of number of values
    
    IF(acctype==i_mean.OR.acctype==i_wmean)THEN
      WHERE(ts(:)>0)
        rs(:) = rs(:) / ts(:)
        cs(:) = cs(:) / ts(:)
      ENDWHERE
    ENDIF
    
  END SUBROUTINE reform_accumulated_periodvalues

  !>Reset the period accumulation variables
  !--------------------------------------------------------------------
  SUBROUTINE reset_periodvalues(rs,cs,ts,n,m)

    !Argument declaration
    INTEGER, INTENT(IN)    :: n              !<Dimension of variables
    INTEGER, INTENT(IN)    :: m              !<Dimension of variables
    REAL, INTENT(INOUT)    :: rs(n,m)        !<Accumulation of recorded value
    REAL, INTENT(INOUT)    :: cs(n,m)        !<Accumulation of computed value
    INTEGER, INTENT(INOUT) :: ts(n,m)        !<Accumulation of number of values
    
    rs(:,:) = 0.
    cs(:,:) = 0.
    ts(:,:) = 0
    
  END SUBROUTINE reset_periodvalues

  !>Calculates Kendalls Tau, rank correlation, between two data sets
  !>with adjustment for ties
  !--------------------------------------------------------------------
  SUBROUTINE kendallstau(n,data1,data2,tau)
    
    !Argument declaration
    INTEGER, INTENT(IN) :: n          !<Dimension of data set
    REAL, INTENT(IN)    :: data1(n)   !<Data set 1
    REAL, INTENT(IN)    :: data2(n)   !<Data set 2
    REAL, INTENT(OUT)   :: tau        !<Calculated Kendalls tau
    
    !Variable declarations
    INTEGER j,k       !loop variables
    INTEGER is,n1,n2  !number of pairs
    REAL a1,a2,aa     !help variables
    
    IF(n==0)RETURN    !no data
    n1=0
    n2=0
    is=0
    DO j=1,n-1
      DO k=j+1,n
        a1=data1(j)-data1(k)
        a2=data2(j)-data2(k)
        aa=a1*a2
        IF(aa.NE.0.)THEN
          n1=n1+1
          n2=n2+1
          IF(aa.GT.0.)THEN
            is=is+1
          ELSE
            is=is-1
          ENDIF
        ELSE
          IF(a1.NE.0.)n1=n1+1
          IF(a2.NE.0.)n2=n2+1
        ENDIF
      ENDDO
    ENDDO
    tau = REAL(is)/SQRT(REAL(n1)*REAL(n2))
    
  END SUBROUTINE kendallstau

  !>Calculate three fitness criteria between two arrays of length n
  !--------------------------------------------------------------------
  SUBROUTINE calculate_fit(n,c,r,m,lim,x)
    
    !Argument declaration
    INTEGER, INTENT(IN) :: n          !<Number of elements
    REAL, INTENT(IN)    :: c(n)       !<Array with computed values
    REAL, INTENT(IN)    :: r(n)       !<Array with recorded values
    REAL, INTENT(IN)    :: m          !<Missing value
    INTEGER, INTENT(IN) :: lim        !<Minimum number of elements
    REAL, INTENT(OUT)   :: x(3)       !<Fitness criteria; 1=R2 (Nash-Sutcliff), 2=RMSE, 3=bias
    
    !Variable declaration
    INTEGER i,nobs
    DOUBLE PRECISION mse,varr,sumc,sumr,sumr2,sume2,s
    
    !> \b Algorithm \n
    !> Initialize 
    nobs = 0
    sumc = 0.D0
    sumr = 0.D0
    sumr2 = 0.D0
    sume2 = 0.D0
    x = m
    !> Calculate sum of values and errors (for not missing values)
    DO i = 1,n
      IF(r(i) == m .OR. c(i) == m) THEN
      ELSE
        nobs = nobs+1
        sumc = sumc+DBLE(c(i))
        sumr = sumr+DBLE(r(i))
        sumr2 = sumr2+DBLE(r(i))**2
        sume2 = sume2+(DBLE(c(i)-r(i)))**2   !Errors^2
      ENDIF
    ENDDO
    !> If enough values not missing; calculate fitness criteria
    IF(nobs>=lim)THEN
      s = sumr2-(sumr**2)/DBLE(nobs)       !Sum of squares around mean value
      varr = s/DBLE(nobs)                  !Variance ("/n")
      mse = sume2/DBLE(nobs)               !Mean square error
      x(1) = REAL(1.D0-mse/varr)                        !R2
      IF(mse>0.D0) x(2) = REAL(DSQRT(mse))              !RMSE
      IF(sumr /= 0.D0) x(3) = REAL((sumc-sumr)/sumr)    !Relative bias
    ENDIF
    
  END SUBROUTINE calculate_fit
  
  !>Calculates fitness criterion, "RA", for two columns of an array of length n.
  !>Criterion RA is similar to NSE but with exponent p instead of 2.
  !--------------------------------------------------------------------
  SUBROUTINE calculate_ra(n,c,r,m,p,lim,x)

    !Argument declaration
    INTEGER, INTENT(IN) :: n          !<Number of elements
    REAL, INTENT(IN)    :: c(n)       !<Array with computed values
    REAL, INTENT(IN)    :: r(n)       !<Array with recorded values
    REAL, INTENT(IN)    :: m          !<Missing value
    REAL, INTENT(IN)    :: p          !<Parameter value of criterion
    INTEGER, INTENT(IN) :: lim        !<Minimum number of elements
    REAL, INTENT(OUT)   :: x          !<Fitness criterion 1=RA
    
    !Variable declaration
    INTEGER i,nobs
    DOUBLE PRECISION sumr,sumea,meanr
    
    !> \b Algorithm \n
    !> Initialize 
    nobs = 0
    sumr = 0.D0
    sumea = 0.D0
    x = m
    !> Calculate sum of values and errors (for not missing values)
    DO i = 1,n
      IF(r(i) == m .OR. c(i) == m) THEN
      ELSE
        nobs = nobs+1
        sumr = sumr+DBLE(r(i))
      ENDIF
    ENDDO
    IF(nobs>0)THEN
      meanr = sumr / DBLE(nobs)
      sumr = 0.D0
      DO i = 1,n
        IF(r(i) == m .OR. c(i) == m) THEN
        ELSE
          IF(ABS(r(i)-meanr)>0.) sumr = sumr+(ABS(DBLE(r(i))-meanr))**p
          IF(ABS(c(i)-r(i))>0.) sumea = sumea+(ABS(DBLE(c(i)-r(i))))**p   !Errors^p
        ENDIF
      ENDDO
    ENDIF
    !> If enough values not missing; calculate fitness criteria
    IF(nobs>=lim)THEN
      x = REAL(1.D0 - sumea / sumr)                    !RA
    ENDIF
    
  END SUBROUTINE calculate_ra

  !>Calculate median of non-missing values in array 
  !--------------------------------------------------------------------
  SUBROUTINE calculate_median(n,x,miss,median)

    !Argument declaration
    INTEGER, INTENT(IN) :: n        !<Number of elements
    REAL, INTENT(IN)    :: x(n)     !<Array with values
    REAL, INTENT(IN)    :: miss     !<Missing value
    REAL, INTENT(OUT)   :: median   !<Median
    
    !Variable declaration
    INTEGER i,num,j,medi
    INTEGER il,ig,narray,ipivot
    REAL pivot
    REAL, ALLOCATABLE :: y(:),yini(:)
    REAL, ALLOCATABLE :: less(:)
    REAL, ALLOCATABLE :: great(:)

    num=COUNT(x/=miss)   !number of elements to be sorted
    
    !> \b Algorithm \n
    !> Calculate the median simle for small arrays, the return
    IF(num==0)THEN
      median = miss
      RETURN
    ELSEIF(num==1)THEN
      DO i=1,n
        IF(x(i)/=miss)THEN
          median = x(i)
          EXIT
        ENDIF
      ENDDO
      RETURN
    ELSEIF(num==2)THEN
      median = 0.
      DO i=1,n
        IF(x(i)/=miss)THEN
          median = median + x(i)
        ENDIF
      ENDDO
      median = median/2.
      RETURN
    ENDIF
    
    !> Calculate the place of the median element
    IF(MOD(num,2)==0)THEN
      medi = num/2
    ELSE
      medi = (num+1)/2
    ENDIF
    
    !Create an array to calculate median for
    ALLOCATE(y(num))
    j = 1
    DO i=1,n
       IF(x(i)/=miss)THEN
          y(j)=x(i)
          j = j + 1
       ENDIF
    ENDDO
    ALLOCATE(yini(num))
    yini=y
    
    !> Calculate the median
    ipivot = 1
    narray = num
    ALLOCATE(less(narray))
    ALLOCATE(great(narray))
    less=0.;great=0.
    
    DO
      pivot = y(ipivot)
      il = ipivot
      ig = narray
      DO i=ipivot+1,narray
        IF(y(i)<pivot)THEN
          less(il)=y(i)
          il=il+1
        ELSE
          great(ig)=y(i)
          ig=ig-1
        ENDIF
      ENDDO
      IF(narray-ig<il-1)THEN  !Add pivot to smaller array
        great(ig)=pivot
        il=il-1
      ELSE
        less(il)=pivot
        ig=ig+1
      ENDIF
      IF(ig==medi)THEN        !The median is the minimum of great
        median = MINVAL(great(ig:narray))
        EXIT
      ELSEIF(il<medi)THEN     !The median is larger than pivot
        y=0.
        y(1:narray-ig+1)=great(ig:narray)
        medi=medi-il
        ipivot=1
        narray=narray-ig+1
        great=0.;less=0.
      ELSEIF(il==medi)THEN    !The median is the maximun of less
        median = MAXVAL(less(1:il))
        EXIT
      ELSE                    !The median is less than pivot
        y=0.
        y(1:il)=less(1:il)
        ipivot=1
        narray=il
        great=0.;less=0.
      ENDIF
    ENDDO
    
    DEALLOCATE(less);DEALLOCATE(great)
    DEALLOCATE(y);DEALLOCATE(yini)

  END SUBROUTINE calculate_median

  !>Calculate (subbasin) mean value of criterion
  !--------------------------------------------------------------------
  SUBROUTINE calculate_mean_criterion(dim,crit,ave)
    
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: dim        !<Dimension of array, subbasins
    REAL, INTENT(IN)    :: crit(dim)  !<Array with criteria
    REAL, INTENT(OUT)   :: ave        !<Average of criteria
    
    ! Local variables
    INTEGER a,isb
    REAL    b
    
    ave = missing_value
    b = 0.
    a = 0
    DO isb = 1,dim
      IF(crit(isb)/=missing_value) THEN 
        b = b + crit(isb)
        a = a + 1
      ENDIF
    ENDDO
    IF(a>0) ave = b / REAL(a)                  !Aritmetic mean of criterion

  END SUBROUTINE calculate_mean_criterion

  !>Calculate subbasin mean absolute value of criterion
  !--------------------------------------------------------------------
  SUBROUTINE calculate_meanabs_criterion(dim,crit,ave)
      
    USE MODVAR, ONLY: missing_value
      
    !Argument declarations
    INTEGER, INTENT(IN) :: dim        !<Dimension of array, subbasins
    REAL, INTENT(IN)    :: crit(dim)  !<Array with criteria
    REAL, INTENT(OUT)   :: ave        !<Average of criteria
      
    !Local variables
    INTEGER a,isb
    REAL    b

    ave = missing_value
    b = 0.
    a = 0
    DO isb = 1,dim
      IF(crit(isb)/=missing_value) THEN 
        b = b + abs(crit(isb))
        a = a + 1
      ENDIF
    ENDDO
    IF(a>0) ave = b / REAL(a)                  !Aritmetic mean 

  END SUBROUTINE calculate_meanabs_criterion

  !>Calculate average value for each subbasin of simulated and
  !>recorded variable
  !--------------------------------------------------------------------
  SUBROUTINE calculate_subbasin_mean_var(i,iac,limityear,dim,rec,sim)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarper, &
                         calvar
    USE MODVAR, ONLY: missing_value
      
    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<Index of criteria
    INTEGER, INTENT(IN) :: iac        !<Index of accumularion criteria
    INTEGER, INTENT(IN) :: limityear  !<Minimum years with data (5)
    INTEGER, INTENT(IN) :: dim        !<Dimension of data arrays, subbasins
    REAL, INTENT(OUT)   :: rec(dim)   !<Recorded variable subbasin yearmean or yearsum
    REAL, INTENT(OUT)   :: sim(dim)   !<Computed variable subbasin yearmeanor yearsum
    
    !Local variables
    INTEGER sr2limitper
    INTEGER isb
    REAL    yearfact
    
    rec = missing_value
    sim = missing_value
    yearfact = compute_yearfact(calvar(i)%vartype,calvarper)
    sr2limitper = NINT(compute_sr2limit(limityear,calvarper))
    DO isb = 1,dim
      IF(critvec(1,isb,iac)>=sr2limitper)THEN
        rec(isb) = yearfact * REAL(critvec(2,isb,iac)/critvec(1,isb,iac))                       
        sim(isb) = yearfact * REAL((critvec(4,isb,iac) + critvec(2,isb,iac))/critvec(1,isb,iac))
      ENDIF
    ENDDO

  END SUBROUTINE calculate_subbasin_mean_var

  !> Calculate Nash-Sutcliffe effciency ("R2") for all subbasins
  !-------------------------------------------------------------
  SUBROUTINE calculate_nash_sutcliffe(iac,dim,r2)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac      !<Index of accumulation criteria
    INTEGER, INTENT(IN) :: dim      !<Dimension of array, number of subbasins
    REAL, INTENT(OUT)   :: r2(dim)  !<Nash-Sutcliffe efficiency
    
    !Local variables
    INTEGER isb
    DOUBLE PRECISION s2, t2
    
    r2 = missing_value
    DO isb=1,dim
      s2  = 0.D0
      t2  = 0.D0
      IF(critvec(1,isb,iac)>=calvarlim)THEN
        t2 = critvec(3,isb,iac) + critvec(7,isb,iac) - 2.D0*critvec(5,isb,iac)
        s2 = critvec(3,isb,iac) - critvec(2,isb,iac) * critvec(2,isb,iac) / critvec(1,isb,iac)
        IF(s2/=0) r2(isb) = REAL(1.D0 - t2 / s2)                                       !NSE
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_nash_sutcliffe

  !>Calculate Nash-Sutcliffe effciency for all time series as one
  !>time series, "regional R2"
  !--------------------------------------------------------------
  SUBROUTINE calculate_regional_nash_sutcliffe(iac,rr2)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value,  &
                      nsub
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac    !<Index of accumulation criteria
    REAL, INTENT(OUT)   :: rr2    !<Calculated criterion, regional R2
    
    !Local variables
    INTEGER m,isb
    DOUBLE PRECISION v2, v3, sd, td
    
    rr2  = missing_value
    m = 0
    td = 0.D0
    v2 = 0.D0 ; v3 = 0.D0
    DO isb = 1,nsub
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        m = m + NINT(critvec(1,isb,iac))
        td = td + critvec(3,isb,iac) + critvec(7,isb,iac) - 2.D0*critvec(5,isb,iac)
        v3 = v3 + critvec(3,isb,iac)
        v2 = v2 + critvec(2,isb,iac)
      ENDIF
    ENDDO
    IF(m>0)THEN
      sd = v3 - v2 * v2 / DBLE(m)
      IF(sd/=0)  rr2 = REAL(1.D0 - td / sd)
    ENDIF
    
  END SUBROUTINE calculate_regional_nash_sutcliffe

  !>Calculate Nash-Sutcliffe efficiency with subbasins as time
  !>series, "spatial R2"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_spatial_nash_sutcliffe(i,iac,limityear,limitbasin,sr2)

    USE MODVAR, ONLY: missing_value,  &
                      nsub
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i            !<Index of criteria
    INTEGER, INTENT(IN) :: iac          !<Index of accumulation criteria
    INTEGER, INTENT(IN) :: limityear    !<Minimum years with data (5)
    INTEGER, INTENT(IN) :: limitbasin   !<Minimum basins with data (5)
    REAL, INTENT(OUT)   :: sr2          !<Criteria, spatial R2
    
    !Local variables
    REAL   :: locqr(nsub)           !Recorded variable subbasin yearmean
    REAL   :: locqc(nsub)           !Computed variable subbasin yearmean
    REAL   :: spatialcrit(3)        !Spatial criteria: R2, RMSE, bias
    
    !> \b Algorithm \n
    !> Calculate average of variables for each subbasin
    CALL calculate_subbasin_mean_var(i,iac,limityear,nsub,locqr,locqc)
    !> Calculate Nash-Sutcliffe efficiency from the subbasin averages
    CALL calculate_fit(nsub,locqc,locqr,missing_value,limitbasin,spatialcrit)
    sr2 = spatialcrit(1)
    
  END SUBROUTINE calculate_spatial_nash_sutcliffe

  !>Calculate Kendalls tau for all subbasins
  !--------------------------------------------------------------------
  SUBROUTINE calculate_all_subbasins_tau(iac,dim,tau)

    USE WORLDVAR, ONLY : ktcomp, ktrec, ktnum,     &
                         ktcomp2,ktrec2,ktnum2,    &
                         ktcomp3,ktrec3,ktnum3,    &
                         ktcomp4,ktrec4,ktnum4
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac      !<Index of accumulated criteria
    INTEGER, INTENT(IN) :: dim      !<Dimension of array, number of subbasins
    REAL, INTENT(OUT)   :: tau(dim) !<Calculated criteria, Kendalls tau
    
    !Local variables
    INTEGER isb
    
    !> \b Algorithm \n
    tau = missing_value
    !> For each subbasin: Calculate kendalls tau
    DO isb = 1,dim
       IF(iac==1) CALL kendallstau(ktnum(isb)-1, ktcomp(isb,1:ktnum(isb)-1),  ktrec(isb,1:ktnum(isb)-1),  tau(isb))
       IF(iac==2) CALL kendallstau(ktnum2(isb)-1,ktcomp2(isb,1:ktnum2(isb)-1),ktrec2(isb,1:ktnum2(isb)-1),tau(isb))
       IF(iac==3) CALL kendallstau(ktnum3(isb)-1,ktcomp3(isb,1:ktnum3(isb)-1),ktrec3(isb,1:ktnum3(isb)-1),tau(isb))
       IF(iac==4) CALL kendallstau(ktnum4(isb)-1,ktcomp4(isb,1:ktnum4(isb)-1),ktrec4(isb,1:ktnum4(isb)-1),tau(isb))
    ENDDO
    
  END SUBROUTINE calculate_all_subbasins_tau

  !>Calculate Nash-Sutcliffe effciency with parameter a instead of 2
  !>for all time series, "RA"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_all_subbasins_ra(i,iac,dim,ra)

    USE WORLDVAR, ONLY : calvar,    &
                         calvarlim, &
                         ktcomp, ktrec, ktnum,     &
                         ktcomp2,ktrec2,ktnum2,    &
                         ktcomp3,ktrec3,ktnum3,    &
                         ktcomp4,ktrec4,ktnum4
    USE MODVAR, ONLY: missing_value
      
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<Index of criteria
    INTEGER, INTENT(IN) :: iac      !<Index of accumulation criteria
    INTEGER, INTENT(IN) :: dim      !<Dimension of array, number of subbasins
    REAL, INTENT(OUT)   :: ra(dim)  !<Calculated criteria, RA
      
    !Local variables
    INTEGER isb
    REAL b
    
    !> \b Algorithm \n
    ra = missing_value
    !> For each subbasin: Calculate criteria RA
    DO isb = 1,dim
       IF(iac==1) CALL calculate_ra(ktnum(isb)-1,ktcomp(isb,1:ktnum(isb)-1),ktrec(isb,1:ktnum(isb)-1),missing_value,calvar(i)%coeff,calvarlim,b)
       IF(iac==2) CALL calculate_ra(ktnum2(isb)-1,ktcomp2(isb,1:ktnum2(isb)-1),ktrec2(isb,1:ktnum2(isb)-1),missing_value,calvar(i)%coeff,calvarlim,b)
       IF(iac==3) CALL calculate_ra(ktnum3(isb)-1,ktcomp3(isb,1:ktnum3(isb)-1),ktrec3(isb,1:ktnum3(isb)-1),missing_value,calvar(i)%coeff,calvarlim,b)
       IF(iac==4) CALL calculate_ra(ktnum4(isb)-1,ktcomp4(isb,1:ktnum4(isb)-1),ktrec4(isb,1:ktnum4(isb)-1),missing_value,calvar(i)%coeff,calvarlim,b)
       ra(isb) = b
    ENDDO
    
  END SUBROUTINE calculate_all_subbasins_ra

  !>Calculate Nash-Sutcliffe effciency with parameter a for all time
  !>series as one time series, "regional RA"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_regional_ra(i,iac,rra)

    USE WORLDVAR, ONLY : calvar,    &
                         calvarlim, &
                         ktcomp, ktrec, ktnum,     &
                         ktcomp2,ktrec2,ktnum2,    &
                         ktcomp3,ktrec3,ktnum3,    &
                         ktcomp4,ktrec4,ktnum4
    USE MODVAR, ONLY: missing_value,  &
                      nsub
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i      !<Index of criteria
    INTEGER, INTENT(IN) :: iac    !<Index of accumulation criteria
    REAL, INTENT(OUT)   :: rra    !<Calculated criteria, regional RA

    !Local variables
    INTEGER isb
    INTEGER dim
    INTEGER kmin
    REAL, ALLOCATABLE :: allcom(:), allrec(:)
    
    rra = missing_value
    
    !> \b Algorithm \n
    !> Prepare arrays to hold all data in the time series
    dim = 0
    DO isb = 1,nsub
      IF(iac==1) dim = dim + ktnum(isb)-1
      IF(iac==2) dim = dim + ktnum2(isb)-1
      IF(iac==3) dim = dim + ktnum3(isb)-1
      IF(iac==4) dim = dim + ktnum4(isb)-1
    ENDDO
    IF(.NOT.ALLOCATED(allcom)) ALLOCATE(allcom(dim))
    IF(.NOT.ALLOCATED(allrec)) ALLOCATE(allrec(dim))
        
    !> Set the arrays with all data in the time series
    kmin = 1
    DO isb = 1,nsub
      IF(iac==1)THEN
        allcom(kmin:kmin+ktnum(isb)-2)   = ktcomp(isb,1:ktnum(isb)-1)
        allrec(kmin:kmin+ktnum(isb)-2)   = ktrec(isb,1:ktnum(isb)-1)
        kmin = kmin + ktnum(isb)-1
      ENDIF
      IF(iac==2)THEN
        allcom(kmin:kmin+ktnum2(isb)-2) = ktcomp2(isb,1:ktnum2(isb)-1)
        allrec(kmin:kmin+ktnum2(isb)-2) = ktrec2(isb,1:ktnum2(isb)-1)
        kmin = kmin + ktnum2(isb)-1
      ENDIF
      IF(iac==3)THEN
        allcom(kmin:kmin+ktnum3(isb)-2) = ktcomp3(isb,1:ktnum3(isb)-1)
        allrec(kmin:kmin+ktnum3(isb)-2) = ktrec3(isb,1:ktnum3(isb)-1)
        kmin = kmin + ktnum3(isb)-1
      ENDIF
      IF(iac==4)THEN
        allcom(kmin:kmin+ktnum4(isb)-2) = ktcomp4(isb,1:ktnum4(isb)-1)
        allrec(kmin:kmin+ktnum4(isb)-2) = ktrec4(isb,1:ktnum4(isb)-1)
        kmin = kmin + ktnum4(isb)-1
      ENDIF
    ENDDO
    
    !> Calculate RA for these large data arrays
    CALL calculate_ra(dim,allcom,allrec,missing_value,calvar(i)%coeff,calvarlim,rra)
    IF(ALLOCATED(allcom)) DEALLOCATE(allcom)
    IF(ALLOCATED(allrec)) DEALLOCATE(allrec)
    
  END SUBROUTINE calculate_regional_ra

  !>Calculate error relative to observed (absolute) mean value for all subbasins, "RE"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_relative_error(iac,dim,re)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac        !<Index of accumulation criteria
    INTEGER, INTENT(IN) :: dim        !<Dimension of array, number of subbasins
    REAL, INTENT(OUT)   :: re(dim)    !<Relative error
    
    !Local variables
    INTEGER isb
    
    re  = missing_value
    DO isb = 1,dim
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        IF(critvec(2,isb,iac)>0.D0)THEN
          re(isb) = REAL(critvec(4,isb,iac) / critvec(2,isb,iac))
        ELSEIF(critvec(2,isb,iac)<0.D0)THEN
          re(isb) = -1.* REAL(critvec(4,isb,iac) / critvec(2,isb,iac))
        ELSEIF(critvec(4,isb,iac)==0.D0)THEN
          re(isb) = 0.
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE calculate_relative_error

  !>Calculate relative volume error for all time series as one time
  !>series, "regional RE"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_regional_relative_error(iac,rve)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value,  &
                      nsub
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac    !<Index of accumulation criteria
    REAL, INTENT(OUT)   :: rve    !<Regional relative error
    
    !Local variables
    INTEGER m
    INTEGER isb
    DOUBLE PRECISION v2, v6
    
    !> \b Algorithm \n
    !> Accumulate data over all subbasins
    rve  = missing_value
    m = 0
    v2 = 0.D0 ; v6 = 0.D0
    DO isb = 1,nsub
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        m = m + NINT(critvec(1,isb,iac))
        v2 = v2 + critvec(2,isb,iac)
        v6 = v6 + critvec(4,isb,iac)
      ENDIF
    ENDDO
    !> Calculate regional relative volume error
    IF(m>0)THEN
! DG: why not do as with the individual relative errors if v2<0??
!      IF(v2>0) rve = REAL(v6 / v2)
      IF(v2>0.D0)THEN
        rve = REAL(v6 / v2)
      ELSEIF(v2<0.D0)THEN
        rve = -REAL(v6 / v2)
      ELSEIF(v6==0.D0)THEN
        rve = 0.
      ENDIF
    ENDIF

  END SUBROUTINE calculate_regional_relative_error
  
  !>Calculate mean absolute error and root mean square error for all
  !>time series, "MAE and RMSE". Calculate mean values for observations
  !>and simulation
  !--------------------------------------------------------------------
  SUBROUTINE calculate_mean_and_errors(iac,dim,rec,sim,rmae,rmse,nrmse)

    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac          !<index of accumulation criteria
    INTEGER, INTENT(IN) :: dim          !<dimension of criteria
    REAL, INTENT(OUT)   :: rec(dim)     !<mean recorded value
    REAL, INTENT(OUT)   :: sim(dim)     !<mean simulated value
    REAL, INTENT(OUT)   :: rmae(dim)    !<mean absolute error
    REAL, INTENT(OUT)   :: rmse(dim)    !<root mean square error
    REAL, INTENT(OUT)   :: nrmse(dim)   !<normalized root mean square error, rmse/(maxRec-MinRec)

    !Local variables
    INTEGER isb
    DOUBLE PRECISION :: mse
    
    rmae = missing_value
    rmse = missing_value
    nrmse = missing_value
    
    DO isb = 1,dim
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        rec(isb) = REAL(critvec(2,isb,iac) / critvec(1,isb,iac))                        !Rec per subbasin
        sim(isb) = REAL((critvec(4,isb,iac)+critvec(2,isb,iac)) / critvec(1,isb,iac))   !Sim per subbasin
        rmae(isb) = REAL(critvec(6,isb,iac) / critvec(1,isb,iac))                       !MAE per subbasin
        mse = (critvec(3,isb,iac)+critvec(7,isb,iac)-2.D0*critvec(5,isb,iac)) / critvec(1,isb,iac)          !MSE
        IF(mse>0.D0)THEN
          rmse(isb) = REAL(DSQRT(mse))                                            !RMSE per subbasin
          IF(critvec(9,isb,iac)>0.D0) nrmse(isb) = rmse(isb)/REAL(critvec(9,isb,iac))   !Normalized RMSE (RMSE/max(OBS)) per subbasin
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_mean_and_errors
  
  !>Calculate mean absolute error and mean square error for all time
  !>series as one time series, "regional MAE and MSE"
  !--------------------------------------------------------------------
  SUBROUTINE calculate_regional_errors(iac,rmae)
    
    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value,  &
                      nsub
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac    !<Index of accumulation criteria
    REAL, INTENT(OUT)   :: rmae   !<Regional mean absolute error
    
    !Local variables
    INTEGER m
    INTEGER isb
    DOUBLE PRECISION v6
    
    rmae = missing_value
    
    !> \b Algorithm \n
    !> Accumulate data over all subbasins
    m = 0
    v6 = 0.D0
    DO isb = 1,nsub
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        m = m + NINT(critvec(1,isb,iac))
        v6 = v6 + critvec(6,isb,iac)
      ENDIF
    ENDDO
    !>Calculate regional MAE of all observations
    IF(m>0)THEN
      rmae = REAL(v6 / DBLE(m))
    ENDIF
    
  END SUBROUTINE calculate_regional_errors

  !>Calculate standard deviation of observations and simulated values
  !>for all time series. Also calculate error, relative error and correlation of
  !>the standard deviation
  !--------------------------------------------------------------------
  SUBROUTINE calculate_variance_based_criteria(iac,dim,recstd,simstd,stderr,stdrelerr,corr)
    
    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac              !<index of accumulation criteria
    INTEGER, INTENT(IN) :: dim              !<dimension of criteria, number of subbasins
    REAL, INTENT(OUT)   :: recstd(dim)      !<recorded standard deviation
    REAL, INTENT(OUT)   :: simstd(dim)      !<simulated standard deviation
    REAL, INTENT(OUT)   :: stderr(dim)      !<error of standard deviation
    REAL, INTENT(OUT)   :: stdrelerr(dim)   !<relative error of standard deviation
    REAL, INTENT(OUT)   :: corr(dim)        !<correlation coefficient
    
    !Local variables
    INTEGER isb
    DOUBLE PRECISION :: n,mrec,msim,vrec,vsim,stdrec,stdsim,cov
    
    !> \b Algorithm \n
    !Default values, missing
    recstd    = missing_value
    simstd    = missing_value
    stderr    = missing_value
    stdrelerr = missing_value
    corr      = missing_value
    
    !>For every subbasin with enough data: calculate criteria
    DO isb = 1,dim
      IF(critvec(1,isb,iac)>=calvarlim) THEN
        n = critvec(1,isb,iac)
        mrec = critvec(2,isb,iac)/n
        msim = (critvec(2,isb,iac)  + critvec(4,isb,iac))/n
        vrec = critvec(3,isb,iac)/n - mrec**2
        vsim = critvec(7,isb,iac)/n - msim**2
        cov  = critvec(5,isb,iac)/n - msim * mrec
        !> Calculate standard deviation, error and correlation
        IF(vrec>0.D0)THEN
          stdrec = DSQRT(vrec)
          recstd(isb) = REAL(stdrec)
        ENDIF
        IF(vsim>0.D0)THEN
          stdsim = DSQRT(vsim)
          simstd(isb) = REAL(stdsim)
        ENDIF
        IF(vrec>0.D0 .AND. vsim>0.D0)THEN
          stderr(isb) = REAL(stdsim - stdrec)
          IF(stdrec>0.D0)THEN
            stdrelerr(isb) = REAL((stdsim - stdrec)/stdrec)
          ELSEIF(stderr(isb) == 0.)THEN
            stdrelerr(isb) = 0.
          ENDIF
          corr(isb) = REAL(cov/stdrec/stdsim)
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_variance_based_criteria
  
  !>\brief Calculate Kling-Gupta effciency for all time series, "KGE"
  !>
  !> \b Reference Gupta et al. 2009
  !> KGE = 1 - sqrt{ (cc-1)^2 + (a-1)^2 + (b-1)^2 }
  !> cc = linear correlation coefficient 
  !> a  = std(SIM)/std(OBS) 
  !> b  = mean(SIM)/mean(OBS) 
  !--------------------------------------------------------------------
  SUBROUTINE calculate_kling_gupta(iac,dim,kge,kge1,kge2,kge3)
    
    USE WORLDVAR, ONLY : critvec,   &
                         calvarlim
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iac      !<Index of accumulation criteria
    INTEGER, INTENT(IN) :: dim      !<Dimension of array, number of subbasins
    REAL, INTENT(OUT)   :: kge(dim) !<Kling-Gupta efficiency
    REAL, INTENT(OUT)   :: kge1(dim) !<Kling-Gupta efficiency part 1, cc
    REAL, INTENT(OUT)   :: kge2(dim) !<Kling-Gupta efficiency part 2, std-quotient
    REAL, INTENT(OUT)   :: kge3(dim) !<Kling-Gupta efficiency part 3, mean-quotient
    
    !Local variables
    INTEGER isb
    DOUBLE PRECISION cc, a, b
    DOUBLE PRECISION n,mrec,msim,vrec,vsim,stdrec,stdsim,cov
    
    !> \b Algorithm \n
    kge = missing_value
    kge1 = missing_value
    kge2 = missing_value
    kge3 = missing_value
    
    !>For every subbasin with enough data: calculate criteria
    DO isb=1,dim
      cc = 0.D0
      a  = 0.D0
      b  = 0.D0
      IF(critvec(1,isb,iac)>=calvarlim)THEN ! if n > minimum
        n = critvec(1,isb,iac)
        mrec = critvec(2,isb,iac)/n
        msim = (critvec(2,isb,iac) + critvec(4,isb,iac))/n
        vrec = critvec(3,isb,iac)/n - mrec**2
        vsim = critvec(7,isb,iac)/n - msim**2
        cov = critvec(5,isb,iac)/n - msim * mrec
        IF(vrec>0.D0) stdrec = DSQRT(vrec)
        IF(vsim>0.D0) stdsim = DSQRT(vsim)
        IF(vrec>0.D0 .AND. vsim>0.D0)THEN
          cc = cov/stdrec/stdsim
          a = stdsim/stdrec
        ENDIF
!        IF(stdsim>0.D0 .AND. stdrec>0.D0)  a  = stdsim/stdrec
        IF(msim>0.D0 .AND. mrec>0.D0) b = msim/mrec
        kge(isb) = REAL(1.D0 -DSQRT((cc-1.D0)**2+(a-1.D0)**2+(b-1.D0)**2))
        kge1(isb) = REAL(cc)
        kge2(isb) = REAL(a)
        kge3(isb) = REAL(b)
      ENDIF
    ENDDO
      
  END SUBROUTINE calculate_kling_gupta
    
END MODULE COMPOUT
