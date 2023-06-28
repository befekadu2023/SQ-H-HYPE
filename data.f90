!> \file data.f90
!> Contains module datamodule.

!> \brief Read and process model input from files, both static data
!> and forcing.
!>
!> Procedures for observations, for reading and preparing input data, 
!> for saving result files, and for loading and saving initial states
MODULE DATAMODULE
  !Copyright 2011-2016 SMHI
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
  !
  !--------------------------------------------------------------------
  USE COMPOUT, ONLY : find_variable_index  
  USE CONVERT, ONLY : lower_case, &
                      string_convert_to_DateType    
  USE LIBDATE
  USE READWRITE_ROUTINES, ONLY : check_obs_timeperiod, &
                                 check_pt_stn, &
                                 check_q_stn,  &
                                 check_xobs,   &
                                 convert_string_to_integer, &
                                 count_data_cols, &
                                 count_data_rows, &
                                 prepare_read_matrix, &
                                 read_basindata5, &
                                 read_basindata6, &
                                 read_column_headings, &
                                 read_geoid, &
                                 read_matrix,  &
                                 read_matrix_line, &
                                 read_next_codestr_on_line, &
                                 read_next_date_on_line, &
                                 write_dataline,  &
                                 read_parameterline
  USE TIMEROUTINES, ONLY : get_dayno_from_monthday, &
                           period_length, &
                           set_timestep_variables,  &
                           add_date_to_array
  !Subroutines also uses modvar, worldvar, modelmodule

  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  ! calculate_special_optpar_parameters
  ! load_lakedata
  ! load_damdata
  ! load_extradata
  ! load_dailydataXXXX
  ! start_lakedata_table
  ! finish_lakedata_table
  ! read_and_calc_basindata
  ! read_pointsourcedata
  ! load_irrigation_data
  ! read_update_data
  ! get_subid_index
  ! open_subbasinfiles
  ! close_subbasinfiles
  ! open_timefiles
  ! close_timefiles
  ! read_parameterline
  ! get_parametervalues
  ! write_ensemble_simulations_heading
  PUBLIC ::  load_observations,&
             load_ascii_pt_observations,&
             load_ascii_ss_observations,&
             load_ascii_qx_observations,&
             reset_observations,&
             prepare_for_update,&
             get_current_date_prec,&
             get_current_temp,&
             get_current_flow,&
             get_current_otherobs,&
             get_current_snowfallfrac,&
             get_current_shortwave,&
             get_current_windspeed,&
             get_current_humidity,&
             get_current_tmin,&
             get_current_tmax,&
             close_observations,&
             load_basindata,&
             set_model_configuration,&
             load_aquiferdata,&
             load_glacierdata,&
             load_branchdata,&
             load_pointsourcedata,&
             get_current_pointsources, &
             calculate_path,&
             reform_inputdata_for_submodel,&
             load_cropdata,&
             get_hyss_arguments,&
             load_coded_info,&
             load_submodel_info,&
             prepare_subbasin_output,&
             save_mapfiles,&
             write_subbasin_assessment,&
             write_simulation_assessment,&
             prepare_outputfiles, &
             close_outputfiles, &
             write_subbasinfiles, &
             write_timefiles, &
             load_parameters,&
             load_optpar,&
             save_respar,&
             initiate_output_routines,&
             initiate_outvar,&
             revise_outvar,&
             save_ensemble_simulations,&
             prepare_save_all_simulations,&
             write_simulation_results,&
             save_loadfiles, &
             checkindata_part1, &
             checkindata_part2, &
             checkindata_stop

CONTAINS
  
  !>Load forcing data files and other observations into the program
  !---------------------------------------------------------------------------------
  SUBROUTINE load_observations(dir,ns,bdate,edate,ndt,n_Result) 
    
    USE WORLDVAR, ONLY : readformat
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir     !<File directory
    INTEGER, INTENT(IN)  :: ns              !<Number of subbasins, basemodel
    TYPE(datetype), INTENT(IN)  :: bdate    !<Begin time for simulation
    TYPE(datetype), INTENT(IN)  :: edate    !<End time for simulation
    INTEGER, INTENT(IN)  :: ndt             !<Number of timesteps
    INTEGER, INTENT(OUT) :: n_Result        !<Status of subroutine
    
    !> \b Algorithm \n
    !>Read observations characteristics
    CALL load_pt_observation_characteristics(dir,ns,n_result)
    IF(n_Result/=0) RETURN

    !>Load observations from files
    IF(readformat==0)THEN   !ASCII-files
      CALL load_ascii_pt_observations(dir,ns,bdate,edate,ndt,n_Result) 
      IF(n_Result/=0) RETURN
      CALL load_ascii_ss_observations(dir,ns,bdate,edate,ndt,n_Result) 
      IF(n_Result/=0) RETURN
      CALL load_ascii_qx_observations(dir,ns,ndt,bdate,edate,n_Result) 
      IF(n_Result/=0) RETURN
    ENDIF
    

  END SUBROUTINE load_observations
  
  !>Reads PT data values in file
  !!
  !>\b Consequences Module worldvar variables dates, tobs, pobs, pobsindex,
  !!tobsindex, pobscol, tobscol, filepath_Tobs and filepath_Pobs will change.
  !!Module modvar variables preci and tempi will be allocated.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_ascii_pt_observations(dir,ns,bdate,edate,ndt,n_Result) 
    
    USE WORLDVAR, ONLY : dates,   &       !OUT (ev.)
                         tobs,    &       !OUT (ev.)
                         pobs,    &       !OUT (ev.)
                         pobsindex,    &  !OUT
                         tobsindex,    &  !OUT
                         pobsid,       &
                         tobsid,       &
                         pobscol,      &  !OUT
                         tobscol,      &  !OUT
                         readdaily,       &
                         readmatlab,      &
                         simsequence,     &
                         set_timeformat,  &
                         fileunit_pobs,   &
                         fileunit_tobs,   &
                         fileunit_temp,   &
                         filepath_Tobs,   &   !OUT
                         filepath_Pobs        !OUT
    USE MODVAR, ONLY : preci,     &   !allocated
                       tempi,     &   !allocated
                       missing_value

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    TYPE(DateType), INTENT(IN)  :: bdate !<Begin simulation date=bdate
    TYPE(DateType), INTENT(IN)  :: edate !<End simulation date=sdate
    INTEGER, INTENT(IN)  :: ndt          !<Number of timesteps in simulation
    INTEGER, INTENT(OUT) :: n_Result     !<Error number
    
    !Local variables
    LOGICAL notimefound
    INTEGER nt
    TYPE(DateType) fbdate,fedate       !Begin and end date of file
    REAL,ALLOCATABLE :: x(:,:)         !Data, only used when reading all obs to memory
    CHARACTER(LEN=3) :: seqnr          !Filename suffix

    !Initial allocation of model variables
    IF(.NOT.ALLOCATED(tempi)) ALLOCATE(tempi(ns))
    IF(.NOT.ALLOCATED(preci)) ALLOCATE(preci(ns))
    IF(.NOT. readdaily)THEN
      IF(.NOT.ALLOCATED(dates)) ALLOCATE(dates(ndt))
    ENDIF

    !Set filename of files to be used
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      filepath_Tobs = TRIM(dir)//'Tobs_'//seqnr//'.txt'
      filepath_Pobs = TRIM(dir)//'Pobs_'//seqnr//'.txt'
    ELSE
      filepath_Tobs = TRIM(dir)//'Tobs.txt'
      filepath_Pobs = TRIM(dir)//'Pobs.txt'
    ENDIF
    
    !Check if forcing data stations are correct and present, and get their index table.
    IF(.NOT.ALLOCATED(pobsindex)) ALLOCATE(pobsindex(ns))
    IF(.NOT.ALLOCATED(tobsindex)) ALLOCATE(tobsindex(ns))
    CALL check_pt_stn(fileunit_temp,filepath_Tobs,ns,tobsid,tobsindex,tobscol,n_Result)
    IF(n_Result.NE.0) RETURN
    CALL check_pt_stn(fileunit_temp,filepath_Pobs,ns,pobsid,pobsindex,pobscol,n_Result)
    IF(n_Result.NE.0) RETURN
    IF(ALLOCATED(pobsid)) DEALLOCATE(pobsid)
    IF(ALLOCATED(tobsid)) DEALLOCATE(tobsid)
    !Check time periods of observations
    CALL check_obs_timeperiod(fileunit_tobs,filepath_Tobs,1,bdate,   &
         edate,fbdate,fedate,notimefound,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR: Forcing data missing'
      RETURN
    ENDIF
    CALL check_obs_timeperiod(fileunit_pobs,filepath_Pobs,1,bdate,   &
         edate,fbdate,fedate,notimefound,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR: Forcing data missing'
      RETURN
    ENDIF
    CALL set_timeformat(notimefound)

    !Load temperature forcing data
    CALL prepare_read_matrix(fileunit_tobs,filepath_Tobs,1,bdate,n_Result)
    IF(n_Result.NE.0) RETURN
    WRITE(6,*) 'File ready: ', TRIM(filepath_Tobs)
    IF(readdaily)THEN
    ELSE
      ALLOCATE(x(ndt,tobscol))
      CALL READ_MATRIX(fileunit_tobs,ndt,tobscol,edate,nt,   &
            dates,x,missing_value,readmatlab) 
      CLOSE(fileunit_tobs)
      IF(.NOT.ALLOCATED(tobs)) ALLOCATE(tobs(ndt,tobscol))
      tobs = NINT(10.*x,2)   !tobs is INTEGER(2)
      DEALLOCATE(x)
    ENDIF

    !Load precipitation forcing data
    CALL prepare_read_matrix(fileunit_pobs,filepath_Pobs,1,bdate,n_Result)
    IF(n_Result.NE.0) RETURN
    WRITE(6,*) 'File ready: ', TRIM(filepath_Pobs)
    IF(readdaily)THEN
    ELSE
      ALLOCATE(x(ndt,pobscol))
      CALL READ_MATRIX(fileunit_pobs,ndt,pobscol,edate,nt,   &
            dates,x,missing_value,readmatlab) 
      CLOSE(fileunit_pobs)
      IF(.NOT.ALLOCATED(pobs)) ALLOCATE(pobs(ndt,pobscol))
      pobs = NINT(10.*x,2)   !pobs is INTEGER(2)
      DEALLOCATE(x)
    ENDIF

    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(x)) DEALLOCATE(x)
    
  END SUBROUTINE load_ascii_pt_observations

  !-----------------------------------------------------------------------------------------
  !>Reads optional forcing data values from file
  !!
  !>\b Consequences Module worldvar variables dates, sfobsindex, swobsindex, uobsindex
  !!rhobsindex, tminobsindex, tmaxobsindex, sfobscol, swobscol, uobscol, rhobscol, tminobscol, tmaxobscol, 
  !!sfobs, swobs, uobs, rhobs, tminobs, tmaxobs, filepath_SFobs, filepath_SWobs, filepath_Uobs, 
  !!filepath_RHobs, filepath_TMINobs and filepath_TMAXobs may change.
  !!Module modvar variables snowfraci, shortwavei, windi, humidi, tmini, tmaxi may be allocated.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_ascii_ss_observations(dir,ns,bdate,edate,ndt,n_Result) 
    
    USE WORLDVAR, ONLY : dates,   &       !OUT (ev.)
                         readdaily,       &
                         readmatlab,      &
                         simsequence,     &
                         set_timeformat,  &
                         fileunit_temp,   &
                         fileunit_sfobs,   &
                         fileunit_swobs,   &
                         fileunit_uobs,   &
                         fileunit_rhobs,   &
                         fileunit_tminobs,   &
                         fileunit_tmaxobs,   &
                         filepath_SFobs,   &  !OUT
                         filepath_SWobs,   &  !OUT
                         filepath_Uobs,   &  !OUT
                         filepath_RHobs,   &  !OUT
                         filepath_TMINobs,   &  !OUT
                         filepath_TMAXobs,   &  !OUT
                         sfobsindex,      &   !OUT 
                         swobsindex,      &   !OUT
                         uobsindex,      &   !OUT
                         rhobsindex,      &   !OUT
                         tminobsindex,      &   !OUT 
                         tmaxobsindex,      &   !OUT
                         sfobsid,       &        
                         swobsid,       &        
                         uobsid,       &        
                         rhobsid,       &        
                         tminobsid,       &        
                         tmaxobsid,       &        
                         sfobscol,      &  !OUT
                         swobscol,      &  !OUT
                         uobscol,      &  !OUT
                         rhobscol,      &  !OUT
                         tminobscol,      &  !OUT
                         tmaxobscol,      &  !OUT
                         sfobs,      &  !OUT
                         swobs,      &  !OUT
                         uobs,      &  !OUT
                         rhobs,      &  !OUT
                         tminobs,      &  !OUT
                         tmaxobs,      &  !OUT
                         readsfobs,  &   
                         readswobs, &
                         readwind, &
                         readhumid, &
                         readtminmaxobs
    USE MODVAR, ONLY : missing_value, &
                       snowfraci, &   !allocated optional forcing
                       shortwavei, &     !allocated optional forcing
                       windi, &
                       humidi, &
                       tmini, &
                       tmaxi
                       
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    TYPE(DateType), INTENT(IN)  :: bdate !<Begin simulation date=bdate
    TYPE(DateType), INTENT(IN)  :: edate !<End simulation date=sdate
    INTEGER, INTENT(IN)  :: ndt          !<Number of timesteps in simulation
    INTEGER, INTENT(OUT) :: n_Result     !<Error number
    
    !Local variables
    LOGICAL notimefound
    LOGICAL fexist
    INTEGER nt                         !Number of simulation time steps
    TYPE(DateType) fbdate,fedate       !Begin and end date of file
    CHARACTER(LEN=3) :: seqnr          !Filename suffix

    IF(.NOT.(readsfobs.OR.readswobs.OR.readwind.OR.readhumid.OR.readtminmaxobs)) RETURN   !No optional forcing used

    !Initial allocation of model variables
    IF(readsfobs) ALLOCATE(snowfraci(ns))
    IF(readswobs) ALLOCATE(shortwavei(ns))
    IF(readwind)  ALLOCATE(windi(ns))
    IF(readhumid) ALLOCATE(humidi(ns))
    IF(readtminmaxobs) ALLOCATE(tmini(ns),tmaxi(ns))
    IF(.NOT. readdaily)THEN
      IF(.NOT.ALLOCATED(dates)) ALLOCATE(dates(ndt))
    ENDIF
    
    !Set filename of files to be used
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      IF(readsfobs) filepath_SFobs = TRIM(dir)//'SFobs_'//seqnr//'.txt'
      IF(readswobs) filepath_SWobs = TRIM(dir)//'SWobs_'//seqnr//'.txt'
      IF(readwind)  filepath_Uobs = TRIM(dir)//'Uobs_'//seqnr//'.txt'
      IF(readhumid) filepath_RHobs = TRIM(dir)//'RHobs_'//seqnr//'.txt'
      IF(readtminmaxobs) filepath_TMINobs = TRIM(dir)//'TMINobs_'//seqnr//'.txt'
      IF(readtminmaxobs) filepath_TMAXobs = TRIM(dir)//'TMAXobs_'//seqnr//'.txt'
    ELSE
      IF(readsfobs) filepath_SFobs = TRIM(dir)//'SFobs.txt'
      IF(readswobs) filepath_SWobs = TRIM(dir)//'SWobs.txt'
      IF(readwind)  filepath_Uobs = TRIM(dir)//'Uobs.txt'
      IF(readhumid) filepath_RHobs = TRIM(dir)//'RHobs.txt'
      IF(readtminmaxobs) filepath_TMINobs = TRIM(dir)//'TMINobs.txt'
      IF(readtminmaxobs) filepath_TMAXobs = TRIM(dir)//'TMAXobs.txt'
    ENDIF
    
    !Snowfall fraction forcing data
    !------------------------------
    IF(readsfobs)THEN
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath_SFobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Snowfall forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_SFobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if stations are correct and present, and get their index table
      IF(.NOT.ALLOCATED(sfobsindex))ALLOCATE(sfobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_SFobs,ns,sfobsid,sfobsindex,sfobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(sfobsid)) DEALLOCATE(sfobsid)

      !Check time period of observations
      CALL check_obs_timeperiod(fileunit_sfobs,filepath_SFobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load snowfall fraction forcing data
      CALL prepare_read_matrix(fileunit_sfobs,filepath_SFobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_SFobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(sfobs(ndt,sfobscol))
        CALL READ_MATRIX(fileunit_sfobs,ndt,sfobscol,edate,nt,   &
             dates,sfobs,missing_value,readmatlab) 
        CLOSE(fileunit_sfobs)
      ENDIF
    ENDIF  

    !Shortwave radiation forcing data
    !------------------------------
    IF(readswobs)THEN
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath_SWobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Shortwave radiation forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_SWobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if optional variables are correct and present, and get their index table
      IF(.NOT.ALLOCATED(swobsindex))ALLOCATE(swobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_SWobs,ns,swobsid,swobsindex,swobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(swobsid)) DEALLOCATE(swobsid)

      !Check time periods of observations
      CALL check_obs_timeperiod(fileunit_swobs,filepath_SWobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load shortwave radiation forcing data
      CALL prepare_read_matrix(fileunit_swobs,filepath_SWobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_SWobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(swobs(ndt,swobscol))
        CALL READ_MATRIX(fileunit_swobs,ndt,swobscol,edate,nt,   &
             dates,swobs,missing_value,readmatlab) 
        CLOSE(fileunit_swobs)
      ENDIF
    ENDIF

    !Wind speed forcing data
    !------------------------------
    IF(readwind)THEN
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath_Uobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Wind speed forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_Uobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if optional variables are correct and present, and get their index table
      IF(.NOT.ALLOCATED(uobsindex))ALLOCATE(uobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_Uobs,ns,uobsid,uobsindex,uobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(uobsid)) DEALLOCATE(uobsid)

      !Check time periods of observations
      CALL check_obs_timeperiod(fileunit_uobs,filepath_Uobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load wind speed forcing data
      CALL prepare_read_matrix(fileunit_uobs,filepath_Uobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_Uobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(uobs(ndt,uobscol))
        CALL READ_MATRIX(fileunit_uobs,ndt,uobscol,edate,nt,   &
             dates,uobs,missing_value,readmatlab) 
        CLOSE(fileunit_uobs)
      ENDIF
    ENDIF

    !Relative humuidity forcing data
    !------------------------------
    IF(readhumid)THEN
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath_RHobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Relative humidity forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_RHobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if optional variables are correct and present, and get their index table
      IF(.NOT.ALLOCATED(rhobsindex))ALLOCATE(rhobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_RHobs,ns,rhobsid,rhobsindex,rhobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(rhobsid)) DEALLOCATE(rhobsid)

      !Check time periods of observations
      CALL check_obs_timeperiod(fileunit_rhobs,filepath_RHobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load relative humidity forcing data
      CALL prepare_read_matrix(fileunit_rhobs,filepath_RHobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_RHobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(rhobs(ndt,rhobscol))
        CALL READ_MATRIX(fileunit_rhobs,ndt,rhobscol,edate,nt,   &
             dates,rhobs,missing_value,readmatlab) 
        CLOSE(fileunit_rhobs)
      ENDIF
    ENDIF

    !Tmin forcing data
    !------------------------------
    IF(readtminmaxobs)THEN
      !Check that tmin-file exist
      INQUIRE(FILE=TRIM(filepath_TMINobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Tmin forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_TMINobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if optional variables are correct and present, and get their index table
      IF(.NOT.ALLOCATED(tminobsindex))ALLOCATE(tminobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_TMINobs,ns,tminobsid,tminobsindex,tminobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(tminobsid)) DEALLOCATE(tminobsid)

      !Check time periods of observations
      CALL check_obs_timeperiod(fileunit_tminobs,filepath_TMINobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load tmin forcing data
      CALL prepare_read_matrix(fileunit_tminobs,filepath_TMINobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_TMINobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(tminobs(ndt,tminobscol))
        CALL READ_MATRIX(fileunit_tminobs,ndt,tminobscol,edate,nt,   &
             dates,tminobs,missing_value,readmatlab) 
        CLOSE(fileunit_tminobs)
      ENDIF
    ENDIF
    
    !Tmax forcing data
    !------------------------------
    IF(readtminmaxobs)THEN
      !Check that tmin-file exist
      INQUIRE(FILE=TRIM(filepath_TMAXobs),EXIST=fexist)
      IF(.NOT.fexist)THEN
        WRITE(6,*) 'ERROR: Tmax forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(filepath_TMAXobs)
        n_Result = 1
        RETURN
      ENDIF
      !Check if optional variables are correct and present, and get their index table
      IF(.NOT.ALLOCATED(tmaxobsindex))ALLOCATE(tmaxobsindex(ns))
      CALL check_pt_stn(fileunit_temp,filepath_TMAXobs,ns,tmaxobsid,tmaxobsindex,tmaxobscol,n_Result)
      IF(n_Result.NE.0) RETURN
      IF(ALLOCATED(tmaxobsid)) DEALLOCATE(tmaxobsid)

      !Check time periods of observations
      CALL check_obs_timeperiod(fileunit_tmaxobs,filepath_TMAXobs,1,bdate,   &
             edate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data missing'
        RETURN
      ENDIF  
    
      !Load Tmin forcing data
      CALL prepare_read_matrix(fileunit_tmaxobs,filepath_TMAXobs,1,bdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(filepath_TMAXobs)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(tmaxobs(ndt,tmaxobscol))
        CALL READ_MATRIX(fileunit_tmaxobs,ndt,tmaxobscol,edate,nt,   &
             dates,tmaxobs,missing_value,readmatlab) 
        CLOSE(fileunit_tmaxobs)
      ENDIF
    ENDIF
  END SUBROUTINE load_ascii_ss_observations

  !>Reads Qobs and Xobs data values in file
  !!
  !>\b Consequences Module worldvar variables qobs, xobs, xcol, readqobs,
  !!qobsindex, numqobsstn, bqdate, eqdatem bxdate and exdate may change.
  !!Module modvar variables qobsi and xobsi may be allocated.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_ascii_qx_observations(dir,ns,ndt,bdate,edate,n_Result) 

    USE WORLDVAR, ONLY : qobs,    & !OUT
         readqobs,     &  !OUT
         qobsindex,    &  !OUT
         numqobsstn,   &  !OUT
         xobs,    &   !OUT
         xcol,    &   !OUT
         readdaily,       &
         readmatlab,      &
         fileunit_qobs,   &
         fileunit_xobs,   &
         fileunit_temp,   &
         filename_Xobs,   &
         filename_Qobs,   &
         bqdate,  & !OUT
         eqdate,  & !OUT
         bxdate,  & !OUT
         exdate     !OUT
    USE MODVAR, ONLY : basin,     &
         qobsi,     & !OUT (allocated)
         xobsi,     & !OUT (allocated)
         outvarid,  &
         max_outvar,   &
         missing_value,   &
         initiate_xobsinformation,  &
         save_xobsinformation

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    INTEGER, INTENT(IN)  :: ns                !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: ndt               !<Number of timesteps in simulation
    TYPE(DateType), INTENT(IN)   :: bdate     !<Begin simulation date
    TYPE(DateType), INTENT(IN)   :: edate     !<End simulation date
    INTEGER, INTENT(OUT) :: n_Result          !<Error number
    
    !Local variables
    INTEGER nt,nskip
    TYPE(DateType) fbdate,fedate          !Begin and end date of file
    TYPE(DateType) tempdate
    LOGICAL qex
    LOGICAL xex
    LOGICAL notimefound
    REAL,ALLOCATABLE :: x2(:,:)           !Data, only used when reading all obs to memory
    REAL,ALLOCATABLE :: x3(:,:)           !Data, only used when reading all obs to memory
    INTEGER,ALLOCATABLE :: subnr(:)       !Help variable
    TYPE(DateType),ALLOCATABLE :: localdates(:) !Help variable
    INTEGER,ALLOCATABLE :: xvar(:,:)      !Help variable

    !Initial allocation of model variables and local variables
    IF(.NOT. readdaily)THEN
      IF(.NOT.ALLOCATED(localdates)) ALLOCATE(localdates(ndt))
    ENDIF
    IF(.NOT.ALLOCATED(subnr)) ALLOCATE(subnr(ns))
    subnr = basin(:)%subid

    !Check if other observation stations are correct and present (only for Q), and get their index tables.
    INQUIRE(FILE=TRIM(dir)//filename_Qobs,EXIST=qex)
    IF(qex)THEN
      IF(.NOT.ALLOCATED(qobsindex)) ALLOCATE(qobsindex(ns))
      CALL check_q_stn(fileunit_temp,TRIM(dir)//filename_Qobs,ns,subnr,   &
            numqobsstn,qobsindex,n_Result)
      IF(n_Result==3)THEN
        qex = .FALSE.
      ELSEIF(n_Result.NE.0)THEN
        RETURN
      ENDIF
    ENDIF
    INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)
    IF(xex)THEN
      CALL count_data_cols(fileunit_temp,TRIM(dir)//filename_Xobs,1, xcol,n_Result)
      IF(n_Result.NE.0) RETURN
      xcol = xcol-1
      IF(xcol>0)THEN
        IF(.NOT.ALLOCATED(xvar)) ALLOCATE(xvar(xcol,2))
        CALL check_xobs(fileunit_temp,TRIM(dir)//filename_Xobs,xcol,    &
             max_outvar,outvarid(:)%shortname,xvar,n_Result) 
        IF(n_Result.NE.0) RETURN
      ENDIF
    ENDIF

    !Check time periods of observations
    IF(qex)THEN
      CALL check_obs_timeperiod(fileunit_qobs,TRIM(dir)//filename_Qobs,1,   &
            bdate,edate,fbdate,fedate,notimefound,n_Result)
      bqdate = bdate
      eqdate = edate
      IF(n_Result==2)THEN  !Shorter time period
        n_Result = 0
        IF(bdate.GT.fedate .OR. edate.LT.fbdate)THEN    !no overlap
          qex = .FALSE.
        ELSE                                      !some overlap
          tempdate = MaxDate(bdate,fbdate)
          bqdate = tempdate
          tempdate = MinDate(edate,fedate)
          eqdate = tempdate
        ENDIF
      ELSEIF(n_Result.NE.0)THEN
        RETURN
      ENDIF
    ENDIF
    IF(xex .AND. xcol>0)THEN
      CALL check_obs_timeperiod(fileunit_xobs,TRIM(dir)//filename_Xobs,3,   &
            bdate,edate,fbdate,fedate,notimefound,n_Result)
      bxdate = bdate
      exdate = edate
      IF(n_Result==2)THEN  !Shorter time period
        n_Result = 0
        IF(bdate.GT.fedate .OR. edate.LT.fbdate)THEN    !no overlap
          xex = .FALSE.
          xcol = 0
        ELSE                                      !some overlap
          bxdate = MaxDate(bdate,fbdate)
          exdate = MinDate(edate,fedate)
        ENDIF
      ELSEIF(n_Result.NE.0)THEN
        RETURN
      ENDIF
    ENDIF

    !Load discharge data
    IF(qex)THEN
      CALL prepare_read_matrix(fileunit_qobs,TRIM(dir)//filename_Qobs,1,bqdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(dir)//filename_Qobs
      IF(.NOT.ALLOCATED(qobsi)) ALLOCATE(qobsi(ns))
      IF(readdaily)THEN
        readqobs = .TRUE.
      ELSE
        IF(.NOT.ALLOCATED(x2)) ALLOCATE(x2(ndt,numqobsstn))
        CALL READ_MATRIX(fileunit_qobs,ndt,numqobsstn,eqdate,nt,   &
             localdates,x2(:,:),missing_value,readmatlab) 
        CLOSE(fileunit_qobs)
        nskip = period_length(bdate,bqdate)
        ALLOCATE (qobs(nskip+1:nskip+nt,numqobsstn))
        qobs(:,:)=x2(1:nt,1:numqobsstn)
      ENDIF
    ENDIF

    !Load other observation data
    IF(xex .AND. xcol>0)THEN
      CALL prepare_read_matrix(fileunit_xobs,TRIM(dir)//filename_Xobs,3,bxdate,n_Result)
      IF(n_Result.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(dir)//filename_Xobs
      IF(.NOT.ALLOCATED(xobsi)) ALLOCATE(xobsi(xcol))
      xobsi = missing_value
      IF(readdaily)THEN
      ELSE
        IF(.NOT.ALLOCATED(x3)) ALLOCATE(x3(ndt,xcol))
        CALL READ_MATRIX(fileunit_xobs,ndt,xcol,exdate,nt,localdates,    &
             x3(:,:),missing_value,readmatlab) 
        CLOSE(fileunit_xobs)
        nskip = period_length(bdate,bxdate)
        IF(.NOT.ALLOCATED(xobs)) ALLOCATE(xobs(nskip+1:nskip+nt,xcol))
        xobs(:,:) = x3(1:nt,:)
      ENDIF
    ENDIF
    CALL initiate_xobsinformation(ns)
    IF(xcol>0) CALL save_xobsinformation(xcol,xvar,ns)

    IF(ALLOCATED(xvar))   DEALLOCATE(xvar)
    IF(ALLOCATED(x2))     DEALLOCATE(x2)
    IF(ALLOCATED(x3))     DEALLOCATE(x3)
    IF(ALLOCATED(subnr))  DEALLOCATE(subnr)
    IF(ALLOCATED(localdates))  DEALLOCATE(localdates)

  END SUBROUTINE load_ascii_qx_observations

  !>Reads PT data characteristics from file
  !
  !> \b Consequences Module modvar variable forcing may be allocated and set
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_pt_observation_characteristics(dir,ns,n_Result) 
    
    USE WORLDVAR, ONLY : simsequence,     &
                         tobsid,          &
                         fileunit_temp,   &
                         maxcharpath,     &
                         i_str,i_real,i_intg
    USE MODVAR, ONLY : forcing    !OUT

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    INTEGER, INTENT(OUT) :: n_Result     !<Error status of subroutine

    !Local parameters
    INTEGER, PARAMETER :: maxcol = 3      !maximum number of column of ForcData.txt

    !Local variables
    INTEGER i ,j
    INTEGER n                             !number of data rows in file (number of tobsid)
    INTEGER mcols                         !number of columns in file
    INTEGER itobsid                       !column for tobsid
    INTEGER itobselev                     !column for tobselev
    CHARACTER(LEN=maxcharpath) filepath
    CHARACTER(LEN=10)  colstr(maxcol)
    CHARACTER(LEN=3)   seqnr
    LOGICAL fileexist
    INTEGER :: code(maxcol)        !Code for column variable
    INTEGER :: rindex(maxcol)      !Index for column real variables
    INTEGER :: iindex(maxcol)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    
    !Initial preparations
    n_Result = 0
    IF(simsequence>0)THEN                     !Set filename of files to be used
      WRITE(seqnr,'(I3.3)') simsequence
      filepath = TRIM(dir)//'ForcData_'//seqnr//'.txt'
    ELSE
      filepath = TRIM(dir)//'ForcData.txt'
    ENDIF
    
    INQUIRE(FILE=filepath,EXIST=fileexist)    !Check for file
    
    IF(fileexist)THEN
      !Initial allocation of model variables and local variables
      IF(.NOT.ALLOCATED(forcing)) ALLOCATE(forcing(ns))
      CALL count_data_rows(fileunit_temp,filepath,1,n,n_Result)
      IF(n_Result/=0) RETURN
      ALLOCATE(xi(n,maxcol))
      ALLOCATE(xr(n,maxcol))

      !Read file and put data in data structure
      OPEN(UNIT=fileunit_temp,FILE=filepath,STATUS='old',ACTION='read',ERR=900)
!      CALL read_and_calc_forcdata(fileunit_temp,ns,forcing,n_Result)
      CALL read_column_headings(fileunit_temp,maxcol,colstr,mcols,n_Result)
      IF(n_Result/=0) RETURN
      code = i_str
      itobsid = 0; itobselev = 0
      DO i = 1,mcols
        IF(colstr(i)(1:10)=='tobsid    ')THEN
          code(i) = i_intg
          itobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='tobselev  ')THEN
          code(i) = i_real
          itobselev = i
        ENDIF  
      ENDDO
      CALL read_basindata5(fileunit_temp,filepath,maxcol,n,mcols,code,rindex,iindex,xi,xr) 
      CLOSE(fileunit_temp)    
      DO i = 1,ns
        DO j = 1,n
          IF(tobsid(i)==xi(j,iindex(itobsid))) forcing(i)%tobselev = xr(j,rindex(itobselev))
        ENDDO
      ENDDO
      DEALLOCATE(xi)
      DEALLOCATE(xr)
    ENDIF
    RETURN
    
900 WRITE(6,*) 'Error open file: ', TRIM(filepath)
    n_Result = 1
    RETURN    
    
  END SUBROUTINE load_pt_observation_characteristics

  !>Prepare forcingdata and other observation for another simulation
  !>in case they are read each time step
  !-----------------------------------------------------------------
  SUBROUTINE reset_observations(dir,n_Result) 

    USE WORLDVAR, ONLY : xcol,    &
         readdaily,       &
         readqobs,        &
         readformat,      &
         bdate,           &
         bqdate,          &
         bxdate,          &
         fileunit_pobs,   &
         fileunit_tobs,   &
         fileunit_qobs,   &
         fileunit_xobs,   &
         filepath_Tobs,   &
         filepath_Pobs,   &
         filename_Xobs,   &
         filename_Qobs,   &
         readsfobs,       &
         readswobs,       &
         readwind,        &
         readhumid,        &
         readtminmaxobs,  &
         fileunit_sfobs,  &
         fileunit_swobs,  &
         fileunit_uobs,  &
         fileunit_rhobs,  &
         fileunit_tminobs,  &
         fileunit_tmaxobs,  &
         filepath_SFobs,    &
         filepath_SWobs,    &
         filepath_Uobs,    &
         filepath_RHobs,    &
         filepath_TMINobs,  &
         filepath_TMAXobs
    USE READWRITE_ROUTINES, ONLY : prepare_read_matrix

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir         !<File directory
    INTEGER, INTENT(OUT) :: n_Result            !<Error number
    
    !Local variables
    LOGICAL xex

    IF(readdaily)THEN
      !Close files
      CALL close_observations(dir)
      IF(readformat==0)THEN
        !Prepare to read forcing data and flow each day, no check
        CALL prepare_read_matrix(fileunit_tobs,filepath_Tobs,1,bdate,n_Result)
        IF(n_Result.NE.0) RETURN
        CALL prepare_read_matrix(fileunit_pobs,filepath_Pobs,1,bdate,n_Result)
        IF(n_Result.NE.0) RETURN
        IF(readqobs)THEN
          CALL prepare_read_matrix(fileunit_qobs,TRIM(dir)//filename_Qobs,1,bqdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)    !Other observations
        IF(xex .and. xcol>0)THEN
          CALL prepare_read_matrix(fileunit_xobs,TRIM(dir)//filename_Xobs,3,bxdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        IF(readsfobs)THEN
          CALL prepare_read_matrix(fileunit_sfobs,filepath_SFobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        IF(readswobs)THEN
          CALL prepare_read_matrix(fileunit_swobs,filepath_SWobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        IF(readwind)THEN
          CALL prepare_read_matrix(fileunit_uobs,filepath_Uobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        IF(readhumid)THEN
          CALL prepare_read_matrix(fileunit_rhobs,filepath_RHobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF
        IF(readtminmaxobs)THEN
          CALL prepare_read_matrix(fileunit_tminobs,filepath_TMINobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
          CALL prepare_read_matrix(fileunit_tmaxobs,filepath_TMAXobs,1,bdate,n_Result)
          IF(n_Result.NE.0) RETURN
        ENDIF

      ENDIF   !readformat
    ENDIF   !readdaily

  END SUBROUTINE reset_observations

  !>Read information for basins to update
  !!
  !>\b Consequences Module modvar variables doupdate, updatestations, 
  !!updatestationsqar, updatestationwar,updatestationsarfact, updatetpcorr, updatetncorr, 
  !!updatetploccorr, updatetnloccorr, and wendupdstations may change.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE prepare_for_update(dir,wobsvarname,locupall,locupnone,locupnone_qar,locupnone_war,wupall,wupnone,ns) 

    USE MODVAR, ONLY : doupdate,  & !OUT
         updatename,            &
         i_quseobs,             &
         i_qar,                 &
         i_war,                 &
         i_tpcorr,              &
         i_tncorr,              &
         i_tploccorr,           &
         i_tnloccorr,           &
         i_wendupd,             &
         updatestations,        & !OUT
         updatestationsqar,     & !OUT
         updatestationswar,     & !OUT
         updatestationsarfact,  & !OUT
         updatetpcorr,          & !OUT
         updatetncorr,          & !OUT
         updatetploccorr,       & !OUT
         updatetnloccorr,       & !OUT
         set_update_nutrientcorr, &
         basin,                 &
         xobsindex,             &
         wendupdstations,       & !OUT
         wobsvar,               &
         i_sp,i_pp,             &
         i_in,i_on,             &
         outvarid,              &
         max_outvar,            &
         nsub_basemodel
    USE WORLDVAR, ONLY : fileunit_temp,         &
         filename_upd,          &
         qobsindex

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER(LEN=4), INTENT(IN) :: wobsvarname !<Short-name of output variable that waterstage is updated with
    LOGICAL, INTENT(IN) :: locupall             !<Flag for quseobs update on all stations
    LOGICAL, INTENT(IN) :: locupnone            !<Flag for no quseobs update on any stations
    LOGICAL, INTENT(IN) :: locupnone_qar        !<Flag for no q_AR update on any stations !JD2012_AR	  
    LOGICAL, INTENT(IN) :: locupnone_war        !<Flag for no Q update with wAR on any stations
    LOGICAL, INTENT(IN) :: wupall               !<Flag for wendupd update on all stations
    LOGICAL, INTENT(IN) :: wupnone              !<Flag for no wendupd update on any stations
    INTEGER, INTENT(IN) :: ns                   !<Number of subbasins
    
    !Local variables
    INTEGER i, j            !loop-index
    INTEGER qusenum
    INTEGER qarnum
    INTEGER warnum
    INTEGER wendnum 
    INTEGER filenum 
    INTEGER qusetemp(nsub_basemodel)
    INTEGER qartemp(nsub_basemodel)
    INTEGER wartemp(nsub_basemodel)
    REAL    qarfac(nsub_basemodel)
    REAL    warfac(nsub_basemodel)
    INTEGER wendtemp(nsub_basemodel)
    REAL    tptemp(nsub_basemodel,2)
    REAL    tntemp(nsub_basemodel,2)
    REAL    tploctemp(nsub_basemodel,2)
    REAL    tnloctemp(nsub_basemodel,2)
    LOGICAL fileexist

    !Initiation
    qusenum = 0
    wendnum = 0
    filenum = 0
    qarnum  = 0
    warnum  = 0
    wobsvar = 0

    !Read update-file if available
    INQUIRE(FILE=TRIM(dir)//filename_upd,EXIST=fileexist) 
    IF(fileexist) CALL read_update_data(fileunit_temp,TRIM(dir)//filename_upd,    & 
         nsub_basemodel,qusenum,qusetemp,qarnum,warnum,qartemp,wartemp,qarfac,warfac,wendnum,  &
         wendtemp,filenum,tptemp,tntemp,tploctemp,tnloctemp)    

    !Dicharge data to use for update  (quseobs)
    IF(doupdate(i_quseobs))THEN
      IF(locupnone)THEN
        doupdate(i_quseobs) = .FALSE. 
      ELSEIF(locupall .OR. qusenum>0) THEN
        IF(.NOT.ALLOCATED(updatestations)) ALLOCATE(updatestations(ns))
        updatestations = .FALSE.             
        IF(locupall) THEN
          WHERE(qobsindex(1:ns)>0)
            updatestations(:) = .TRUE. 
          ENDWHERE
        ELSEIF(qusenum>0) THEN
          DO i = 1,ns
            DO j = 1, qusenum
              IF(basin(i)%subid == qusetemp(j)) THEN
                updatestations(i)=.TRUE.
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    !Dicharge data to use for AR-update  (qar)
    IF(doupdate(i_qar))THEN
      IF(locupnone_qar)THEN
        doupdate(i_qar) = .FALSE. 
      ELSEIF(qarnum>0) THEN
        IF(.NOT.ALLOCATED(updatestationsqar)) ALLOCATE(updatestationsqar(ns))
        IF(.NOT.ALLOCATED(updatestationsarfact)) ALLOCATE(updatestationsarfact(ns))
        updatestationsqar = .FALSE.             
        updatestationsarfact = 0
        DO i = 1,ns
          DO j = 1, qarnum
            IF(basin(i)%subid == qartemp(j)) THEN
              updatestationsqar(i)=.TRUE.
              updatestationsarfact(i)=qarfac(j)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    !Waterstage data to use for AR-update  (war)
    IF(doupdate(i_war))THEN
      IF(locupnone_war)THEN
        doupdate(i_war) = .FALSE. 
      ELSEIF(warnum>0) THEN
        IF(.NOT.ALLOCATED(updatestationswar)) ALLOCATE(updatestationswar(ns))
        IF(.NOT.ALLOCATED(updatestationsarfact))THEN
          ALLOCATE(updatestationsarfact(ns))
          updatestationsarfact = 0
        ENDIF
        updatestationswar = .FALSE.             
        DO i = 1,ns
          DO j = 1, warnum
            IF(basin(i)%subid == wartemp(j)) THEN
              updatestationswar(i)=.TRUE.
              updatestationsarfact(i)=warfac(j)
              EXIT
            ENDIF
          ENDDO
        ENDDO

        !Check for missing input data
        DO wobsvar=1,max_outvar
          IF(outvarid(wobsvar)%shortname==wobsvarname)EXIT
        ENDDO
        IF(wobsvar>max_outvar)THEN
          WRITE(6,*) 'WARNING: observation variable for wAR-updating not found in info.txt file, no wAR updating'
          IF(ALLOCATED(updatestationswar)) DEALLOCATE(updatestationswar)
          doupdate(i_war) = .FALSE. 
        ENDIF
      ENDIF

    ENDIF
    
    !Correction of phosphorus concentration (tpcorr,tploccorr)
    IF(doupdate(i_tpcorr).OR.doupdate(i_tploccorr))THEN
      IF(i_sp==0.OR.i_pp==0)THEN
        WRITE(6,*) 'ERROR: Updating of phsosphorus not possible unless P is simulated'
        STOP 1
      ENDIF
      IF(filenum>0)THEN
        IF(doupdate(i_tpcorr))    CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tpcorr,tptemp,updatetpcorr)
        IF(doupdate(i_tploccorr)) CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tploccorr,tploctemp,updatetploccorr)
      ELSEIF(.NOT.fileexist)THEN
        WRITE(6,*) 'WARNING: missing update.txt file, necessary for update of phosphourus'
        doupdate(i_tpcorr) = .FALSE.
        doupdate(i_tploccorr) = .FALSE.
      ELSE
        WRITE(6,*) 'WARNING: missing data in update.txt file, no phosphourus updating'
        doupdate(i_tpcorr) = .FALSE.
        doupdate(i_tploccorr) = .FALSE.
      ENDIF
    ENDIF

    !Correction of nitrogen concentration out of subbasin (tncorr,tnloccorr)
    IF(doupdate(i_tncorr).OR.doupdate(i_tnloccorr))THEN
      IF(i_in==0.OR.i_on==0)THEN
        WRITE(6,*) 'ERROR: Updating of nitrogen not possible unless N is simulated'
        STOP 1
      ENDIF
      IF(filenum>0)THEN
        IF(doupdate(i_tncorr))    CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tncorr,tntemp,updatetncorr)
        IF(doupdate(i_tnloccorr)) CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tnloccorr,tnloctemp,updatetnloccorr)
      ELSEIF(.NOT.fileexist)THEN
        WRITE(6,*) 'WARNING: missing update.txt file, necessary for update nitrogen'
        doupdate(i_tncorr) = .FALSE.
        doupdate(i_tnloccorr) = .FALSE.
      ELSE
        WRITE(6,*) 'WARNING: missing data in update.txt file, no nitrogen updating'
        doupdate(i_tncorr) = .FALSE.
        doupdate(i_tnloccorr) = .FALSE.
      ENDIF
    ENDIF

    !Update lake waterstage at end of time step (wendupd)
    IF(doupdate(i_wendupd))THEN
      IF(wupnone)THEN      
        doupdate(i_wendupd) = .FALSE. 
      ELSEIF(wupall .OR. wendnum>0) THEN
        IF(.NOT.ALLOCATED(wendupdstations)) ALLOCATE(wendupdstations(ns))
        wendupdstations = .FALSE.             
        IF(wupall) THEN
          wendupdstations(:) = .TRUE. 
        ELSEIF(wendnum>0) THEN
          DO i = 1,ns
            DO j = 1, wendnum
              IF(basin(i)%subid == wendtemp(j)) THEN
                wendupdstations(i)=.TRUE.
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        !Check for missmatches
        DO wobsvar=1,max_outvar
          IF(outvarid(wobsvar)%shortname==wobsvarname)EXIT
        ENDDO
        IF(wobsvar>max_outvar)THEN
          WRITE(6,*) 'WARNING: missing variable in info.txt file, no wendupd updating'
          IF(ALLOCATED(wendupdstations)) DEALLOCATE(wendupdstations)
          doupdate(i_wendupd) = .FALSE. 
        ELSE
          DO i = 1,ns
            IF(wendupdstations(i).AND. xobsindex(wobsvar,i)==0) wendupdstations(i)=.FALSE.
          ENDDO
        ENDIF
      ELSE
        WRITE(6,*) 'WARNING: missing stations in update.txt, no wendupd updating'
        doupdate(i_wendupd) = .FALSE. 
      ENDIF
    ENDIF

    !Check for incompatible update methods (war not to be used with qar or wendupd)
    IF(doupdate(i_war).AND.doupdate(i_qar))THEN
      DO i = 1,ns
        IF(updatestationswar(i).AND.updatestationsqar(i))THEN
          updatestationswar(i)=.FALSE.
          WRITE(6,*)'WARNING: Update with q-AR and w-AR at the same subbasin is not allowed. w-AR turned off for subbasin',basin(i)%subid
          warnum = warnum - 1
        ENDIF
      ENDDO
      IF(warnum == 0) doupdate(i_war)=.FALSE.
    ENDIF
    IF(doupdate(i_war).AND.doupdate(i_wendupd))THEN
      DO i = 1,ns
        IF(updatestationswar(i).AND.wendupdstations(i))THEN
          updatestationswar(i)=.FALSE.
          WRITE(6,*)'WARNING: Update with wendupd and w-AR at the same subbasin is not allowed. w-AR turned off for subbasin',basin(i)%subid
          warnum = warnum - 1
        ENDIF
      ENDDO
      IF(warnum == 0) doupdate(i_war)=.FALSE.
    ENDIF

    IF(doupdate(i_quseobs).OR.doupdate(i_qar).OR.doupdate(i_war).OR.doupdate(i_wendupd)  &
       .OR.doupdate(i_tpcorr).OR.doupdate(i_tncorr).OR.doupdate(i_tploccorr).OR.doupdate(i_tnloccorr))THEN
      WRITE(6,*)
      WRITE(6,*)'-----Information about the simulation (cont.)----'
      IF(doupdate(i_quseobs))   WRITE(6,*)'Simulation with updated discharge, ',TRIM(updatename(i_quseobs))
      IF(doupdate(i_qar))       WRITE(6,*)'Simulation with qAR-updated discharge, ',TRIM(updatename(i_qar))
      IF(doupdate(i_war))       WRITE(6,*)'Simulation with wAR-updated discharge, ',TRIM(updatename(i_war))
      IF(doupdate(i_wendupd))   WRITE(6,*)'Simulation with updated lake water stage, ',TRIM(updatename(i_wendupd))
      IF(doupdate(i_tpcorr))    WRITE(6,*)'Simulation with updated TP-concentrations, ',TRIM(updatename(i_tpcorr))
      IF(doupdate(i_tncorr))    WRITE(6,*)'Simulation with updated TN-concentrations, ',TRIM(updatename(i_tncorr))
      IF(doupdate(i_tploccorr)) WRITE(6,*)'Simulation with updated TP-concentrations, ',TRIM(updatename(i_tploccorr))
      IF(doupdate(i_tnloccorr)) WRITE(6,*)'Simulation with updated TN-concentrations, ',TRIM(updatename(i_tnloccorr))
      WRITE(6,*)'-------------------------------------------------'
    ENDIF
  END SUBROUTINE prepare_for_update

  !>Get current date from file or memory
  !--------------------------------------------------------------------
  SUBROUTINE get_current_date_prec(i,n,current_date,current_prec)

    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         fileunit_pobs,   &
                         pobsindex,       &
                         readmatlab,      &
                         pobscol,         &
                         get_current_date_memory, &
                         get_current_pobs
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    INTEGER, INTENT(IN)  :: i         !<current timestep index
    INTEGER, INTENT(IN)  :: n         !<number of subbasins
    TYPE(DateType), INTENT(OUT) :: current_date !<current date
    REAL, INTENT(OUT)    :: current_prec(n)     !<precipitation for current timestep [mm]
    
    ! Local variables
    TYPE(DateType) d 
    REAL    yp(pobscol)

    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_pobs,pobscol,d,yp,missing_value,readmatlab) 
        current_date = d
        current_prec = yp(pobsindex)
      ELSE
        current_date = get_current_date_memory(i)
        yp = get_current_pobs(i,pobscol)
        current_prec = yp(pobsindex)
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE get_current_date_prec

  !>Get current air temperature for selected subbasins from file or
  !>memory 
  !----------------------------------------------------------------
  FUNCTION get_current_temp(i,n) RESULT(current_temp)

    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         tobsindex,       &
                         fileunit_tobs,   &
                         readmatlab,      &
                         tobscol,         &
                         get_current_tobs
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    INTEGER, INTENT(IN) :: i  !<current timestep
    INTEGER, INTENT(IN) :: n  !<number of subbasins
    REAL current_temp(n)      !<temperature for current timestep
    
    !Local variables
    TYPE(DateType) d    !date
    REAL    yt(tobscol)

    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL read_matrix_line(fileunit_tobs,tobscol,d,yt,missing_value,readmatlab) 
        current_temp = yt(tobsindex)
      ELSE
        yt = get_current_tobs(i,tobscol)
        current_temp = yt(tobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_temp

  !>Get current discharge from file or memory
  !---------------------------------------------------------
  FUNCTION get_current_flow(i,cd,n,no) RESULT(current_flow)

    USE WORLDVAR, ONLY : readdaily,        &
                         readformat,       &
                         fileunit_qobs,    &
                         readmatlab,       &
                         qobsindex,        &
                         get_current_qobs, &
                         bqdate,           &
                         eqdate
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current time step
    TYPE(DateType), INTENT(IN) :: cd  !<current date
    INTEGER, INTENT(IN) :: n          !<number of subbasins
    INTEGER, INTENT(IN) :: no         !<number of observation stations
    REAL current_flow(n)              !<current flow
    
    !Local variables
    TYPE(DateType) d
    REAL    y(n),x(no)

    IF(cd.LT.bqdate .OR. cd.GT.eqdate)THEN
      current_flow = missing_value
      RETURN
    ENDIF

    IF(readformat==0.OR.readformat==4)THEN   !ASCII
      IF(readdaily)THEN
        y = missing_value
        CALL READ_MATRIX_LINE(fileunit_qobs,no,d,x,missing_value,readmatlab) 
        WHERE(qobsindex(1:n)>0)
          y=x(qobsindex(1:n))
        ENDWHERE
        current_flow = y
        IF(d.NE.cd) WRITE(*,*) 'ERROR read Qobs, date not equal to Pobs date'   !borde ej kunna inträffa
      ELSE
        current_flow = get_current_qobs(i,n)
      ENDIF
    ENDIF

  END FUNCTION get_current_flow

  !>Get current value of other observations from file or memory
  !------------------------------------------------------------
  FUNCTION get_current_otherobs(i,cd) RESULT(current_obs)

    USE WORLDVAR, ONLY : readdaily,        &
                         readformat,       &
                         fileunit_xobs,    &
                         get_current_xobs, &
                         readmatlab,       &
                         xcol,             &
                         bxdate,           &
                         exdate
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current time step
    TYPE(DateType), INTENT(IN) :: cd  !<current date 
    REAL current_obs(xcol)            !<current observations
    
    ! Local variables
    TYPE(DateType) d   !date
    REAL y(xcol)

    IF(xcol==0) RETURN        !no observations

    IF(readformat==0.OR.readformat==4)THEN   !ASCII
      current_obs = missing_value
      IF(cd.LT.bxdate .OR. cd.GT.exdate) RETURN
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_xobs,xcol,d,y,missing_value,readmatlab) 
        current_obs = y
        IF(d.NE.cd) WRITE(*,*) 'ERROR read Xobs, date not equal to Pobs date'   !borde ej kunna inträffa
      ELSE
        current_obs = get_current_xobs(i)
      ENDIF
    ENDIF

  END FUNCTION get_current_otherobs

  !>Get current fraction of snowfall in precipitation for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_snowfallfrac(i,n) RESULT(current_snowfallfrac)
  
    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         sfobsindex,       &
                         fileunit_sfobs,   &
                         readmatlab,      &
                         sfobscol,         &
                         get_current_sfobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_snowfallfrac(n)  !<current snowfall fraction
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL yt(sfobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_sfobs,sfobscol,d,yt,missing_value,readmatlab) 
        current_snowfallfrac = yt(sfobsindex)
      ELSE
        yt = get_current_sfobs(i,sfobscol)
        current_snowfallfrac = yt(sfobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_snowfallfrac
  
  !>Get current downward shortwave radiation for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_shortwave(i,n) RESULT(current_shortwave)
  
    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         swobsindex,      &
                         fileunit_swobs,  &
                         readmatlab,      &
                         swobscol,        &
                         get_current_swobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_shortwave(n)     !<current shortwave radiation
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL    yt(swobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_swobs,swobscol,d,yt,missing_value,readmatlab) 
        current_shortwave = yt(swobsindex)
      ELSE
        yt = get_current_swobs(i,swobscol)
        current_shortwave = yt(swobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_shortwave

  !>Get current wind speed for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_windspeed(i,n) RESULT(current_windspeed)
  
    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         uobsindex,      &
                         fileunit_uobs,  &
                         readmatlab,      &
                         uobscol,        &
                         get_current_uobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_windspeed(n)     !<current wind speed
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL    yt(uobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_uobs,uobscol,d,yt,missing_value,readmatlab) 
        current_windspeed = yt(uobsindex)
      ELSE
        yt = get_current_uobs(i,uobscol)
        current_windspeed = yt(uobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_windspeed

  !>Get current relative humidity for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_humidity(i,n) RESULT(current_humidity)
  
    USE WORLDVAR, ONLY : readdaily,       &
                         readformat,      &
                         rhobsindex,      &
                         fileunit_rhobs,  &
                         readmatlab,      &
                         rhobscol,        &
                         get_current_rhobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_humidity(n)     !<current relative humidity
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL    yt(rhobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_rhobs,rhobscol,d,yt,missing_value,readmatlab) 
        current_humidity = yt(rhobsindex)
      ELSE
        yt = get_current_rhobs(i,rhobscol)
        current_humidity = yt(rhobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_humidity

  !>Get current daily minimum air temperature for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_tmin(i,n) RESULT(current_tmin)
  
    USE WORLDVAR, ONLY : readdaily,         &
                         readformat,        &
                         tminobsindex,      &
                         fileunit_tminobs,  &
                         readmatlab,        &
                         tminobscol,        &
                         get_current_tminobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_tmin(n)          !<current observation
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL    yt(tminobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_tminobs,tminobscol,d,yt,missing_value,readmatlab) 
        current_tmin = yt(tminobsindex)
      ELSE
        yt = get_current_tminobs(i,tminobscol)
        current_tmin = yt(tminobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_tmin

  !>Get current daily maximum air temperature for selected subbasins from file or memory
  !---------------------------------------
  FUNCTION get_current_tmax(i,n) RESULT(current_tmax)
  
    USE WORLDVAR, ONLY : readdaily,       &
         readformat,      &
         tmaxobsindex,       &
         fileunit_tmaxobs,   &
         readmatlab,      &
         tmaxobscol,         &
         get_current_tmaxobs
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current time step
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    REAL current_tmax(n)          !<current observation
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL yt(tmaxobscol)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_tmaxobs,tmaxobscol,d,yt,missing_value,readmatlab) 
        current_tmax = yt(tmaxobsindex)
      ELSE
        yt = get_current_tmaxobs(i,tmaxobscol)
        current_tmax = yt(tmaxobsindex)
      ENDIF
    ENDIF
  END FUNCTION get_current_tmax

  !>Close files with observations used for forcing and calibration
  !--------------------------------------------------------------------
  SUBROUTINE close_observations(dir)

    USE WORLDVAR, ONLY : fileunit_tobs,  &
                         fileunit_pobs,  &
                         fileunit_qobs,  &
                         fileunit_xobs,  &
                         filename_xobs,  &
                         xcol,           &
                         readqobs,       &
                         readdaily,      &
                         readsfobs,      &
                         readswobs,      &
                         readwind,       &
                         readhumid,      &
                         readtminmaxobs, &
                         fileunit_sfobs, &
                         fileunit_swobs, &
                         fileunit_uobs, &
                         fileunit_rhobs, &
                         fileunit_tminobs, &
                         fileunit_tmaxobs

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory

    !Local variables
    LOGICAL xex

    IF(readdaily)THEN  
      CLOSE(fileunit_tobs)
      CLOSE(fileunit_pobs)
      IF(readqobs) CLOSE(fileunit_qobs)
      INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)
      IF(xex .AND. xcol>0) CLOSE(fileunit_xobs)
      IF(readsfobs)CLOSE(fileunit_sfobs)
      IF(readswobs)CLOSE(fileunit_swobs)
      IF(readwind)CLOSE(fileunit_uobs)
      IF(readhumid)CLOSE(fileunit_rhobs)
      IF(readtminmaxobs)THEN
        CLOSE(fileunit_tminobs)
        CLOSE(fileunit_tmaxobs)
      ENDIF
    ENDIF

  END SUBROUTINE close_observations

  !>Gets the information about soil-landuse-classes and subbasins
  !>from GeoClass and GeoData.
  !!
  !!\b Consequences Module worldvar variables pobsid, tobsid, sfobsid, 
  !!swobsid, uobsid, rhobsid, tminobsid, tmsxobsid may be set.
  !!Module modvar variables basin, classdata, vegtype, pathsubid, load, 
  !!wetland, nclass, nluse, nsoil, tiledepth, streamdepth, wetlandexist 
  !!and rotation is set. 
  !----------------------------------------------------------------
  SUBROUTINE load_basindata(dir,lfile,infile,n,n_Result,lakestate) 

    USE WORLDVAR, ONLY : fileunit_temp,   &
                         pobsid,          & !OUT
                         tobsid,          & !OUT
                         sfobsid,         & !OUT
                         swobsid,         & !OUT
                         uobsid,          & !OUT
                         rhobsid,         & !OUT
                         tminobsid,       & !OUT
                         tmaxobsid,       & !OUT
                         readsfobs,       &
                         readswobs,       &
                         readwind,        &
                         readhumid,       &
                         readtminmaxobs,  &
                         vegtypereset       !OUT
    USE MODVAR, ONLY : basin,             & !OUT
                       classdata,         & !OUT 
                       vegtype,           & !OUT
                       wetland,           & !OUT
                       nluse,nsoil,       & !OUT
                       nclass,            & !OUT
                       numsubstances,     &
                       tiledepth,         & !OUT
                       streamdepth,       & !OUT
                       set_soildepth,           &
                       set_soillayerdepth,      &
                       set_coded_classes,       &
                       allocate_basinvariables, &
                       wetlandexist,      & !OUT
                       rotation,          & !OUT
                       dam,               &
                       damindex
   USE STATETYPE_MODULE, ONLY : lakestatetype

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<File directory
    CHARACTER(LEN=*), INTENT(IN) :: lfile    !<Name of GeoClass file to be read
    CHARACTER(LEN=*), INTENT(IN) :: infile   !<Name of GeoData file to be read
    INTEGER, INTENT(OUT) :: n                !<Number of subbasins
    INTEGER, INTENT(OUT) :: n_Result         !<Error status
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate !MH2017: used to initialize conditioning
    
    !Local parameters
    INTEGER,PARAMETER :: maxclass = 200     !Maximum number classes to be read
    INTEGER,PARAMETER :: rcols = 50         !Number of data columns
    INTEGER,PARAMETER :: maxgccol = 13      !Maximum number of data columns in GeoClass
    INTEGER :: ncols                        !Total number of columns (rcols+number of s-lu-combinations (slc,dhslc,cr2) 
    
    !Local variables
    INTEGER i
    INTEGER mcols                   !Actual number of land-soil combinations in GeoData.txt
    REAL id(maxgccol,maxclass)      !id-data from Geoclass-file
    INTEGER idnum                   !Number of land-soil use combination in GeoClass.txt
    INTEGER temp                    !help variable
    INTEGER,ALLOCATABLE :: lakedataid(:)    !lakedataid from GeoData

    !> \b Algorithm \n
    !Start of subroutine
    n_Result = 0

    !>Read the key to basin characteristics file (GeoClass.txt)
    CALL read_geoid(fileunit_temp,TRIM(dir)//lfile,maxgccol,maxclass,id,idnum)

    nclass = idnum
    IF(.NOT.ALLOCATED(classdata)) ALLOCATE(classdata(nclass))
    IF(.NOT.ALLOCATED(vegtype)) ALLOCATE(vegtype(nclass))
    IF(.NOT.ALLOCATED(tiledepth)) ALLOCATE(tiledepth(nclass))
    IF(.NOT.ALLOCATED(streamdepth)) ALLOCATE(streamdepth(nclass))
    IF(.NOT.ALLOCATED(rotation)) ALLOCATE(rotation(nclass))
    nluse = NINT(MAXVAL(id(1,1:idnum)))
    nsoil = NINT(MAXVAL(id(2,1:idnum)))
    classdata(1:idnum)%luse  = NINT(id(1,1:idnum))
    classdata(1:idnum)%soil  = NINT(id(2,1:idnum))
    classdata(1:idnum)%crop  = NINT(id(3,1:idnum))
    classdata(1:idnum)%crop2 = NINT(id(4,1:idnum))
    rotation(1:idnum) = NINT(id(5,1:idnum))
    vegtype(1:idnum)     = NINT(id(6,1:idnum))
    tiledepth(1:idnum)   = id(8,1:idnum)
    streamdepth(1:idnum) = id(9,1:idnum)
    DO i = 1, idnum 
      temp = set_soildepth(NINT(id(10,i)),i,id(11:10+NINT(id(10,i)),i))
    ENDDO
    temp = set_soillayerdepth(NINT(id(10,1)),id(11:10+NINT(id(10,1)),1))
    !Check vegtype>0
    vegtypereset = .FALSE.
    DO i = 1,nclass
      IF(vegtype(i)==0)THEN
        vegtype(i)=1
        vegtypereset = .TRUE.
        WRITE(6,*) 'Warning: vegetation type not given in GeoClass.txt for slc-class',i,' vegtype=1 (open) is used'
      ENDIF
    ENDDO
    
    !Find coded classes
    CALL set_coded_classes(idnum,id(7,1:idnum))

    !>Count number of columns in GeoData
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,0,ncols,n_Result)
    IF(n_Result/=0)RETURN

    !>Count number of subbasins in GeoData
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,1,n,n_Result)
    IF(n_Result/=0)RETURN
    
    !>Initiate subbasin data
    CALL allocate_basinvariables(n,numsubstances)
    
!    basin(1:n)%ilakecatch = 1.
    basin(1:n)%ilakecatch = 0.
    wetland(1:n,1:2)%area = 0.
    
    !MH2017: crashing from stack overflow (170)
    !symptoms characteristic of oversize local array (crashing on CALL read_and_calc...)
    !lakedataid(n) below is the likely suspect; at this point n(rows) = 6668
    !this is probably far too large for an array on the stack
    IF(.NOT.ALLOCATED(lakedataid)) ALLOCATE(lakedataid(n))
    IF(.NOT.ALLOCATED(pobsid)) ALLOCATE(pobsid(n))
    IF(.NOT.ALLOCATED(tobsid)) ALLOCATE(tobsid(n))
    lakedataid = 0
    pobsid = 0; tobsid = 0
    IF(readsfobs)THEN
      IF(.NOT.ALLOCATED(sfobsid))ALLOCATE(sfobsid(n))
      sfobsid = 0
    ENDIF
    IF(readswobs)THEN
      IF(.NOT.ALLOCATED(swobsid))ALLOCATE(swobsid(n))
      swobsid = 0
    ENDIF
    IF(readwind)THEN
      IF(.NOT.ALLOCATED(uobsid))ALLOCATE(uobsid(n))
      uobsid = 0
    ENDIF
    IF(readhumid)THEN
      IF(.NOT.ALLOCATED(rhobsid))ALLOCATE(rhobsid(n))
      rhobsid = 0
    ENDIF
    IF(readtminmaxobs)THEN
      IF(.NOT.ALLOCATED(tminobsid))ALLOCATE(tminobsid(n))
      tminobsid = 0
      IF(.NOT.ALLOCATED(tmaxobsid))ALLOCATE(tmaxobsid(n))
      tmaxobsid = 0
    ENDIF   

    !>Read the description of the basins (GeoData.txt)
    CALL read_and_calc_basindata(fileunit_temp,TRIM(dir)//infile,n,lakedataid,ncols,mcols)

    !Check if wetland is present
    IF(SUM(wetland(1:n,1:2)%area)==0)THEN
       wetlandexist = .FALSE.
       IF(ALLOCATED(wetland)) DEALLOCATE(wetland)
    ENDIF

    !>Read irrigation information of the basins (MgmtData.txt)
    CALL load_irrigation_data(fileunit_temp,TRIM(dir),n_Result)
    IF(n_Result/=0)RETURN

    !>Read observation coupling and elevation from separate file? (ForcData.txt)
    CALL load_forcing_information_data(fileunit_temp,TRIM(dir),n,n_Result)
    IF(n_Result/=0)RETURN

    !>Read additional description of important lakes (LakeData.txt) (possibly overwrite GeoData information)
    CALL load_lakedata(fileunit_temp,dir,n,lakedataid,n_Result)
    IF(n_Result/=0)RETURN
    IF(ALLOCATED(lakedataid)) DEALLOCATE(lakedataid)

    !>Read additional description of important dams (DamData.txt)
    CALL load_damdata(fileunit_temp,dir,n,n_Result)
    IF(n_Result/=0)RETURN
    
    !************************************************************
    !MH2017: Read additional description of special dams (ExtraData.txt)
    CALL load_extradata(fileunit_temp,dir,n,n_result,lakestate)
    IF(n_Result/=0)RETURN
    !************************************************************
    !MH2017: Read daily array for each special dam (DailyDataXXXX.txt)
    DO i=1,SIZE(dam)
      IF (dam(i)%purpose>=5)THEN
        CALL load_dailydataXXXX(fileunit_temp,dir,n,n_result,i)
        IF(n_result/=0)RETURN
      ENDIF
    ENDDO
    !************************************************************

  END SUBROUTINE load_basindata

  !>Set model choice configuration
  !>
  !>/b Consequences Module modvar variables dimriverlag, glacierexist 
  !>and doirrigation is set
  !------------------------------------------------------------------
  SUBROUTINE set_model_configuration(nsubb,idir)
  !Vet ej än hur mycket denna ska täcka. Bara modell, även simulering?
  
  USE MODVAR, ONLY : nclass,classmodel, &
                     glacierexist,doirrigation, &  !OUT
                     deallocate_glacier,  &
                     dimriverlag  !OUT
  USE MODELMODULE, ONLY : set_special_models, &
                          set_special_parameters_from_model,  &
                          calculate_special_model_parameters
  
  !Argument declaration
  INTEGER, INTENT(IN) :: nsubb  !<number of subbasin, submodel
  CHARACTER(LEN=*), INTENT(IN) :: idir    !<file directory
  
  !Local variables
  REAL optrivvel,optdamp
  INTEGER rvindex,dpindex   !index of rivvel and damp in modparid (a hypevariable)
  
  CALL set_special_models(nclass,classmodel,glacierexist,doirrigation)
  CALL set_special_parameters_from_model(rvindex,dpindex)

  !Find maximum dimension of river queue variables
  CALL calculate_special_optpar_parameters(idir,rvindex,dpindex,optrivvel,optdamp)  !optpar.txt
  CALL calculate_special_model_parameters(nsubb,optrivvel,optdamp,dimriverlag)  !par.txt,reg_par.txt
  
  !Deallocate glacier-variables
  IF(.NOT.glacierexist) CALL deallocate_glacier()
  
  END SUBROUTINE set_model_configuration

  !>Get the information about branches from BranchData.txt
  !
  !>\b Consequences Module modvar variable branchsubid will be set
  !------------------------------------------------------------------
  SUBROUTINE load_branchdata(dir,n_Result) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp, &
                         maxcharpath
    USE MODVAR, ONLY : branchsubid            !OUT

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir    !<File directory
    INTEGER, INTENT(OUT) :: n_Result        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of branches in BranchData.txt
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    n_Result = 0
    filename = TRIM(dir)//'BranchData.txt'

    !> \b Algorithm \n
    !>Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN   

    !Count number of columns in BranchData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,n_Result)
    IF(n_Result/=0)RETURN

    !Count number of rows in BranchData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,n_Result)
    IF(n_Result/=0)RETURN

    !>Initiate branch data
    IF(.NOT.ALLOCATED(branchsubid)) ALLOCATE(branchsubid(nrows))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    branchsubid%source = 0
    branchsubid%branch = 0
    branchsubid%mainpart  = 1.
    branchsubid%maxQ = 0.
    branchsubid%minQ = 0.
    branchsubid%maxQbranch = 0.

    !>Read BranchData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,n_Result)
    IF(n_Result.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='sourceid  ') code(i) = i_intg
      IF(str(i)(1:10)=='branchid  ') code(i) = i_intg
      IF(str(i)(1:10)=='mainpart  ') code(i) = i_real
      IF(str(i)(1:10)=='maxqmain  ') code(i) = i_real
      IF(str(i)(1:10)=='maxqbranch') code(i) = i_real
      IF(str(i)(1:10)=='minqmain  ') code(i) = i_real
    ENDDO

    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !>Save loaded branch data temporary
    DO i = 1,mcols
      IF(str(i)(1:10)=='sourceid  ')   branchsubid(1:nrows)%source     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='branchid  ')   branchsubid(1:nrows)%branch     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='mainpart  ')   branchsubid(1:nrows)%mainpart   = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='maxqmain  ')   branchsubid(1:nrows)%maxQ       = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='minqmain  ')   branchsubid(1:nrows)%minQ       = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='maxqbranch')   branchsubid(1:nrows)%maxQbranch = xr(1:nrows,rindex(i))
    ENDDO
    
    !>Check for illegal combinations of branch data
    DO i = 1,nrows
      IF(branchsubid(i)%maxQbranch>0)THEN
        IF(branchsubid(i)%maxQ>0)THEN
          branchsubid(i)%maxQbranch = 0.
          WRITE(6,*) 'WARNING: Maximum flow in both main channel and branch is given for source: ',branchsubid(i)%source
          WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of branch will be removed.'
          WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
        ENDIF
      ENDIF
      IF(branchsubid(i)%maxQbranch<0)THEN
        branchsubid(i)%maxQbranch = 0.
        WRITE(6,*) 'WARNING: Maximum flow in branch is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of branch will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
      IF(branchsubid(i)%maxQ<0)THEN
        branchsubid(i)%maxQ = 0.
        WRITE(6,*) 'WARNING: Maximum flow in main channel is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of main channel will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
      IF(branchsubid(i)%minQ<0)THEN
        branchsubid(i)%minQ = 0.
        WRITE(6,*) 'WARNING: Minimum flow in main channel is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Minimum flow of main channel will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
    ENDDO

    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_branchdata

  !>Get the information about aquifers from AquiferData.txt
  !!
  !>\b Consequences Module modvar variable aquifer and pathsubid will be set. 
  !!Module modvar varible nregions may be increased.
  !------------------------------------------------------------------
  SUBROUTINE load_aquiferdata(dir,nsub,numaqu,n_Result) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp, &
                         maxcharpath
    USE MODVAR, ONLY : basin,     &
                       pathsubid, &  !OUT
                       aquifer,   &  !OUT
                       nregions,  &  !OUT
                       modeloption, &
                       p_deepgroundwater

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir    !<File directory
    INTEGER, INTENT(IN)  :: nsub            !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: numaqu          !<Number of aquifers
    INTEGER, INTENT(OUT) :: n_Result        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,j
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of rows in AquiferData.txt
    INTEGER subid_i,aquid_i,recha_i,outfr_i   !data column index for these columns
    INTEGER area_i,depth_i,basdp_i,topdp_i   !data column index for these columns
    INTEGER rate_i,delay_i,por_i,reg_i  !data column index for these columns
    INTEGER temp_i,coIN_i,coSP_i  !data column index for these columns
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    n_Result = 0
    numaqu = 0
    filename = TRIM(dir)//'AquiferData.txt'

    !Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(fex)THEN
      IF(modeloption(p_deepgroundwater)/=2)THEN
        WRITE(6,*) 'Warning: Existing AquiferData.txt will not be used, since modeloption not set.'
        RETURN
      ENDIF
    ELSE
      IF(modeloption(p_deepgroundwater)==2)THEN
        WRITE(6,*) 'Warning: Use modeloption with aquifer without AquiferData.txt not possible.'
        WRITE(6,*) 'Warning: Deep groundwater modeloption turned off.'
        modeloption(p_deepgroundwater)=0
        RETURN
      ENDIF
      RETURN   
    ENDIF

    !Count number of columns in AquiferData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,n_Result)
    IF(n_Result/=0)RETURN

    !Count number of rows in AquiferData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,n_Result)
    IF(n_Result/=0)RETURN

    !Initiate branch data
    IF(.NOT.ALLOCATED(pathsubid)) ALLOCATE(pathsubid(nsub))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    pathsubid%aquid = 0
    pathsubid%rechargebasin = .FALSE.
    pathsubid%recievefraction = 0.

    !Read AquiferData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,n_Result)
    IF(n_Result.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='aquid     ') code(i) = i_intg
      IF(str(i)(1:10)=='parreg    ') code(i) = i_intg
      IF(str(i)(1:10)=='recharge  ') code(i) = i_intg
      IF(str(i)(1:10)=='area      ') code(i) = i_real
      IF(str(i)(1:10)=='inidepth  ') code(i) = i_real
      IF(str(i)(1:10)=='topdepth  ') code(i) = i_real
      IF(str(i)(1:10)=='basedepth ') code(i) = i_real
      IF(str(i)(1:10)=='porosity  ') code(i) = i_real
      IF(str(i)(1:10)=='retfrac   ') code(i) = i_real
      IF(str(i)(1:10)=='retrate   ') code(i) = i_real
      IF(str(i)(1:10)=='delay     ') code(i) = i_real
      IF(str(i)(1:10)=='temp      ') code(i) = i_real
      IF(str(i)(1:10)=='conc_in   ') code(i) = i_real
      IF(str(i)(1:10)=='conc_sp   ') code(i) = i_real
    ENDDO

    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Find and check data columns
    subid_i=0;aquid_i=0;recha_i=0;area_i=0;depth_i=0;basdp_i=0
    topdp_i=0;por_i=0;outfr_i=0;rate_i=0;delay_i=0;temp_i=0;coIN_i=0;coSP_i=0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')   subid_i = iindex(i)
      IF(str(i)(1:10)=='aquid     ')   aquid_i = iindex(i)
      IF(str(i)(1:10)=='parreg    ')   reg_i   = iindex(i)
      IF(str(i)(1:10)=='recharge  ')   recha_i = iindex(i)
      IF(str(i)(1:10)=='area      ')   area_i  = rindex(i)
      IF(str(i)(1:10)=='inidepth  ')   depth_i = rindex(i)
      IF(str(i)(1:10)=='basedepth ')   basdp_i = rindex(i)
      IF(str(i)(1:10)=='topdepth  ')   topdp_i = rindex(i)
      IF(str(i)(1:10)=='porosity  ')   por_i   = rindex(i)
      IF(str(i)(1:10)=='retfrac   ')   outfr_i = rindex(i)
      IF(str(i)(1:10)=='retrate   ')   rate_i  = rindex(i)
      IF(str(i)(1:10)=='delay     ')   delay_i = rindex(i)
      IF(str(i)(1:10)=='temp      ')   temp_i  = rindex(i)
      IF(str(i)(1:10)=='conc_in   ')   coIN_i  = rindex(i)
      IF(str(i)(1:10)=='conc_sp   ')   coSP_i  = rindex(i)
    ENDDO
    IF(subid_i==0 .OR. aquid_i==0 .OR. recha_i==0 .OR. area_i==0 .OR. &
       depth_i==0 .OR. basdp_i==0 .OR. por_i==0  .OR. &
       outfr_i==0 .OR. rate_i==0  .OR. delay_i==0)THEN
      WRITE(6,*) 'Some columns are missing in AquiferData.txt'
      n_Result=1
      RETURN
    ENDIF
    
    !Save loaded data
    numaqu = MAXVAL(xi(:,aquid_i))
    ALLOCATE(aquifer(numaqu))
    aquifer%temperature = 0.
    aquifer%conc_IN = 0.
    aquifer%conc_SP = 0.
    DO i = 1,nrows
      IF(xi(i,subid_i)>0)THEN   !Row with aquifer subbasin coupling information
        DO j = 1,nsub
          IF(basin(j)%subid==xi(i,subid_i)) EXIT
        ENDDO
        IF(j<=nsub)THEN
          pathsubid(j)%aquid = xi(i,aquid_i)
          pathsubid(j)%rechargebasin = (xi(i,recha_i)==1)
          pathsubid(j)%recievefraction = xr(i,outfr_i)
        ELSE
          WRITE(6,*) 'Warning: subid in AquiferData not in GeoData', xi(i,subid_i)
        ENDIF
      ELSEIF(xi(i,subid_i)==0)THEN  !Row with aquifer information
        aquifer(xi(i,aquid_i))%area = xr(i,area_i)
        aquifer(xi(i,aquid_i))%parregion = xi(i,reg_i)
        aquifer(xi(i,aquid_i))%retrate = xr(i,rate_i)
        aquifer(xi(i,aquid_i))%percdelay = xr(i,delay_i)
        aquifer(xi(i,aquid_i))%inivol = (xr(i,depth_i)-xr(i,basdp_i))*xr(i,area_i)*xr(i,por_i)  !m3
        IF(topdp_i>0) aquifer(xi(i,aquid_i))%maxvol = (xr(i,topdp_i)-xr(i,basdp_i))*xr(i,area_i)*xr(i,por_i)  !m3
        aquifer(xi(i,aquid_i))%basedepth = xr(i,basdp_i)
        aquifer(xi(i,aquid_i))%porosity = xr(i,por_i)
        IF(temp_i>0) aquifer(xi(i,aquid_i))%temperature = xr(i,temp_i)
        IF(coIN_i>0) aquifer(xi(i,aquid_i))%conc_IN = xr(i,coIN_i)*1.E-3
        IF(coSP_i>0) aquifer(xi(i,aquid_i))%conc_SP = xr(i,coSP_i)*1.E-3
      ENDIF
    ENDDO
    
    !Check/Increase number of parameter regions
    nregions = MAX(nregions,MAXVAL(aquifer(1:numaqu)%parregion))

    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_aquiferdata

  !>Get the information about glaciers from GlacierData.txt
  !!
  !>\b Consequences Module modvar variables glacier and glacierindex will be set. 
  !------------------------------------------------------------------
  SUBROUTINE load_glacierdata(dir,nsb,n_Result) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp, &
                         maxcharpath
    USE MODVAR, ONLY : basin,     &
                       glacier,   &     !OUT
                       nglaciers, &     !OUT
                       glacierindex, &  !OUT
                       missing_value

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: nsb           !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: n_Result      !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,j
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of rows in AquiferData.txt
    INTEGER subid_i         !data column index for these columns
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: glaciersubid(:)  !subid from GlacierData
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    n_Result = 0
    filename = TRIM(dir)//'GlacierData.txt'

    !Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN   

    !Count number of columns in GlacierData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,n_Result)
    IF(n_Result/=0)RETURN

    !Count number of rows in GlacierData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,n_Result)
    IF(n_Result/=0)RETURN

    !Initiate
    IF(.NOT.ALLOCATED(glaciersubid)) ALLOCATE(glaciersubid(nrows))
    IF(.NOT.ALLOCATED(glacier))      ALLOCATE(glacier(nrows))
    IF(.NOT.ALLOCATED(glacierindex)) ALLOCATE(glacierindex(nsb))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    glacier(1:nrows)%gtype = INT(missing_value)     !glacier type (0 glacier, 1 ice cap)
    glacier(1:nrows)%volcorr = 0.                   !log vol correction
    glacier(1:nrows)%glacdate = INT(missing_value)  !date of SLC area (used to rescale initial volume)
    glacier(1:nrows)%glacinimb = 0.                 !annual mass balance for initial glacier volume correction

    !Read GlacierData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,n_Result)
    IF(n_Result.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    subid_i=0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        subid_i = i
      ENDIF
      IF(str(i)(1:10)=='glactype  ') code(i) = i_intg
      IF(str(i)(1:10)=='logvolcor ') code(i) = i_real
      IF(str(i)(1:10)=='glacdate  ') code(i) = i_intg
      IF(str(i)(1:10)=='glacinimb ') code(i) = i_real
    ENDDO

    !Check data columns
    IF(subid_i==0)THEN
      WRITE(6,*) 'SUBID columns are missing in GlacierData.txt'
      n_Result=1
      RETURN
    ENDIF

    !Read data
    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Set glacier characheristics array
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') glaciersubid(1:nrows)    = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='glactype  ') glacier(1:nrows)%gtype   = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='logvolcor ') glacier(1:nrows)%volcorr = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='glacdate  ') glacier(1:nrows)%glacdate= xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='glacinimb ') glacier(1:nrows)%glacinimb = xr(1:nrows,rindex(i))
    ENDDO
       
    !Calculate glacier finding index array
    DO j = 1,nrows
      DO i = 1,nsb
        IF(basin(i)%subid==glaciersubid(j))THEN
          glacierindex(i) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO
    
    nglaciers = nrows
        
    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(glaciersubid)) DEALLOCATE(glaciersubid)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_glacierdata

  !>Collects the information about special lakes from LakeData.txt
  !>
  !>\b Consequences Module modvar variables lake, lakebasin, lakeindex and lakebasinindex
  !>may be allocated and set. Module modvar variables basin, nbasinlakes may be changed.
  !------------------------------------------------------------------
  SUBROUTINE load_lakedata(funit,dir,n,lakedataid,n_Result) 

    USE WORLDVAR, ONLY : maxcharpath,  &
                         checkdata, &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         lake,         &    !OUT
                         lakebasin,    &    !OUT
                         lakeindex,    &    !INOUT
                         lakebasinindex, &  !OUT
                         max_par,      &    
                         modparid,     &    
                         m_ldpar,      &
                         nbasinlakes        !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                   !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER, INTENT(IN)  :: n                       !<Number of subbasins (base model)
    INTEGER, INTENT(IN)  :: lakedataid(n)           !<lakedataid from GeoData
    INTEGER, INTENT(OUT) :: n_Result                !<Error status
    
    !Local parameters
    INTEGER,PARAMETER :: maxcol = 50             !Maximum number of columns in LakeData
    INTEGER,PARAMETER :: i_simple = 1            !Lakedata type; simple olake
    INTEGER,PARAMETER :: i_combil = 2            !Lakedata type; lake composed of several lake basins
    INTEGER,PARAMETER :: i_lbasin = 3            !Lakedata type; lake basin
    INTEGER,PARAMETER :: i_llbasin = 4           !Lakedata type; last lake basin
    
    !Local variables
    LOGICAL fex
    INTEGER i,j,ilake,ibasin
    INTEGER ictype                    !Column with lakedata type ('ldtype')
    INTEGER iclakeid                  !Column with lakeid
    INTEGER iclakedataid              !Column with id (lakedataid in GeoData)
    INTEGER ica                       !Column with area
    INTEGER icdepth                   !Column with depth
    INTEGER nrows                     !Number of lakes in file
    INTEGER ncols                     !Number of columns in file
    INTEGER nlake,nbasin,nlakefb,nlakedim
    INTEGER code(maxcol)              !Code for column variable
    INTEGER rindex(maxcol)            !Index for column real variables
    INTEGER iindex(maxcol)            !Index for column integer variables
    INTEGER MonthDayDate              !Non-formated date from lakeData.txt
    INTEGER dayNumber                 !Corresponding date in day number format
    REAL    lakearea                  !lakearea calculated from GeoData
    REAL    parea
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10) str(maxcol)     !Content string
    INTEGER,ALLOCATABLE :: lakeindex2(:)      !row->isub
    INTEGER,ALLOCATABLE :: lakeindex3(:)      !row->ilake
    INTEGER,ALLOCATABLE :: lastlakeindex(:)   !ilake->row of last lakebasin for this lakebasinlake
    INTEGER,ALLOCATABLE :: qprod1date(:)
    INTEGER,ALLOCATABLE :: qprod2date(:)
    INTEGER,ALLOCATABLE :: lakeidindex(:,:)
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL,ALLOCATABLE :: regvolume(:)          !Regulation volume (M m3)
    REAL,ALLOCATABLE :: wamp(:)               !Regulation amplitude (real) (m)

    n_Result = 0
    nbasinlakes = 0
    filename = TRIM(dir)//'LakeData.txt'

    !Check if file exist and if lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with LakeData.txt'
      n_Result = 1
      RETURN
    ENDIF

    !Count number of data rows, and allocate variables for reading LakeData
    CALL count_data_rows(funit,filename,1,nrows,n_Result)
    IF(n_Result/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No lakes in LakeData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))

    !Open LakeData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    ica = 0
    ictype = 0
    iclakedataid = 0
    iclakeid = 0
    icdepth = 0

    DO i = 1,ncols
      IF(str(i)(1:10)=='area      ')THEN
        code(i) = i_real
        ica = i
      ENDIF
      IF(str(i)(1:10)=='ldtype    ')THEN
        code(i) = i_intg   !type of info: lake, basin or whole lake (that is divided in basins)
        ictype = i
      ENDIF
      IF(str(i)(1:10)=='lakedataid')THEN
        code(i) = i_intg   !lakedataid in GeoData
        iclakedataid = i
      ENDIF
      IF(str(i)(1:10)=='lake_depth')THEN
        code(i) = i_real   !lake_depth
        icdepth = i
      ENDIF
      IF(str(i)(1:10)=='w0ref     ') code(i) = i_real   !w0
      IF(str(i)(1:10)=='deltaw0   ') code(i) = i_real   !deltaw0
      IF(str(i)(1:10)=='rate      ') code(i) = i_real   !lake_rate
      IF(str(i)(1:10)=='exp       ') code(i) = i_real   !lake_exp
      IF(str(i)(1:10)=='regvol    ') code(i) = i_real   !regvol
      IF(str(i)(1:10)=='wamp      ') code(i) = i_real   !reg. amplitude
      IF(str(i)(1:10)=='qprod1    ') code(i) = i_real   !two date-based production rates, qprod1 and qprod2, now replacing the unique rate qprod
      IF(str(i)(1:10)=='qprod2    ') code(i) = i_real
      IF(str(i)(1:10)=='datum1    ') code(i) = i_intg   !two dates for the corresponding two dam production flows (date of start)               
      IF(str(i)(1:10)=='datum2    ') code(i) = i_intg
      IF(str(i)(1:10)=='qamp      ') code(i) = i_real   !lake_qamp
      IF(str(i)(1:10)=='qpha      ') code(i) = i_real   !lake_qpha
      IF(str(i)(1:10)=='lakeid    ')THEN
        code(i) = i_intg
        iclakeid = i
      ENDIF
      !All model parameters that can be given in LakeData is real
      DO j=1,max_par
        IF(str(i)(1:10)==modparid(j)%shortname .AND. modparid(j)%deptype==m_ldpar) code(i) = i_real
      ENDDO
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    !Count number of laketypes: lakes, lakes formed by basins and lakebasins.
    IF(ictype==0)THEN
      WRITE(6,*) 'ERROR: ldtype not found in LakeData. ',TRIM(filename)
      n_Result = 1
      RETURN
    ELSE
      nlake = 0
      nbasin = 0
      nlakefb = 0
      DO j = 1,nrows
        IF(xi(j,iindex(ictype))==i_simple) nlake = nlake + 1
        IF(xi(j,iindex(ictype))==i_combil) nlakefb = nlakefb + 1
        IF(xi(j,iindex(ictype))==i_lbasin.OR.xi(j,iindex(ictype))==i_llbasin) nbasin = nbasin + 1
      ENDDO
    ENDIF
    nlakedim = nlake+nlakefb

    !Allocate variables for lake information
    IF(.NOT.ALLOCATED(lakebasin)) ALLOCATE(lakebasin(nbasin))
    IF(.NOT.ALLOCATED(lake)) ALLOCATE(lake(nlakedim))
    IF(.NOT.ALLOCATED(regvolume)) ALLOCATE(regvolume(nlakedim))
    IF(.NOT.ALLOCATED(wamp)) ALLOCATE(wamp(nlakedim))
    IF(.NOT.ALLOCATED(qprod1date)) ALLOCATE(qprod1date(nlakedim)) 
    IF(.NOT.ALLOCATED(qprod2date)) ALLOCATE(qprod2date(nlakedim))
    regvolume = 0.; wamp = missing_value
    qprod1date = 0; qprod2date = 0  

    !Find subbasin-lake coupling: lakeindex(isb),lakebasinindex(isb), lakeindex2(ldrow) and lakeindex3(ldrow)
    IF(iclakedataid==0)THEN
      WRITE(6,*) 'ERROR: id not found in LakeData. ',TRIM(filename)
      n_Result = 1
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(lakeindex)) ALLOCATE(lakeindex(n))
    IF(nbasin>0 .AND.(.NOT.ALLOCATED(lakebasinindex))) ALLOCATE(lakebasinindex(n))
    IF(.NOT.ALLOCATED(lakeindex2)) ALLOCATE(lakeindex2(nrows))
    IF(.NOT.ALLOCATED(lakeindex3)) ALLOCATE(lakeindex3(nrows))
    IF(.NOT.ALLOCATED(lastlakeindex)) ALLOCATE(lastlakeindex(nrows))
    lakeindex = 0
    IF(nbasin>0) lakebasinindex = 0
    lakeindex2 = 0
    lakeindex3 = 0
    lastlakeindex = 0
    DO i = 1,n
      IF(lakedataid(i)==0)CYCLE
      DO j = 1,nrows
        IF(lakedataid(i)==xi(j,iindex(iclakedataid)))THEN
          lakeindex2(j)=i
          EXIT
        ENDIF
      ENDDO
      IF(j==nrows+1)THEN
        WRITE(6,*) 'Error: A lake in GeoData is missing in LakeData. lakedataid: ',lakedataid(i)
        STOP 1
      ENDIF  
    ENDDO

    !First take care of the lakes that are formed by lakebasins
    IF(nlakefb>0)THEN
      IF(ica==0.OR.iclakeid==0)THEN             !Check for necessary data
        WRITE(6,*) 'ERROR: area or lakeid not found in LakeData. ',TRIM(filename)
        WRITE(6,*) 'ERROR: These are necessary for lakes formed by basins (ldtype 3).'
        n_Result = 1
        RETURN
      ENDIF
      IF(.NOT.ALLOCATED(lakeidindex)) ALLOCATE(lakeidindex(nlakefb,2)) !lakeid-jrow coupling for later use in subroutine
      ilake = 0
      DO j = 1,nrows
        IF(xi(j,iindex(ictype))==i_combil)THEN
          ilake = ilake + 1
          lakeidindex(ilake,1) = xi(j,iindex(iclakeid))
          lakeidindex(ilake,2) = j
        ENDIF
      ENDDO
      DO ilake = 1,nlakefb                      !Save information permanent or temporary
        DO i = 1,ncols
          IF(str(i)=='area      ')   lake(ilake)%area      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref     = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0   = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1    = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2    = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='datum1    ')   qprod1date(ilake)     = xi(lakeidindex(ilake,2),iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)     = xi(lakeidindex(ilake,2),iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp       = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)           = xr(lakeidindex(ilake,2),rindex(i))
        ENDDO
        lakeindex3(lakeidindex(ilake,2))=ilake
      ENDDO
    ENDIF

    !Simple lakes and lakebasins
    ilake = nlakefb
    ibasin = 0
    DO j = 1,nrows
      IF(lakeindex2(j)==0) CYCLE    !skip rows with lakes formed by basins and lakes not included in simulation
      IF(xi(j,iindex(ictype))==i_simple)THEN     !single lake
        ilake = ilake + 1
        lake(ilake)%area = xr(j,rindex(ica))    !is this used? remove?
        IF(icdepth>0)THEN
          IF(xr(j,rindex(icdepth)).NE.missing_value) basin(lakeindex2(j))%lakedepth(2) = xr(j,rindex(icdepth))
        ENDIF
        DO i = 1,ncols
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref    = xr(j,rindex(i))
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0  = xr(j,rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1   = xr(j,rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2   = xr(j,rindex(i))
          IF(str(i)=='datum1    ')   qprod1date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate     = xr(j,rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp      = xr(j,rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp     = xr(j,rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha     = xr(j,rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)     = xr(j,rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)          = xr(j,rindex(i))
        ENDDO
        lakeindex(lakeindex2(j))=ilake
        lakeindex3(j)=ilake
      ELSEIF(xi(j,iindex(ictype))==i_lbasin.OR.xi(j,iindex(ictype))==i_llbasin)THEN     !lake basin
        ibasin = ibasin + 1
        IF(xr(j,rindex(icdepth)).NE.missing_value) basin(lakeindex2(j))%lakedepth(2) = xr(j,rindex(icdepth))
        DO i = 1,nlakefb
          IF(lakeidindex(i,1)==xi(j,iindex(iclakeid))) EXIT
        ENDDO
        IF(i>nlakefb)THEN 
          WRITE(6,*) 'ERROR: lakebasin misses lake. lakeid',xi(j,iindex(iclakeid))
          WRITE(6,*) 'ERROR: can not find row with info on whole lake.'
          n_Result = 1
          RETURN
        ENDIF
        lakebasin(ibasin)%ilk = i
        IF(xi(j,iindex(ictype))==i_llbasin) lastlakeindex(lakebasin(ibasin)%ilk) = j
        IF(xi(j,iindex(ictype))==i_llbasin) lakebasin(ibasin)%last = .TRUE.
        lakebasinindex(lakeindex2(j))=ibasin
        !These have all w0ref=0, since not set
      ENDIF
    ENDDO
    IF(ilake<nlakedim) WRITE(6,*) ilake,nlakedim
    IF(ibasin>0) nbasinlakes = lakebasin(ibasin)%ilk

    !Calculate other lake variables, w0ref,wmin, datum, w01, w02,wampcoeff
    lake(:)%wmin = missing_value
    lake(:)%wampcoeff = missing_value
    DO j = 1,nrows
      IF(lakeindex3(j)>0)THEN
        IF(regvolume(lakeindex3(j))>0.)THEN
          IF(xi(j,iindex(ictype))==i_simple)THEN
            lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          ELSEIF(xi(j,iindex(ictype))==i_combil)THEN
            lakearea = lake(lakeindex3(j))%area
          ENDIF
          IF(lakearea>0)THEN
            lake(lakeindex3(j))%wmin = 0. - regvolume(lakeindex3(j)) * 1000000. / lakearea
            IF(wamp(lakeindex3(j))/=missing_value) &
              lake(lakeindex3(j))%wampcoeff = wamp(lakeindex3(j))/(-1.*lake(lakeindex3(j))%wmin)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    DO ilake = 1,nlakedim 
      MonthDayDate = qprod1date(ilake)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      lake(ilake)%datum1 = dayNumber
      MonthDayDate = qprod2date(ilake)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      lake(ilake)%datum2 = dayNumber
    ENDDO

    !Check lake area against slc-area
    IF(checkdata(1,1))THEN
      fex = .FALSE.
      DO j=1,nrows
        lakearea = 0
        parea = 0.
        IF(lakeindex2(j)>0 .AND. lakeindex3(j)>0)THEN
          lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          IF(lakearea>0)THEN
            parea = ABS(lake(lakeindex3(j))%area-lakearea)/lakearea
          ELSEIF(lake(lakeindex3(j))%area>0)THEN
            parea = 0.
          ENDIF
        ENDIF
        IF(parea > 0.001)THEN
          WRITE(6,'(a29,es14.7,a28)') ' INFO: lakearea in LakeData (',lake(lakeindex3(j))%area,') differ more than 0.1% from'
          WRITE(6,'(a29,es14.7,a1)')  ' INFO: lakearea in GeoData  (',lakearea,')'
          WRITE(6,'(a18,I8)') ' INFO: for subid =',basin(lakeindex2(j))%subid 
          fex = .TRUE.
        ENDIF
      ENDDO
      IF(fex) WRITE(6,*)
    ENDIF

    !Build table with lakedata parameter values
    CALL start_lakedata_table(maxcol,ncols,nrows,str,xi,xr,iindex,rindex,n,lakedataid)

    !Deallocate local variables
    IF(ALLOCATED(lakeindex2)) DEALLOCATE(lakeindex2)
    IF(ALLOCATED(lakeindex3)) DEALLOCATE(lakeindex3)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(regvolume)) DEALLOCATE(regvolume)
    IF(ALLOCATED(qprod1date)) DEALLOCATE(qprod1date) 
    IF(ALLOCATED(qprod2date)) DEALLOCATE(qprod2date)
    IF(ALLOCATED(lakeidindex)) DEALLOCATE(lakeidindex)
    IF(ALLOCATED(lastlakeindex)) DEALLOCATE(lastlakeindex)

    WRITE(6,*) 'File read: ', TRIM(filename)

  END SUBROUTINE load_lakedata

  !>Collects the information about dams from DamData.txt
  !>
  !>\b Consequences Module modvar variables dam, damindex may be allocated and set. 
  !>Module modvar variable basin may be changed (lake_depth).
  !------------------------------------------------------------------
  SUBROUTINE load_damdata(funit,dir,n,n_Result)

    USE WORLDVAR, ONLY : maxcharpath,  &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         dam,         &    !OUT
                         damindex        !INOUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                   !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER, INTENT(IN)  :: n                       !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: n_Result                !<Error status
    
    !Local parameters
    INTEGER,PARAMETER :: maxcol = 50             !Maximum number of columns in DamData
    
    !Local variables
    LOGICAL fex
    INTEGER i,idam
    INTEGER idtype                    !Column with subid
    INTEGER nrows                     !Number of dams in file
    INTEGER ncols                     !Number of columns in file
    INTEGER code(maxcol)              !Code for column variable
    INTEGER rindex(maxcol)            !Index for column real variables
    INTEGER iindex(maxcol)            !Index for column integer variables
    INTEGER MonthDayDate              !Non-formated date from DamData.txt
    INTEGER dayNumber                 !Corresponding date in day number format
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10) str(maxcol)     !Content string
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL damarea                            ! area of dam reservoir (same as lakearea)
    INTEGER,ALLOCATABLE :: qprod1date(:)
    INTEGER,ALLOCATABLE :: qprod2date(:)
    INTEGER, ALLOCATABLE :: damindex2(:)   ! row where  
    REAL,ALLOCATABLE :: wamp(:)               !Regulation amplitude (real) (m)
    REAL,ALLOCATABLE :: qinflow(:,:)          !Monthly natural flow (m3/s)
    REAL,ALLOCATABLE :: minbranchQ(:,:)     !Monthly minimum dam branch flow parameter (m3/s)

    n_Result = 0
    ALLOCATE(damindex(n))
    damindex = 0
    filename = TRIM(dir)//'DamData.txt'

    !Check if file exist and if outlet lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with DamData.txt'
      n_Result = 1
      RETURN
    ENDIF

    !Count number of data rows, and allocate variables for reading DamData
    CALL count_data_rows(funit,filename,1,nrows,n_Result)
    IF(n_Result/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No dam in DamData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))
    IF(.NOT.ALLOCATED(qprod1date)) ALLOCATE(qprod1date(nrows)) 
    IF(.NOT.ALLOCATED(qprod2date)) ALLOCATE(qprod2date(nrows))
    IF(.NOT.ALLOCATED(wamp)) ALLOCATE(wamp(nrows))
    IF(.NOT.ALLOCATED(qinflow)) ALLOCATE(qinflow(12,nrows))
    IF(.NOT.ALLOCATED(minbranchQ)) ALLOCATE(minbranchQ(12,nrows))
    qprod1date = 0; qprod2date = 0  
    wamp = missing_value
    qinflow = missing_value

    !Open DamData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    idtype = 0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        idtype = i                                        !! idtype is column where SUBID is
      ENDIF
      
      IF(str(i)(1:10)=='lake_depth') code(i) = i_real
      IF(str(i)(1:10)=='purpose   ') code(i) = i_intg   !Purpose of dam (1=irrigation,2=water supply, 3=Flood control, 4 =hydroelectricity)
      IF(str(i)(1:10)=='regvol    ') code(i) = i_real   !dam_regvol
      IF(str(i)(1:10)=='datum1    ') code(i) = i_intg   !date for two date based prod rate
      IF(str(i)(1:10)=='datum2    ') code(i) = i_intg   !date for two date based prod rate
      IF(str(i)(1:10)=='rate      ') code(i) = i_real   !dam_rate
      IF(str(i)(1:10)=='exp       ') code(i) = i_real   !dam_exp
      IF(str(i)(1:10)=='w0ref     ') code(i) = i_real   !dam_w0
      IF(str(i)(1:10)=='wamp      ') code(i) = i_real   !wamp (m)
      IF(str(i)(1:10)=='qprod1    ') code(i) = i_real   !two date-based production rates, qprod1 and qprod2, now replacing the unique rate qprod
      IF(str(i)(1:10)=='qprod2    ') code(i) = i_real
      IF(str(i)(1:10)=='qamp      ') code(i) = i_real   !dam_qamp
      IF(str(i)(1:10)=='qpha      ') code(i) = i_real   !dam_qpha
      IF(str(i)(1:10)=='qinfjan   ') code(i) = i_real   !dam_qinfjan !Mean monthly inflow to dam (simulated using natural run)
      IF(str(i)(1:10)=='qinffeb   ') code(i) = i_real   !dam_qinffeb
      IF(str(i)(1:10)=='qinfmar   ') code(i) = i_real   !dam_qinfmar
      IF(str(i)(1:10)=='qinfapr   ') code(i) = i_real   !dam_qinfapr
      IF(str(i)(1:10)=='qinfmay   ') code(i) = i_real   !dam_qinfmay
      IF(str(i)(1:10)=='qinfjun   ') code(i) = i_real   !dam_qinfjun
      IF(str(i)(1:10)=='qinfjul   ') code(i) = i_real   !dam_qinfjul
      IF(str(i)(1:10)=='qinfaug   ') code(i) = i_real   !dam_qinfaug
      IF(str(i)(1:10)=='qinfsep   ') code(i) = i_real   !dam_qinfsep
      IF(str(i)(1:10)=='qinfoct   ') code(i) = i_real   !dam_qinfoct
      IF(str(i)(1:10)=='qinfnov   ') code(i) = i_real   !dam_qinfnov
      IF(str(i)(1:10)=='qinfdec   ') code(i) = i_real   !dam_qinfdec
      IF(str(i)(1:10)=='snowfrac  ') code(i) = i_real   !dam_snowfrac 
     
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    !Allocate variables for dam information
    IF(.NOT.ALLOCATED(dam)) ALLOCATE(dam(nrows))
    ALLOCATE(damindex2(nrows))

    !Find subbasin-dam coupling: damindex(isb)
    IF(idtype==0)THEN
      WRITE(6,*) 'ERROR: subid not found in DamData. ',TRIM(filename)
      n_Result = 1
      RETURN
    ENDIF
    DO i = 1,n
      DO idam = 1,nrows
        IF(basin(i)%subid==xi(idam,iindex(idtype)))THEN
          damindex(i)=idam
          damindex2(idam)=i    !! find row in geodata for dam row
          EXIT
        ENDIF
      ENDDO
!!      IF(j==nrows+1)THEN      !!CD20141009 commented out
!!      WRITE(6,*) 'Error: A dam in DamData is missing in GeoData. subid: ',xi(j,iindex(idtype))
!!        STOP 1
!!      ENDIF  
    ENDDO

    !Save dam information
    DO idam = 1,nrows
      DO i = 1,ncols
        IF(str(i)=='regvol    ')   dam(idam)%regvol    = xr(idam,rindex(i))
        IF(str(i)=='rate      ')   dam(idam)%rate      = xr(idam,rindex(i))
        IF(str(i)=='exp       ')   dam(idam)%exp       = xr(idam,rindex(i))
        IF(str(i)=='w0ref     ')   dam(idam)%w0ref     = xr(idam,rindex(i))
        IF(str(i)=='wamp      ')   wamp(idam)          = xr(idam,rindex(i))
        IF(str(i)=='qprod1    ')   dam(idam)%qprod1    = xr(idam,rindex(i))
        IF(str(i)=='qprod2    ')   dam(idam)%qprod2    = xr(idam,rindex(i))
        IF(str(i)=='qamp      ')   dam(idam)%qamp      = xr(idam,rindex(i))
        IF(str(i)=='qpha      ')   dam(idam)%qpha      = xr(idam,rindex(i))
        IF(str(i)=='qinfjan   ')   qinflow(1,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinffeb   ')   qinflow(2,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfmar   ')   qinflow(3,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfapr   ')   qinflow(4,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfmay   ')   qinflow(5,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfjun   ')   qinflow(6,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfjul   ')   qinflow(7,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfaug   ')   qinflow(8,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfsep   ')   qinflow(9,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfoct   ')   qinflow(10,idam)    = xr(idam,rindex(i))
        IF(str(i)=='qinfnov   ')   qinflow(11,idam)    = xr(idam,rindex(i))
        IF(str(i)=='qinfdec   ')   qinflow(12,idam)    = xr(idam,rindex(i))
        IF(str(i)=='snowfrac  ')   dam(idam)%snowfrac  = xr(idam,rindex(i))
        IF(str(i)=='lake_depth')THEN
          IF(xr(idam,rindex(i)).NE.missing_value)   basin(damindex2(idam))%lakedepth(2) = xr(idam,rindex(i))
        ENDIF
      
        !Integer input
        IF(str(i)=='datum1    ')   qprod1date(idam)    = xi(idam,iindex(i))
        IF(str(i)=='datum2    ')   qprod2date(idam)    = xi(idam,iindex(i))
        IF(str(i)=='purpose   ')   dam(idam)%purpose   = xi(idam,iindex(i))
       
      ENDDO
    ENDDO
   
    !Calculate other lake/dam variables; w0ref, wmin, datum, wampcoeff, qinflow
    dam(:)%wmin = missing_value
    dam(:)%wampcoeff = missing_value
    DO i = 1,n
      IF(damindex(i)>0)THEN
        damarea = basin(i)%area * classbasin(i,slc_olake)%part
        IF(damarea>0)THEN
          dam(damindex(i))%wmin = 0. - dam(damindex(i))%regvol * 1000000. / damarea
          IF(wamp(damindex(i))/=missing_value) &
              dam(damindex(i))%wampcoeff = wamp(damindex(i))/(-1.*dam(damindex(i))%wmin)
        ENDIF
      ENDIF
    ENDDO
    DO idam = 1,nrows 
      MonthDayDate = qprod1date(idam)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      dam(idam)%datum1 = dayNumber
      MonthDayDate = qprod2date(idam)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      dam(idam)%datum2 = dayNumber
      dam(idam)%qinfmed = SUM(qinflow(:,idam))/12.
      dam(idam)%qinfmin = MINVAL(qinflow(:,idam))
      dam(idam)%qinfmax = MAXVAL(qinflow(:,idam))     
    ENDDO
 
    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(qprod1date)) DEALLOCATE(qprod1date) 
    IF(ALLOCATED(qprod2date)) DEALLOCATE(qprod2date)
    IF(ALLOCATED(wamp)) DEALLOCATE(wamp)
    IF(ALLOCATED(damindex2)) DEALLOCATE(damindex2)

    WRITE(6,*) 'File read: ', TRIM(filename)
    

  END SUBROUTINE load_damdata


  
  
   !------------------------------------------------------------------
  !>\brief MH2017 collects extra data about dams from ExtraData.txt
  !! modified version of load_extradata for style 2 dams, will replace load_extradata
  !! IMPORTANT: dams must be in same order in ExtraData.txt as DamData.txt
  !>\b Consequences: dam() may change
  !------------------------------------------------------------------
  SUBROUTINE load_extradata(funit,dir,n,n_Result,lakestate)
  
    USE WORLDVAR, ONLY : maxcharpath,  &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         dam,         &    !INOUT
                         damindex        !INOUT
    USE STATETYPE_MODULE, ONLY : lakestatetype

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                   !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER, INTENT(IN)  :: n                       !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: n_Result                !<Error status
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate !used to initialize conditioning
    
    !Local parameters
    INTEGER,PARAMETER :: maxcol = 200             !Maximum number of columns in ExtraData
    
    !Local variables
    LOGICAL fex
    INTEGER i,idam
    INTEGER j                    !MH2017: current row in ExtraData
                                      !IMPORTANT: dams must be in same order in DamData and ExtraData
    INTEGER counter,ii,jj             !used to fill monthly data array
    INTEGER idtype                    !Column with subid
    INTEGER nrows                     !Number of dams in file
    INTEGER ncols                     !Number of columns in file
    INTEGER code(maxcol)              !Code for column variable
    INTEGER rindex(maxcol)            !Index for column real variables
    INTEGER iindex(maxcol)            !Index for column integer variables
    INTEGER MonthDayDate              !Non-formated date from DamData.txt
    INTEGER dayNumber                 !Corresponding date in day number format
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10) str(maxcol)     !Content string
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL damarea                            ! area of dam reservoir (same as lakearea)


    n_Result = 0
    filename = TRIM(dir)//'ExtraData.txt'
    j = 0

    !Check if file exist and if outlet lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with ExtraData.txt'
      n_Result = 1
      RETURN
    ENDIF

    !Count number of data rows, and allocate variables for reading ExtraData
    CALL count_data_rows(funit,filename,1,nrows,n_Result)
    IF(n_Result/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No dam in ExtraData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))

   

    !Open ExtraData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    idtype = 0
    
    DO i = 1,ncols !don't read all days if possible
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        idtype = i                                        !! idtype is column where SUBID is
      ENDIF
      !special dam data
      !real data
      IF(str(i)(1:10)=='maxqmult  ') code(i) = i_real   !Coefficient to Scale Transition Band Type
      IF(str(i)(1:10)=='wslinit   ') code(i) = i_real   !Lake Water Surface Level used to Initialize Model
      IF(str(i)(1:10)=='stora     ') code(i) = i_real   !where: S = A * (WSL - B) ^ C + D
      IF(str(i)(1:10)=='storb     ') code(i) = i_real   !
      IF(str(i)(1:10)=='storc     ') code(i) = i_real   !
      IF(str(i)(1:10)=='stord     ') code(i) = i_real   !
      IF(str(i)(1:10)=='conda     ') code(i) = i_real   !where: Qo = A * (XX - B) ^ C + D
      IF(str(i)(1:10)=='condb     ') code(i) = i_real   !
      IF(str(i)(1:10)=='condc     ') code(i) = i_real   !
      IF(str(i)(1:10)=='condd     ') code(i) = i_real   !
      IF(str(i)(1:10)=='minqa     ') code(i) = i_real   !where: Qo = A * (WSL - B) ^ C + D
      IF(str(i)(1:10)=='minqb     ') code(i) = i_real   !
      IF(str(i)(1:10)=='minqc     ') code(i) = i_real   !
      IF(str(i)(1:10)=='minqd     ') code(i) = i_real   !
      IF(str(i)(1:10)=='midqa     ') code(i) = i_real   !
      IF(str(i)(1:10)=='midqb     ') code(i) = i_real   !
      IF(str(i)(1:10)=='midqc     ') code(i) = i_real   !
      IF(str(i)(1:10)=='midqd     ') code(i) = i_real   !
      IF(str(i)(1:10)=='maxqa     ') code(i) = i_real   !
      IF(str(i)(1:10)=='maxqb     ') code(i) = i_real   !
      IF(str(i)(1:10)=='maxqc     ') code(i) = i_real   !
      IF(str(i)(1:10)=='maxqd     ') code(i) = i_real   !
      IF(str(i)(1:10)=='natminqa  ') code(i) = i_real   !where: Qo = A * (WSL - B) ^ C + D
      IF(str(i)(1:10)=='natminqb  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natminqc  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natminqd  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmidqa  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmidqb  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmidqc  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmidqd  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmaxqa  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmaxqb  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmaxqc  ') code(i) = i_real   !
      IF(str(i)(1:10)=='natmaxqd  ') code(i) = i_real   !
      IF(str(i)(1:10)=='deltaw    ') code(i) = i_real   !Weekly Flow Change Limit
      IF(str(i)(1:10)=='deltad    ') code(i) = i_real   !Daily Flow Change Limit
      IF(str(i)(1:10)=='natdeltaw ') code(i) = i_real   !Weekly Flow Change Limit
      IF(str(i)(1:10)=='natdeltad ') code(i) = i_real   !Daily Flow Change Limit
      IF(str(i)(1:10)=='phqo      ') code(i) = i_real   !Powerhouse Limit Flow Value
      IF(str(i)(1:10)=='drqo      ') code(i) = i_real   !Minimum Fixed Flow Value
      IF(str(i)(1:10)=='flqo      ') code(i) = i_real   !Maximum Fixed Flow Value
      IF(str(i)(1:10)=='natphqo   ') code(i) = i_real   !Powerhouse Limit Flow Value
      IF(str(i)(1:10)=='natdrqo   ') code(i) = i_real   !Minimum Fixed Flow Value
      IF(str(i)(1:10)=='natflqo   ') code(i) = i_real   !Maximum Fixed Flow Value
      IF(str(i)(1:10)=='mina1     ') code(i:i+107) = i_real   !Used for monthly regression calculations, aka weighted flow, by column
                                                        !column 1:minA 2:midA 3:highA 4:natminA 5:natmidA 6: nathighA 7:minB 8:midB 9:highB
      
      !integer data
      IF(str(i)(1:10)=='identify  ') code(i) = i_intg   !Unique 4-digit code used to open DailyDataXXXX.txt files
      IF(str(i)(1:10)=='purpose2  ') code(i) = i_intg   !purpose of dam prior to being made a special dam;
                                                        !allows turning special dams 'off' so they function as normal dams
      IF(str(i)(1:10)=='natswitch ') code(i) = i_intg   !0 = No Natural Outflow, 1 = Has Natural Outflow
      IF(str(i)(1:10)=='locswitch ') code(i) = i_intg   !0 = No Local Inflow, 1 = Has Local Inflow 
      IF(str(i)(1:10)=='branchswt ') code(i) = i_intg   !0 = reservoir is branch, natural is main, 1 = reservoir is main, natural is branch
      IF(str(i)(1:10)=='condswitch') code(i) = i_intg   !0 = No Flow Condition, 1 = Has Conditioning 
      IF(str(i)(1:10)=='condrestr ') code(i) = i_intg   !0 = Low to High, 1 = Drought to Flood
      IF(str(i)(1:10)=='condqowsl ') code(i) = i_intg   !0 = Conditioned by Qo, 1 = Conditioned by WSL
      IF(str(i)(1:10)=='condsubid ') code(i) = i_intg   !HYPE Basin SubID Used for Conditioning Value                                                      
      IF(str(i)(1:10)=='drtype    ') code(i) = i_intg   !0 = Dr. Percentile, 1 = Fixed Min.
      IF(str(i)(1:10)=='mintype   ') code(i) = i_intg   !0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='midtype   ') code(i) = i_intg   !3 = Ideal Band Q, 0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='maxtype   ') code(i) = i_intg   !0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='fltype    ') code(i) = i_intg   !0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
      IF(str(i)(1:10)=='phtype    ') code(i) = i_intg   !0 = No PH Restriction, 1 = Hold Water
      IF(str(i)(1:10)=='natdrtype ') code(i) = i_intg   !0 = Dr. Percentile, 1 = Fixed Min.
      IF(str(i)(1:10)=='natmintype') code(i) = i_intg   !0 =  Min. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='natmidtype') code(i) = i_intg   !0 =   Avg. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='natmaxtype') code(i) = i_intg   !3 = Storage Change, 0 =  Max. Percentile, 1 = Monthly Regressed, 2 = A-B-C-D
      IF(str(i)(1:10)=='natfltype ') code(i) = i_intg   !0 = Fl. Percentile, 1 = Fixed Max., 2 = Extend Max. Type
      IF(str(i)(1:10)=='natphtype ') code(i) = i_intg   !0 = No PH Restriction, 1 = Hold Water
      
      
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)
 
    !dam should already be allocated
  
    IF(.NOT.ALLOCATED(lakestate%useswitch)) ALLOCATE(lakestate%useswitch(n))
    lakestate%useswitch = 0
    !cycle through each dam 
    damCycle: DO idam = 1,SIZE(dam)
      !assign extra variables
      specialdam: IF (dam(idam)%purpose >= 5) THEN
        j = j + 1
        colCycle: DO i = 1,ncols
          !real data
          IF(str(i)=='maxqmult  ') dam(idam)%maxqmult    = xr(j,rindex(i))
          IF(str(i)=='wslinit   ') dam(idam)%wslinit     = xr(j,rindex(i))
          IF(str(i)=='stora     ') dam(idam)%stora       = xr(j,rindex(i))
          IF(str(i)=='storb     ') dam(idam)%storb       = xr(j,rindex(i))
          IF(str(i)=='storc     ') dam(idam)%storc       = xr(j,rindex(i))
          IF(str(i)=='stord     ') dam(idam)%stord       = xr(j,rindex(i))
          IF(str(i)=='conda     ') dam(idam)%conda       = xr(j,rindex(i))
          IF(str(i)=='condb     ') dam(idam)%condb       = xr(j,rindex(i))
          IF(str(i)=='condc     ') dam(idam)%condc       = xr(j,rindex(i))
          IF(str(i)=='condd     ') dam(idam)%condd       = xr(j,rindex(i))
          IF(str(i)=='minqa     ') dam(idam)%minqa       = xr(j,rindex(i))
          IF(str(i)=='minqb     ') dam(idam)%minqb       = xr(j,rindex(i))
          IF(str(i)=='minqc     ') dam(idam)%minqc       = xr(j,rindex(i))
          IF(str(i)=='minqd     ') dam(idam)%minqd       = xr(j,rindex(i))
          IF(str(i)=='midqa     ') dam(idam)%midqa       = xr(j,rindex(i))
          IF(str(i)=='midqb     ') dam(idam)%midqb       = xr(j,rindex(i))
          IF(str(i)=='midqc     ') dam(idam)%midqc       = xr(j,rindex(i))
          IF(str(i)=='midqd     ') dam(idam)%midqd       = xr(j,rindex(i))
          IF(str(i)=='maxqa     ') dam(idam)%maxqa       = xr(j,rindex(i))
          IF(str(i)=='maxqb     ') dam(idam)%maxqb       = xr(j,rindex(i))
          IF(str(i)=='maxqc     ') dam(idam)%maxqc       = xr(j,rindex(i))
          IF(str(i)=='maxqd     ') dam(idam)%maxqd       = xr(j,rindex(i))
          IF(str(i)=='natminqa  ') dam(idam)%natminqa    = xr(j,rindex(i))
          IF(str(i)=='natminqb  ') dam(idam)%natminqb    = xr(j,rindex(i))
          IF(str(i)=='natminqc  ') dam(idam)%natminqc    = xr(j,rindex(i))
          IF(str(i)=='natminqd  ') dam(idam)%natminqd    = xr(j,rindex(i))
          IF(str(i)=='natmidqa  ') dam(idam)%natmidqa    = xr(j,rindex(i))
          IF(str(i)=='natmidqb  ') dam(idam)%natmidqb    = xr(j,rindex(i))
          IF(str(i)=='natmidqc  ') dam(idam)%natmidqc    = xr(j,rindex(i))
          IF(str(i)=='natmidqd  ') dam(idam)%natmidqd    = xr(j,rindex(i))
          IF(str(i)=='natmaxqa  ') dam(idam)%natmaxqa    = xr(j,rindex(i))
          IF(str(i)=='natmaxqb  ') dam(idam)%natmaxqb    = xr(j,rindex(i))
          IF(str(i)=='natmaxqc  ') dam(idam)%natmaxqc    = xr(j,rindex(i))
          IF(str(i)=='natmaxqd  ') dam(idam)%natmaxqd    = xr(j,rindex(i))
          IF(str(i)=='deltaw    ') dam(idam)%deltaw      = xr(j,rindex(i))
          IF(str(i)=='deltad    ') dam(idam)%deltad      = xr(j,rindex(i))
          IF(str(i)=='natdeltaw ') dam(idam)%natdeltaw   = xr(j,rindex(i))
          IF(str(i)=='natdeltad ') dam(idam)%natdeltad   = xr(j,rindex(i))
          IF(str(i)=='phqo      ') dam(idam)%phqo        = xr(j,rindex(i))
          IF(str(i)=='drqo      ') dam(idam)%drqo        = xr(j,rindex(i))
          IF(str(i)=='flqo      ') dam(idam)%flqo        = xr(j,rindex(i))
          IF(str(i)=='natphqo   ') dam(idam)%natphqo     = xr(j,rindex(i))
          IF(str(i)=='natdrqo   ') dam(idam)%natdrqo     = xr(j,rindex(i))
          IF(str(i)=='natflqo   ') dam(idam)%natflqo     = xr(j,rindex(i))
          IF(str(i)=='mina1     ') THEN
            IF(.NOT.ALLOCATED(dam(idam)%monthly)) ALLOCATE(dam(idam)%monthly(12,9))
            counter = i
            DO ii=1,12
              DO jj=1,9
                dam(idam)%monthly(ii,jj) = xr(j,rindex(counter))
                counter = counter + 1
              ENDDO
            ENDDO
          ENDIF
          
          !integer data
          IF(str(i)=='identify  ') dam(idam)%identify    = xi(j,iindex(i))
          IF(str(i)=='purpose2  ') dam(idam)%purpose2    = xi(j,iindex(i))
          IF(str(i)=='natswitch ') dam(idam)%natswitch   = xi(j,iindex(i))
          IF(str(i)=='locswitch ') dam(idam)%locswitch   = xi(j,iindex(i))
          IF(str(i)=='branchswt ') dam(idam)%branchswt   = xi(j,iindex(i))
          IF(str(i)=='condswitch') dam(idam)%condswitch  = xi(j,iindex(i))
          IF(str(i)=='condrestr ') dam(idam)%condrestr   = xi(j,iindex(i))
          IF(str(i)=='condqowsl ') dam(idam)%condqowsl   = xi(j,iindex(i))
          IF(str(i)=='condsubid ')THEN
            dam(idam)%condsubid   = xi(j,iindex(i))
            IF(dam(idam)%condsubid>0)THEN !conditioning required
              dam(idam)%condindex   = get_subid_index(dam(idam)%condsubid,n)    
              lakestate%useswitch(dam(idam)%condindex) = 1
              IF(dam(idam)%condindex < damindex(idam)) dam(idam)%orderswitch = 1                            
            ENDIF
          ENDIF 
          IF(str(i)=='drtype    ') dam(idam)%drtype      = xi(j,iindex(i))
          IF(str(i)=='mintype   ') dam(idam)%mintype     = xi(j,iindex(i))
          IF(str(i)=='midtype   ') dam(idam)%midtype     = xi(j,iindex(i))
          IF(str(i)=='maxtype   ') dam(idam)%maxtype     = xi(j,iindex(i))
          IF(str(i)=='fltype    ') dam(idam)%fltype      = xi(j,iindex(i))
          IF(str(i)=='phtype    ') dam(idam)%phtype      = xi(j,iindex(i))
          IF(str(i)=='natdrtype ') dam(idam)%natdrtype   = xi(j,iindex(i))
          IF(str(i)=='natmintype') dam(idam)%natmintype  = xi(j,iindex(i))
          IF(str(i)=='natmidtype') dam(idam)%natmidtype  = xi(j,iindex(i))
          IF(str(i)=='natmaxtype') dam(idam)%natmaxtype  = xi(j,iindex(i))
          IF(str(i)=='natfltype ') dam(idam)%natfltype   = xi(j,iindex(i))
          IF(str(i)=='natphtype ') dam(idam)%natphtype   = xi(j,iindex(i))
        ENDDO colCycle
        IF(.NOT.ALLOCATED(dam(idam)%qin14day)) ALLOCATE(dam(idam)%qin14day(14))
        IF(.NOT.ALLOCATED(dam(idam)%sabd14day)) ALLOCATE(dam(idam)%sabd14day(14))
        IF(.NOT.ALLOCATED(dam(idam)%weekly7day)) ALLOCATE(dam(idam)%weekly7day(7))
        IF(.NOT.ALLOCATED(dam(idam)%resqo7day)) ALLOCATE(dam(idam)%resqo7day(7))
        IF(.NOT.ALLOCATED(dam(idam)%qin7day)) ALLOCATE(dam(idam)%qin7day(7))
        IF(.NOT.ALLOCATED(dam(idam)%natweekly7day)) ALLOCATE(dam(idam)%natweekly7day(7))
        IF(.NOT.ALLOCATED(lakestate%wsl14day)) ALLOCATE(lakestate%wsl14day(n,15))
        IF(.NOT.ALLOCATED(lakestate%qo14day)) ALLOCATE(lakestate%qo14day(n,15))
        IF(.NOT.ALLOCATED(dam(idam)%wsl14day)) ALLOCATE(dam(idam)%wsl14day(15))
        IF(.NOT.ALLOCATED(dam(idam)%qo14day)) ALLOCATE(dam(idam)%qo14day(15))
        !size 15 is intentional             
        dam(idam)%qin14day = 0.
        dam(idam)%sabd14day = 0.
        dam(idam)%weekly7day = 0.
        dam(idam)%resqo7day = 0.
        dam(idam)%qin7day = 0.
        dam(idam)%natweekly7day = 0.
        lakestate%wsl14day = 0.
        lakestate%qo14day = 0.        
        dam(idam)%wsl14day = 0.
        dam(idam)%qo14day = 0.
        !turns a special dam 'off' so it functions as a normal dam
        IF(dam(idam)%identify==9999) dam(idam)%purpose = dam(idam)%purpose2
      ENDIF specialdam
    ENDDO damCycle
  
    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)


    WRITE(6,*) 'File read: ', TRIM(filename)
  
  
  END SUBROUTINE load_extradata
  
     !------------------------------------------------------------------
  !>\brief MH2017 collects daily data for special dams from DailyDataXXXX.txt
  !>\b 'XXXX' is a unique 4-digit code read from column 'IDENTIFY' in ExtraData.txt
  !>\b in order to open file 'DailyData[code].txt'
  !>\b Consequences: dam()%daily will be allocated and change
  !------------------------------------------------------------------
  SUBROUTINE load_dailydataXXXX(funit,dir,n,n_Result,idam)
  
    USE WORLDVAR, ONLY : maxcharpath,  &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         dam,         &    !INOUT
                         damindex        !INOUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                   !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER, INTENT(IN)  :: n                       !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: n_Result                !<Error status
    INTEGER, INTENT(IN)  :: idam                    !<index of current dam
    
    !Local parameters
    INTEGER,PARAMETER :: maxcol = 30             !Maximum number of columns in ExtraData
    
    !Local variables
    LOGICAL fex
    INTEGER i,j                       !current column, row (opposite standard nomenclature)           
    INTEGER nrows                     !Number of dams in file
    INTEGER ncols                     !Number of columns in file
    INTEGER code(maxcol)              !Code for column variable
    INTEGER rindex(maxcol)            !Index for column real variables
    INTEGER iindex(maxcol)            !Index for column integer variables
    INTEGER dayNumber                 !Corresponding date in day number format
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10) str(maxcol)     !Content string
    CHARACTER(LEN=4) identify             !4-digit code used to open DailyData[identify].txt
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
   

    WRITE(identify,'(i4)') dam(idam)%identify
    n_Result = 0
    filename = TRIM(dir)//'DailyData'//identify//'.txt'
    j = 0
    
    !WRITE(6,*) 'Attempting to open ',filename !uncomment to check for dailydata/extradata mismatch

    !Check if file exist and if outlet lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with ',filename
      n_Result = 1
      RETURN
    ENDIF

    !Count number of data rows, and allocate variables for reading DailyData
    CALL count_data_rows(funit,filename,1,nrows,n_Result)
    IF(n_Result/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No days in DailyData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))

   

    !Open DailyData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    
    DO i = 1,ncols 
      !daily special dam data
      !real data
      IF(str(i)(1:10)=='drlevel   ') code(i) = i_real   !column 1: drought WSL
      IF(str(i)(1:10)=='lowops    ') code(i) = i_real   !column 2: lower operations WSL
      IF(str(i)(1:10)=='lowband   ') code(i) = i_real   !column 3: lower band WSL
      IF(str(i)(1:10)=='avglevel  ') code(i) = i_real   !column 4: average WSL
      IF(str(i)(1:10)=='upband    ') code(i) = i_real   !column 5: upper band WSL
      IF(str(i)(1:10)=='upops     ') code(i) = i_real   !column 6: upper operations WSL
      IF(str(i)(1:10)=='fllevel   ') code(i) = i_real   !column 7: flood WSL
      IF(str(i)(1:10)=='drprcntl  ') code(i) = i_real   !column 8: smoothed drought percentile flow
      IF(str(i)(1:10)=='minprcntl ') code(i) = i_real   !column 9: ...minimum...
      IF(str(i)(1:10)=='avgprcntl ') code(i) = i_real   !column 10: ...average...
      IF(str(i)(1:10)=='maxprcntl ') code(i) = i_real   !column 11: ...maximum...
      IF(str(i)(1:10)=='flprcntl  ') code(i) = i_real   !column 12: ...flood...
      IF(str(i)(1:10)=='natdrprc  ') code(i) = i_real   !column 13: natural drought percentile flow
      IF(str(i)(1:10)=='natminprc ') code(i) = i_real   !column 14: ...minimum...
      IF(str(i)(1:10)=='natavgprc ') code(i) = i_real   !column 15: ...average...
      IF(str(i)(1:10)=='natmaxprc ') code(i) = i_real   !column 16: ...maximum...
      IF(str(i)(1:10)=='natflprc  ') code(i) = i_real   !column 17: ...flood...
      IF(str(i)(1:10)=='mweight   ') code(i) = i_real   !column 18: monthly weighting
    ENDDO
     
    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)
 
    !cannot use (nrows,ncols) because only some columns have data
    IF(.NOT.ALLOCATED(dam(idam)%daily)) ALLOCATE(dam(idam)%daily(nrows,18))
    
    !fill dam(idam)%daily with data from file
    rowCycle: DO j = 1,nrows
      colCycle: DO i = 1,ncols
        !real data
        IF(str(i)=='drlevel   ') dam(idam)%daily(j,1)  = xr(j,rindex(i))
        IF(str(i)=='lowops    ') dam(idam)%daily(j,2)  = xr(j,rindex(i))
        IF(str(i)=='lowband   ') dam(idam)%daily(j,3)  = xr(j,rindex(i))
        IF(str(i)=='avglevel  ') dam(idam)%daily(j,4)  = xr(j,rindex(i))
        IF(str(i)=='upband    ') dam(idam)%daily(j,5)  = xr(j,rindex(i))
        IF(str(i)=='upops     ') dam(idam)%daily(j,6)  = xr(j,rindex(i))
        IF(str(i)=='fllevel   ') dam(idam)%daily(j,7)  = xr(j,rindex(i))
        IF(str(i)=='drprcntl  ') dam(idam)%daily(j,8)  = xr(j,rindex(i))
        IF(str(i)=='minprcntl ') dam(idam)%daily(j,9)  = xr(j,rindex(i))
        IF(str(i)=='avgprcntl ') dam(idam)%daily(j,10) = xr(j,rindex(i))
        IF(str(i)=='maxprcntl ') dam(idam)%daily(j,11) = xr(j,rindex(i))
        IF(str(i)=='flprcntl  ') dam(idam)%daily(j,12) = xr(j,rindex(i))
        IF(str(i)=='natdrprc  ') dam(idam)%daily(j,13) = xr(j,rindex(i))
        IF(str(i)=='natminprc ') dam(idam)%daily(j,14) = xr(j,rindex(i))
        IF(str(i)=='natavgprc ') dam(idam)%daily(j,15) = xr(j,rindex(i))
        IF(str(i)=='natmaxprc ') dam(idam)%daily(j,16) = xr(j,rindex(i))
        IF(str(i)=='natflprc  ') dam(idam)%daily(j,17) = xr(j,rindex(i))
        IF(str(i)=='mweight   ') dam(idam)%daily(j,18) = xr(j,rindex(i))
      ENDDO colCycle
    ENDDO rowCycle
  
    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)


    WRITE(6,*) 'File read: ', TRIM(filename)
  
  
  END SUBROUTINE load_dailydataXXXX
  
  !------------------------------------------------------------------
  !>\brief Set lakedatapar values to those model variable values specified in LakeData.txt. 
  !!
  !>\b Consequences Module modvar variables lakedatapar, lakedataparindex may be allocated and set. 
  !------------------------------------------------------------------
  SUBROUTINE start_lakedata_table(maxcol,ncols,nrows,str,xi,xr,iindex,rindex,ns,lakedataid)

    USE MODVAR, ONLY : modparid,          &
                       lakedatapar,       &   !OUT
                       lakedataparindex,  &   !OUT
                       max_par,           &
                       m_ldpar,           &
                       m_lrpar,           &
                       nlakeregions,      &
                       missing_value,     &
                       basin

    !Argument declaration
    INTEGER, INTENT(IN)  :: maxcol                    !<Maximum possible amount of columns on LakeData.txt
    INTEGER, INTENT(IN)  :: ncols, nrows              !<Amount of columns and rows effectively used in LakeData.txt
    CHARACTER(LEN=10), INTENT(IN) :: str(maxcol)      !<Character array with labels of ALL columns given in LakeData.txt (column headers)
    INTEGER, INTENT(IN)  :: xi(nrows,maxcol)          !<Integer values aquired by reading LakeData.txt
    REAL, INTENT(IN)     :: xr(nrows,maxcol)          !<Float          "                 "
    INTEGER, INTENT(IN)  :: iindex(maxcol)            !<Index for column real variables
    INTEGER, INTENT(IN)  :: rindex(maxcol)            !<Index for column integer variables
    INTEGER, INTENT(IN)  :: ns                        !<Number of subbasins (base model)
    INTEGER, INTENT(IN)  :: lakedataid(ns)            !<lakedataid from GeoData.txt
    
    !Local variables
    INTEGER rowCounter, colCounter, modparCounter, basinCounter
    INTEGER nlakedatapar, tempVal1

    nlakedatapar = 0          !Amount of lakedatapar model variables; set to 0, then find the largest one
    DO modparCounter = 1,max_par
      IF(modparid(modparCounter)%deptype==m_ldpar)THEN
        nlakedatapar = MAX(nlakedatapar,modparid(modparCounter)%parno)
      ENDIF
    ENDDO

    IF(.NOT.ALLOCATED(lakedatapar)) ALLOCATE(lakedatapar(nlakeregions+nrows, nlakedatapar)) !This is more rows than needed
    lakedatapar = missing_value       !Initialize the table of parameter values with the missing value (-9999)
    IF(.NOT.ALLOCATED(lakedataparindex)) ALLOCATE(lakedataparindex(ns,2))
    lakedataparindex = 0              !Initialize the line index table with 0 (useful for 2nd colum, with the row number)

    !Parameter by parameter (= column by column), copy those parameter values given in LakeData.txt in the table "lakedatapar" (it's a modvar)
    !Values are put at the bottom of the table. Eventual parameter values coming from par.txt are appended later in the table top rows, in another subroutine.
    !The order in which lakes are given in LakeData.txt is respected linewise
    DO colCounter = 1,ncols

       !Loop through all model variables to look for label match with the currently considered parameter, labelled str(colCounter) in lakeData
       DO modparCounter = 1,max_par
          !When a matching variable shortname is found, determine if the model variable has the parameter value read from lakedata.txt (deptype = m_ldpar), or that from par.txt (deptype = others)
          IF(str(colCounter)==modparid(modparCounter)%shortname .AND. modparid(modparCounter)%deptype==m_ldpar)THEN

             !Parameter values coming from lakedata are stored in either xi for integer, or xr for float values; take them from there
             ![Note: discriminating between integer and float is likely completely useless since the lakedata parameter values are all floats, I think. No?]
             IF(iindex(colCounter).NE.0)THEN
                DO rowCounter = 1,nrows
                   lakedatapar(nlakeregions+rowCounter, modparid(modparCounter)%parno) = xi(rowCounter, iindex(colCounter))
                ENDDO
             ELSEIF(rindex(colCounter).NE.0)THEN
                DO rowCounter = 1,nrows
                   lakedatapar(nlakeregions+rowCounter, modparid(modparCounter)%parno) = xr(rowCounter, rindex(colCounter))
                ENDDO
             ENDIF   !End of integer vs float discrimination test
             EXIT    !CP110523  no need to loop further when found
          ENDIF   !End of the check on label match and on the origin of the parameter value (lakeData.txt or par.txt)
       ENDDO   !End of model variable name loop

       IF(str(colCounter)=='lakedataid')THEN
         DO rowCounter = 1,nrows
           tempVal1 = xi(rowCounter, iindex(colCounter))
           IF(tempVal1>0)THEN
             DO basinCounter = 1,ns
               IF(lakedataid(basinCounter) == tempVal1)THEN    
!                 tempVal2 = basinCounter
                 lakedataparindex(basinCounter, 2) = nlakeregions+rowCounter  !CP130524
                 EXIT
               ENDIF
             ENDDO
!             lakedataparindex(tempVal2, 2) = nlakeregions+rowCounter   !Error when tempVal1 not found, then old tempVal2 reset!
           ENDIF
         ENDDO
       ENDIF

    ENDDO   !End of lakeData column labels loop

    !Point ilake and non-LakeData lakes to lakeregion-row in lakedatapar
    DO basinCounter = 1,ns
       lakedataparindex(basinCounter,1) = basin(basinCounter)%lakeregion
       IF(lakedataparindex(basinCounter,2) == 0)THEN
          lakedataparindex(basinCounter,2) = basin(basinCounter)%lakeregion
       ENDIF
    ENDDO

  END SUBROUTINE start_lakedata_table

  !>Set lakedatapar values to those model variable values specified
  !>in par.txt. FC, 20.03.11, changed CP110525
  !------------------------------------------------------------------
  SUBROUTINE finish_lakedata_table(lakedatafile,ns)

    USE MODVAR, ONLY : modparid,    &
                       lakedatapar, &
                       lakedataparindex,  &
                       basin,       &  
                       genpar,      &
                       lregpar,     &    
                       max_par,     &
                       m_gpar,      &
                       m_lrpar,     &
                       m_ldpar,     &
                       nlakeregions,  &  
                       missing_value

    !Argument declaration
    LOGICAL, INTENT(IN) :: lakedatafile   !<status of LakeData.txt
    INTEGER, INTENT(IN) :: ns             !<number of subbasins (nsub_basemodel)
    
    !Local variables
    LOGICAL found         !corresponding parameter found
    INTEGER parCounter    !index for loop to find corresponding par
    INTEGER ldparCounter  !index for loop to find ldpar
    INTEGER lakeCounter   !lake/row in lakedata
    INTEGER isb           !subbasin index/number
    INTEGER dimlakepar    !size of lakedatapar
    INTEGER dimlakepar2   !second dimension of lakedatapar
    INTEGER,ALLOCATABLE :: subbasinindex(:) !subbasin index for lakedata (rows)

    !Start of subroutine
    IF(lakedatafile)THEN

       !Set index table for finding subbasin index for lakedata (row)
       dimlakepar = SIZE(lakedatapar,1)
       IF(.NOT.ALLOCATED(subbasinindex)) ALLOCATE(subbasinindex(dimlakepar))
       subbasinindex = 0
       DO lakeCounter = nlakeregions+1,dimlakepar
         DO isb = 1,ns
           IF(lakedataparindex(isb,2)==lakeCounter)THEN
             EXIT
           ENDIF
         ENDDO
         IF(isb<=ns) subbasinindex(lakeCounter) = isb
       ENDDO

       DO ldparCounter = 1,max_par
          !For every lakedata-parameter find the corresponding parameter in par.txt
          IF(modparid(ldparCounter)%deptype==m_ldpar)THEN

             !Loop through all model parameters to look for label match with the 
             !currently considered parameter, labelled str(colCounter) in LakeData.txt
             found = .FALSE.
             DO parCounter = 1,max_par
                IF(modparid(parCounter)%deptype==m_gpar .OR. modparid(parCounter)%deptype==m_lrpar)THEN
                   IF(modparid(parCounter)%shortname==modparid(ldparCounter)%shortname)THEN  !Check label match with model parameter of ldpar-type
                      !Corresponding parameter found!
                      found = .TRUE.
                      !Set general/lakeregion parameter values in lakedatapar
                      IF(modparid(parCounter)%deptype==m_gpar)THEN 
                         lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                      ELSEIF(modparid(parCounter)%deptype==m_lrpar)THEN
                         lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = lregpar(modparid(parCounter)%parno,1:nlakeregions)
                      ENDIF

                      !Set missing parameters from LakeData file to general/lakeregion parameter value
                      DO lakeCounter = nlakeregions+1,dimlakepar
                         IF(lakedatapar(lakeCounter, modparid(ldparCounter)%parno)==missing_value)THEN
                            IF(modparid(parCounter)%deptype==m_gpar)THEN
                               lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                            ENDIF
                            IF(modparid(parCounter)%deptype==m_lrpar)THEN
                               IF(subbasinindex(lakeCounter)>0)THEN
                                  lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = lregpar(modparid(parCounter)%parno,basin(subbasinindex(lakeCounter))%lakeregion)
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDDO

                      !Next ldpar
                      EXIT
                   ENDIF !shortname is equal
                ENDIF !gpar or lrpar
             ENDDO !parCounter

             IF(.NOT.found)THEN  !This variable does not have a corresponing par.txt parameter !CP110526
                !Set non-LakeData set lake to zero
                lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = 0.
                DO lakeCounter = nlakeregions+1,dimlakepar
                   IF(lakedatapar(lakeCounter, modparid(ldparCounter)%parno)==missing_value)THEN
                      lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = 0.
                   ENDIF
                ENDDO
             ENDIF

          ENDIF !ldpar
       ENDDO !ldparCounter

       IF(ALLOCATED(subbasinindex)) DEALLOCATE(subbasinindex)

    ELSE
       !No LakeData.txt: set lakedatapar with only the general/lakeregion parameters 
       dimlakepar = nlakeregions
       dimlakepar2 = 0          !Amount of lakedatapar model variables; set to 0, then find the largest one
       DO parCounter = 1,max_par
         IF(modparid(parCounter)%deptype==m_ldpar)THEN
           dimlakepar2 = MAX(dimlakepar2,modparid(parCounter)%parno)
         ENDIF
       ENDDO
       IF(.NOT.ALLOCATED(lakedatapar)) ALLOCATE(lakedatapar(nlakeregions,dimlakepar2))
       IF(.NOT.ALLOCATED(lakedataparindex)) ALLOCATE(lakedataparindex(ns,2))

       !Point ilake and non-LakeData olakes to lakeregion-row in lakedatapar
       lakedataparindex(:,1) = basin(:)%lakeregion
       lakedataparindex(:,2) = lakedataparindex(:,1)

       DO ldparCounter = 1,max_par
          !For every lakedata-parameter find the corresponding parameter in par.txt
          IF(modparid(ldparCounter)%deptype==m_ldpar)THEN

             !Loop through all model parameters to look for label match with the currently considered parameter, labelled str(colCounter) in LakeData.txt
             DO parCounter = 1,max_par
                IF(modparid(parCounter)%deptype==m_gpar .OR. modparid(parCounter)%deptype==m_lrpar)THEN
                   IF(modparid(parCounter)%shortname==modparid(ldparCounter)%shortname)THEN  !Check label match with model parameter of ldpar-type
                      !Corresponding parameter found!

                      !Set general/lakeregion parameter values in lakedatapar
                      IF(modparid(parCounter)%deptype==m_gpar)THEN 
                         lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                      ELSEIF(modparid(parCounter)%deptype==m_lrpar)THEN
                         lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = lregpar(modparid(parCounter)%parno,1:nlakeregions)
                      ENDIF

                      !Next ldpar
                      EXIT
                   ENDIF !shortname is equal
                ENDIF !gpar or lrpar
             ENDDO !parCounter
          ENDIF !ldpar
       ENDDO !ldparCounter
    ENDIF

  END SUBROUTINE finish_lakedata_table

  !>\brief Gets the information about subbasins from GeoData.
  !!Reads the matrix of basin data values from the file vith mcols columns
  !!
  !>\b Consequences Module variables will be allocated and set
  !------------------------------------------------------------------------------
  SUBROUTINE read_and_calc_basindata(funit,infile,n,lakedataid,maxcol,mcols) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         pobsid,      &   !OUT
                         tobsid,      &   !OUT
                         sfobsid,     &   !OUT (ev.)
                         swobsid,     &   !OUT (ev.)
                         uobsid,      &   !OUT (ev.)
                         rhobsid,     &   !OUT (ev.)
                         tminobsid,     &   !OUT (ev.)
                         tmaxobsid,     &   !OUT (ev.)
                         readsfobs,  &
                         readswobs, &
                         readwind,  &
                         readhumid,  &
                         readtminmaxobs,  &
                         readobsid            
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp,  &
                       max_pstype, &
                       basin,        &    !OUT
                       classbasin,   &    !OUT
                       pathsubid,    &    !OUT
                       load,         &    !OUT
                       wetland,      &    !OUT
                       nlakeregions, &    !OUT  !Number of lakeregions
                       nregions,     &    !OUT  !Number of parameter regions
                       nwqregions,   &    !OUT  !Number of water quality parameter regions
                       nilakeregions,&    !OUT  !Number of ilake parameter regions
                       nolakeregions      !OUT  !Number of olake parameter regions

    !Argumetn declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile        !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: n                      !<Number of subbasins
    INTEGER, INTENT(OUT) :: lakedataid(n)          !<lakedataid from GeoData
    INTEGER, INTENT(IN)  :: maxcol                 !<Maximum number of data columns
    INTEGER, INTENT(OUT) :: mcols                  !<Actual number of columns
    
    !Local variables
    !MH2017: making arrays with (n) allocatable to overcome stack overflow
    INTEGER i,k
    INTEGER idslc
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER code(maxcol)               !Code for column variable
    INTEGER rindex(maxcol)             !Index for column real variables
    INTEGER iindex(maxcol)             !Index for column integer variables
    INTEGER status                     !Error status
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL, ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL, ALLOCATABLE :: pointsource(:,:,:) !Read pointsource info (subbasins,column,pstype)
    CHARACTER(LEN=10) str(maxcol)      !Content string

    !Start of subroutine
    status = 0
    onetspday = 1./REAL(timesteps_per_day)
    pointsource = 0.
    
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(n,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(n,maxcol))
    IF(.NOT.ALLOCATED(pointsource)) ALLOCATE(pointsource(n,5,max_pstype))
    
    
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     !Open GeoData-file

    !Reads the column headings from file
    CALL read_column_headings(funit,maxcol,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:10)=='area      ') code(i) = i_real
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='xcoord    ') code(i) = i_real
      IF(str(i)(1:10)=='ycoord    ') code(i) = i_real
      IF(str(i)(1:10)=='longitude ') code(i) = i_real
      IF(str(i)(1:10)=='latitude  ') code(i) = i_real
      IF(str(i)(1:10)=='elev_mean ') code(i) = i_real
      IF(str(i)(1:10)=='elev_std  ') code(i) = i_real
      IF(str(i)(1:10)=='slope_mean') code(i) = i_real
      IF(str(i)(1:10)=='slope_std ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='loc_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_in    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='region    ') code(i) = i_intg
      IF(str(i)(1:10)=='lakeregion') code(i) = i_intg
      IF(str(i)(1:10)=='lake_depth') code(i) = i_real
      IF(str(i)(1:10)=='lake_whigh'.OR.str(i)(1:10)=='lake_qavg '.OR.   &
         str(i)(1:10)=='lake_qamp '.OR.str(i)(1:10)=='lake_rate '.OR.   &
         str(i)(1:10)=='lake_qhigh'.OR.str(i)(1:10)=='lake_exp  '.OR.   &
         str(i)(1:10)=='mq        '.OR.str(i)(1:10)=='regvol    '.OR.   &
         str(i)(1:10)=='lake_qpha '.OR.str(i)(1:10)=='lake_wref ')THEN
        WRITE(6,*) 'WARNING: specific lake parameter no longer read from GeoData'
        WRITE(6,*)
      ENDIF
      IF(str(i)(1:10)=='icatch    ') code(i) = i_real
      IF(str(i)(1:10)=='mainfl    ') THEN
        WRITE(6,*) 'WARNING: mainfl in GeoData is an column name no longer used' 
        WRITE(6,*)
      ENDIF
      IF(str(i)(1:10)=='branch    ') WRITE(6,*) 'WARNING: branch in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='grwflow1  ') WRITE(6,*) 'WARNING: grwflow1 in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='branchdown') WRITE(6,*) 'WARNING: branchdown in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='maindown  ') code(i) = i_intg
      IF(str(i)(1:10)=='grwdown   ') code(i) = i_intg
      IF(str(i)(1:10)=='grwolake  ') code(i) = i_real
      IF(str(i)(1:3)=='slc')         code(i) = i_real
      IF(str(i)(1:5)=='dhslc')       code(i) = i_real
      IF(str(i)(1:3)=='scr')         code(i) = i_real
      IF(str(i)(1:10)=='parreg    ') code(i) = i_intg
      IF(str(i)(1:10)=='wqparreg  ') code(i) = i_intg
      IF(str(i)(1:10)=='rivlen    ') code(i) = i_real
      IF(str(i)(1:10)=='wetdep_n  ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n1 ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n2 ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n3 ') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_area') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_area') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_dep ') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_dep ') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_part') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_part') code(i) = i_real
      IF(str(i)(1:10)=='buffer    ') code(i) = i_real
      IF(str(i)(1:10)=='close_w   ') code(i) = i_real
      IF(str(i)(1:10)=='lakedataid') code(i) = i_intg   !id i LakeData
      IF(str(i)(1:10)=='pobsid    ') code(i) = i_intg   !id i Pobs
      IF(str(i)(1:10)=='tobsid    ') code(i) = i_intg   !id i Tobs
      IF(str(i)(1:10)=='sfobsid   ') code(i) = i_intg   !id in SFobs
      IF(str(i)(1:10)=='swobsid   ') code(i) = i_intg   !id in SWobs
      IF(str(i)(1:10)=='uobsid    ') code(i) = i_intg   !id in Uobs
      IF(str(i)(1:10)=='rhobsid   ') code(i) = i_intg   !id in RHobs
      IF(str(i)(1:10)=='tminobsid ') code(i) = i_intg   !id in TMINobs
      IF(str(i)(1:10)=='tmaxobsid ') code(i) = i_intg   !id in TNAXons
      IF(str(i)(1:10)=='ilakereg  ') code(i) = i_intg
      IF(str(i)(1:10)=='olakereg  ') code(i) = i_intg
      IF(str(i)(1:10)=='noncontr') code(i) = i_real     !MM2016: non-contributing fraction
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,maxcol,n,mcols,code,rindex,iindex,xi,xr)

    !600   FORMAT(A18000)
    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    DO i = 1,mcols
      IF(str(i)(1:10)=='area      ')   basin(1:n)%area        = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='subid     ')   basin(1:n)%subid       = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='xcoord    ')   basin(1:n)%xcoord      = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ycoord    ')   basin(1:n)%ycoord      = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='longitude ')   basin(1:n)%longitude   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='latitude  ')   basin(1:n)%latitude    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='elev_mean ')   basin(1:n)%elev        = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='elev_std  ')   basin(1:n)%selev       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='slope_mean')   basin(1:n)%slope       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='slope_std ')   basin(1:n)%sslope      = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_tp    ')   pointsource(1:n,1,1)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_sp    ')   pointsource(1:n,2,1)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_tn    ')   pointsource(1:n,3,1)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_in    ')   pointsource(1:n,4,1)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_vol   ')   pointsource(1:n,5,1)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_tp    ')   pointsource(1:n,1,2)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_sp    ')   pointsource(1:n,2,2)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_tn    ')   pointsource(1:n,3,2)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_in    ')   pointsource(1:n,4,2)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_vol   ')   pointsource(1:n,5,2)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_tp    ')   pointsource(1:n,1,3)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_sp    ')   pointsource(1:n,2,3)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_tn    ')   pointsource(1:n,3,3)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_in    ')   pointsource(1:n,4,3)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_vol   ')   pointsource(1:n,5,3)   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_tp    ')  load(1:n)%tpconcloc    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_sp    ')  load(1:n)%sppartloc    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_tn    ')  load(1:n)%tnconcloc    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_in    ')  load(1:n)%inpartloc    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_vol   ')  load(1:n)%volloc       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='region    ')  basin(1:n)%region      = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='lakeregion')  basin(1:n)%lakeregion  = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='lake_depth')  basin(1:n)%lakedepth(2)= xr(1:n,rindex(i))
      IF(str(i)(1:10)=='icatch    ')  basin(1:n)%ilakecatch  = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='maindown  ')  pathsubid(1:n)%main    = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='grwdown   ')  pathsubid(1:n)%grw1      = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='grwolake  ')  pathsubid(1:n)%grwtolake = xr(1:n,rindex(i))
      IF(str(i)(1:3)=='slc') THEN
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        classbasin(1:n,idslc)%part=xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:5)=='dhslc') THEN
        idslc = 0
        IF(ICHAR(str(i)(7:7))>=49 .AND. ICHAR(str(i)(7:7))<=57)THEN
          idslc = ICHAR(str(i)(7:7)) - 48
        ENDIF
        IF(ICHAR(str(i)(8:8))>=48 .AND. ICHAR(str(i)(8:8))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(8:8)) - 48
        ENDIF
        classbasin(1:n,idslc)%deltah = xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:3)=='scr') THEN    !area part with secondary crop
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        classbasin(1:n,idslc)%part2cr = xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:10)=='parreg    ')THEN
        basin(1:n)%parregion    = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='wqparreg  ')  basin(1:n)%wqparregion  = xi(1:n,iindex(i))
      ! Start CryoProcesses: regions for ilake and olake parameters
      IF(str(i)(1:10)=='ilakereg  ')THEN
        basin(1:n)%ilakeregion  = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='olakereg  ')  basin(1:n)%olakeregion  = xi(1:n,iindex(i))
       ! ...and non-contributing fraction
      IF(str(i)(1:10)=='noncontr  ')  basin(1:n)%noncontr     = xr(1:n,rindex(i))
      ! End CryoProcesses    
      IF(str(i)(1:10)=='rivlen    ')  basin(1:n)%rivlen       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='wetdep_n  ')  load(1:n)%inwetdep      = xr(1:n,rindex(i))*1.E-3    !ug/L -> mg/L
      IF(str(i)(1:10)=='drydep_n1 ')  load(1:n)%indrydep(1)   = xr(1:n,rindex(i))*onetspday  !kg/km2/d->kg/km2/ts
      IF(str(i)(1:10)=='drydep_n2 ')  load(1:n)%indrydep(2)   = xr(1:n,rindex(i))*onetspday/timesteps_per_day
      IF(str(i)(1:10)=='drydep_n3 ')  load(1:n)%indrydep(3)   = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:10)=='lrwet_area')  wetland(1:n,1)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_area')  wetland(1:n,2)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lrwet_dep ')  wetland(1:n,1)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_dep ')  wetland(1:n,2)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lrwet_part')  wetland(1:n,1)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_part')  wetland(1:n,2)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='buffer    ')  basin(1:n)%buffer       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='close_w   ')  basin(1:n)%closewater   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lakedataid')  lakedataid(1:n)         = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='pobsid    ')  pobsid(1:n)             = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='tobsid    ')  tobsid(1:n)             = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='sfobsid   ')THEN
        IF(readsfobs)  sfobsid(1:n)        = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='swobsid   ')THEN
        IF(readswobs) swobsid(1:n)         = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='uobsid    ')THEN
        IF(readwind) uobsid(1:n)         = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='rhobsid   ')THEN
        IF(readhumid) rhobsid(1:n)         = xi(1:n,iindex(i))
      ENDIF
      IF(str(i)(1:10)=='tminobsid   ')THEN
        IF(readtminmaxobs)  tminobsid(1:n) = xi(1:n,iindex(i))  
      ENDIF
      IF(str(i)(1:10)=='tmaxobsid   ')THEN
        IF(readtminmaxobs)  tmaxobsid(1:n) = xi(1:n,iindex(i))  
      ENDIF
    ENDDO

    !Check groundwater flow direction      
    IF(SUM(pathsubid(1:n)%grw1)==0) pathsubid(1:n)%grw1 = pathsubid(1:n)%main
    !Calculate number of parameter region and lakeregion
    nregions = MAXVAL(basin(1:n)%parregion)   !temporary value, aquifers may have more regions
    nwqregions = MAXVAL(basin(1:n)%wqparregion)
    
    !ilake and olake parameter regions (please note difference from lakeregions - the ilake and olake regions are more similar to wqparregion)
    nilakeregions = MAXVAL(basin(1:n)%ilakeregion)
    nolakeregions = MAXVAL(basin(1:n)%olakeregion)
    
    nlakeregions = MAXVAL(basin(1:n)%lakeregion)
    !Check if id for forcing is present and to be used
    IF((SUM(pobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
      pobsid(1:n) = basin(1:n)%subid
    ENDIF
    IF((SUM(tobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
      tobsid(1:n) = basin(1:n)%subid
    ENDIF
    IF(readsfobs)THEN
      IF((SUM(sfobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        sfobsid(1:n) = basin(1:n)%subid      
      ENDIF
    ENDIF    
    IF(readswobs)THEN
      IF((SUM(swobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        swobsid(1:n) = basin(1:n)%subid      
      ENDIF
    ENDIF
    IF(readwind)THEN
      IF((SUM(uobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        uobsid(1:n) = basin(1:n)%subid      
      ENDIF
    ENDIF
    IF(readhumid)THEN
      IF((SUM(rhobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        rhobsid(1:n) = basin(1:n)%subid      
      ENDIF
    ENDIF
    IF(readtminmaxobs)THEN
      IF((SUM(tminobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        tminobsid(1:n) = basin(1:n)%subid      
      ENDIF
      IF((SUM(tmaxobsid(1:n))==0) .OR. (.NOT.readobsid))THEN
        tmaxobsid(1:n) = basin(1:n)%subid      
      ENDIF
    ENDIF
    !Calculate point source loads (PS in GeoData) or zeroing the loads.
    DO i = 1,n
      DO k = 1,max_pstype
        load(i)%psvol(k) = pointsource(i,5,k)/86400.  !m3/d->m3/s
      ENDDO
      IF(i_in>0)THEN
        DO k = 1,max_pstype
          help = pointsource(i,3,k)*pointsource(i,5,k)/timesteps_per_day*1.E-3    !kg TN/timestep
          load(i)%psload(k,i_in) = help*pointsource(i,4,k)
          load(i)%psload(k,i_on) = help*(1.-pointsource(i,4,k))
        ENDDO
      ENDIF
      IF(i_sp>0)THEN
        DO k = 1,max_pstype
          help = pointsource(i,1,k)*pointsource(i,5,k)/timesteps_per_day*1.E-3    !kg TP/timestep
          load(i)%psload(k,i_sp) = help*pointsource(i,2,k)
          load(i)%psload(k,i_pp) = help*(1.-pointsource(i,2,k))
        ENDDO
      ENDIF
    ENDDO
    
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(pointsource)) DEALLOCATE(pointsource)
    

  END SUBROUTINE read_and_calc_basindata

  !>\brief Gets the information about nutrient point sources from PointSourceData.txt
  !>Recalculate into form suitable for model. Set permanent loads.
  !!
  !>\b Consequences Module modvar variable load is set. Module worldvar variable 
  !>psdates may be allocated and set.
  !>
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE load_pointsourcedata(dir,infile,n,status) 

    USE WORLDVAR, ONLY : steplen, &
                         psdates    !OUT
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp, &
                       i_t1,i_t2, &
                       basin,       &
                       load             !OUT
    USE READWRITE_ROUTINES, ONLY : find_reorder_index
 
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(IN)  :: n                    !<Number of subbasins
    INTEGER, INTENT(OUT) :: status               !<Error status
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    LOGICAL fexist                     !Status of file exist
    INTEGER i,j,k
    INTEGER nrows
    INTEGER maxdates
    INTEGER,ALLOCATABLE :: aindex(:)                  !Index connectingPS-row to basinindex
    INTEGER,ALLOCATABLE :: readsubid(:)
    INTEGER,ALLOCATABLE :: readpstype(:)
    INTEGER,ALLOCATABLE :: readsource(:)
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL,ALLOCATABLE :: pointsource(:,:) !Read pointsource info (subbasins,column)
    CHARACTER(LEN=50),ALLOCATABLE :: readdate1(:)  !String date read from file
    CHARACTER(LEN=50),ALLOCATABLE :: readdate2(:)  !String date read from file
    TYPE(DateType) :: onedate                 !One date
    TYPE(DateType),ALLOCATABLE :: psdateslocal(:)  !Array of dates from file

    !Start of subroutine
    status = 0
    !Note. If(when) old format is removed from GeoData-reading, zeroing of load%psvol and %psload must be done here.

    !Read PointsourceData.txt, subroutine will allocate readsubid etc.
    CALL read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 
    IF(.NOT.fexist) RETURN
    IF(status/=0)RETURN
    
    !Allocate and initiate variables 
    IF(.NOT.ALLOCATED(aindex)) ALLOCATE(aindex(nrows))
    DO i=1,n    !Zeroing the point source loads for accumulation
      load(i)%psvol(:) = 0.
      load(i)%psload(:,:) = 0.
    ENDDO
 
    !Collect dates where point source load changes; fromdata and todate+1
    maxdates = 0
    IF(.NOT.ALLOCATED(psdateslocal)) ALLOCATE(psdateslocal(2*nrows))
    DO i = 1,nrows
      IF(readdate1(i)(1:1)=='0')THEN
      ELSE
        CALL string_convert_to_datetype(readdate1(i),onedate)
        CALL add_date_to_array(onedate,2*nrows,maxdates,psdateslocal)
      ENDIF
      IF(readdate2(i)(1:1)=='0')THEN
      ELSE
        CALL string_convert_to_datetype(readdate2(i),onedate)
        onedate = AddDates(onedate,steplen)
        CALL add_date_to_array(onedate,2*nrows,maxdates,psdateslocal)
      ENDIF
    ENDDO
    IF(maxdates>0)THEN
      IF(.NOT.ALLOCATED(psdates)) ALLOCATE(psdates(maxdates))
      psdates=psdateslocal(1:maxdates)
    ENDIF
    
    IF(.NOT.ALLOCATED(psdates))THEN
      !Calculate permanent point source loads
      onetspday = 1./REAL(timesteps_per_day)
      CALL find_reorder_index(infile,nrows,readsubid,n,basin%subid,.FALSE.,aindex,status) !status=1 is not an error
      DO j=1,nrows
        IF(aindex(j)>0)THEN !Point sources not in model are ignored
          IF(readdate1(j)=='0'.AND.readdate2(j)=='0')THEN
            i=aindex(j)
            IF(pointsource(j,5)<0.)THEN !Negative point source, abstraction
              load(i)%abstrvol(readsource(j)) = load(i)%abstrvol(readsource(j)) - pointsource(j,5)/86400.  !m3/d->m3/s
            ELSEIF(pointsource(j,5)>0.)THEN
              k=readpstype(j)
              IF(k>0)THEN
                load(i)%psvol(k) = load(i)%psvol(k) + pointsource(j,5)/86400.  !m3/d->m3/s
                IF(i_in>0)THEN
                  help = pointsource(j,3)*pointsource(j,5)*onetspday*1.E-3    !kg TN/timestep
                  load(i)%psload(k,i_in) = load(i)%psload(k,i_in) + help*pointsource(j,4)
                  load(i)%psload(k,i_on) = load(i)%psload(k,i_on) + help*(1.-pointsource(j,4))
                ENDIF
                IF(i_sp>0)THEN
                  help = pointsource(j,1)*pointsource(j,5)*onetspday*1.E-3    !kg TP/timestep
                  load(i)%psload(k,i_sp) = load(i)%psload(k,i_sp) + help*pointsource(j,2)
                  load(i)%psload(k,i_pp) = load(i)%psload(k,i_pp) + help*(1.-pointsource(j,2))
                ENDIF
                IF(i_t1>0)THEN
                  help = pointsource(j,6)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                  load(i)%psload(k,i_t1) = load(i)%psload(k,i_t1) + help
                ENDIF
                IF(i_t2>0)THEN
                  help = pointsource(j,7)*pointsource(j,5)*onetspday*1.E-3    !amount T2/timestep
                  load(i)%psload(k,i_t2) = load(i)%psload(k,i_t2) + help
                ENDIF
              ELSE
                WRITE(6,*) 'Warning: An permanent point source without ps_type found'
                WRITE(6,*) 'Warning: This is not a recommended use of PointSourceData.txt'
                WRITE(6,*) 'Warning: Point source may be ignored!'
              ENDIF
            ENDIF
          ENDIF
        ELSE
          WRITE(6,*) 'Point source in subid',readsubid(j),'ignored, not in model domain.'
        ENDIF
      ENDDO
    ENDIF
    
    !Deallocate local arrays
    IF(ALLOCATED(aindex)) DEALLOCATE(aindex)
    IF(ALLOCATED(readsubid)) DEALLOCATE(readsubid)
    IF(ALLOCATED(readpstype)) DEALLOCATE(readpstype)
    IF(ALLOCATED(readsource)) DEALLOCATE(readsource)
    IF(ALLOCATED(readdate1)) DEALLOCATE(readdate1)
    IF(ALLOCATED(readdate2)) DEALLOCATE(readdate2)
    IF(ALLOCATED(pointsource)) DEALLOCATE(pointsource)
    IF(ALLOCATED(psdateslocal)) DEALLOCATE(psdateslocal)

  END SUBROUTINE load_pointsourcedata

  !>\brief Reads the matrix of point source data from the file. 
  !!
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(OUT) :: nrows                !<Number of rows in file
    INTEGER, INTENT(OUT) :: status               !<Error status
    LOGICAL, INTENT(OUT) :: fexist               !<Status of file exist
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readsubid(:)    !<subid if point sources
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readpstype(:)   !<type of point source
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readsource(:)   !<abstraction source 
    REAL,ALLOCATABLE, INTENT(INOUT) :: pointsource(:,:)   !<Read pointsource info (subbasins,column)
    CHARACTER(LEN=50),ALLOCATABLE, INTENT(INOUT) :: readdate1(:)  !<From date of point source
    CHARACTER(LEN=50),ALLOCATABLE, INTENT(INOUT) :: readdate2(:)  !<To date of point source
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    INTEGER i
    INTEGER ncols,mcols
    INTEGER,ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER,ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)             !Index for column integer variables
    INTEGER,ALLOCATABLE :: sindex(:)             !Index for column string variables
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=50),ALLOCATABLE :: xs(:,:)  !String data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)   !Content string

    !Start of subroutine
    status = 0
    
    !Check existence of file
    INQUIRE(FILE=TRIM(dir)//infile,EXIST=fexist)
    IF(.NOT.fexist) RETURN
    
    !Count number of columns and rows in PointSourceData
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,0,ncols,status)
    IF(status/=0)RETURN
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,1,nrows,status)
    IF(status/=0)RETURN

    !Allocate variables for holding file data
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(sindex)) ALLOCATE(sindex(ncols))
    IF(.NOT.ALLOCATED(readsubid)) ALLOCATE(readsubid(nrows))
    IF(.NOT.ALLOCATED(readpstype)) ALLOCATE(readpstype(nrows))
    IF(.NOT.ALLOCATED(readsource)) ALLOCATE(readsource(nrows))
    IF(.NOT.ALLOCATED(readdate1)) ALLOCATE(readdate1(nrows))
    IF(.NOT.ALLOCATED(readdate2)) ALLOCATE(readdate2(nrows))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    IF(.NOT.ALLOCATED(xs)) ALLOCATE(xs(nrows,ncols))
    IF(.NOT.ALLOCATED(pointsource)) ALLOCATE(pointsource(nrows,7))
    readpstype = 0  !default point source is first
    readsource = 1  !default abstraction source is main river
    pointsource = 0.
    readdate1 = '0'
    readdate2 = '0'
    
    !Read PointSourceData.txt file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(dir)//infile, STATUS = 'old')     !Open PointSourceData-file

    !Reads the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(TRIM(dir)//infile)
      RETURN
    ENDIF

    !Code column for variable type
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_source ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_type   ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_tpconc ') code(i) = i_real
      IF(str(i)(1:10)=='ps_spfrac ') code(i) = i_real
      IF(str(i)(1:10)=='ps_tnconc ') code(i) = i_real
      IF(str(i)(1:10)=='ps_infrac ') code(i) = i_real
      IF(str(i)(1:10)=='ps_vol    ') code(i) = i_real
      IF(str(i)(1:10)=='ps_t1     ') code(i) = i_real
      IF(str(i)(1:10)=='ps_t2     ') code(i) = i_real
    ENDDO

    !Read all data values
    CALL read_basindata6(fileunit_temp,TRIM(dir)//infile,ncols,nrows,ncols,50,code,rindex,iindex,sindex,xi,xr,xs)
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(TRIM(dir)//infile)

    !Move data from read matrix to pointsource variables
    DO i = 1,ncols
      IF(str(i)(1:10)=='fromdate  ') readdate1(1:nrows)     = xs(1:nrows,sindex(i)) !String date
      IF(str(i)(1:10)=='todate    ') readdate2(1:nrows)     = xs(1:nrows,sindex(i))
      IF(str(i)(1:10)=='subid     ') readsubid(1:nrows)     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_type   ') readpstype(1:nrows)    = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_source ') readsource(1:nrows)    = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_tpconc ') pointsource(1:nrows,1) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_spfrac ') pointsource(1:nrows,2) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_tnconc ') pointsource(1:nrows,3) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_infrac ') pointsource(1:nrows,4) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_vol    ') pointsource(1:nrows,5) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_t1     ') pointsource(1:nrows,6) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_t2     ') pointsource(1:nrows,7) = xr(1:nrows,rindex(i))
    ENDDO
    
    !Deallocate local arrays
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(sindex)) DEALLOCATE(sindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(xs)) DEALLOCATE(xs)

  END SUBROUTINE read_pointsourcedata

  !>\brief Gets the information about temporary point sources from PointSourceData.txt
  !>Recalculate into form suitable for model.
  !!
  !>\b Consequences Module modvar variable load and t1load are changed
  !>
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_pointsources(dir,infile,ns,time,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         bdate,  &
                         psdates,     &
                         fileunit_temp
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp, &
                       i_t1,i_t2,   &
                       max_pstype,  &
                       basin,       &
                       load,        & !OUT
                       tload,       & !OUT
                       tloadexist     !OUT
    USE READWRITE_ROUTINES, ONLY : find_reorder_index
 
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(IN)  :: ns                   !<Number of subbasins (submodel)
     TYPE(DateType), INTENT(IN) :: time           !<Current time (e.g. 2005-01-26 18:00)
    INTEGER, INTENT(OUT) :: status               !<Error status
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    LOGICAL fexist                     !Status of file exist
    INTEGER i,j,k
    INTEGER nrows
    INTEGER maxdates
    INTEGER numt1
    INTEGER,ALLOCATABLE :: aindex(:)                  !Index connectingPS-row to basinindex
    INTEGER,ALLOCATABLE :: readsubid(:)
    INTEGER,ALLOCATABLE :: readpstype(:)
    INTEGER,ALLOCATABLE :: readsource(:)
    INTEGER,ALLOCATABLE :: tloadindexlocal(:,:)
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL,ALLOCATABLE :: pointsource(:,:) !Read pointsource info (subbasins,column)
    REAL,ALLOCATABLE :: tloadlocal(:,:)  !Read pointsource info on T1 sources
    CHARACTER(LEN=50),ALLOCATABLE :: readdate1(:)  !String date read from file
    CHARACTER(LEN=50),ALLOCATABLE :: readdate2(:)  !String date read from file
    TYPE(DateType),ALLOCATABLE :: psdateslocal(:,:)  !Array of dates from file

    !Start of subroutine
    status = 0
    
    !Check if beginning of simulation or current time in change array
    IF(time>bdate)THEN
      maxdates = SIZE(psdates)
      DO i = 1,maxdates
        IF(time==psdates(i)) EXIT
      ENDDO
      IF(i>maxdates) RETURN
    ENDIF
    
    !Read PointsourceData.txt, subroutine will allocate readsubid etc.
    CALL read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 
    IF(status/=0)RETURN

    !Allocate and initiate variables 
    IF(.NOT.ALLOCATED(aindex)) ALLOCATE(aindex(nrows))
    IF(.NOT.ALLOCATED(psdateslocal)) ALLOCATE(psdateslocal(nrows,2))
    IF(.NOT.ALLOCATED(tloadindexlocal)) ALLOCATE(tloadindexlocal(2,nrows))
    IF(.NOT.ALLOCATED(tloadlocal)) ALLOCATE(tloadlocal(3,nrows))
    DO i=1,ns    !Zeroing the point source loads for accumulation
      load(i)%psvol(:) = 0.
      load(i)%psload(:,:) = 0.
    ENDDO

    !Transform dates to DateType
    DO j=1,nrows
      IF(readdate1(j)=='0')THEN
        psdateslocal(j,1)=datetype(3000,0,0,0,0)
      ELSE
        CALL string_convert_to_datetype(readdate1(j),psdateslocal(j,1))
      ENDIF
      IF(readdate2(j)=='0')THEN
        psdateslocal(j,2)=datetype(0,0,0,0,0)
      ELSE
        CALL string_convert_to_datetype(readdate2(j),psdateslocal(j,2))
      ENDIF
    ENDDO

    !Calculate and add permanent point source loads and 
    !temporary point source load discharging today
    numt1 = 0
    onetspday = 1./REAL(timesteps_per_day)
    CALL find_reorder_index(infile,nrows,readsubid,ns,basin%subid,.FALSE.,aindex,status) !status=1 is not an error
    DO j=1,nrows
      IF(aindex(j)>0)THEN !Point sources not in (sub)model are ignored
        i=aindex(j)
        IF(readdate1(j)=='0'.AND.readdate2(j)=='0'.OR. &              !permanent point source
          (readdate1(j)=='0'.AND.psdateslocal(j,2)>=time) .OR. &      !ending point source
          (psdateslocal(j,1)<=time.AND.readdate2(j)=='0') .OR. &      !starting point source
          (psdateslocal(j,1)<=time.AND.psdateslocal(j,2)>=time))THEN  !temporary point source
          IF(pointsource(j,5)<0.)THEN   !Negative point source, abstraction
            load(i)%abstrvol(readsource(j)) = load(i)%abstrvol(readsource(j)) - pointsource(j,5)/86400.  !m3/d->m3/s
          ELSEIF(pointsource(j,5)>0.)THEN
            k=readpstype(j)
            IF(k>0)THEN
              load(i)%psvol(k) = load(i)%psvol(k) + pointsource(j,5)/86400.  !m3/d->m3/s
              IF(i_in>0)THEN
                help = pointsource(j,3)*pointsource(j,5)*onetspday*1.E-3    !kg TN/timestep
                load(i)%psload(k,i_in) = load(i)%psload(k,i_in) + help*pointsource(j,4)
                load(i)%psload(k,i_on) = load(i)%psload(k,i_on) + help*(1.-pointsource(j,4))
              ENDIF
              IF(i_sp>0)THEN
                help = pointsource(j,1)*pointsource(j,5)*onetspday*1.E-3    !kg TP/timestep
                load(i)%psload(k,i_sp) = load(i)%psload(k,i_sp) + help*pointsource(j,2)
                load(i)%psload(k,i_pp) = load(i)%psload(k,i_pp) + help*(1.-pointsource(j,2))
              ENDIF
              IF(i_t1>0)THEN
                help = pointsource(j,6)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                load(i)%psload(k,i_t1) = load(i)%psload(k,i_t1) + help
              ENDIF
              IF(i_t2>0)THEN
                help = pointsource(j,7)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                load(i)%psload(k,i_t2) = load(i)%psload(k,i_t2) + help
              ENDIF
            ELSEIF(readsource(j)>0)THEN
              IF(i_t1>0)THEN
                numt1 = numt1 + 1
                tloadlocal(1,numt1) = pointsource(j,5)/86400.  !m3/d->m3/s
                tloadlocal(2,numt1) = pointsource(j,6)
                tloadlocal(3,numt1) = pointsource(j,7)
                tloadindexlocal(1,numt1) = readsource(j)
                tloadindexlocal(2,numt1) = i
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ELSE
        WRITE(6,*) 'Point source in subid',readsubid(j),'ignored, not in model domain.'
      ENDIF
    ENDDO

    !Set T1 loads if present
    IF(numt1>0)THEN
      IF(.NOT.ALLOCATED(tload))THEN
        ALLOCATE(tload(numt1))
        ALLOCATE(tloadexist(ns))
        tloadexist = .FALSE.
      ENDIF
      tload%psvol = tloadlocal(1,1:numt1)
      tload%psconc = tloadlocal(2,1:numt1)
      tload%pstemp = tloadlocal(3,1:numt1)
      tload%sw_code = tloadindexlocal(1,1:numt1)
      tload%subindex = tloadindexlocal(2,1:numt1)
      tloadexist(tloadindexlocal(2,1:numt1)) = .TRUE.
    ELSE
      IF(ALLOCATED(tload))THEN
        DEALLOCATE(tload)
        DEALLOCATE(tloadexist)
      ENDIF
    ENDIF
    
    !Deallocate local variables
    IF(ALLOCATED(readsubid)) DEALLOCATE(readsubid)
    IF(ALLOCATED(readpstype)) DEALLOCATE(readpstype)
    IF(ALLOCATED(readsource)) DEALLOCATE(readsource)
    IF(ALLOCATED(readdate1)) DEALLOCATE(readdate1)
    IF(ALLOCATED(readdate2)) DEALLOCATE(readdate2)
    IF(ALLOCATED(pointsource)) DEALLOCATE(pointsource)
    IF(ALLOCATED(aindex)) DEALLOCATE(aindex)
    IF(ALLOCATED(psdateslocal)) DEALLOCATE(psdateslocal)
    IF(ALLOCATED(tloadindexlocal)) DEALLOCATE(tloadindexlocal)
    IF(ALLOCATED(tloadlocal)) DEALLOCATE(tloadlocal)
    
  END SUBROUTINE get_current_pointsources

  !>\brief Gets the information about irrigation from MgmtData.
  !!Reads the matrix of basin data values from the file vith mcols
  !!columns. Output: data variables, heading variable and actual
  !!number of data columns
  !--------------------------------------------------------------------
  SUBROUTINE load_irrigation_data(funit,dir,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         maxcharpath
    USE MODVAR, ONLY : irrigationsystem    !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(OUT) :: status                 !<Error status
    
    !Local variables
    INTEGER ncols,nrows                !size of data in file
    INTEGER mcol
    INTEGER i,j
    LOGICAL fileex                     !Status of file
    LOGICAL swcol,gwcol                !Status of sw/gw_part columns
    CHARACTER(LEN=maxcharpath) infile             !Name of irrigation characteristics file 
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER, ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER, ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER, ALLOCATABLE :: iindex(:)             !Index for column integer variables
    REAL, ALLOCATABLE    :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=10), ALLOCATABLE :: str(:)      !File column content string

    status = 0
    infile = TRIM(dir)//'MgmtData.txt'

    !Check if file exist
    INQUIRE(FILE=TRIM(infile),EXIST=fileex)
    IF(.NOT.fileex)THEN
       WRITE(6,*) 'No MgmtData.txt file found.'
       RETURN
    ENDIF

    !Count number of columns in MgmtData
    CALL count_data_cols(funit,TRIM(infile),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of subbasins in MgmtData
    CALL count_data_rows(funit,TRIM(infile),1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
       WRITE(6,*) 'No data in MgmtData.txt file found.'
       RETURN
    ENDIF

    !Allocate and initiate local variables
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))  
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))  
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols)) 
    swcol = .FALSE.; gwcol = .FALSE.

    !Open file
    OPEN(UNIT = funit,FILE = TRIM(infile), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(funit,ncols,str,mcol,status)
    IF(status.NE.0) THEN
       WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
       RETURN
    ENDIF

    !Allocate irrigation variable 
    IF(.NOT.ALLOCATED(irrigationsystem)) ALLOCATE(irrigationsystem(nrows))

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,ncols
       IF(str(i)(1:10)=='subid     ') code(i) = i_intg
       IF(str(i)(1:10)=='regsrcid  ') code(i) = i_intg
       IF(str(i)(1:10)=='local_eff ') code(i) = i_real
       IF(str(i)(1:10)=='region_eff') code(i) = i_real
       IF(str(i)(1:10)=='gw_part   ') code(i) = i_real
       IF(str(i)(1:10)=='demandtype') code(i) = i_intg
       IF(str(i)(1:10)=='irrdam    ') code(i) = i_intg
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,ncols,nrows,ncols,code,rindex,iindex,xi,xr)

    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    DO i = 1,ncols
       IF(str(i)(1:10)=='subid     ') irrigationsystem(1:nrows)%subid       = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='regsrcid  ') irrigationsystem(1:nrows)%regsourceid = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='local_eff ') irrigationsystem(1:nrows)%local_eff   = xr(1:nrows,rindex(i))
       IF(str(i)(1:10)=='region_eff') irrigationsystem(1:nrows)%reg_eff     = xr(1:nrows,rindex(i))
       IF(str(i)(1:10)=='demandtype') irrigationsystem(1:nrows)%demandtype  = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='gw_part   ')THEN
          irrigationsystem(1:nrows)%gw_part     = xr(1:nrows,rindex(i))
          gwcol = .TRUE.
       ENDIF
       IF(str(i)(1:10)=='irrdam    ')THEN
          DO j = 1,nrows
             irrigationsystem(j)%dam     = xi(j,iindex(i))==1
          ENDDO
       ENDIF
    ENDDO

    !Check irrigation used (if file used for other things in the future)
    IF(gwcol)THEN
       irrigationsystem(1:nrows)%sw_part = 1. - irrigationsystem(1:nrows)%gw_part
    ELSE
       WRITE(6,*) 'No irrigation source division found. Irrigation not used.'
       IF(ALLOCATED(irrigationsystem)) DEALLOCATE(irrigationsystem)
    ENDIF

    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)  
    IF(ALLOCATED(code)) DEALLOCATE(code)  
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(str)) DEALLOCATE(str) 

  END SUBROUTINE load_irrigation_data

!--------------------------------------------------------------------
  !>\brief Gets information about forcing data coupling to subbasin 
  !!and other information.
  !!Reads the matrix of basin data values from the file with mcols
  !!columns. Output: data variables, heading variable and actual
  !!number of data columns
  !!
  !>/b Consequences Module worldvar variables tobsid, pobsid, sfobsid, 
  !>swobsid,uobsid,rhobsid,tminobsid,tmaxobsid are set. 
  !>Module modvar variable forcing is allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_forcing_information_data(funit,dir,ns,status) 

    USE WORLDVAR, ONLY : simsequence, &
                         tobsid,      &   !OUT
                         pobsid,      &   !OUT
                         sfobsid,     &   !OUT
                         swobsid,     &   !OUT
                         uobsid,     &   !OUT
                         rhobsid,     &   !OUT
                         tminobsid,     &   !OUT
                         tmaxobsid,     &   !OUT
                         readsfobs,   &
                         readswobs,   &
                         readwind,    &
                         readhumid,    &
                         readtminmaxobs,   &
                         readobsid,   &
                         i_str,       &
                         i_intg,      &
                         i_real,      &
                         maxcharpath
    USE MODVAR, ONLY : basin,   &
                       forcing    !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit number for file
    CHARACTER (LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(IN) :: ns                      !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: status                 !<Error status
    
    !Local variables
    INTEGER ncols                      !size of data in file
    INTEGER mcols
    INTEGER i,j
    INTEGER isubid
    INTEGER itobsid,ipobsid,isfobsid,iswobsid,iuobsid,irhobsid,itminobsid,itmaxobsid,itobselev  !Column in file
    LOGICAL fileex                     !Existence of file
    CHARACTER(LEN=3) seqnr             !Sequence number
    CHARACTER(LEN=maxcharpath) infile             !Name and path of file 
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER, ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER, ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER, ALLOCATABLE :: iindex(:)             !Index for column integer variables
    REAL, ALLOCATABLE    :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=10), ALLOCATABLE :: colstr(:)   !File column content string

    status = 0
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      infile = TRIM(dir)//'ForcKey_'//seqnr//'.txt'
    ELSE
      infile = TRIM(dir)//'ForcKey.txt'
    ENDIF
    
    !Check if file exist
    INQUIRE(FILE=TRIM(infile),EXIST=fileex)
    IF(.NOT.fileex) RETURN

    !Count number of columns in file
    CALL count_data_cols(funit,TRIM(infile),0,ncols,status)
    IF(status/=0)RETURN

    !Allocation of local variables
    ALLOCATE(xi(ns,ncols))  
    ALLOCATE(code(ncols))  
    ALLOCATE(rindex(ncols))
    ALLOCATE(iindex(ncols))
    ALLOCATE(xr(ns,ncols))
    ALLOCATE(colstr(ncols)) 

    !Read file
    OPEN(UNIT=funit,FILE=TRIM(infile),STATUS='old',ACTION='read',ERR=900)
      CALL read_column_headings(funit,ncols,colstr,mcols,status)
      IF(status/=0) RETURN
      code = i_str
      isubid = 0; itobsid = 0; ipobsid = 0; isfobsid = 0; iswobsid = 0; iuobsid = 0; irhobsid = 0; itminobsid = 0; itmaxobsid = 0; itobselev = 0
      DO i = 1,mcols
        IF(colstr(i)(1:10)=='subid     ')THEN
          code(i) = i_intg
          isubid = i
        ENDIF  
        IF(colstr(i)(1:10)=='sfobsid   ')THEN
          code(i) = i_intg
          IF(readsfobs) isfobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='swobsid   ')THEN
          code(i) = i_intg
          IF(readswobs) iswobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='uobsid    ')THEN
          code(i) = i_intg
          IF(readwind) iuobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='rhobsid   ')THEN
          code(i) = i_intg
          IF(readhumid) irhobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='tminobsid ')THEN
          code(i) = i_intg
          IF(readtminmaxobs) itminobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='tmaxobsid ')THEN
          code(i) = i_intg
          IF(readtminmaxobs) itmaxobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='tobsid    ')THEN
          code(i) = i_intg
          itobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='pobsid    ')THEN
          code(i) = i_intg
          ipobsid = i
        ENDIF  
        IF(colstr(i)(1:10)=='tobselev  ')THEN
          code(i) = i_real
          itobselev = i
        ENDIF  
      ENDDO
      IF(isubid==0)THEN
        WRITE(6,*) 'ERROR: subid missing in file:', TRIM(infile)
        status = 1
        RETURN
      ENDIF
      CALL read_basindata5(funit,infile,ncols,ns,mcols,code,rindex,iindex,xi,xr) 
      CLOSE(funit) 

      !Set forcing data information to subbasin array
      IF(itobselev>0)THEN
        IF(.NOT.ALLOCATED(forcing)) ALLOCATE(forcing(ns))
        DO i = 1,ns
          IF(basin(i)%subid==xi(i,iindex(isubid)))THEN
            forcing(i)%tobselev = xr(i,rindex(itobselev))
          ELSE
            DO j = 1,ns
              IF(basin(i)%subid==xi(j,iindex(isubid)))THEN
                forcing(i)%tobselev = xr(j,rindex(itobselev))
                EXIT
              ENDIF  
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      !Set subbasin - observation coupling variables
      IF(readobsid)THEN
        DO i = 1,ns
          IF(basin(i)%subid==xi(i,iindex(isubid)))THEN
            IF(ipobsid>0)  pobsid(i)  = xi(i,iindex(ipobsid))
            IF(itobsid>0)  tobsid(i)  = xi(i,iindex(itobsid))
            IF(isfobsid>0) sfobsid(i) = xi(i,iindex(isfobsid))
            IF(iswobsid>0) swobsid(i) = xi(i,iindex(iswobsid))
            IF(iuobsid>0) uobsid(i) = xi(i,iindex(iuobsid))
            IF(irhobsid>0) rhobsid(i) = xi(i,iindex(irhobsid))
            IF(itminobsid>0)  tminobsid(i)  = xi(i,iindex(itminobsid))
            IF(itmaxobsid>0)  tmaxobsid(i)  = xi(i,iindex(itmaxobsid))
          ELSE
            DO j = 1,ns
              IF(basin(i)%subid==xi(j,iindex(isubid)))THEN
                IF(ipobsid>0)  pobsid(i)  = xi(j,iindex(ipobsid))
                IF(itobsid>0)  tobsid(i)  = xi(j,iindex(itobsid))
                IF(isfobsid>0) sfobsid(i) = xi(j,iindex(isfobsid))
                IF(iswobsid>0) swobsid(i) = xi(j,iindex(iswobsid))
                IF(iuobsid>0) uobsid(i) = xi(j,iindex(iuobsid))
                IF(irhobsid>0) rhobsid(i) = xi(j,iindex(irhobsid))
                IF(itminobsid>0)  tminobsid(i)  = xi(j,iindex(itminobsid))
                IF(itmaxobsid>0)  tmaxobsid(i)  = xi(j,iindex(itmaxobsid))
                EXIT
              ENDIF  
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      !Check if id for forcing is present and to be used
      IF((SUM(pobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
        pobsid(1:ns) = basin(1:ns)%subid
      ENDIF
      IF((SUM(tobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
        tobsid(1:ns) = basin(1:ns)%subid
      ENDIF
      IF(readsfobs)THEN
        IF((SUM(sfobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          sfobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
      ENDIF    
      IF(readswobs)THEN
        IF((SUM(swobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          swobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
      ENDIF
      IF(readwind)THEN
        IF((SUM(uobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          uobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
      ENDIF
      IF(readhumid)THEN
        IF((SUM(rhobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          rhobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
      ENDIF
      IF(readtminmaxobs)THEN
        IF((SUM(tminobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          tminobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
        IF((SUM(tmaxobsid(1:ns))==0) .OR. (.NOT.readobsid))THEN
          tmaxobsid(1:ns) = basin(1:ns)%subid      
        ENDIF
      ENDIF    

    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)  
    IF(ALLOCATED(code)) DEALLOCATE(code)  
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(colstr)) DEALLOCATE(colstr) 

    RETURN
    
900 WRITE(6,*) 'Error open file: ', TRIM(infile)
    status = 1
    RETURN

  END SUBROUTINE load_forcing_information_data

  !\brief Gets the information about update from update.txt.
  !-----------------------------------------------------------------
  SUBROUTINE read_update_data(funit,infile,ns,nqsub,qupdsub,nqarsub,nwarsub, &
       qarupdsub,warupdsub,qarupdar,warupdar,nwsub,wendsub,nrow,tpcorrarr,tncorrarr,tploccorrarr,tnloccorrarr)

    USE WORLDVAR, ONLY : i_str,i_intg,i_real

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: ns                   !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: nqsub                !<Number of subbasins with quseobs update
    INTEGER, INTENT(OUT) :: qupdsub(ns)          !<Subid of subbasins with quseobs update
    INTEGER, INTENT(OUT) :: nqarsub              !<Number of subbasins with qAR date
    INTEGER, INTENT(OUT) :: nwarsub              !<Number of subbasins with wAR date
    INTEGER, INTENT(OUT) :: qarupdsub(ns)        !<Subid of subbasins with qAR update
    INTEGER, INTENT(OUT) :: warupdsub(ns)        !<Subid of subbasins with wAR update
    REAL, INTENT(OUT)    :: qarupdar(ns)         !<AR-factor for qAR-update
    REAL, INTENT(OUT)    :: warupdar(ns)         !<AR-factor for wAR-update
    INTEGER, INTENT(OUT) :: nwsub                !<Number of subbasins with wendupd update
    INTEGER, INTENT(OUT) :: wendsub(ns)          !<Subid of subbasins with wendupd update
    INTEGER, INTENT(OUT) :: nrow                 !<Number of data rows in file
    REAL,    INTENT(OUT) :: tpcorrarr(ns,2)      !<Subid,value of subbasins with tpcorr update
    REAL,    INTENT(OUT) :: tncorrarr(ns,2)      !<Subid,value of subbasins with tncorr update
    REAL,    INTENT(OUT) :: tploccorrarr(ns,2)   !<Subid,value of subbasins with tploccorr update
    REAL,    INTENT(OUT) :: tnloccorrarr(ns,2)   !<Subid,value of subbasins with tnloccorr update

    !Local parameters
    INTEGER, PARAMETER :: maxcol = 30     !Maximum number of columns (read) in file
    INTEGER, PARAMETER :: nskip  = 1      !Number of initial heading rows 

    !Local variables declarations
    INTEGER i,irow                     !Loop counters
    INTEGER isub                       !Data column with subid
    INTEGER iquseobs                   !Data column with status for quseobs
    INTEGER iqar,iwar                  !Data column with status for q-AR and w-AR
    INTEGER iqarar                     !Data column with AR-factor for q-AR
    INTEGER iwendupd                   !Data column with status for wendupd
    INTEGER itpcorr                    !Data column with value for tp correction
    INTEGER itncorr                    !Data column with value for tn correction
    INTEGER itploccorr                 !Data column with value for local tp correction
    INTEGER itnloccorr                 !Data column with value for local tn correction
    INTEGER mcols                      !Actual number of columns
    INTEGER status                     !Error status
    INTEGER code(maxcol)               !Code for column variable
    INTEGER rindex(maxcol)             !Index for column real variables
    INTEGER iindex(maxcol)             !Index for column integer variables
    INTEGER, ALLOCATABLE :: xi(:,:)    !Integer data read from file
    INTEGER, ALLOCATABLE :: temparray(:) !Temporary data array
    REAL, ALLOCATABLE :: temparray2(:) !Temporary data array2 !JD2012_AR	  
    REAL, ALLOCATABLE    :: xr(:,:)    !Real data read from file
    CHARACTER(LEN=10) str(maxcol)      !Content string

    !Initiations
    status = 0
    nqsub  = 0
    qupdsub = 0
    nqarsub  = 0 
    nwarsub  = 0 
    qarupdsub = 0
    warupdsub = 0
    qarupdar = 0 
    warupdar = 0 
    nrow = 0
    tpcorrarr = 0.
    tncorrarr = 0.
    tploccorrarr = 0.
    tnloccorrarr = 0.

    !Count the number of rows and allocate variables accordingly
    CALL count_data_rows(funit,infile,nskip,nrow,status)
    IF(status/=0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrow,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrow,maxcol))

    !Open file for reading data
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     

    !Reads the column headings
    CALL read_column_headings(funit,maxcol,str,mcols,status)
    IF(status/=0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='quseobs   ') code(i) = i_intg
      IF(str(i)(1:10)=='wendupd   ') code(i) = i_intg
      IF(str(i)(1:10)=='tpcorr    ') code(i) = i_real
      IF(str(i)(1:10)=='tncorr    ') code(i) = i_real
      IF(str(i)(1:10)=='tploccorr ') code(i) = i_real
      IF(str(i)(1:10)=='tnloccorr ') code(i) = i_real
      IF(str(i)(1:10)=='qarupd    ') code(i) = i_intg !Q-AR on(1), off(0)
      IF(str(i)(1:10)=='warupd    ') code(i) = i_intg !W-AR on(1), off(0)
      IF(str(i)(1:10)=='arfact    ') code(i) = i_real !factor for AR
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,maxcol,nrow,mcols,code,rindex,iindex,xi,xr)

    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    !Find variables for update function
    isub = 0; iquseobs = 0; iwendupd=0; itpcorr=0; itncorr=0; itploccorr=0; itnloccorr=0; iqar=0; iqarar=0; iwar=0
    DO i = 1,mcols
      IF(str(i)(1:10)=='subid     ')   isub = i
      IF(str(i)(1:10)=='quseobs   ')   iquseobs = i
      IF(str(i)(1:10)=='wendupd   ')   iwendupd = i
      IF(str(i)(1:10)=='tpcorr    ')   itpcorr = i
      IF(str(i)(1:10)=='tncorr    ')   itncorr = i
      IF(str(i)(1:10)=='tploccorr ')   itploccorr = i
      IF(str(i)(1:10)=='tnloccorr ')   itnloccorr = i
      IF(str(i)(1:10)=='qarupd    ')   iqar = i  
      IF(str(i)(1:10)=='warupd    ')   iwar = i  
      IF(str(i)(1:10)=='arfact    ')   iqarar = i
    ENDDO
    IF(isub==0)THEN
      WRITE(6,*) 'No subid information in update.txt'
      RETURN
    ENDIF

    !Set variables for update function quseobs      
    IF(iquseobs>0)THEN      !Save subid for stations to be updated
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      temparray = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iquseobs))==1)THEN
          nqsub = nqsub + 1
          temparray(nqsub) = xi(irow,iindex(isub))
        ENDIF
      ENDDO
      IF(nqsub>0)THEN
        qupdsub(1:nqsub)=temparray(1:nqsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
    ENDIF

    !Set variables for update function qAR  
    IF(iqar>0)THEN
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      IF(.NOT.ALLOCATED(temparray2)) ALLOCATE(temparray2(nrow))
      temparray = 0
      temparray2 = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iqar))==1)THEN
          nqarsub = nqarsub + 1
          temparray(nqarsub) = xi(irow,iindex(isub))
          temparray2(nqarsub) = xr(irow,rindex(iqarar))
        ENDIF
      ENDDO
      IF(nqarsub>0)THEN
        qarupdsub(1:nqarsub) = temparray(1:nqarsub)
        qarupdar(1:nqarsub) = temparray2(1:nqarsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
      IF(ALLOCATED(temparray2)) DEALLOCATE(temparray2)
    ENDIF

    !Set variables for update function wAR  
    IF(iwar>0)THEN
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      IF(.NOT.ALLOCATED(temparray2)) ALLOCATE(temparray2(nrow))
      temparray = 0
      temparray2 = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iwar))==1)THEN
          nwarsub = nwarsub + 1
          temparray(nwarsub) = xi(irow,iindex(isub))
          temparray2(nwarsub) = xr(irow,rindex(iqarar))
        ENDIF
      ENDDO
      IF(nwarsub>0)THEN
        warupdsub(1:nwarsub) = temparray(1:nwarsub)
        warupdar(1:nwarsub) = temparray2(1:nwarsub) !Same as for qAR
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
      IF(ALLOCATED(temparray2)) DEALLOCATE(temparray2)
    ENDIF

    !Set variables for update function tpcorr
    IF(itpcorr>0)THEN      !Save subid/value for all stations to be updated
      tpcorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tpcorrarr(1:nrow,2) = xr(1:nrow,rindex(itpcorr))
    ENDIF
    !Set variables for update function tncorr
    IF(itncorr>0)THEN      !Save subid/value for all stations to be updated
      tncorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tncorrarr(1:nrow,2) = xr(1:nrow,rindex(itncorr))
    ENDIF
    !Set variables for update function tploccorr
    IF(itploccorr>0)THEN      !Save subid/value for all stations to be updated
      tploccorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tploccorrarr(1:nrow,2) = xr(1:nrow,rindex(itploccorr))
    ENDIF
    !Set variables for update function tnloccorr
    IF(itnloccorr>0)THEN      !Save subid/value for all stations to be updated
      tnloccorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tnloccorrarr(1:nrow,2) = xr(1:nrow,rindex(itnloccorr))
    ENDIF

    !Set variables for update function wendupd      
    IF(iwendupd>0)THEN      !Save subid for stations to be updated
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      temparray = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iwendupd))==1)THEN
          nwsub = nwsub + 1
          temparray(nwsub) = xi(irow,iindex(isub))
        ENDIF
      ENDDO
      IF(nwsub>0)THEN
        wendsub(1:nwsub)=temparray(1:nwsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
    ENDIF

    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE read_update_data

  !>Calculate the coupling between subbasins by index from the existing
  !>coupling by subid
  !>
  !>\b Consequences Module modvar variables path, branchdata and branchindex 
  !>may be allocated and set. Module modvar variables pathsubid and branchsubid 
  !>is deallocated. 
  !--------------------------------------------------------------------
  SUBROUTINE calculate_path(ns) 

    USE MODVAR, ONLY : path,        &   !OUT
                       pathsubid,   &
                       branchdata,  &   !OUT
                       branchsubid, &
                       branchindex      !OUT

    !Argument declarations
    INTEGER, INTENT(IN) :: ns                   !<Number of subbasins (of submodel)
    
    !Local variables
    INTEGER i
    INTEGER dim   !size of branchdata table
    INTEGER mainflow(ns)
    INTEGER grwflow(ns)
    INTEGER,ALLOCATABLE :: source(:),branch(:)

    !>\b Algorithm \n
    !>Set main path and groundwater flow path
    IF(.NOT.ALLOCATED(path)) ALLOCATE(path(ns))
    mainflow = pathsubid%main
    grwflow = pathsubid%grw1
    path%grwtolake = pathsubid%grwtolake
    path%aquid     = pathsubid%aquid
    path%rechargebasin = pathsubid%rechargebasin
    path%recievefraction = pathsubid%recievefraction
    DO i = 1,ns
      path(i)%main   = get_subid_index(mainflow(i),ns)      !uses basin%subid!
      path(i)%grw1   = get_subid_index(grwflow(i),ns)
    ENDDO
    IF(ALLOCATED(pathsubid)) DEALLOCATE(pathsubid)

    !>Set branch path and data, if present
    dim = 0
    IF(ALLOCATED(branchsubid))THEN
      dim = SIZE(branchsubid)
      IF(.NOT.ALLOCATED(branchdata)) ALLOCATE(branchdata(dim))
      IF(.NOT.ALLOCATED(branchindex)) ALLOCATE(branchindex(ns))
      IF(.NOT.ALLOCATED(source)) ALLOCATE(source(dim))
      IF(.NOT.ALLOCATED(branch)) ALLOCATE(branch(dim))
      branchindex = 0
      source = branchsubid%source
      branch = branchsubid%branch
      
      !Copy branch flow data
      branchdata%mainpart   = branchsubid%mainpart
      branchdata%maxQ       = branchsubid%maxQ
      branchdata%minQ       = branchsubid%minQ
      branchdata%maxQbranch = branchsubid%maxQbranch
      
      !Calculate corresponding index
      DO i = 1,dim
        branchdata(i)%source = get_subid_index(source(i),ns)
        IF(branchdata(i)%source>0) branchindex(branchdata(i)%source) = i
        branchdata(i)%branch = get_subid_index(branch(i),ns)
      ENDDO
      IF(ALLOCATED(branchsubid)) DEALLOCATE(branchsubid)
      IF(ALLOCATED(source)) DEALLOCATE(source)
      IF(ALLOCATED(branch)) DEALLOCATE(branch)
      IF(SUM(branchindex)==0)THEN   !No branches in model
        IF(ALLOCATED(branchindex)) DEALLOCATE(branchindex)
        IF(ALLOCATED(branchdata)) DEALLOCATE(branchdata)
        WRITE(6,*) 'No branches in (sub-)model found'
      ENDIF
    ENDIF

  END SUBROUTINE calculate_path

  !>Reform the information about subbasins to arrays for submodel
  !-------------------------------------------------------------------
  SUBROUTINE reform_inputdata_for_submodel(nsfile,ns,indexarray)

    USE MODVAR, ONLY : basin,        &
                       classbasin,   &
                       pathsubid,    &
                       load,         &
                       glacierexist, &
                       wetland,      &
                       wetlandexist,     &
                       forcing,          &
                       basintype,        &
                       classbasintype,   &
                       pathtype,         &
                       npcloadtype,      &
                       wetlandtype,      &
                       forcingtype,      &
                       basinpar,         &
                       damindex,         & 
                       lakeindex,        & 
                       lakebasinindex,   & 
                       lakedataparindex, & 
                       glacierindex,     &
                       nclass,           &
                       numsubstances,    &
                       modparid,         &
                       m_bpar,           &
                       max_par,          &
                       allocate_basinvariables, &
                       xobsindex,    &
                       max_outvar,   &
                       max_pstype
    USE WORLDVAR, ONLY : qobsindex,    &
                         pobsindex,    &
                         tobsindex,    &
                         sfobsindex,   &
                         swobsindex,   &
                         uobsindex,   &
                         rhobsindex,   &
                         tminobsindex, &
                         tmaxobsindex

    !Argument declarations
    INTEGER, INTENT(IN) :: nsfile           !<number of subbasins in file
    INTEGER, INTENT(IN) :: ns               !<number of subbasins to be simulated
    INTEGER, INTENT(IN) :: indexarray(ns)   !<index for basemodel
    
    !Local variables 
    INTEGER i,k,s,nbasinpar
    TYPE(BASINTYPE),ALLOCATABLE :: temp_basin(:)               !help variable for basin
    TYPE(PATHTYPE),ALLOCATABLE  :: temp_pathsubid(:)           !help variable for pathsubid
    TYPE(NPCLOADTYPE),ALLOCATABLE  :: temp_load(:)             !help variable for load
    TYPE(CLASSBASINTYPE),ALLOCATABLE :: temp_classbasin(:,:)   !help variable for classbasin
    TYPE(WETLANDTYPE), ALLOCATABLE :: temp_wetland(:,:)        !help variable for wetland
    TYPE(FORCINGTYPE),ALLOCATABLE :: temp_forcing(:)           !help variable for forcing
    REAL,ALLOCATABLE    :: temp_basinpar(:,:)                  !help variable for basinpar
    INTEGER,ALLOCATABLE :: temp_lakeindex(:)                   !help variable for lakeindex/damindex
    INTEGER,ALLOCATABLE :: temp_glacierindex(:)                !help variable for glacerindex
    INTEGER,ALLOCATABLE :: temp_lakebasinindex(:)              !help variable for lakebasinindex 
    INTEGER,ALLOCATABLE :: temp_lakedataparindex(:,:)          !help variable for lakedataparindex
    INTEGER,ALLOCATABLE :: tmp(:),xtmp(:,:)

    !Copy subbasin information from GeoData to temporary arrays
    IF(.NOT.ALLOCATED(temp_basin)) ALLOCATE(temp_basin(nsfile))
    temp_basin      = basin
    DEALLOCATE(basin)
    IF(.NOT.ALLOCATED(temp_pathsubid)) ALLOCATE(temp_pathsubid(nsfile))
    temp_pathsubid  = pathsubid
    DEALLOCATE(pathsubid)
    IF(.NOT.ALLOCATED(temp_classbasin)) ALLOCATE(temp_classbasin(nsfile,nclass))
    temp_classbasin = classbasin
    DEALLOCATE(classbasin)
    ALLOCATE(temp_load(nsfile))
    DO i = 1,nsfile
      ALLOCATE(temp_load(i)%psvol(max_pstype))
      ALLOCATE(temp_load(i)%psload(max_pstype,numsubstances))
      DO k=1,max_pstype
        temp_load(i)%psvol(k) = load(i)%psvol(k)
        DO s=1,numsubstances
          temp_load(i)%psload(k,s) = load(i)%psload(k,s)
        ENDDO
      ENDDO  
    ENDDO
    temp_load%inwetdep    = load%inwetdep
    temp_load%indrydep(1) = load%indrydep(1)
    temp_load%indrydep(2) = load%indrydep(2)
    temp_load%indrydep(3) = load%indrydep(3)
    temp_load%volloc      = load%volloc
    temp_load%tpconcloc   = load%tpconcloc
    temp_load%tnconcloc   = load%tnconcloc
    temp_load%inpartloc   = load%inpartloc
    temp_load%sppartloc   = load%sppartloc
    temp_load%abstrvol(1) = load%abstrvol(1)
    temp_load%abstrvol(2) = load%abstrvol(2)
    DO i = 1,nsfile
      DEALLOCATE(load(i)%psvol)
      DEALLOCATE(load(i)%psload)
    ENDDO
    DEALLOCATE(load)
    IF(.NOT.ALLOCATED(temp_wetland)) ALLOCATE(temp_wetland(nsfile,2))
    IF(wetlandexist)THEN  
      temp_wetland = wetland
      DEALLOCATE(wetland)
    ENDIF

    !Allocate for and copy reform subbasin information
    CALL allocate_basinvariables(ns,numsubstances)
    basin(:)        = temp_basin(indexarray(:))
    DEALLOCATE(temp_basin)
    pathsubid(:)    = temp_pathsubid(indexarray(:))
    DEALLOCATE(temp_pathsubid)
    classbasin(:,:) = temp_classbasin(indexarray(:),:)
    DEALLOCATE(temp_classbasin)
    load(:)%inwetdep    = temp_load(indexarray(:))%inwetdep
    load(:)%indrydep(1) = temp_load(indexarray(:))%indrydep(1)
    load(:)%indrydep(2) = temp_load(indexarray(:))%indrydep(2)
    load(:)%indrydep(3) = temp_load(indexarray(:))%indrydep(3)
    load(:)%volloc      = temp_load(indexarray(:))%volloc
    load(:)%tpconcloc   = temp_load(indexarray(:))%tpconcloc
    load(:)%tnconcloc   = temp_load(indexarray(:))%tnconcloc
    load(:)%inpartloc   = temp_load(indexarray(:))%inpartloc
    load(:)%sppartloc   = temp_load(indexarray(:))%sppartloc
    load(:)%abstrvol(1) = temp_load(indexarray(:))%abstrvol(1)
    load(:)%abstrvol(2) = temp_load(indexarray(:))%abstrvol(2)
    DO i = 1,ns
      DO k=1,max_pstype
        load(i)%psvol(k) = temp_load(indexarray(i))%psvol(k)
        DO s=1,numsubstances
          load(i)%psload(k,s) = temp_load(indexarray(i))%psload(k,s)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,nsfile
      DEALLOCATE(temp_load(i)%psvol)
      DEALLOCATE(temp_load(i)%psload)
    ENDDO
    DEALLOCATE(temp_load)
    IF(wetlandexist)THEN  
      wetland(:,:)    = temp_wetland(indexarray(:),:)
      IF(SUM(wetland%area)==0)THEN
        wetlandexist = .FALSE.
        DEALLOCATE(wetland)
      ENDIF
    ELSE
      DEALLOCATE(wetland)
    ENDIF
    DEALLOCATE(temp_wetland)
    ! nregions shall not be changed!

    !Reform observation indexarray for submodel
    IF(ALLOCATED(pobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = pobsindex(indexarray(:))
      DEALLOCATE(pobsindex)
      CALL MOVE_ALLOC(tmp,pobsindex)
    ENDIF
    IF(ALLOCATED(tobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = tobsindex(indexarray(:))
      DEALLOCATE(tobsindex)
      CALL MOVE_ALLOC(tmp,tobsindex)
    ENDIF
    IF(ALLOCATED(qobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = qobsindex(indexarray(:))
      DEALLOCATE(qobsindex)
      CALL MOVE_ALLOC(tmp,qobsindex)
    ENDIF
    IF(ALLOCATED(xobsindex))THEN 
      IF(.NOT.ALLOCATED(xtmp)) ALLOCATE(xtmp(max_outvar,ns))
      xtmp = xobsindex(:,indexarray(:))
      DEALLOCATE(xobsindex)
      CALL MOVE_ALLOC(xtmp,xobsindex)
    ENDIF
    IF(ALLOCATED(sfobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = sfobsindex(indexarray(:))
      DEALLOCATE(sfobsindex)
      CALL MOVE_ALLOC(tmp,sfobsindex)
    ENDIF
    IF(ALLOCATED(swobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = swobsindex(indexarray(:))
      DEALLOCATE(swobsindex)
      CALL MOVE_ALLOC(tmp,swobsindex)
    ENDIF
    IF(ALLOCATED(uobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = uobsindex(indexarray(:))
      DEALLOCATE(uobsindex)
      CALL MOVE_ALLOC(tmp,uobsindex)
    ENDIF
    IF(ALLOCATED(rhobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = rhobsindex(indexarray(:))
      DEALLOCATE(rhobsindex)
      CALL MOVE_ALLOC(tmp,rhobsindex)
    ENDIF
    IF(ALLOCATED(tminobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = tminobsindex(indexarray(:))
      DEALLOCATE(tminobsindex)
      CALL MOVE_ALLOC(tmp,tminobsindex)
    ENDIF
    IF(ALLOCATED(tmaxobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = tmaxobsindex(indexarray(:))
      DEALLOCATE(tmaxobsindex)
      CALL MOVE_ALLOC(tmp,tmaxobsindex)
    ENDIF

    !Reform subbasin information about forcing
    IF(ALLOCATED(forcing))THEN
      ALLOCATE(temp_forcing(nsfile))
      temp_forcing = forcing
      DEALLOCATE(forcing)
      ALLOCATE(forcing(ns))
      forcing(:)   = temp_forcing(indexarray(:))
      DEALLOCATE(temp_forcing)
    ENDIF
    
    !Reform subbasin information from LakeData and DamData
    !Only index for lakes handled the lakedata is kept as is
    IF(ALLOCATED(lakeindex))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = lakeindex
      DEALLOCATE(lakeindex)
      ALLOCATE(lakeindex(ns))
      lakeindex(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF

    IF(ALLOCATED(lakebasinindex))THEN
      IF(.NOT.ALLOCATED(temp_lakebasinindex)) ALLOCATE(temp_lakebasinindex(nsfile))
      temp_lakebasinindex = lakebasinindex
      DEALLOCATE(lakebasinindex)
      ALLOCATE(lakebasinindex(ns))
      lakebasinindex(:) = temp_lakebasinindex(indexarray(:))
      DEALLOCATE(temp_lakebasinindex)
    ENDIF

    IF(ALLOCATED(lakedataparindex))THEN
      IF(.NOT.ALLOCATED(temp_lakedataparindex)) ALLOCATE(temp_lakedataparindex(nsfile,2))
      temp_lakedataparindex = lakedataparindex
      DEALLOCATE(lakedataparindex)
      ALLOCATE(lakedataparindex(ns,2))
      lakedataparindex(:,:) = temp_lakedataparindex(indexarray(:),:)
      DEALLOCATE(temp_lakedataparindex)
    ENDIF

    IF(ALLOCATED(damindex))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = damindex
      DEALLOCATE(damindex)
      ALLOCATE(damindex(ns))
      damindex(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF

    !IF(glacierexist)THEN
      IF(ALLOCATED(glacierindex))THEN
        IF(.NOT.ALLOCATED(temp_glacierindex)) ALLOCATE(temp_glacierindex(nsfile))
        temp_glacierindex    = glacierindex
        DEALLOCATE(glacierindex)
        ALLOCATE(glacierindex(ns))
        glacierindex(:)      = temp_glacierindex(indexarray(:))
        DEALLOCATE(temp_glacierindex)
      ENDIF
    !ENDIF

    !Model parameters reformation, only subbasin dependent parameters need to be changed
    nbasinpar = 0
    DO i = 1,max_par
      IF(modparid(i)%deptype==m_bpar)  nbasinpar = MAX(nbasinpar,modparid(i)%parno)
    ENDDO
    IF(.NOT.ALLOCATED(temp_basinpar)) ALLOCATE(temp_basinpar(nbasinpar,nsfile))
    temp_basinpar = basinpar
    IF(ALLOCATED(basinpar)) DEALLOCATE(basinpar)
    ALLOCATE(basinpar(nbasinpar,ns))
    basinpar = temp_basinpar(:,indexarray)
    IF(ALLOCATED(temp_basinpar)) DEALLOCATE(temp_basinpar)

  END SUBROUTINE reform_inputdata_for_submodel
  
  !>Finds the corresponding index for a subid
  !--------------------------------------------
  INTEGER FUNCTION get_subid_index(id,ns)

    USE MODVAR, ONLY : basin

    !Argument declarations
    INTEGER, INTENT(IN)  :: id    !<subid
    INTEGER, INTENT(IN)  :: ns    !<number of subbasins
    
    !Local variables
    INTEGER j

    !>\b Algorithm \n
    !>Loop through subbasin id: to find the index of the current one.
    get_subid_index = 0
    DO j=1,ns
      IF(basin(j)%subid==id)THEN
        get_subid_index = j
        EXIT
      ENDIF
    ENDDO

  END FUNCTION get_subid_index

  !>Reads information about crops from file
  !>
  !>\b Consequences Module modvar variables cropdata, cropirrdata 
  !>and cropdataindex is allocated and set.
  !----------------------------------------------------------
  SUBROUTINE load_cropdata(dir,infile,n,n_Result) 

    USE WORLDVAR, ONLY : fileunit_temp, &
                         i_str,i_intg,i_real
    USE MODVAR, ONLY : cropindex,      &  !OUT
                       cropdata,       &  !OUT
                       cropirrdata,    &  !OUT
                       numsubstances,  &
                       set_cropdataindex, &
                       initiate_cropdata, &
                       i_in,i_on,i_sp,i_pp
    
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir           !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile        !<Name of file to be read (CropData.txt)
    INTEGER, INTENT(OUT) :: n                      !<Number of rows in file
    INTEGER, INTENT(OUT) :: n_Result               !<Error status
    
    !Local parameters
    INTEGER, PARAMETER :: nskip = 1   !heading on row 1
    
    !Local variables
    INTEGER i
    INTEGER mcols                      !Actual number of coulmns in file
    INTEGER sumirr
    INTEGER ncols                      !Number of columns in file
    INTEGER, ALLOCATABLE :: code(:)   !Data type code
    INTEGER, ALLOCATABLE :: iindex(:) !Index integer data from file
    INTEGER, ALLOCATABLE :: rindex(:) !Index real data from file
    INTEGER, ALLOCATABLE :: xi(:,:)   !integer data from file
    INTEGER, ALLOCATABLE :: cropregion(:) !CropData region read from file
    INTEGER, ALLOCATABLE :: cropid(:) !CropData cropid read from file
    REAL, ALLOCATABLE    :: xr(:,:)   !real data from file
    CHARACTER (LEN=10),ALLOCATABLE :: colstr(:)   !Content string of data from file
    LOGICAL calcNP                     !Tells if NP is modelled
    LOGICAL fileexist                  !CropData.txt file exist?  !CP111201

    !>\b Algorithm \n
    n_Result = 0

    !NP simulation?
    calcNP = .FALSE.
    DO i = 1, numsubstances
      IF(i==i_in .OR. i==i_on .OR. i==i_sp .OR. i==i_pp)    calcNP = .TRUE.
    ENDDO

    !>Check if file exist and if it is required
    INQUIRE(FILE=TRIM(dir)//infile,EXIST=fileexist)
    IF(.NOT.fileexist)THEN
      !N / P modelling require the file
      IF(.NOT.calcNP)THEN
        RETURN
      ELSE
        WRITE(6,*) ' ERROR: Missing CropData.txt for NP simulation'  
        n_Result = 1
        RETURN
      ENDIF
    ENDIF

    !>Count number of crops and data columns
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,nskip,n,n_Result)
    IF(n_Result/=0)RETURN
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,nskip,ncols,n_Result)  !CP120522
    IF(n_Result/=0)RETURN
    IF(.NOT.ALLOCATED(colstr)) ALLOCATE(colstr(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(n,ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(n,ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))

    !>Allocate and initiate cropdata variables
    IF(.NOT.ALLOCATED(cropregion)) ALLOCATE(cropregion(n))   
    IF(.NOT.ALLOCATED(cropid)) ALLOCATE(cropid(n))   
    IF(.NOT.ALLOCATED(cropdata)) ALLOCATE(cropdata(n))   
    IF(.NOT.ALLOCATED(cropirrdata)) ALLOCATE(cropirrdata(n))
    CALL initiate_cropdata()

    !>Read column headings from file
    OPEN(UNIT=fileunit_temp,FILE=TRIM(dir)//infile,STATUS='old',ACTION='read')
    CALL read_column_headings(fileunit_temp,ncols,colstr,mcols,n_Result)
    IF(n_Result.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !>Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
       IF(colstr(i)(1:10)=='cropid    ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='reg       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fn1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fn2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fp1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fp2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mn1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mn2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mp1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mp2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fday1     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fday2     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='mday1     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='mday2     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fdown1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fdown2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mdown1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mdown2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resn      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resp      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resc      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resday    ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='resdown   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resfast   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up3       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='bd1       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd2       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd3       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd4       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd5       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='ccmax1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='ccmax2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gcmax1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gcmax2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbini    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbmid    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbend    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='dlref     ')   code(i) = i_real
       IF(colstr(i)(1:10)=='plantday  ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthini ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthdev ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthmid ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthlate')   code(i) = i_intg
       IF(colstr(i)(1:10)=='imm_start ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='imm_end   ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='upupper   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='pnupr     ')   code(i) = i_real
       IF(colstr(i)(1:10)=='daylength ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gddsow    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='basetemp  ')   code(i) = i_real
       IF(colstr(i)(1:10)=='firstday  ')   code(i) = i_intg
    ENDDO

    !>Read the data of the crops to matrix, first column is string which is skipped
    CALL read_basindata5(fileunit_temp,infile,ncols,n,ncols,code,rindex,iindex,xi,xr) 
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(infile)

    !>Set cropdata variables according to matrix with read data
    DO i = 1,mcols
      IF(colstr(i)(1:10)=='cropid    ')   cropid(1:n)               = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='reg       ')   cropregion(1:n)           = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fn1       ')   cropdata(1:n)%fertnamount1 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fn2       ')   cropdata(1:n)%fertnamount2 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fp1       ')   cropdata(1:n)%fertpamount1 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fp2       ')   cropdata(1:n)%fertpamount2 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mn1       ')   cropdata(1:n)%mannamount1  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mn2       ')   cropdata(1:n)%mannamount2  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mp1       ')   cropdata(1:n)%manpamount1  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mp2       ')   cropdata(1:n)%manpamount2  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fday1     ')   cropdata(1:n)%fertday1     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fday2     ')   cropdata(1:n)%fertday2     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='mday1     ')   cropdata(1:n)%manday1      = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='mday2     ')   cropdata(1:n)%manday2      = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fdown1    ')   cropdata(1:n)%fertdown1    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='fdown2    ')   cropdata(1:n)%fertdown2    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='mdown1    ')   cropdata(1:n)%mandown1     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='mdown2    ')   cropdata(1:n)%mandown2     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='resn      ')   cropdata(1:n)%resnamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resp      ')   cropdata(1:n)%respamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resc      ')   cropdata(1:n)%rescamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resday    ')   cropdata(1:n)%resdayno     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='resdown   ')   cropdata(1:n)%resdown      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='resfast   ')   cropdata(1:n)%resfast      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up1       ')   cropdata(1:n)%uptake1      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up2       ')   cropdata(1:n)%uptake2      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up3       ')   cropdata(1:n)%uptake3      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='bd1       ')   cropdata(1:n)%baredayno1   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd2       ')   cropdata(1:n)%baredayno2   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd3       ')   cropdata(1:n)%baredayno3   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd4       ')   cropdata(1:n)%baredayno4   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd5       ')   cropdata(1:n)%baredayno5   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='ccmax1    ')   cropdata(1:n)%ccmax1       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='ccmax2    ')   cropdata(1:n)%ccmax2       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gcmax1    ')   cropdata(1:n)%gcmax1       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gcmax2    ')   cropdata(1:n)%gcmax2       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbini    ')   cropirrdata(1:n)%kcbini    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbmid    ')   cropirrdata(1:n)%kcbmid    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbend    ')   cropirrdata(1:n)%kcbend    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='dlref     ')   cropirrdata(1:n)%dlref     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='plantday  ')   cropirrdata(1:n)%plantingdayno   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthini ')   cropirrdata(1:n)%lengthini       = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthdev ')   cropirrdata(1:n)%lengthdev       = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthmid ')   cropirrdata(1:n)%lengthmid       = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthlate')   cropirrdata(1:n)%lengthlate      = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='imm_start ')   cropirrdata(1:n)%imm_start       = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='imm_end   ')   cropirrdata(1:n)%imm_end         = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='upupper   ')   cropdata(1:n)%uptakeupper        = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='pnupr     ')   cropdata(1:n)%PNuptakeRatio      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='daylength ')   cropdata(1:n)%daylength    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gddsow    ')   cropdata(1:n)%gddsow       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='basetemp  ')   cropdata(1:n)%basetemp     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='firstday  ')   cropdata(1:n)%firstday     = xi(1:n,iindex(i))
    ENDDO

    !>Check if cropdata is needed?
    IF(.NOT.calcNP) DEALLOCATE(cropdata)
    sumirr = SUM(cropirrdata(1:n)%plantingdayno)
    IF(sumirr==0) DEALLOCATE(cropirrdata)

    IF(calcNP .OR. (sumirr>0))THEN
      !>Calculate index variable to find crops
      CALL set_cropdataindex(n,cropid,cropregion,cropindex)
    ENDIF

    !Deallocate local arrays
    DEALLOCATE(xr,xi,iindex,rindex,code,colstr,cropid,cropregion)
    
  END SUBROUTINE load_cropdata
  
  !>Get command line argument: infodir and simsequence
  !>Or read them from filedir.txt (old variant)
  !-----------------------------------------------------------------------
  SUBROUTINE get_hyss_arguments(dir,iseq)
    
    USE WORLDVAR, ONLY : fileunit_temp,   &
         maxcharpath

    !Argument declarations
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: dir  !<Directory for information about simulation (infodir)
    INTEGER, INTENT(OUT) :: iseq                    !<Sequence number (simsequence)
    
    !Local variables
    INTEGER i       !Argument number
    INTEGER narg    !Number of command line argument
    INTEGER st      !Status
    INTEGER pos,oldpos  !Position on line
    LOGICAL nostr
    CHARACTER(LEN=maxcharpath) argument,argi  !Command line argument i
    CHARACTER(LEN=500) line

    !>\b Algorithm \n
    iseq = 0
    pos = 1
    dir=''
    narg = COMMAND_ARGUMENT_COUNT()
    !>If no command line argument, read information from filedir.txt
    IF(narg==0)THEN
      OPEN(fileunit_temp,FILE = 'filedir.txt',STATUS = 'old',ACTION='read',ERR=200)
      READ(fileunit_temp,*,END=300,ERR=300) line
      DO
        nostr = .FALSE.
        oldpos = pos
        CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err',.TRUE.)
        IF(nostr) EXIT
        argi=TRIM(argument)
        SELECT CASE(argi)
        CASE('-sequence','-s')
          CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err')
          READ(argument,*,END=101,ERR=101) iseq
        CASE('-infodir','-i')
          CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err',.TRUE.)
          CALL GET_COMMAND_ARGUMENT(i,argument,STATUS=st)
          IF(.NOT.nostr)THEN
            dir=TRIM(argument)
          ELSE
            WRITE(*,*) 'ERROR reading infodir from filedir.txt'
          ENDIF
        CASE DEFAULT
          IF(oldpos==1)THEN
            dir=argi
          ELSE
            WRITE(*,*) 'ERROR: No file path found in filedir.txt.'
            STOP 1
          ENDIF
        END SELECT
      END DO
       CLOSE(fileunit_temp)
       RETURN
    ENDIF
    
    !>Get command line argument
    i = 1
    DO WHILE(i<=narg)
       CALL GET_COMMAND_ARGUMENT(i,argument)
       argi=TRIM(argument)
       SELECT CASE(argi)
       CASE('-sequence','-s')
          i = i + 1
          CALL GET_COMMAND_ARGUMENT(i,argument)
          READ(argument,*,END=100,ERR=100) iseq
       CASE('-infodir','-i')
          i = i + 1
          CALL GET_COMMAND_ARGUMENT(i,argument,STATUS=st)
          IF(st==0)THEN
             dir=TRIM(argument)
          ELSE
             WRITE(*,*) 'ERROR reading infodir command line argument'
          ENDIF
       CASE DEFAULT
          IF(i==1)THEN
             dir=argi
          ELSE
             WRITE(*,*) 'ERROR in command line argument'
             STOP 1
          ENDIF
       END SELECT
       i = i + 1
    END DO
    RETURN

    !Error handling
100 WRITE(*,*) 'ERROR reading sequence command line argument'      
    STOP 1
101 WRITE(*,*) 'ERROR reading sequence from filedir.txt'      
    STOP 1
200 WRITE(*,*) 'ERROR no command line argument and no filedir.txt found' 
    STOP 1
300 WRITE(*,*) 'ERROR reading filedir.txt' 
    STOP 1

  END SUBROUTINE get_hyss_arguments

  !\brief Reads file with information about the current simulation 
  !!
  !!Format of file: first on line is a coded string followed by the value 
  !!The string "!!" is used for comment rows and are not read
  !!
  !>\b Consequences The following worldvar module variables is set: readdaily,
  !> readobsid, outvarallbasin, writeload, readmatlab, writematlab, thirtydaymonth,
  !> ncrit, readformat, simsubmodel, resultseq, temp_outbasins, steplen, 
  !> readsfobs, readswobs, readwind, readhumid, readtminmaxobs, and checkdata. The following modvar 
  !> module variables is set: i_in,i_on,i_sp,i_pp,i_t1,i_t2,i_t3,i_oc, 
  !> conductN,conductP,conductC,conductT,conductload,conductregest and modeloption.
  !--------------------------------------------------------------------
  SUBROUTINE load_coded_info(dir,n_Result,date1,date2,date3,skip,step,nsubst,     &
                             stateinput,somax,sonum,bstateoutdate,pstateoutdate,  &
                             mdir,rdir, nsubout,                          &
                             locupall,locupnone,locupnone_qar,locupnone_war,wupall,wupnone,wupdname) 

    USE WORLDVAR, ONLY : fileunit_temp,   &
                         maxcharpath,     &
                         maxoutbasins,    &
                         maxcrit,         &
                         comment_str,     &
                         set_calibration, &
                         set_outvar,      &
                         set_maxmap,      &
                         readdaily,       &   !OUT
                         readobsid,       &   !OUT
                         outvarallbasin,  &   !OUT
                         writeload,       &   !OUT
                         readmatlab,      &   !OUT
                         writematlab,     &   !OUT
                         thirtydaymonth,  &   !OUT
                         ncrit,           &   !OUT
                         readformat,      &   !OUT
                         simsubmodel,     &   !OUT
                         simsequence,     &
                         resultseq,       &   !OUT
                         temp_outbasins,  &   !OUT
                         i_d,             &
                         steplen,         &   !(OUT)
                         readsfobs,       &   !OUT
                         readswobs,       &   !OUT
                         readwind,        &   !OUT
                         readhumid,       &   !OUT
                         readtminmaxobs,  &   !OUT
                         checkdata            !OUT
    USE MODVAR, ONLY : i_in,i_on,i_sp,i_pp, &   !OUT
                       i_t1,i_t2,i_t3,i_oc, &   !OUT
                       i_sm,i_gm,i_rn,i_li,i_ri,i_si,i_sr,i_dr, & !OUT
                       doupdate,       &
                       conductN,conductP,conductC,conductT,conductWS, & !OUT
                       conductxoms, & !OUT
                       conductwb,   & !OUT
                       conductload,    &  !OUT
                       conductregest,  &  !OUT
                       soiliniwet,     &
                       max_outvar,     &   
                       i_quseobs,      &
                       i_qar,i_war,    &
                       i_wendupd,      &
                       i_tpcorr,       &
                       i_tncorr,       &
                       i_tploccorr,    &
                       i_tnloccorr,    &
                       irrunlimited,   &
                       modeloption,    &    !OUT
                       p_snowfall,     &   
                       p_snowmelt,     &
                       p_snowevap,     &
                       p_lakeriverice, & 
                       p_petmodel,     &
                       p_deepgroundwater, &
                       p_swtemperature, &
                       p_growthstart, &
                       num_modelprocess, &
                       p_infiltration, &  !for CryoProcesses, infilt in frozen soil                           
                       p_frozenrunoff

    !Argument declarations
    CHARACTER(LEN=maxcharpath), INTENT(IN) :: dir             !<File directory
    INTEGER, INTENT(OUT) :: n_Result                !<Status of subroutine
    TYPE(DateType), INTENT(OUT) :: date1            !<Begin date of simulation
    TYPE(DateType), INTENT(OUT) :: date2            !<End date of simulation
    TYPE(DateType), INTENT(OUT) :: date3            !<Begin date for criteria/output date
    INTEGER, INTENT(OUT) :: skip                    !<Length of warm-up period
    INTEGER, INTENT(OUT) :: step                    !<Length of simulation period
    INTEGER, INTENT(OUT) :: nsubst                  !<Number of substances simulated
    LOGICAL, INTENT(OUT) :: stateinput              !<Code for reading input state
    INTEGER, INTENT(IN)  :: somax                   !<Dimension dates for output state
    INTEGER, INTENT(OUT) :: sonum                   !<Number of dates for output state
    TYPE(DateType), INTENT(OUT) :: bstateoutdate(somax) !<Dates for writing state variables
    TYPE(DateType), INTENT(OUT) :: pstateoutdate(somax) !<Dates previous to date for writing 
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: mdir           !<File directory model
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: rdir           !<File directory result
    INTEGER, INTENT(OUT) :: nsubout                 !<Number of subbasins output in basinfiles
    LOGICAL, INTENT(OUT) :: locupall                !<Flag for Q-update for all stations
    LOGICAL, INTENT(OUT) :: locupnone               !<Flag for Q-update on no stations
    LOGICAL, INTENT(OUT) :: locupnone_qar           !<Flag for Q-update on no stations
    LOGICAL, INTENT(OUT) :: locupnone_war           !<Flag for WAR-update on no stations
    LOGICAL, INTENT(OUT) :: wupall                  !<Flag for W-update all stations
    LOGICAL, INTENT(OUT) :: wupnone                 !<Flag for W-update nostations
    CHARACTER(LEN=4), INTENT(OUT) :: wupdname       !<Variable name for wendupd-function
    
    !Local parameters
    INTEGER, PARAMETER :: linelen = 18000      
    CHARACTER (LEN=2), PARAMETER :: errstr = '##'  
    
    !Local variables
    LOGICAL error
    LOGICAL onoff
    LOGICAL calibration
    LOGICAL nostrfound,notimefound
    LOGICAL exitbigloop
    LOGICAL cyclebigloop
    INTEGER io
    INTEGER temp
    INTEGER maxdim
    INTEGER intvalue
    CHARACTER(LEN=3), ALLOCATABLE :: c(:)
    INTEGER, ALLOCATABLE :: critc(:)
    INTEGER, ALLOCATABLE :: critr(:)
    INTEGER, ALLOCATABLE :: flowc(:)
    INTEGER, ALLOCATABLE :: flowr(:)
    REAL, ALLOCATABLE :: critweight(:)
    REAL, ALLOCATABLE :: critpar(:)
    INTEGER istate
    INTEGER ioutbasin
    INTEGER readbasins(maxoutbasins)
    INTEGER ivarbasinoutput
    INTEGER ivarmapoutput
    INTEGER ivartimeoutput
    INTEGER ivarindex
    INTEGER varindex,flowtype
    INTEGER basinoutputperiod 
    INTEGER mapoutputperiod   
    INTEGER timeoutputperiod  
    INTEGER basinoutputdecimal
    INTEGER mapoutputdecimal 
    INTEGER timeoutputdecimal
    INTEGER varbasinoutput(max_outvar,2)
    INTEGER varmapoutput(max_outvar,2)
    INTEGER vartimeoutput(max_outvar,2)
    INTEGER icrit,critperiod,critlimit
    INTEGER tstep,indx !Time step length [stepunit] and string index (used with stepstr)
    INTEGER linepos
    INTEGER icheck
    CHARACTER (LEN=linelen) :: line
    CHARACTER (LEN=16) :: code,code2, code3
    CHARACTER (LEN=10) :: strvalue
    CHARACTER (LEN=16) :: strdate,strdate1,strdate2,strdate3
    CHARACTER (LEN=16), ALLOCATABLE :: strstates(:)
    CHARACTER (LEN=4)  :: varstr
    CHARACTER (LEN=3)  :: yesnostr
    CHARACTER (LEN=3)  :: critstr
    CHARACTER (LEN=maxcharpath) :: filedir
    CHARACTER (LEN=maxcharpath+8) :: filename
    CHARACTER (LEN=4)  :: stepunit      !Time step unit; mo,d,h or min
    CHARACTER (LEN=10) :: stepstr       !Temporary storage for time step string
    LOGICAL, ALLOCATABLE :: critcond(:) !flag that this criteria should only be used conditionally /DG
    REAL, ALLOCATABLE :: critthres(:)   !acceptance criteria for conditional criteria /DG
    character (LEN=3) str3
    
    n_Result = 0
    error = .FALSE.
    ncrit = maxcrit
    ALLOCATE(c(ncrit),critweight(ncrit),critc(ncrit),flowc(ncrit),    &
             critr(ncrit),flowr(ncrit),critpar(ncrit), &
             critcond(ncrit),critthres(ncrit)) !additional info needed for conditional acceptance
    ncrit = 0
    !initialize the conditional criteria and threshold
    critcond = .FALSE.
    critthres = 0.

    !Open file      
    filename=TRIM(dir)//'info.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')

    !Default values
    calibration    = .FALSE.
    stateinput     = .FALSE.
    writematlab    = .FALSE.
    readmatlab     = .FALSE.
    thirtydaymonth = .FALSE.
    outvarallbasin = .FALSE.
    exitbigloop    = .FALSE.
    simsubmodel    = .FALSE.
    doupdate       = .FALSE.
    locupall       = .FALSE.
    locupnone      = .FALSE.
    locupnone_qar  = .FALSE.
    locupnone_war  = .FALSE.
    wupall         = .FALSE.
    wupnone        = .FALSE.
    readdaily      = .FALSE.
    writeload      = .FALSE.
    readobsid      = .TRUE. 
    irrunlimited   = .FALSE.                !Unlimited irrigation, defaults to False
    soiliniwet     = .FALSE.
    conductregest  = .FALSE.
    resultseq      = .FALSE.  
    readsfobs      = .FALSE.
    readswobs      = .FALSE.
    readwind       = .FALSE.
    readhumid      = .FALSE.
    readtminmaxobs = .FALSE.
    conductxoms    = .FALSE.
    conductwb      = .FALSE.
    ALLOCATE(modeloption(num_modelprocess))
    modeloption = 0
    checkdata      = .FALSE.
    conductN       = .FALSE.
    conductP       = .FALSE.
    conductC       = .FALSE.
    conductT       = .FALSE.
    conductload    = .FALSE.
    IF(simsequence>0) resultseq = .TRUE.  
    date3          = DateType(0,0,0,0,0)
    bstateoutdate  = DateType(0,0,0,0,0)
    ncrit          = 0
    readformat     = 0         !ascii
    nsubst         = 0
    i_in=0;i_on=0;i_sp=0;i_pp=0;i_t1=0;i_t2=0;i_t3=0;i_oc=0
    i_sm=0;i_gm=0;i_rn=0;i_li=0;i_ri=0
    mdir           = dir
    rdir           = dir
    ioutbasin          = 0
    ivarbasinoutput    = 0
    ivarmapoutput      = 0
    ivartimeoutput     = 0
    varbasinoutput     = 0
    varmapoutput       = 0
    vartimeoutput      = 0
    basinoutputperiod  = 1   !daily
    mapoutputperiod    = 5   !mean over simulation period
    timeoutputperiod   = 1   !daily
    basinoutputdecimal = 1
    mapoutputdecimal   = 1
    timeoutputdecimal  = 3
    critperiod         = i_d  !for calvarper
    critlimit          = 3    !for calvarlim
    wupdname=''
    sonum              = 0
    tstep              = 1    ! Default time step length is 1 ..
    stepunit           = 'd'  ! .. day

    !Read content of file      
    DO
      READ(fileunit_temp,600,END=200,ERR=800,IOSTAT=io) line
      linepos = 1
      CALL read_next_codestr_on_line(linelen,16,linepos,line,code,nostrfound,errstr)    !Read code
      IF(code(1:2)==errstr)THEN
        error = .TRUE.
        code=line(1:16)
        EXIT
      ENDIF
      IF(code(1:2)==comment_str) CYCLE
      IF(nostrfound) CYCLE  !empty row
      cyclebigloop = .FALSE.

      !Read values:
      !Read date   
      IF(code(1:5)=='bdate' .OR. &
         code(1:5)=='cdate' .OR. &
         code(1:5)=='edate' .OR. &
         code(1:12)=='outstatedate')THEN
        CALL read_next_date_on_line(linelen,16,linepos,line,strdate,nostrfound,notimefound,errstr)
        IF(nostrfound.OR.strdate(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        SELECT CASE(code)
        CASE ('bdate')
          CALL string_convert_to_DateType(strdate,date1)
        CASE ('cdate')
          CALL string_convert_to_DateType(strdate,date3)
        CASE ('edate')
          CALL string_convert_to_DateType(strdate,date2)
        CASE ('outstatedate')
          DO
            sonum = sonum + 1
            IF(sonum>somax) THEN
              WRITE(6,*) 'ERROR: Too many outstatedate i info.txt'
              n_Result = 1
              RETURN
            ENDIF
            CALL string_convert_to_DateType(strdate,bstateoutdate(sonum))
            CALL read_next_date_on_line(linelen,16,linepos,line,strdate,nostrfound,notimefound,errstr)
            IF(nostrfound.OR.strdate(1:2)==errstr)THEN
              EXIT
            ENDIF
          ENDDO             
        END SELECT
      ENDIF

      !Read time step, if present
      IF(code(1:10)=='steplength')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,stepstr,nostrfound,errstr,.TRUE.)
        IF(nostrfound)THEN   ! No error here if not defined, only using default values
          WRITE(6,'(A60,I2,A3)') &
                  'WARNING: No value for steplength found, using default values',tstep,stepunit
        ELSE
          indx=SCAN(stepstr,'dhm') ! Looking for unit: d,h or min
          IF(indx.GT.0)THEN
            READ(stepstr(1:indx-1),'(I4)') tstep
            stepunit=stepstr(indx:LEN(stepstr))
!                WRITE(6, '(A17,I4,A3)') 'Found steplength ',tstep,stepunit
          ELSE
            WRITE(6,'(A65,I2,A3)') &
                   'WARNING: No valid unit for steplength found, using default values',tstep,stepunit
          ENDIF
        ENDIF
      ENDIF

      !Read yes or no coded variables
      IF(code(1:11)=='calibration' .OR. &
         code(1:9) =='readdaily'   .OR. &
         code(1:9) =='printload'   .OR. &
         code(1:8) =='submodel'    .OR. &
         code(1:7) =='instate'     .OR. &
         code(1:8) =='resseqnr'    .OR. &
         code(1:12)=='irrunlimited'.OR. &
         code(1:9) =='readsfobs'   .OR. &
         code(1:9) =='readswobs'   .OR. &
         code(1:8) =='readwind'    .OR. &
         code(1:9) =='readhumid'   .OR. &
         code(1:9) =='readobsid'   .OR. &
         code(1:14)=='readtminmaxobs' .OR. &
         code(1:13)=='printwaterbal' .OR. &
         code(1:10)=='soiliniwet' .OR. &
         code(1:11)=='regestimate' .OR. &
         code(1:13)=='readxomsfiles')THEN
        CALL read_next_codestr_on_line(linelen,3,linepos,line,yesnostr,nostrfound,errstr)
        IF(nostrfound.OR.yesnostr(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(yesnostr(1:1)=='Y' .OR. yesnostr(1:1)=='y' .OR.    &
           yesnostr(1:1)=='J' .OR. yesnostr(1:1)=='j')THEN
          onoff = .TRUE.
        ELSEIF(yesnostr(1:1)=='N' .OR. yesnostr(1:1)=='n')THEN
          onoff = .FALSE.
        ELSE
          error = .TRUE.
          EXIT  
        ENDIF
        SELECT CASE(code)
        CASE ('calibration')
          calibration = onoff
        CASE ('readdaily')
          readdaily = onoff
        CASE ('printload')
          writeload = onoff
        CASE ('submodel')
          simsubmodel = onoff
        CASE ('instate')
          stateinput = onoff
        CASE ('resseqnr')
          resultseq = onoff
        CASE ('readobsid')
          readobsid = onoff
        CASE ('readsfobs')
          readsfobs = onoff
        CASE ('readswobs')
          readswobs = onoff
        CASE ('readwind')
          readwind = onoff
        CASE ('readhumid')
          readhumid = onoff
        CASE ('readtminmaxobs')
          readtminmaxobs = onoff
        CASE ('readxomsfiles')
          conductxoms = onoff
        CASE ('printwaterbal')
          conductwb = onoff
        CASE ('irrunlimited')
          irrunlimited = onoff
        CASE ('soiliniwet')
          soiliniwet = onoff
        CASE ('regestimate')
          conductregest = onoff
        END SELECT
      ENDIF

      !Read integer coded variables
      IF(code(1:10)=='readformat' .OR. &
         code(1:11)=='writeformat')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        READ(strvalue,*) intvalue
        IF(code(1:10)=='readformat')THEN
          IF(intvalue==1) readmatlab = .TRUE.     
          IF(intvalue==2) thirtydaymonth = .TRUE.
        ELSEIF(code(1:11)=='writeformat')THEN
          IF(intvalue==1) writematlab = .TRUE.
        ENDIF
      ENDIF

      !Read optional model structures coded with integers
      IF(code(1:13)=='snowfallmodel' .OR. &
         code(1:13)=='snowmeltmodel' .OR. &
         code(1:15)=='snowevaporation' .OR. &
         code(1:12)=='lakeriverice' .OR. &
         code(1:8)=='petmodel' .OR. &
         code(1:10)=='deepground' .OR.  &
         code(1:13)=='swtemperature' .OR. &
         code(1:16)=='growthstartmodel' .OR.  &
         code(1:12)=='infiltration')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        READ(strvalue,*) intvalue
        IF(code(1:13)=='snowfallmodel')THEN
          modeloption(p_snowfall) = intvalue
        ELSEIF(code(1:13)=='snowmeltmodel')THEN
          modeloption(p_snowmelt) = intvalue
        ELSEIF(code(1:15)=='snowevaporation')THEN
          modeloption(p_snowevap) = intvalue
        ELSEIF(code(1:8)=='petmodel')THEN
          modeloption(p_petmodel) = intvalue
        ELSEIF(code(1:12)=='lakeriverice')THEN
          modeloption(p_lakeriverice) = intvalue
          !Make sure T2 temperature model is switched on if lake and river ice model is used
          IF(i_t2.LE.0 .AND. modeloption(p_lakeriverice).GE.1)THEN
            nsubst=nsubst+1
            i_t2 = nsubst
            conductT = .TRUE.
          ENDIF
        ELSEIF(code(1:10)=='deepground')THEN
          modeloption(p_deepgroundwater) = intvalue
        ELSEIF(code(1:13)=='swtemperature')THEN
          modeloption(p_swtemperature) = intvalue
        ELSEIF(code(1:16)=='growthstartmodel')THEN
          modeloption(p_growthstart) = intvalue
        ELSEIF(code(1:12)=='infiltration')THEN
          modeloption(p_infiltration) = intvalue
        ELSEIF(code(1:12)=='frozenrunoff')THEN
          modeloption(p_frozenrunoff) = intvalue
        ENDIF
      ENDIF

      !Read substances to be simulated
      IF(code(1:9) =='substance')THEN
        DO
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(strvalue(1:2)==errstr)THEN
            error = .TRUE.
            exitbigloop = .TRUE.
            EXIT
          ENDIF
          IF(nostrfound)THEN
            cyclebigloop = .TRUE.
            EXIT
          ENDIF
          SELECT CASE(strvalue(1:1))
          CASE('N','n')
            i_in = nsubst+1
            i_on = nsubst+2
            nsubst=nsubst+2
            conductN = .TRUE.
          CASE('P','p')
            i_sp = nsubst+1
            i_pp = nsubst+2
            nsubst=nsubst+2
            conductP = .TRUE.
          CASE('C','c')
            nsubst=nsubst+1
            i_oc = nsubst
            conductC = .TRUE.
          CASE('T','t')
            SELECT CASE(strvalue(1:2))
            CASE('T1','t1')
              nsubst=nsubst+1
              i_t1 = nsubst
            CASE('T2','t2')
              !Make sure T2 is only initiated once (see lake and river ice setting above)
              IF(i_t2.LE.0)THEN
                nsubst=nsubst+1
                i_t2 = nsubst
                !Is T2 model relevant without lake and river ice model?
              ENDIF
            CASE('T3','t3')
              nsubst=nsubst+1
              i_t3 = nsubst
            END SELECT
            conductT = .TRUE.
          CASE('W','w')
            str3 = strvalue(1:3)
            CALL lower_case(str3,1)
            SELECT CASE(str3)
            CASE('wsm') !water from Snow
              conductWS = .TRUE.
            CASE('wgm') !water from Glacier
              conductWS = .TRUE.
            CASE('wrn') !water from Rain
              conductWS = .TRUE.
            CASE('wli') !water from Lake(ini)
              conductWS = .TRUE.
            CASE('wri') !water from River(ini)
              conductWS = .TRUE.
            CASE('wsi') !water from Soil(ini)
              conductWS = .TRUE.
            CASE('wsr') !water from Surface runoff
              conductWS = .TRUE.
            CASE('wdr') !water from Drainage
              conductWS = .TRUE.
            END SELECT
            IF(conductWS)THEN
              IF(i_sm<=0)THEN
                nsubst=nsubst+1
                i_sm = nsubst
              ENDIF
              IF(i_gm<=0)THEN
                nsubst=nsubst+1
                i_gm = nsubst
              ENDIF
              IF(i_rn<=0)THEN
                nsubst=nsubst+1
                i_rn = nsubst
              ENDIF
              IF(i_li<=0)THEN
                nsubst=nsubst+1
                i_li = nsubst
              ENDIF
              IF(i_ri<=0)THEN
                nsubst=nsubst+1
                i_ri = nsubst
              ENDIF
              IF(i_si<=0)THEN
                nsubst=nsubst+1
                i_si = nsubst  
              ENDIF
              IF(i_sr<=0)THEN
                nsubst=nsubst+1
                i_sr = nsubst  
              ENDIF
              IF(i_dr<=0)THEN
                nsubst=nsubst+1
                i_dr = nsubst  
              ENDIF
            ENDIF
          END SELECT
        ENDDO
        IF(exitbigloop)EXIT
        IF(cyclebigloop)CYCLE
      ENDIF

      !Read directory paths
      IF(code(1:8)=='modeldir' .OR. &
         code(1:9)=='resultdir')THEN
        CALL read_next_codestr_on_line(linelen,maxcharpath,linepos,line,filedir,nostrfound,errstr,.TRUE.)
        IF(nostrfound.OR.filedir(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(code(1:8)=='modeldir')  mdir = filedir
        IF(code(1:9)=='resultdir') rdir = filedir
      ENDIF

      !Read output information
      IF(code(1:11)=='basinoutput' .OR. &
         code(1:9) =='mapoutput'   .OR. &
         code(1:10)=='timeoutput') THEN
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(code2(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(nostrfound) CYCLE
        IF(code2(1:8)=='allbasin')THEN
          outvarallbasin = .TRUE.
        ENDIF
        IF(code2(1:10)=='meanperiod')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) intvalue
          IF(code(1:11)=='basinoutput')THEN
            basinoutputperiod = intvalue
          ELSEIF(code(1:9)=='mapoutput')THEN
            mapoutputperiod = intvalue
          ELSEIF(code(1:10)=='timeoutput')THEN
            timeoutputperiod = intvalue
          ENDIF
        ENDIF
        IF(code2(1:8)=='decimals')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) intvalue
          IF(code(1:11)=='basinoutput')THEN
            basinoutputdecimal = intvalue
          ELSEIF(code(1:9)=='mapoutput')THEN
            mapoutputdecimal = intvalue
          ELSEIF(code(1:10)=='timeoutput')THEN
            timeoutputdecimal = intvalue
          ENDIF
        ENDIF
        IF(code2(1:8)=='variable')THEN
          DO
            CALL read_next_codestr_on_line(linelen,4,linepos,line,varstr,nostrfound,errstr)
            IF(varstr(1:2)==errstr)THEN
              error = .TRUE.
              WRITE(6,*)  'ERROR: In variable for ',TRIM(code)
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(nostrfound)THEN
              cyclebigloop = .TRUE.
              EXIT
            ENDIF
            n_Result = find_variable_index(varstr,varindex,flowtype)
            IF(n_Result>0)THEN
              WRITE(6,*)  'ERROR: In variable index calculation: ', varstr,' for ',TRIM(code)
              error = .TRUE.
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(code(1:11)=='basinoutput')THEN
              ivarbasinoutput = ivarbasinoutput + 1
              IF(ivarbasinoutput>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many basinoutput variables, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              varbasinoutput(ivarbasinoutput,1) = varindex
              varbasinoutput(ivarbasinoutput,2) = flowtype
            ELSEIF(code(1:9)=='mapoutput')THEN
              ivarmapoutput = ivarmapoutput + 1
              IF(ivarmapoutput>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many mapoutput variable, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              varmapoutput(ivarmapoutput,1) = varindex
              varmapoutput(ivarmapoutput,2) = flowtype
            ELSEIF(code(1:10)=='timeoutput')THEN
              ivartimeoutput = ivartimeoutput + 1
              IF(ivartimeoutput>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many timeoutput variable, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              vartimeoutput(ivartimeoutput,1) = varindex
              vartimeoutput(ivartimeoutput,2) = flowtype
            ENDIF
          ENDDO
          IF(exitbigloop)EXIT
          IF(cyclebigloop)CYCLE
        ENDIF
        IF(code2(1:8)=='subbasin')THEN
          DO
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(strvalue(1:2)==errstr)THEN
              error = .TRUE.
              WRITE(6,*)  'ERROR: In subbasin for ',TRIM(code)
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(nostrfound)THEN
              cyclebigloop = .TRUE.
              EXIT
            ENDIF
            CALL convert_string_to_integer(10,strvalue,intvalue,n_Result) 
            IF(n_result/=0)THEN
              error = .TRUE.
              exitbigloop = .TRUE.
            ENDIF  
            !READ(strvalue,*) intvalue
            ioutbasin = ioutbasin + 1
            IF(ioutbasin>maxoutbasins)THEN
              WRITE(6,*) 'Warning: Number of "basinoutput" "subbasin"s larger than max',maxoutbasins
              WRITE(6,*) 'Warning: Excess subbasins will be ignored'
            ELSE
              readbasins(ioutbasin) = intvalue
            ENDIF
          ENDDO
          IF(exitbigloop)EXIT
          IF(cyclebigloop)CYCLE
        ENDIF
      ENDIF

      !Read update information
      IF(code(1:6)=='update') THEN
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(nostrfound.OR.code2(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF

        !Read information for updating with qobs
        IF(code2(1:7)=='quseobs') THEN
          doupdate(i_quseobs) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(code3(1:10)=='allstation')THEN
            locupall = .TRUE.
          ENDIF
          IF(code3(1:9)=='nostation')THEN
            locupnone = .TRUE.
          ENDIF

        ELSEIF(code2(1:3)=='qar') THEN
          doupdate(i_qar) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)          
          IF(code3(1:9)=='nostation')THEN
            locupnone_qar = .TRUE.
          ENDIF

          !Read information for updating with wstr
        ELSEIF(code2(1:7)=='wendupd') THEN
          doupdate(i_wendupd) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(nostrfound.OR.code3(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          wupdname = code3(1:4)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(code3(1:10)=='allstation')THEN
            wupall = .TRUE.
          ENDIF
          IF(code3(1:9)=='nostation')THEN
            wupnone = .TRUE.
          ENDIF
        ELSEIF(code2(1:3)=='war') THEN
          doupdate(i_war) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(nostrfound.OR.code3(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          wupdname = code3(1:4)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)          
          IF(code3(1:9)=='nostation')THEN
            locupnone_war = .TRUE.
          ENDIF

        !Read information for updating nitrogen and phosphorus concentrations
        ELSEIF(code2(1:6)=='tpcorr') THEN
          doupdate(i_tpcorr) = .TRUE. 
        ELSEIF(code2(1:6)=='tncorr') THEN
          doupdate(i_tncorr) = .TRUE. 
        ELSEIF(code2(1:9)=='tploccorr') THEN
          doupdate(i_tploccorr) = .TRUE. 
        ELSEIF(code2(1:9)=='tnloccorr') THEN
          doupdate(i_tnloccorr) = .TRUE. 
        ENDIF
      ENDIF

      !Read criteria information
      IF(code(1:4)=='crit')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(strvalue(1:10)=='meanperiod')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) critperiod
        ELSEIF(strvalue(1:9)=='datalimit')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) critlimit
        ELSE
          READ(strvalue,*) icrit
          IF(icrit>maxcrit)THEN
            WRITE(6,*)  'ERROR: Too many criteria, max: ',maxcrit
            error = .TRUE.
            EXIT
          ENDIF
          ncrit = MAX(ncrit,icrit)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
          IF(nostrfound.OR.code2(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          IF(code2(1:9)=='criterion')THEN
            CALL read_next_codestr_on_line(linelen,3,linepos,line,critstr,nostrfound,errstr)
            IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            c(icrit) = critstr
          ENDIF
          IF(code2(1:9)=='cvariable')THEN
            CALL read_next_codestr_on_line(linelen,4,linepos,line,varstr,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            n_Result = find_variable_index(varstr,varindex,flowtype)
            IF(n_Result>0)THEN
              WRITE(6,*)  'ERROR in variable index calculation: ', varstr,' for ',TRIM(code2)
              error = .TRUE.
              EXIT
            ENDIF
            critc(icrit) = varindex
            flowc(icrit) = flowtype
          ENDIF
          IF(code2(1:9)=='rvariable')THEN
            CALL read_next_codestr_on_line(linelen,4,linepos,line,varstr,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            n_Result = find_variable_index(varstr,varindex,flowtype)
            IF(n_Result>0)THEN
              WRITE(6,*)  'ERROR in variable index calculation: ', varstr,' for ',TRIM(code2)
              error = .TRUE.
              EXIT
            ENDIF
            critr(icrit) = varindex
            flowr(icrit) = flowtype
          ENDIF
          IF(code2(1:6)=='weight')THEN
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critweight(icrit)
          ENDIF
          IF(code2(1:9)=='parameter')THEN
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critpar(icrit)
          ENDIF
          IF(code2(1:11)=='conditional')THEN
            critcond(icrit) = .TRUE.
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critthres(icrit)
          ENDIF
        ENDIF
      ENDIF

      !Read check flag
      IF(code(1:5)=='check')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(strvalue(1:6)=='indata')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.varstr(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) icheck
          checkdata(icheck,1) = .TRUE.
        ENDIF
      ENDIF

    ENDDO
                
    IF(error)THEN
      WRITE(6,*) 'ERROR: reading info.txt (',TRIM(filename),'). Code=',code
      n_Result = 1
      CLOSE(fileunit_temp)
      RETURN
    ENDIF

200 CONTINUE
    ! Check consistency between optional model structure and optional forcing data input
    IF(.NOT.readsfobs .AND. modeloption(p_snowfall).EQ.1)THEN
      WRITE(6,*) 'ERROR: (info.txt) modeloption SNOWFALLMODEL can not be 1 if READSFOBS is N - please provide input file SFobs.txt and change option READSFOBS to Y!'
      n_Result = 1
      RETURN
    ENDIF
    IF(.NOT.readswobs .AND. .NOT.readtminmaxobs .AND. modeloption(p_petmodel).GE.3)THEN
      WRITE(6,*)  'WARNING: model option PETMODEL', modeloption(p_petmodel), 'requires shortwave radiation data (SWobs.txt) and/or Tmin and Tmax data. Model will run, but derives diurnal temperature range from clear sky turbidity, which might cause strange results....' 
    ENDIF
    IF(.NOT.readtminmaxobs .AND. modeloption(p_petmodel).GE.4)THEN
      WRITE(6,*)  'WARNING: model option PETMODEL', modeloption(p_petmodel), 'requires Tmin data to estimate relative humidity (assuming actual vapour pressure = saturation pressure at Tmin). Model will run, but derives diurnal temperature range from clear sky turbidity, which might cause strange results.... (btw RHUM data will be optional in future releases)' 
    ENDIF
    !Check compatible simulation settings
    IF((modeloption(p_lakeriverice)>=1.AND.i_t2==0) .OR. (modeloption(p_lakeriverice)==0.AND.i_t2>0))THEN
      WRITE(6,*)  'ERROR: lakerivericemodel requires simulation of water temperature (T2), and vise versa'
      n_Result = 1
      RETURN
    ENDIF
       
    !Some calculations
    temp = set_timestep_variables(tstep,stepunit,date1)    
    IF(date3%Year==0.AND.date3%Month==0.AND.date3%Day==0.AND.date3%Hour==0.AND.date3%Minute==0) date3=date1      
    skip = period_length(date1,date3)         !Calculate length of warmup period
    step = period_length(date1,AddDates(date2,steplen)) !Calculate length of simulation period
    temp = set_calibration(calibration,ncrit,c(1:ncrit),critperiod,critlimit,   &
         critweight(1:ncrit),critc(1:ncrit),critr(1:ncrit),   &
         flowc(1:ncrit),critpar(1:ncrit), &   !NOTE only flowc is saved, combination of different type of variables may go wrong
         critcond(1:ncrit),critthres(1:ncrit)) ! used for conditional acceptance in DEMC calibration
    maxdim = MAX(ivarbasinoutput,ivarmapoutput,ivartimeoutput)
    temp = set_outvar(1,max_outvar,maxdim,varbasinoutput(:,1),    &
         varbasinoutput(:,2),ivarbasinoutput,basinoutputperiod,    &
         basinoutputdecimal)
    temp = set_outvar(2,max_outvar,maxdim,varmapoutput(:,1),varmapoutput(:,2),  &
         ivarmapoutput,mapoutputperiod,mapoutputdecimal)
    temp = set_maxmap(step-skip+1,mapoutputperiod)
    temp = set_outvar(3,max_outvar,maxdim,vartimeoutput(:,1),   &
         vartimeoutput(:,2),ivartimeoutput,timeoutputperiod,       &
         timeoutputdecimal)
    nsubout = ioutbasin
    IF(ioutbasin>0)THEN      
      ALLOCATE(temp_outbasins(ioutbasin))
      temp_outbasins=readbasins(1:ioutbasin)
    ENDIF
    DO ivarindex = 1,sonum     !Previous time step for output of state
      pstateoutdate(ivarindex)=SubtractDates(bstateoutdate(ivarindex),steplen)
    ENDDO

    !Dates to formatted strings for hyss.log
    CALL format_date(date1,'yyyy-mm-dd HH:MM',strdate1)
    CALL format_date(date2,'yyyy-mm-dd HH:MM',strdate2)
    CALL format_date(date3,'yyyy-mm-dd HH:MM',strdate3)
    IF(sonum>0) ALLOCATE(strstates(sonum))
    DO istate=1,sonum
      CALL format_date(bstateoutdate(istate),'yyyy-mm-dd HH:MM',strstates(istate))
    END DO

    !Check consistency between optional model structure and time step (temporary not working)
    IF(.NOT.(steplen==datetype(0,0,1,0,0)))THEN
      IF(modeloption(p_petmodel).GE.1)THEN
        WRITE(6,*)  'ERROR: model option PETMODEL', modeloption(p_petmodel), 'does not work with other time step than day (to come in future release)' 
        n_Result = 1
        RETURN
      ENDIF  
      IF(i_t2.GE.1 .OR. modeloption(p_lakeriverice).GE.1)THEN
        WRITE(6,*) 'ERROR: modeloption lakeriverice and water temperature model does not work with other time step than day (to come in future release)'
        n_Result = 1
        RETURN
      ENDIF
      IF(.NOT.readswobs.AND.modeloption(p_snowmelt)==2)THEN
        WRITE(6,*) 'ERROR: modeloption snowmeltmodel=2 require SWobs.txt to work with other time step than day (for now)'
        n_Result = 1
        RETURN
      ENDIF
    ENDIF    

    !Write information about simulation to hyss.log
    WRITE(6,*)
    WRITE(6,*)'--------Information about the simulation---------'
    WRITE(6,*)'Current model set-up:  ',TRIM(mdir)
    WRITE(6,*)'Simulation sequence:   ',simsequence
    WRITE(6,*)
    WRITE(6,*)'Simulation begin date  ',strdate1
    WRITE(6,*)'Crit/Output start date ',strdate3
    WRITE(6,*)'Simulation end date    ',strdate2
    WRITE(6,'(A23,I2,A3)') 'Simulation steplength ',tstep,stepunit
    IF(calibration) WRITE(6,*)'Calibration run'
    WRITE(6,*)'Number of substances modelled',nsubst
    WRITE(6,*)'Model options:'
    WRITE(6,*)'               PET model',modeloption(p_petmodel)
    WRITE(6,*)'          Snowmelt model',modeloption(p_snowmelt)
    WRITE(6,*)'          Snowfall model',modeloption(p_snowfall)
    WRITE(6,*)'          Snowevap model',modeloption(p_snowevap)
    WRITE(6,*)'     Lakeriver ice model',modeloption(p_lakeriverice)
    WRITE(6,*)'  Deep groundwater model',modeloption(p_deepgroundwater)
    WRITE(6,*)'Surface water temp model',modeloption(p_swtemperature)
    WRITE(6,*)'     Growth season model',modeloption(p_growthstart)
    WRITE(6,*)'    Infiltration model',modeloption(p_infiltration)
    WRITE(6,*)'Groundwater runoff model',modeloption(p_frozenrunoff)
    WRITE(6,*)
    WRITE(6,*)'Output information:'
    IF(ioutbasin>0)THEN
      WRITE(6,*)ivarbasinoutput,'output variables will be saved for ',ioutbasin,'subbasins'
    ELSEIF(outvarallbasin)THEN
      WRITE(6,*)ivarbasinoutput,'output variables will be saved for all subbasins'
    ENDIF
    WRITE(6,*)ivarmapoutput,'output variables will be saved as files for mapping'
    WRITE(6,*)ivartimeoutput,'output variables will be saved as timeseries output'
    IF(writeload.AND.(conductN.OR.conductP))THEN
      WRITE(6,*) 'Yearly loads will be saved'
      conductload = .TRUE.
    ELSEIF(writeload)THEN
      WRITE(6,*) 'Yearly loads will not be saved. Need to simulate N or P.'
      writeload = .FALSE.
    ENDIF
    IF(sonum>0)THEN
      WRITE(6,*)'State variables will be saved for: ',strstates(1:sonum)
    ENDIF
    WRITE(6,*)'-------------------------------------------------'
    WRITE(6,*)

    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)
    IF(ALLOCATED(c)) DEALLOCATE(c)
    IF(ALLOCATED(critweight)) DEALLOCATE(critweight)
    IF(ALLOCATED(critc)) DEALLOCATE(critc)
    IF(ALLOCATED(critr)) DEALLOCATE(critr)
    IF(ALLOCATED(flowc)) DEALLOCATE(flowc)
    IF(ALLOCATED(flowr)) DEALLOCATE(flowr)
    IF(ALLOCATED(critpar)) DEALLOCATE(critpar)
    IF(ALLOCATED(strstates)) DEALLOCATE(strstates)
    IF(ALLOCATED(critcond)) DEALLOCATE(critcond)
    IF(ALLOCATED(critthres)) DEALLOCATE(critthres)

    !Checks for error in model set-up
    IF(.NOT.readdaily)THEN   
      IF(step>26*366)THEN
        WRITE(6,*)
        WRITE(6,*) 'WARNING: Simulation period >26 year'
        WRITE(6,*) 'WARNING: Use readdaily if problem with saving PT to memory'
        WRITE(6,*)
      ENDIF
    ENDIF

600 FORMAT(A18000)    !linelen!  kanske kan kortas ner?
    RETURN

800 CONTINUE
    n_Result = 1
    WRITE(6,*) 'ERROR: reading info.txt. io=',io
    CLOSE(fileunit_temp)
    RETURN

  END SUBROUTINE load_coded_info

   !>\brief Reads file with information about submodel for the current
  !run and >gets the coupling to the model set-up
  !!
  !!Format of file: first row is number of subbasins, and the
  !!following holds the subid of these subbasins
  !--------------------------------------------------------------------
  SUBROUTINE load_submodel_info(dir,submodelrun,msub,n_Result)

    USE WORLDVAR, ONLY : fileunit_temp,   &
                         ibasemodel
    USE MODVAR, ONLY : nsub_basemodel, &
                       basin

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory for subinfo-file
    LOGICAL, INTENT(IN)  :: submodelrun  !<Submodel is supposed to be run
    INTEGER, INTENT(OUT) :: msub         !<Number of subbasins of submodel
    INTEGER, INTENT(OUT) :: n_Result     !<Error status
    
    !Local variables
    LOGICAL fileexist
    INTEGER i,j
    CHARACTER (LEN=100) filename
    INTEGER, ALLOCATABLE :: selectedaro(:)

    n_Result = 0

    !File exist for simulation with submodel?
    filename=TRIM(dir)//'pmsf.txt'
    INQUIRE(FILE=filename, EXIST=fileexist)
    IF(submodelrun.AND.(.NOT.fileexist))THEN
       WRITE(6,*) 'ERROR: File is missing:',TRIM(filename)
       n_Result = 1
       RETURN
    ENDIF

    !Read subbasins to be simulated
    IF(submodelrun)THEN
       OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
       READ(fileunit_temp,*) msub
       IF(.NOT.ALLOCATED(selectedaro)) ALLOCATE(selectedaro(msub))
       READ(fileunit_temp,*) selectedaro
       CLOSE(fileunit_temp)
       WRITE(6,*) 'File read: ',TRIM(filename)
    ELSE
       msub = nsub_basemodel
    ENDIF

    !Set reordering of subbasins if necessary
    IF(.NOT.ALLOCATED(ibasemodel)) ALLOCATE(ibasemodel(msub))
    ibasemodel(:)=-9999
    IF(submodelrun)THEN
       DO i = 1,nsub_basemodel
          DO j = 1,msub
             IF(basin(i)%subid==selectedaro(j))THEN
                ibasemodel(j)=i
                EXIT
             ENDIF
          ENDDO
       ENDDO
       !Check if alla subbasins found
       j=0 
       DO i=1,msub 
          IF(ibasemodel(i)==-9999) THEN
             IF (j==0) THEN
                WRITE(6,*) ' '
                WRITE(6,*) '---------------------------------------------------'
                WRITE(6,*) 'ERROR Subbasins in pmsf.txt not found i geodata.txt'
                WRITE(6,*) '    SUBBASIN:       PMSF#'
                j=1
             ENDIF
             WRITE(6,*) selectedaro(i), i
          ENDIF
       ENDDO
       IF(j==1) THEN
          WRITE(6,*)
          WRITE(6,*) 'PROGRAM TERMINATED!'
          WRITE(6,*) '---------------------------------------------------'
          n_Result = 1
          RETURN
       ENDIF

       IF(ALLOCATED(selectedaro)) DEALLOCATE(selectedaro)   

    ELSE

       !Set index straight if no submodel
       DO i = 1,nsub_basemodel
          ibasemodel(i)=i
       ENDDO
    ENDIF

    !Write information to hyss.log
    IF(submodelrun)THEN
       WRITE(6,*)
       WRITE(6,*)'-----Information about the simulation (cont.)----'
       WRITE(6,*)'Submodel defined by file: ',TRIM(filename)
       WRITE(6,*)'Number of chosen subbasins to simulate',msub
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)
    ENDIF

  END SUBROUTINE load_submodel_info

  !>Preparation for subbasin output
  !--------------------------------------------------------------------
  SUBROUTINE prepare_subbasin_output()

    USE MODVAR, ONLY : basin,            &
                       nsub
    USE WORLDVAR, ONLY : temp_outbasins, &
                         outvarbasinindex,    &
                         outvarallbasin                  
    !Local variables
    INTEGER i,isub

    IF(.NOT.ALLOCATED(outvarbasinindex)) ALLOCATE(outvarbasinindex(nsub))
    outvarbasinindex = 0

    IF(outvarallbasin)THEN
       DO isub = 1,nsub
          outvarbasinindex(isub) = isub
       ENDDO
       RETURN
    ELSE
       IF(.NOT.ALLOCATED(temp_outbasins)) RETURN
       DO i = 1,SIZE(temp_outbasins)
          DO isub = 1,nsub
             IF(temp_outbasins(i)==basin(isub)%subid)THEN
                outvarbasinindex(isub) = i
                EXIT
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(temp_outbasins)
    ENDIF
  END SUBROUTINE prepare_subbasin_output

  !>Writes the files with data in format suitable for mapping.
  !--------------------------------------------------------------------
  SUBROUTINE save_mapfiles(dir,numsub,n,iens,runens,allens)

    USE WORLDVAR, ONLY : mapdata,     &
                         maptime,     &
                         writematlab, &
                         resultseq,   &
                         simsequence, &
                         outvarindex, &
                         outvarinfo,  &
                         o_nout,      &
                         o_ndecout,   &
                         o_wperiod,   &
                         i_t,i_d,i_w,i_m,i_s
    USE MODVAR, ONLY : outvarid,      &
                       basin,         &
                       i_sum

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN) :: numsub        !<Number of subbasins
    INTEGER, INTENT(IN) :: n             !<Number of values in vector = number of timesteps
    INTEGER, INTENT(IN) :: iens          !<Current simulation
    LOGICAL, INTENT(IN) :: runens        !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allens        !<Flag for writing all ensemble results

    !Local variables
    INTEGER i          !Loop variables
    INTEGER j          !Timestep
    INTEGER isub       !Subbasin
    INTEGER nn         !File prefix number
    CHARACTER(LEN=20) filename  !Current filename
    CHARACTER(LEN=20) var
    CHARACTER(LEN=30) unit
    REAL x(n)          !Values
    INTEGER lt,lout
    CHARACTER (LEN=16)   t
    CHARACTER (LEN=32000) outtxt           !32000=maxmap*16
    LOGICAL writeperiod

    filename=''
    DO i = 1,outvarinfo(2,o_nout)

       !Find variable and file
       IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. mapCOUT_005.txt
          IF(resultseq.AND.simsequence>0) nn=simsequence
          IF(runens) nn=iens
          filename(1:7) = outvarid(outvarindex(2,i,1))%mapfilename(1:7)
          filename(8:8) = '_'
          IF(allens)THEN
            WRITE(filename(9:15),602) nn
            filename(16:19) = '.txt'
          ELSE
            WRITE(filename(9:11),601) nn
            filename(12:15) = '.txt'
          ENDIF
       ELSE  !Ex. mapCOUT.txt
          filename=outvarid(outvarindex(2,i,1))%mapfilename
       ENDIF
601    FORMAT(I3.3)
602    FORMAT(I7.7)
       var=outvarid(outvarindex(2,i,1))%longname
       unit=''
       unit=outvarid(outvarindex(2,i,1))%longunit
       IF(outvarindex(2,i,2)==i_sum)THEN
          SELECT CASE(outvarinfo(2,o_wperiod))
          CASE(i_t)
             unit=TRIM(unit)//' per timestep'
          CASE(i_d)
             unit=TRIM(unit)//' per day'
          CASE(i_w)
             unit=TRIM(unit)//' per week'
          CASE(i_m)
             unit=TRIM(unit)//' per month'
          CASE DEFAULT
             unit=TRIM(unit)//' per year'
          END SELECT
       ENDIF
       OPEN(file=TRIM(dir)//filename,unit=10,status='unknown',form='formatted')

       !Write heading
       IF(writematlab)THEN
          WRITE(10,*) '%This is a file with '//TRIM(var)//' in '//TRIM(unit)//' for GIS mapping'
       ELSE
          WRITE(10,*) 'This is a file with '//TRIM(var)//' in '//TRIM(unit)//' for GIS mapping'
       ENDIF
       outtxt = ''
       IF(writematlab)THEN
          outtxt(1:7) = '%SUBID'//CHAR(44)
          lout = 7
       ELSE
          outtxt(1:6) = 'SUBID'//CHAR(44)
          lout = 6
       ENDIF
       IF(outvarinfo(2,o_wperiod)==i_s) THEN
          writeperiod=.TRUE.
       ELSE
          writeperiod=.FALSE.
       ENDIF
       DO j = 1,n
          t = ''
          t = TRIM(maptime(j))
          t = ADJUSTL(t)    !behövs?
          lt = LEN_TRIM(t)
          outtxt(lout+1:lout+lt) = t(1:lt)
          lout = lout+lt
          IF(j < n) THEN
             outtxt(lout+1:lout+1) = CHAR(44)   !Also the last one necessary
             lout = lout+1                      !for reading in "free format"
          ENDIF
       ENDDO
       WRITE(10,'(a)') outtxt(1:lout) 

       !Write data
       DO isub=1,numsub
          x=mapdata(isub,i,1:n)
          CALL write_dataline(10,n,x,outvarinfo(2,o_ndecout),&
               0,',',0,writematlab,id=basin(isub)%subid)
       ENDDO
       CLOSE(10)
    ENDDO

    RETURN
  END SUBROUTINE save_mapfiles

  !>Writes subbasin assessment to file subass.txt
  !--------------------------------------------------------------------
  SUBROUTINE write_subbasin_assessment(dir,numsub,sas2,iens,runens)

    USE WORLDVAR, ONLY : acccalvar,     &
                         calvarper,     &
                         fileunit_temp, &
                         writematlab,   &
                         resultseq,     &
                         simsequence,   &
                         maxsubass,       &
                         nacrit
    USE MODVAR, ONLY :   outvarid,      &
                         basin,         &
                         missing_value      

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(IN) :: numsub                 !<Number of subbasins
    REAL, INTENT(IN)    :: sas2(numsub,maxsubass,nacrit)  !<Subbasin assessment
    INTEGER, INTENT(IN) :: iens                   !<Current simulation
    LOGICAL, INTENT(IN) :: runens                 !<Flag for ensemble simulation
    
    !Local variables
    INTEGER i          !Criteria number
    INTEGER isub       !Subbasin
    INTEGER nn         !File suffix number
    CHARACTER(LEN=20) filename
    REAL x(maxsubass)  !Values
    INTEGER ndec       !Number of decimals to be written (rounding to ndec decimals will occurr)
    CHARACTER(LEN=1) tal

    tal = '0'
    DO i = 1, nacrit
      tal = CHAR(ICHAR(tal)+1)    !tal=i
      IF(tal==':') tal = 'A'      !för fler kriterier än 9 används bokstäver
      IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. subass1_005.txt
        IF(resultseq.AND.simsequence>0) nn=simsequence
        IF(runens) nn=iens
        filename(1:7) = 'subass'//tal
        filename(8:8) = '_'
        WRITE(filename(9:11),601) nn
        filename(12:15) = '.txt'
      ELSE  !Ex. subass1.txt
        filename='subass'//tal//'.txt'
      ENDIF
601   FORMAT(I3.3)
      OPEN(file=TRIM(dir)//filename,unit=fileunit_temp,status='unknown',form='formatted')
      !Write heading
      IF(writematlab)THEN
        WRITE(fileunit_temp,'(a58,I2,a13,a4,a2,a4,a8,a6)')    &
              '%Subbasin assessment. Criteria is calculated for period',calvarper,'. Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname,' Unit: ',outvarid(acccalvar(i)%comp)%shortunit
        WRITE(fileunit_temp,'(a)') '%SUBID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'
      ELSE
        WRITE(fileunit_temp,'(a57,I2,a13,a4,a2,a4,a8,a6)')    &
              'Subbasin assessment. Criteria is calculated for period',calvarper,'. Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname,' Unit: ',outvarid(acccalvar(i)%comp)%shortunit
        WRITE(fileunit_temp,'(a)') 'SUBID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'
      ENDIF

      !Write data
      ndec=3
      DO isub=1,numsub
        x=sas2(isub,1:maxsubass,i)
        IF(x(1)/=missing_value) CALL write_dataline(fileunit_temp,maxsubass,x,ndec,0,CHAR(9),0,writematlab,id=basin(isub)%subid)
      ENDDO
      CLOSE(fileunit_temp)
    ENDDO

  END SUBROUTINE write_subbasin_assessment

  !--------------------------------------------------------------------
  !>Writes the files with subbasin assessment.
  !! 
  !! @param[in] dir         File directory
  !! @param[in] iens        Current simulation
  !! @param[in] n           Dimension of performance measures
  !! @param[in] optcrit     ?
  !! @param[in] performance ?
  !! @param[in] runens      Flag for ensemble simulation
  !! @param[in] ccrit       conditional criteria
  !! @param[in] cthres      acceptance threshold for conditional criteria
  !--------------------------------------------------------------------
  SUBROUTINE write_simulation_assessment(dir,iens,n,optcrit,performance,runens,ccrit,cthres)

    USE WORLDVAR, ONLY : acccalvar,   &
         calvarper,   &
         resultseq,   &
         maxperf,     &
         simsequence, &
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
         fileunit_temp
    USE MODVAR, ONLY :   outvarid

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir           !File directory
    INTEGER, INTENT(IN)          :: iens          !Current simulation
    INTEGER, INTENT(IN)          :: n             !Dimension of performance measures
    REAL, INTENT(IN)             :: optcrit
    REAL, INTENT(IN)             :: performance(maxperf,n)
    LOGICAL, INTENT(IN)          :: runens        !Flag for ensemble simulation
    REAL, OPTIONAL, INTENT(IN)   :: ccrit         !conditional criteria
    REAL, OPTIONAL, INTENT(IN)   :: cthres        !conditional criteria threshold
 
    !Local variables
    INTEGER i
    INTEGER nn    !file suffix number
    CHARACTER(LEN=20) filename

    !Write simulation assessment to log-file
    WRITE(6,*)
    IF(PRESENT(cthres))THEN
       WRITE(6,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: ',ccrit,', threshold: ', cthres
    ELSE
       WRITE(6,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: -9999, threshold: -9999'
    ENDIF
    DO i = 1, n
       WRITE(6,*)
       WRITE(6,*) 'Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname
       WRITE(6,*) 'Period: ',calvarper
       WRITE(6,*) 'Regional NSE:',performance(i_rnse,i)
       WRITE(6,*) 'Regional  RA:',performance(i_rra,i)
       WRITE(6,*) 'Regional  RE:',performance(i_rre,i)
       WRITE(6,*) 'Regional MAE:',performance(i_rmae,i)
       WRITE(6,*) 'Average  NSE:',performance(i_mnse,i)
       WRITE(6,*) 'Average   RA:',performance(i_mra,i)
       WRITE(6,*) 'Average   RE:',performance(i_mre,i)
       WRITE(6,*) 'Average RSDE:',performance(i_mstdre,i)
       WRITE(6,*) 'Average   CC:',performance(i_mcc,i)
       WRITE(6,*) 'Average  ARE:',performance(i_mabsre,i)
       WRITE(6,*) 'Spatial  NSE:',performance(i_snse,i)
       WRITE(6,*) 'Spatial   RA:',performance(i_sra,i)
       WRITE(6,*) 'Spatial   RE:',performance(i_sre,i)
       WRITE(6,*) 'Kendalls Tau:',performance(i_tau,i)
       WRITE(6,*) 'Median   NSE:',performance(i_mdnse,i)
       WRITE(6,*) 'Median    RA:',performance(i_mdra,i)
       WRITE(6,*) 'Median   KGE:',performance(i_mdkg,i)
       WRITE(6,*) 'Median NRMSE:',performance(i_mnrmse,i)
    ENDDO
    WRITE(6,*)

    !Open file
    IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. simass_005.txt
       IF(resultseq.AND.simsequence>0) nn=simsequence
       IF(runens) nn=iens
       filename(1:6) = 'simass'
       filename(7:7) = '_'
       WRITE(filename(8:10),601) nn
       filename(11:14) = '.txt'
    ELSE
       filename='simass.txt'
    ENDIF
601 FORMAT(I3.3)
    OPEN(FILE=TRIM(dir)//filename,UNIT=fileunit_temp,STATUS='unknown',FORM='formatted')
    !Write heading
    WRITE(fileunit_temp,*) 'Simulation assessment all variable used in criterium'

    !Write data to simass-file
    WRITE(fileunit_temp,*) 'Simulation number: ',iens
    IF(PRESENT(cthres))THEN
      WRITE(fileunit_temp,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: ',ccrit,', threshold: ', cthres
    ELSE
      WRITE(fileunit_temp,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: -9999, threshold: -9999'
    ENDIF
    DO i = 1, n
       WRITE(fileunit_temp,*)
       WRITE(fileunit_temp,*) 'Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname
       WRITE(fileunit_temp,*) 'Period: ',calvarper
       WRITE(fileunit_temp,*) 'Regional NSE:',performance(i_rnse,i)
       WRITE(fileunit_temp,*) 'Regional  RA:',performance(i_rra,i)
       WRITE(fileunit_temp,*) 'Regional  RE:',performance(i_rre,i)
       WRITE(fileunit_temp,*) 'Regional MAE:',performance(i_rmae,i)
       WRITE(fileunit_temp,*) 'Average  NSE:',performance(i_mnse,i)
       WRITE(fileunit_temp,*) 'Average   RA:',performance(i_mra,i)
       WRITE(fileunit_temp,*) 'Average   RE:',performance(i_mre,i)
       WRITE(fileunit_temp,*) 'Average RSDE:',performance(i_mstdre,i)
       WRITE(fileunit_temp,*) 'Average   CC:',performance(i_mcc,i)
       WRITE(fileunit_temp,*) 'Average  ARE:',performance(i_mabsre,i)
       WRITE(fileunit_temp,*) 'Spatial  NSE:',performance(i_snse,i)
       WRITE(fileunit_temp,*) 'Spatial   RA:',performance(i_sra,i)
       WRITE(fileunit_temp,*) 'Spatial   RE:',performance(i_sre,i)
       WRITE(fileunit_temp,*) 'Kendalls Tau:',performance(i_tau,i)
       WRITE(fileunit_temp,*) 'Median   NSE:',performance(i_mdnse,i)
       WRITE(fileunit_temp,*) 'Median    RA:',performance(i_mdra,i)
       WRITE(fileunit_temp,*) 'Median   KGE:',performance(i_mdkg,i)
       WRITE(fileunit_temp,*) 'Median NRMSE:',performance(i_mnrmse,i)
    ENDDO
    WRITE(fileunit_temp,*)
    CLOSE(fileunit_temp)
    
  END SUBROUTINE write_simulation_assessment

  !>Prepare output files; opens files and write heading
  !--------------------------------------------------------------------
  SUBROUTINE prepare_outputfiles(dir,n,na,iens,runens,allsim)

    USE MODELMODULE, ONLY : open_modeldefined_outputfiles
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allsim       !<Flag for writing all simulation results

    CALL open_subbasinfiles(dir,n,iens,runens,allsim)
    CALL open_timefiles(dir,n,iens,runens,allsim)
    CALL open_modeldefined_outputfiles(dir,n,na,iens,runens)

  END SUBROUTINE prepare_outputfiles
  
  !>Prepare output files; close files
  !--------------------------------------------------------------------
  SUBROUTINE close_outputfiles(n,na,iens)

    USE MODELMODULE, ONLY : close_modeldefined_outputfiles
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation

    CALL close_subbasinfiles(n)
    CALL close_timefiles(iens)
    CALL close_modeldefined_outputfiles(n,na,iens)

  END SUBROUTINE close_outputfiles
  
  !>Opens files for subbasin printout
  !--------------------------------------------------------------------
  SUBROUTINE open_subbasinfiles(dir,n,iens,runens,allens)

    USE WORLDVAR, ONLY : outvarinfo,      &
                         outvarindex,     &
                         o_nout,          &
                         writematlab,     &
                         resultseq,       &
                         simsequence,     &
                         outvarbasinindex,  &
                         fileunit_add_sub
    USE MODVAR, ONLY : outvarid,        &
                       basin

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allens       !<Flag for writing all ensemble results

    !Local variables
    INTEGER i
    INTEGER nn                  !file suffix number
    INTEGER funit               !fileunit for subbasin output
    CHARACTER(LEN=15) str       !filename (without .txt)
    CHARACTER(LEN=5) var
    CHARACTER(LEN=6) varunit
    INTEGER j,lt,lout,ls,lout2
    CHARACTER (LEN=16)   t,s
    CHARACTER (LEN=2000) outtxt,outtxt2           !max_outvar*7

    !Open files
    DO i = 1,n
      IF(outvarbasinindex(i)>0)THEN
        funit=fileunit_add_sub+i
        str=''
        IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. 0000124_005.txt
          IF(resultseq.AND.simsequence>0) nn=simsequence
          IF(runens) nn=iens
          WRITE(str(1:7),600) basin(i)%subid
          str(8:8) ='_'
          IF(allens)THEN
            WRITE(str(9:15),602) nn   !Ex. 0000124_0001001.txt
          ELSE
            WRITE(str(9:11),601) nn
          ENDIF
        ELSE  !Ex. 0000124.txt
          WRITE(str,600) basin(i)%subid
        ENDIF
        OPEN(UNIT=funit,FILE=TRIM(dir)//TRIM(ADJUSTL(str))//'.txt',status='unknown',form='formatted')
      ENDIF
    ENDDO
600 FORMAT(I7.7)
601 FORMAT(I3.3)
602 FORMAT(I7.7)

    !Preparation of heading
    IF(writematlab)THEN
      outtxt(1:6) = '%DATE'//CHAR(9)
      outtxt2(1:7) = '%UNITS'//CHAR(9)
      lout = 6
      lout2 = 7
    ELSE
      outtxt(1:5) = 'DATE'//CHAR(9)
      outtxt2(1:6) = 'UNITS'//CHAR(9)
      lout = 5
      lout2 = 6
    ENDIF
    DO j = 1,outvarinfo(1,o_nout)
      var  = outvarid(outvarindex(1,j,1))%shortname
      varunit = outvarid(outvarindex(1,j,1))%shortunit
      WRITE(t,'(a16)') var
      WRITE(s,'(a16)') varunit
      t = ADJUSTL(t)
      s = ADJUSTL(s)
      lt = LEN_TRIM(t)
      ls = LEN_TRIM(s)
      outtxt(lout+1:lout+lt) = t(1:lt)
      outtxt2(lout2+1:lout2+ls) = s(1:ls)
      lout = lout+lt
      lout2 = lout2+ls
      IF(j < outvarinfo(1,o_nout)) THEN
        outtxt(lout+1:lout+1) = CHAR(9)
        outtxt2(lout2+1:lout2+1) = CHAR(9)    !Also the last one necessary
        lout = lout+1                      !for reading in "free format"
        lout2 = lout2+1
      ENDIF
    ENDDO

    !Write headings to files
    DO i = 1,n
      IF(outvarbasinindex(i)>0)THEN
        funit=fileunit_add_sub+i
        WRITE(funit,'(a)') outtxt(1:lout)     !Write heading
        WRITE(funit,'(a)') outtxt2(1:lout2)
      ENDIF
    ENDDO
    
  END SUBROUTINE open_subbasinfiles

  !>Calculate and write subbasin files
  !--------------------------------------------------------------------
  SUBROUTINE write_subbasinfiles(idt,ndt,d)

    USE LIBDATE, ONLY : DateType
    USE WORLDVAR, ONLY : outvarinfo,      &
                         o_nout,          &
                         o_ndecout,       &
                         o_wperiod,       &
                         writematlab,     &
                         outvarbasinindex,  &
                         outstartdate,    &
                         fileunit_add_sub
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_outvar
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    TYPE(DateType), INTENT(IN) :: d     !<timestep date-time

    !Local variables
    INTEGER isb   !subbasin index
    LOGICAL pwrite              !Flag for periodend, time to write to file
    REAL,ALLOCATABLE :: x(:)    !Help variable for continous print out

     IF(outvarinfo(1,o_nout)>0)THEN  !Write period mean for selected subbasins
       ALLOCATE(x(outvarinfo(1,o_nout)))
       DO isb = 1,nsub
         IF(outvarbasinindex(isb)>0)THEN
           CALL compute_outvar(d,isb,outvarinfo(1,o_nout),x,pwrite,idt,ndt)
           IF(pwrite) CALL write_dataline(fileunit_add_sub+isb,outvarinfo(1,o_nout),x,outvarinfo(1,o_ndecout),  &
                           outvarinfo(1,o_wperiod),CHAR(9),0,writematlab,d=d,odate=outstartdate)
         ENDIF
       ENDDO
     ENDIF
     IF(ALLOCATED(x)) DEALLOCATE(x)

  END SUBROUTINE write_subbasinfiles

  !>Close files for subbasin printout
  !--------------------------------------------------------------------
  SUBROUTINE close_subbasinfiles(n)
  
    USE WORLDVAR, ONLY : outvarbasinindex,  &
                         fileunit_add_sub

    INTEGER, INTENT(IN) :: n           !<Number of subbasins

    !Local variables
    INTEGER i

    DO i = 1,n 
       IF(outvarbasinindex(i)>0)THEN
          CLOSE(fileunit_add_sub+i)
       ENDIF
    ENDDO

  END SUBROUTINE close_subbasinfiles

  !>Opens files for timeserie printout
  !--------------------------------------------------------------------
  SUBROUTINE open_timefiles(dir,n,iens,runens,allens)
  
    USE WORLDVAR, ONLY : outvarinfo,        &
                         outvarindex,       &
                         writematlab,       &
                         resultseq,         &
                         simsequence,       &
                         i_t,i_h,i_d,i_w,i_m,i_y,       &
                         o_wperiod,         &
                         o_nout,            &
                         fileunit_add_time
    USE MODVAR, ONLY : outvarid,           &
                       i_sum,              &
                       basin

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<File directory
    INTEGER, INTENT(IN) :: n            !<Number of data columns (subbasins)
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble run
    LOGICAL, INTENT(IN) :: allens       !<Flag for writing all ensemble results
    
    !Local variables
    INTEGER i
    INTEGER nn                  !file suffix number
    INTEGER funit
    CHARACTER(LEN=20) filename
    CHARACTER(LEN=20) var
    CHARACTER(LEN=30) unit
    INTEGER j,lt,lout
    CHARACTER (LEN=16)   t
    CHARACTER (LEN=500000) outtxt           !17600=1100*16, !20110224 JD ändrat från 300000 till 500000

    filename=''
    DO i = 1, outvarinfo(3,o_nout)

      !Find variable and file
      IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. timeCOUT_005.txt
        IF(resultseq.AND.simsequence>0) nn=simsequence
        IF(runens) nn=iens
        filename(1:8) = outvarid(outvarindex(3,i,1))%timefilename(1:8)
        filename(9:9) = '_'
        IF(allens)THEN
          WRITE(filename(10:16),602) nn   !Ex. timeCOUT_0001001.txt
          filename(17:20) = '.txt'
        ELSE
          WRITE(filename(10:12),601) nn
          filename(13:16) = '.txt'
        ENDIF
      ELSE  !Ex. timeCOUT.txt
        filename=outvarid(outvarindex(3,i,1))%timefilename
      ENDIF
601    FORMAT(I3.3)
602    FORMAT(I7.7)

      var=outvarid(outvarindex(3,i,1))%longname
      unit=outvarid(outvarindex(3,i,1))%longunit
      IF(outvarindex(3,i,2)==i_sum)THEN
        SELECT CASE(outvarinfo(3,o_wperiod))
        CASE(i_h)
          unit=TRIM(unit)//' per hour'
        CASE(i_t)
          unit=TRIM(unit)//' per timestep'
        CASE(i_d)
          unit=TRIM(unit)//' per day'   
        CASE(i_w)
          unit=TRIM(unit)//' per week'
        CASE(i_m)
          unit=TRIM(unit)//' per month'
        CASE(i_y)
          unit=TRIM(unit)//' per year'
        END SELECT
      ENDIF
      funit = iens*fileunit_add_time+i

      !Open file
      OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=funit,STATUS='unknown',FORM='formatted')

      !Write headings
      IF(writematlab)THEN
        WRITE(funit,*) '%This is a file with timeseries of '//TRIM(var)//' in '//TRIM(unit)
        outtxt(1:6) = '%DATE'//CHAR(9)
        lout = 6
        DO j = 1,n
          WRITE(t,'(i16)') basin(j)%subid
          t = ADJUSTL(t)
          lt = LEN_TRIM(t)
          outtxt(lout+1:lout+lt) = t(1:lt)
          lout = lout+lt
          IF(j < n) THEN
            outtxt(lout+1:lout+1) = CHAR(9)    !Also the last one necessary
            lout = lout+1                      !for reading in "free format"
          ENDIF
        ENDDO
        WRITE(funit,'(a)') outtxt(1:lout)     !Write heading
      ELSE
        WRITE(funit,*) 'This is a file with timeseries of '//TRIM(var)//' in '//TRIM(unit)
        outtxt(1:5) = 'DATE'//CHAR(9)
        lout = 5
        DO j = 1,n
          WRITE(t,'(i16)') basin(j)%subid
          t = ADJUSTL(t)
          lt = LEN_TRIM(t)
          outtxt(lout+1:lout+lt) = t(1:lt)
          lout = lout+lt
          IF(j < n) THEN
            outtxt(lout+1:lout+1) = CHAR(9)    !Also the last one necessary
            lout = lout+1                      !for reading in "free format"
          ENDIF
        ENDDO
        WRITE(funit,'(a)') outtxt(1:lout)     !Write heading
      ENDIF
    ENDDO

  END SUBROUTINE open_timefiles

  !>Opens files for timeserie printout
  !--------------------------------------------------------------------
  SUBROUTINE write_timefiles(idt,ndt,iens,d)
  
    USE WORLDVAR, ONLY : outvarinfo,        &
                         writematlab,       &
                         o_wperiod,         &
                         o_ndecout,         &
                         o_nout,            &
                         outstartdate,      &
                         fileunit_add_time
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_timevar
    USE READWRITE_ROUTINES, ONLY : write_dataline
    USE LIBDATE, ONLY : DateType

    !Argument declarations
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: d     !<timestep date-time
    
    !Local variables
    INTEGER ivar    !index for variable
    LOGICAL pwrite              !Flag for periodend, time to write to file
    REAL,ALLOCATABLE :: y(:)    !Help variable for continous print out

    IF(outvarinfo(3,o_nout)>0) ALLOCATE(y(nsub))
    DO ivar = 1, outvarinfo(3,o_nout)      !Write time serie file
      CALL compute_timevar(d,ivar,nsub,y,pwrite,idt,ndt)
      IF(pwrite) CALL write_dataline(iens*fileunit_add_time+ivar,nsub,y,outvarinfo(3,o_ndecout),  &
                      outvarinfo(3,o_wperiod),CHAR(9),0,writematlab,d=d,odate=outstartdate)
    ENDDO 
    IF(ALLOCATED(y)) DEALLOCATE(y)

  END SUBROUTINE write_timefiles

  !>Close files for timeserie printout
  !--------------------------------------------------------------------
  SUBROUTINE close_timefiles(iens)
  
    USE WORLDVAR, ONLY : outvarinfo,  &
                         o_nout,      &
                         fileunit_add_time

    !Argument declarations
    INTEGER, INTENT(IN) :: iens        !<Current simulation
    
    !Local variables
    INTEGER i
    INTEGER funit

    DO i = 1, outvarinfo(3,o_nout)
      funit = iens * fileunit_add_time + i
      CLOSE(funit)
    ENDDO

  END SUBROUTINE close_timefiles

  !>Read parameter values from file
  !>
  !>/b Consequences The modvar module variables soilpar, landpar, genpar,
  !> basinpar, regpar, wqparreg, lregpar, lakedatapar and monthpar will 
  !> be allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_parameters(dir,ns,name) 

    USE WORLDVAR, ONLY : fileunit_temp,     &
                         comment_str
    USE MODVAR, ONLY : soilpar,    &   !OUT, parameters depending on soil type
                       landpar,    &   !OUT, parameters depending on land use
                       genpar,     &   !OUT, parameters not dependent
                       basinpar,   &   !OUT, parameters depending on subbbasin
                       regpar,     &   !OUT, parameters depending on parregion
                       wqregpar,   &   !OUT, parameters depending on wqparregion
                       lregpar,    &   !OUT, parameters depending on lake region
                       lakedatapar, &  !OUT, parameters from LakeData 
                       monthpar,   &   !OUT, parameters depending on month
                       modparid,   &   !definition of model parameters
                       nregions,  &
                       nwqregions,  &
                       nlakeregions,   &
                       max_par,    &
                       m_gpar,  &
                       m_bpar,  &
                       m_spar,  &
                       m_lpar,  &
                       m_rpar,  &
                       m_wqrpar,  &
                       m_lrpar, &
                       m_ldpar, &
                       m_mpar, &
                       nluse,   &
                       nsoil,   &
                       allocate_modelparameters, & 
                       nilakeregions, &
                       nolakeregions, &
                       m_ilrpar, & 
                       m_olrpar, &
                       ilregpar,   &   !OUT, parameters depending on ilakeregion
                       olregpar   !OUT, parameters depending on olakeregion
                       

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER,          INTENT(IN) :: ns             !<Number of subbasins (nsub_basemodel)
    CHARACTER(LEN=*), INTENT(IN) :: name           !<File name
    
    !Local variables
    CHARACTER (LEN=220) filename
    CHARACTER (LEN=10) varstr       !string with variable name
    CHARACTER(LEN=18000) line
    INTEGER j
    INTEGER dim
    INTEGER dimcheck                !current parameter supposed size
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER varindex                !parameter index in array
    LOGICAL lakedatafile            !Status of LakeData file
    REAL,ALLOCATABLE :: values(:)               !values of parameter read from file

    !>\b Algorithm \n
    !Initialisations
    lakedatafile=.FALSE.
    IF(ALLOCATED(lakedatapar)) lakedatafile=.TRUE.

    !>Open par.txt file and read heading
    filename=TRIM(dir)//TRIM(name)
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')

    nskip = 1
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO

    !>Allocate and initiate parameter variables
    dim = MAX(1,ns,nsoil,nluse,nlakeregions,nregions,nwqregions,12,nilakeregions,nolakeregions)
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dim))
    CALL allocate_modelparameters(ns)

    !>Read parameter values from file
    nline = nskip + 1
    DO 
      READ(fileunit_temp,'(a)',END=100) line
      IF(line(1:2)==comment_str)CYCLE
      CALL read_parameterline(line,dim,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file

      !Find corresponding varindex
      varindex = 0
      DO j = 1,max_par
        IF(varstr==modparid(j)%shortname)THEN
          IF(modparid(j)%deptype/=m_ldpar)THEN 
            varindex = j
            EXIT
          ENDIF
        ENDIF
      ENDDO
      IF(varindex==0)THEN
        WRITE(6,*) 'Warning - unknown variable on line. Parameter ',TRIM(varstr),' ignored.'
!        WRITE(6,*) ' Error - check file',filename
!        STOP 1
!      ENDIF
      ELSE
        !Find and check variable dimension and set parameter values
        SELECT CASE(modparid(varindex)%deptype)
        CASE(m_gpar)
          dimcheck = 1                   !genpar
        CASE(m_bpar)
          dimcheck = ns                  !basinpar
        CASE(m_spar)
          dimcheck = nsoil
        CASE(m_lpar)
          dimcheck = nluse
        CASE(m_lrpar)
          dimcheck = nlakeregions
        CASE(m_rpar)
          dimcheck = nregions
        CASE(m_wqrpar)
          dimcheck = nwqregions
        CASE(m_ilrpar)
          dimcheck = nilakeregions
        CASE(m_olrpar)
          dimcheck = nolakeregions
        CASE(m_mpar)
          dimcheck = 12
        END SELECT
        IF(dimcheck /= nvalues)THEN
          IF(modparid(varindex)%deptype==m_rpar .AND. dimcheck<nvalues)THEN
            nvalues = dimcheck     !Fullösning bugg2650(10716) antal regioner minskas vid selectaro, men ej antalet i par.txt
          ELSEIF(modparid(varindex)%deptype==m_lrpar .AND. dimcheck<nvalues)THEN
            nvalues = dimcheck     !Fullösning bugg2650(10716) antal lake regioner minskas vid selectaro, men ej antalet i par.txt
          ELSEIF(modparid(varindex)%deptype==m_wqrpar .AND. dimcheck<nvalues)THEN
            nvalues = dimcheck     !Fullösning bugg2650(10716) antal water quality regioner minskas vid selectaro, men ej antalet i par.txt
          ELSEIF(modparid(varindex)%deptype==m_ilrpar .AND. dimcheck<nvalues)THEN
            nvalues = dimcheck     !Fullösning - Antal ilakeregion minskas vid selectaro, men ej antalet i par.txt
          ELSEIF(modparid(varindex)%deptype==m_olrpar .AND. dimcheck<nvalues)THEN
            nvalues = dimcheck     !Fullösning - Antal olakeregion minskas vid selectaro, men ej antalet i par.txt
          ELSE
            WRITE(6,*) 'Error - wrong number of values in parameter file'
            WRITE(6,*) 'Error - check parameter ',varstr
            STOP 1
          ENDIF
        ENDIF
        IF(dimcheck.GE.1)THEN ! To avoid problems when dimcheck=0 (for instance, no olregpar defined and allocated)
          SELECT CASE(modparid(varindex)%deptype)
          CASE(m_gpar)
            genpar(modparid(varindex)%parno) = values(1)
          CASE(m_bpar)
            basinpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_spar)
            soilpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_lpar)
            landpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_rpar)
            regpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_wqrpar)
            wqregpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_ilrpar)
            ilregpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_olrpar)
            olregpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_lrpar)
            lregpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          CASE(m_mpar)
            monthpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
          END SELECT
        ENDIF
      ENDIF

      nline = nline + 1
    ENDDO
100 CONTINUE
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)
    IF(ALLOCATED(values)) DEALLOCATE(values)

    !>Set parameter values for missing lakes with general values
    CALL finish_lakedata_table(lakedatafile,ns)     !Fill missing values in lakedatapar with general values  


  END SUBROUTINE load_parameters

  !>Reads parameters to be calibrated from file and saves them in
  !>worldvar-arrays.
  !!
  !>/b Consequences Module worldvar variables optparmin, optparmax, optparprecision,
  !> optparid, parindex, optim, dimpar, numoptimpar will be allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_optpar(dir) 

    USE WORLDVAR, ONLY : optparmin, optparmax, & 
                         optparprecision, optparid, &
                         dimpar,maxoptpar,    &
                         parindex,       &
                         optim,               &
                         numoptimpar,         &
                         fileunit_temp
    USE MODVAR, ONLY : modparid,max_par,  &
                       nsoil,nluse,nsub,nregions,nwqregions,nlakeregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<File directory

    !Local variables
    INTEGER j, caseFlag
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER numpar                  !number of different parameters to be calibrated (parameter names in optpar.txt)
    INTEGER varindex                !index for parameter in varstr
    REAL,ALLOCATABLE :: values(:)   !values of parameter read from file (different soil or landuse)
    REAL,ALLOCATABLE :: par(:,:)    !Tempary storage of parameter value 
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10) varstr       !string with parameter name
    CHARACTER(LEN=18000) line       !line in file
    CHARACTER(LEN=2) taskchar       !code for optimation task
    CHARACTER(LEN=1) flagChar       !character string for one-character flags (Y/N)

    !Local parameters
    INTEGER, PARAMETER :: ninfolines = 20

    !Initiation of variables for optimization parameters
    dimpar = MAX(1,nsub,nsoil,nluse,nregions,nwqregions,nlakeregions)
    IF(.NOT.ALLOCATED(optparmin)) ALLOCATE(optparmin(maxoptpar,dimpar)) !worldvar
    IF(.NOT.ALLOCATED(optparmax)) ALLOCATE(optparmax(maxoptpar,dimpar)) !worldvar
    IF(.NOT.ALLOCATED(optparprecision)) ALLOCATE(optparprecision(maxoptpar,dimpar))   !worldvar
    optparmin = 0.
    optparmax = 0.
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dimpar))           !local
    IF(.NOT.ALLOCATED(par)) ALLOCATE(par(maxoptpar,dimpar))       !local

    filename=TRIM(dir)//'optpar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    nskip = 1
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO
    nline = nskip + 1

    DO j = 1, ninfolines     !Read lines with optimisation information
      READ(fileunit_temp,'(a)',END=100) line

      IF(line(1:4)=='task')THEN
        READ(line(6:16000),*,END=102) taskchar
        SELECT CASE(taskchar)
        CASE('MC')
          optim%task_MC        = .TRUE.
        CASE('BP')
          optim%task_boundps   = .TRUE.
        CASE('SC')                        ! Scanning mode for 2-parameter problems, execute an organised sampling of crit(param)  [Fred, 06.10.11]
          optim%task_Scanning  = .TRUE.
        CASE('SM')                        ! Stage-wise centering Monte Carlo routine  [Fred, 26.08.10]
          optim%task_stageMC   = .TRUE.
        CASE('BN')                        ! Brent optimisation routine, with new line search  [Fred, 27.11.10]
          optim%task_BrentNew  = .TRUE.
        CASE('SD')                        ! Steepest descent method  [Fred, 08.10.11]
          optim%task_stpstDesc = .TRUE.
        CASE('Q1')                        ! Quasi-Newton DFP gradient-based optimisation routine  [Fred, 14.07.10]
          optim%task_DFP       = .TRUE.
        CASE('Q2')                        ! Quasi-Newton BFGS gradient-based optimisation routine  [Fred, 09.12.10]
          optim%task_BFGS      = .TRUE.
        CASE('WA') 
          optim%task_writeall  = .TRUE.
        CASE('WS') 
          optim%task_writesim  = .TRUE.  !Write simulation result for all ensembles (MC-methods)
        CASE('DE')
          optim%task_DEMC      = .TRUE.   !Differential-Evolution Markov Chain method [David, 2013.02.11]
        CASE('DR')
          optim%task_DREAM     = .TRUE.   !Differential-Evolution Adaptive Metropolis-Hastings (Markov Chain)
        END SELECT
      ENDIF

102   IF(line(1:6)=='num_mc')THEN
        READ(line(8:16000),*,END=103) optim%nruns_MC
      ENDIF

103   IF(line(1:7)=='num_ens')THEN
        READ(line(9:16000),*,END=104) optim%nruns_best
      ENDIF

104   IF(line(1:8)=='num_bpmc')THEN
        READ(line(10:16000),*,END=105) optim%nruns_MCI
      ENDIF

105   IF(line(1:9)=='num_bpmax')THEN
        READ(line(11:16000),*,END=106) optim%nruns_MCImax
      ENDIF

106   IF(line(1:10)=='num_stages')THEN                        ! Read 'num_stages' (amount of successive zooming MC stages) for zoom MC routine  [Fred, 26.08.10]
        READ(line(11:16000),*,END=107) optim%nruns_stages
      ENDIF

107   IF(line(1:8)=='num_zoom')THEN                           ! Read numerical parameter 'num_zoom' (zooming factor) for zoom MC routine  [Fred, 26.08.10]
        READ(line(11:16000),*,END=108) optim%nruns_zoom       ! Checked later (in stage MC procedure) that it is not larger than 1
      ENDIF

108   IF(line(1:10)=='num_dbgMod')THEN                        ! Flag for debug scenario (0 = HYPE, see optim.f90 for alternatives)
        READ(line(11:16000),*,END=109) optim%cal_debugCase
      ENDIF

109   IF(line(1:7)=='cal_log')THEN                            ! Y/N flag to write or not the file "calibration.log", that contains all details on the calibration routine
        READ(line(11:16000),*,END=110) flagChar
        IF(flagChar == 'N')THEN
          optim%cal_log = .FALSE.
        ELSEIF(flagChar .NE. 'N' .AND. flagChar .NE. 'Y')THEN
          WRITE(6,*) 'ERROR: flag to enable/disable writing calibration.log must be either Y or N'
          STOP 1
        ENDIF
      ENDIF

      !Numerical parameter for interruption of non-MC methods [Fred; checked for consistency with optim-type on 29.09.11]
110   IF(line(1:10)=='num_maxItr')THEN                        ! Max amount of iterations allowed
        READ(line(11:16000),*,END=111) optim%cal_maxIterat
      ENDIF

111   IF(line(1:10)=='num_maxTim')THEN                        ! Max amout of time (hours) allowed to calibration routine
        READ(line(11:16000),*,END=112) optim%cal_maxTime
      ENDIF

112   IF(line(1:10)=='num_criItr')THEN                        ! Amount of last optimisation iterations taken into account for criteria improvement monitoring
        READ(line(11:16000),*,END=113) optim%cal_improvCritIter
      ENDIF

113   IF(line(1:10)=='num_criTol')THEN                        ! Tolerance to consider criteria as optimised (delta/mean of criteria over the last 'num_critIt' iterations)
        READ(line(11:16000),*,END=114) optim%cal_improvCritTol
      ENDIF

114   IF(line(1:10)=='num_parItr')THEN                        ! Amount of last optimisation iterations taken into account for parameter improvement monitoring
        READ(line(11:16000),*,END=115) optim%cal_improvParamIter
      ENDIF

      !Numerical parameters for Quasi-Newton method [Fred]
115   IF(line(1:10)=='QN_nrmTol')THEN                         ! Tolerance for gradient norm to be considered zero
        READ(line(11:16000),*,END=116) optim%QN_flatTol
      ENDIF

116   IF(line(1:10)=='QN_pctDerv')THEN                        ! Offset (percentage of parameter value) for numerical derivative
        READ(line(11:16000),*,END=117) optim%QN_factorDeriv
      ENDIF

117   IF(line(1:10)=='QN_stencil')THEN                        ! Stencil type
        READ(line(11:16000),*,END=118) optim%QN_stencil       ! QN optimisation algorithm checkes that it is either 2, 4, 6 or 8
      ENDIF

118   IF(line(1:10)=='QN_lambMax')THEN                        ! Factor containing lambda, to avoid taking points for gradient outside allowed parameter space
        READ(line(11:16000),*,END=119) optim%QN_lambdaMaxFac
      ENDIF

119   IF(line(1:10)=='QN_lambAcc')THEN                        ! Factor increasing the step length proposed by QN algorithms (case lambda = 1 to be replaced by lambdaAccel), in order to allow for faster iteration progression
        READ(line(11:16000),*,END=120) optim%QN_lambdaAccel
      ENDIF

120   IF(line(1:10)=='BR_diagStp')THEN                        ! Flag to disable diagonal step at the end of each Brent iteration
        READ(line(11:16000),*,END=122) flagChar
        IF(flagChar == 'N')THEN
          optim%Brent_diagonalStep = .FALSE.
        ELSEIF(flagChar .NE. 'N' .AND. flagChar .NE. 'Y')THEN
          WRITE(6,*) 'ERROR: flag to enable/disable diagonal Brent step must be either Y or N'
          STOP 1
        ENDIF
      ENDIF

      !Numerical parameters for line search [Fred]
122   IF(line(1:10)=='lnS_maxItr')THEN                        ! Maximum amount of allowed line search iterations
        READ(line(11:16000),*,END=123) optim%lineSearch_maxIter
      ENDIF

123   IF(line(1:7)=='lnS_tol')THEN                            ! Tolerance for numerical line search
        READ(line(11:16000),*,END=124) optim%lineSearch_tol
      ENDIF

124   IF(line(1:9)=='scan_numx')THEN                          ! Amount of scanning points in the dimension of the 1st parameter
        READ(line(11:16000),*,END=125) optim%scan_xpoints
      ENDIF

125   IF(line(1:9)=='scan_numy')THEN                          ! Amount of scanning points in the dimension of the 2nd parameter
        READ(line(11:16000),*,END=126) optim%scan_ypoints
      ENDIF

      !Parameters for DE-MC [David, 2013-02-11]
126   IF(line(1:9)=='DEMC_ngen')THEN                          ! Number of generations in DE-MC simulation
        READ(line(11:16000),*,END=127) optim%DEMC_ngen
      ENDIF

127   IF(line(1:9)=='DEMC_npop')THEN                          ! Number of populations in DE-MC simulation
        READ(line(11:16000),*,END=128) optim%DEMC_npop
      ENDIF
128   IF(line(1:15)=='DEMC_gammascale')THEN                   ! Scaling of the default gamma factor (2.38/sqrt(2*npar)) in DE-MC
        READ(line(17:16000),*,END=129) optim%DEMC_gammascale
      ENDIF
129   IF(line(1:14)=='DEMC_crossover')THEN                     !Crossover probability in DE-MC simulation
        READ(line(16:16000),*,END=130) optim%DEMC_crossover
      ENDIF
130   IF(line(1:10)=='DEMC_sigma')THEN                         !Standard deviation of sample perturbations in DE-MC simulation
        READ(line(12:16000),*,END=131) optim%DEMC_sigma
      ENDIF
131   IF(line(1:12)=='DEMC_accprob')THEN                         !Standard deviation of sample perturbations in DE-MC simulation
        READ(line(14:16000),*,END=132) optim%DEMC_accprob
      ENDIF
!Parameters for DREAM [David, 2013-10-08]
132     IF(line(1:10)=='DREAM_ngen')THEN
          READ(line(12:16000),*,END=133) optim%DREAM_ngen
        ENDIF
133    IF(line(1:10)=='DREAM_npop')THEN
          READ(line(12:16000),*,END=134) optim%DREAM_npop
       ENDIF
134    IF(line(1:11)=='DREAM_npair')THEN
          READ(line(13:16000),*,END=135) optim%DREAM_npair
       ENDIF
135    IF(line(1:14)=='DREAM_longjump')THEN
          READ(line(16:16000),*,END=136) optim%DREAM_longjump
       ENDIF
136    IF(line(1:16)=='DREAM_ncrossover')THEN
          READ(line(18:16000),*,END=137) optim%DREAM_ncrossover
       ENDIF
137    IF(line(1:20)=='DREAM_crossover')THEN
          READ(line(22:16000),*,END=138) optim%DREAM_crossover
       ENDIF
138    IF(line(1:16)=='DREAM_gammascale')THEN
          READ(line(18:16000),*,END=139) optim%DREAM_gammascale
       ENDIF
139    IF(line(1:13)=='DREAM_eminmax')THEN
          READ(line(15:16000),*,END=140) optim%DREAM_eminmax
       ENDIF
140    IF(line(1:14)=='DREAM_epssigma')THEN
          READ(line(16:16000),*,END=141) optim%DREAM_epssigma
       ENDIF
141    IF(line(1:13)=='DREAM_accprob')THEN
          READ(line(15:16000),*,END=142) optim%DREAM_accprob
       ENDIF
142    IF(line(1:16)=='DREAM_crtunestep')THEN
          READ(line(18:16000),*,END=143) optim%DREAM_crtunestep
       ENDIF
143    IF(line(1:17)=='DREAM_grthreshold')THEN
          READ(line(19:16000),*,END=144) optim%DREAM_grthreshold
       ENDIF
144    IF(line(1:14)=='DREAM_keepbest')THEN
          READ(line(16:16000),*,END=145) optim%DREAM_keepbest
       ENDIF
145    IF(line(1:13)=='DREAM_accmult')THEN
          READ(line(15:16000),*,END=146) optim%DREAM_accmult
       ENDIF
146    IF(line(1:13)=='DREAM_critopt')THEN
          READ(line(15:16000),*,END=147) optim%DREAM_critopt
       ENDIF
147    IF(line(1:12)=='DREAM_accexp')THEN
          READ(line(14:16000),*,END=148) optim%DREAM_accexp
       ENDIF
148    IF(line(1:15)=='DREAM_accmethod')THEN
          READ(line(17:16000),*,END=149) optim%DREAM_accmethod
       ENDIF
149    IF(line(1:13)=='DREAM_criteps')THEN
          READ(line(15:16000),*,END=150) optim%DREAM_criteps
       ENDIF
150    CONTINUE
    ENDDO

    nline = nline + ninfolines 

    IF(optim%task_MC .OR. optim%task_boundps .OR. optim%task_stageMC) optim%task_runens = .TRUE.
    IF(optim%task_DEMC .OR. optim%task_DREAM) optim%task_runens = .TRUE.
    IF(optim%task_DEMC)                       optim%nruns_best = optim%DEMC_npop+1 !Number of populations corresponds to the number of best, add 1 to store median of the populations as "the best" NO1
    IF(optim%task_DREAM)                      optim%nruns_best = optim%DREAM_npop+1 !Number of populations corresponds to the number of best, add 1 to store median of the populations as "the best" NO1
    IF(.NOT. optim%task_runens) optim%nruns_best = 1  
    IF(optim%task_writesim .AND. .NOT.(optim%task_MC .OR. optim%task_boundps .OR. optim%task_DEMC .OR. optim%task_DREAM))THEN
      WRITE(6,*) 'Warning: It is only allowed to write all simulation results for some of the '
      WRITE(6,*) 'Warning: MonteCarlo methods. Write all simulation results is turned off'
      optim%task_writesim = .FALSE.
    ENDIF
    
    numpar = 0
    caseFlag = 1
    DO 
      READ(fileunit_temp,'(a)',END=240) line
      CALL read_parameterline(line,dimpar,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file
      varindex = 0
      DO j = 1,max_par                !find index of parameter
        IF(varstr==modparid(j)%shortname)THEN
          varindex = j
          EXIT
        ENDIF
      ENDDO
      IF(varindex==0)THEN             !check if parameter index found
        WRITE(6,*) ' Error - unknown variable on line: ', nline
        WRITE(6,*) ' Error - check file',filename
        STOP 1
      ENDIF
      IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
        numpar = numpar + 1
        optparmin(numpar,1:nvalues) = values(1:nvalues)
        optparid(numpar)=varindex
        caseFlag = 2
      ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
        optparmax(numpar,1:nvalues) = values(1:nvalues)
        caseFlag = 3
      ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
        optparprecision(numpar,1:nvalues) = values(1:nvalues)
        caseFlag = 1
      ENDIF

      nline = nline + 1
    ENDDO

240 CONTINUE
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Sort min and max values
    WHERE(optparmin>optparmax)
      par=optparmax
      optparmax=optparmin
      optparmin=par
    ENDWHERE

    !Set and check worldvar numoptimpar = amount of model parameters to be optimized
    CALL count_optim_par(maxoptpar,dimpar)
    IF(numoptimpar == 0) THEN
      WRITE(6,*) 'ERROR: only zero-length model parameters intervals found in optpar.txt!'
      STOP 1 
    ENDIF

    !Allocate parindex
    IF(.NOT.ALLOCATED(parindex)) ALLOCATE(parindex(numoptimpar,2))

    IF(.NOT.(optim%task_MC .OR. optim%task_boundps .OR. optim%task_Scanning .OR. optim%task_stageMC .OR. optim%task_BrentNew .OR. optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS .OR. optim%task_DEMC .OR. optim%task_DREAM)) THEN     ! [Fred 06.10.11] Added task_Scanning to checklist
      WRITE(6,*) 'ERROR: no task set in optpar'
      STOP 1
    ENDIF
    IF(ALLOCATED(values)) DEALLOCATE(values)
    IF(ALLOCATED(par)) DEALLOCATE(par)
    RETURN

100 WRITE(6,*) 'End of file during read, optpar.txt' 
    CLOSE(fileunit_temp)
    STOP 1

  END SUBROUTINE load_optpar

  !>Reads parameters to be calibrated to find dimension of river states
  !--------------------------------------------------------------------
  SUBROUTINE calculate_special_optpar_parameters(dir,velindex,dampindex,rivvel,damp) 

    USE WORLDVAR, ONLY : fileunit_temp
    USE MODVAR, ONLY : modparid,nsoil,nluse,nsub,nregions,nwqregions,nlakeregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<file directory
    INTEGER, INTENT(IN)  :: velindex   !<index of rivvel in modparid
    INTEGER, INTENT(IN)  :: dampindex  !<index of damp in modparid
    REAL, INTENT(OUT)    :: rivvel   !<lower river velocity boundary
    REAL, INTENT(OUT)    :: damp     !<lower damp boundary

    !Local variables
    INTEGER j, caseFlag
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER localdimpar             !maximal dimention of parameters
    INTEGER nvalues                 !number of values read from line
    LOGICAL fexist                  !presence of file
    REAL minvalv,maxvalv,minvald,maxvald
    REAL,ALLOCATABLE :: values(:)   !values of parameter read from file (different soil or landuse)
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10) varstr       !string with parameter name
    CHARACTER(LEN=18000) line       !line in file

    !Local parameters
    INTEGER, PARAMETER :: ninfolines = 20

    !Initiate default value
    rivvel = 9999.; damp = 9999.   !unreasonable large values
    minvald = 0.; maxvald = 0.
    minvalv = 0.; maxvalv = 0.
    
    !Check file presence
    filename=TRIM(dir)//'optpar.txt'
    INQUIRE(FILE=TRIM(filename),EXIST=fexist)
    IF(.NOT.fexist) RETURN
    
    !Initiation of variables for optimization parameters
    localdimpar = MAX(1,nsub,nsoil,nluse,nregions,nwqregions,nlakeregions)
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(localdimpar))           !local

    filename=TRIM(dir)//'optpar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    nskip = 1
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO
    nline = nskip + 1
    DO j = 1, ninfolines     !Skip lines with optimisation information
      READ(fileunit_temp,'(a)',END=100) line
    ENDDO
    nline = nline + ninfolines 

    caseFlag = 1
    DO 
      READ(fileunit_temp,'(a)',END=140) line
      CALL read_parameterline(line,localdimpar,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file

      IF(varstr==modparid(velindex)%shortname)THEN
        IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
          minvalv = values(1)
          caseFlag = 2
        ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
          maxvalv = values(1)
          caseFlag = 3
        ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
          caseFlag = 1
        ENDIF
      ENDIF
      IF(varstr==modparid(dampindex)%shortname)THEN
        IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
          minvald = values(1)
          caseFlag = 2
        ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
          maxvald = values(1)
          caseFlag = 3
        ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
          caseFlag = 1
        ENDIF
      ENDIF

      nline = nline + 1
    ENDDO

140 CONTINUE
    CLOSE(fileunit_temp)

    !Set output
    IF(minvalv>0. .OR. maxvalv>0.) rivvel = MIN(minvalv,maxvalv)
    IF(minvald>0. .OR. maxvald>0.) damp = MIN(minvald,maxvald)
    IF(ALLOCATED(values)) DEALLOCATE(values)
    RETURN
    
100 WRITE(6,*) 'End of file during read optpar.txt for river dimension' 
    CLOSE(fileunit_temp)
    STOP 1

  END SUBROUTINE calculate_special_optpar_parameters

  !>Save optimal values of optimized parameters to file respar.txt
  !-----------------------------------------------------------------
  SUBROUTINE save_respar(dir,numpar,n) 

    USE WORLDVAR, ONLY : optparid,  &
                         parindex,  &
                         fileunit_temp
    USE MODVAR, ONLY   : modparid,  & 
                         nluse,     &
                         nsoil,     &
                         nregions,  &
                         nwqregions,  &
                         nlakeregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<File directory
    INTEGER, INTENT(IN)          :: numpar !<number of parameters that has been calibrated
    INTEGER, INTENT(IN)          :: n      !<number of subbasins
    
    !Local variables 
    INTEGER i,j
    INTEGER dim                     !max number of values for parameter
    INTEGER num                     !number of values for parameter
    INTEGER ndec
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10)  numparstr   !number of parameter values
    CHARACTER (LEN=50)  form_str    !format for writing parameter
    INTEGER varindex                !array-index of parameter
    LOGICAL written                 !status of parameter
    REAL,ALLOCATABLE :: par(:)      !parameter values

    !>\b Algoritm \n
    !>Allocate and initiate local variables
    dim = MAX(1,n,nsoil,nluse,nregions,nwqregions,nlakeregions)
    IF(.NOT.ALLOCATED(par)) ALLOCATE(par(dim))

    !>Open file for writing parameter and write heading
    filename=TRIM(dir)//'respar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'unknown')
    WRITE(fileunit_temp,*) 'Optimal value of parameters found during automatic calibration'    !heading

    !>For every calibrated parameter:
    ndec = 7    !7 decimals for parameter values
    DO i = 1,numpar
      varindex = optparid(parindex(i,1))
      written = .FALSE.
      DO j = 1, i-1
        IF(varindex==optparid(parindex(j,1))) written = .TRUE.       !Check if parameter already written to file
      ENDDO
      IF(.NOT. written)THEN
        !>\li Get current parameter value
        CALL get_parametervalues(varindex,n,dim,par,num)
        !>\li Write parameter to file
        form_str=''
        WRITE(numparstr,'(I10)') num   !Integer to character
        form_str = '(A10,'//TRIM(ADJUSTL(numparstr))//'(1x,F16.'//CHAR(ndec+48)//'))'
        WRITE(fileunit_temp,form_str) modparid(varindex)%shortname,par(1:num)       !write parameter values to file
      ENDIF
    ENDDO

    CLOSE(fileunit_temp)
    IF(ALLOCATED(par)) DEALLOCATE(par)

  END SUBROUTINE save_respar

  !>Collect parameter values from model variables
  !--------------------------------------------------------------
  SUBROUTINE get_parametervalues(varindex,n,dim,values,num)

    USE MODVAR, ONLY : soilpar,  &
                       landpar,  &
                       genpar,   &
                       regpar,   &
                       wqregpar, &
                       lregpar,  &
                       basinpar, &
                       monthpar, &
                       modparid, &
                       m_gpar,   &
                       m_bpar,   &
                       m_spar,   &
                       m_lpar,   &
                       m_rpar,   &
                       m_wqrpar, &
                       m_lrpar,  &
                       m_mpar,   &
                       nluse,    &
                       nsoil,    &
                       nregions,  &
                       nwqregions,  &
                       nlakeregions, &
                       m_ilrpar, &
                       m_olrpar, &
                       nilakeregions, &
                       nolakeregions, &
                       ilregpar, &
                       olregpar

    !Argument declarations
    INTEGER, INTENT(IN)  :: varindex    !<model parameter index
    INTEGER, INTENT(IN)  :: n           !<number of subbasins
    INTEGER, INTENT(IN)  :: dim         !<max number of parameter values (soil types/land uses/subbasins)
    REAL, INTENT(OUT)    :: values(dim) !<parameter values
    INTEGER, INTENT(OUT) :: num         !<number of parameter values (soil types/land uses/subbasins) that are used

    !>\b Algoritm \n
    !>Depending on variable type: Find variable dimension and parameter values
    IF(modparid(varindex)%deptype==m_gpar)THEN
      num = 1                    !genpar
      values(1:num) = genpar(modparid(varindex)%parno)
    ELSEIF(modparid(varindex)%deptype==m_bpar)THEN
      num = n                    !basinpar
      values(1:num) = basinpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_spar)THEN
      num = nsoil
      values(1:num) = soilpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_lpar)THEN
      num = nluse
      values(1:num) = landpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_rpar)THEN
      num = nregions
      values(1:num) = regpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_wqrpar)THEN
      num = nwqregions
      values(1:num) = wqregpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_ilrpar)THEN
      num = nilakeregions
      values(1:num) = ilregpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_olrpar)THEN
      num = nolakeregions
      values(1:num) = olregpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_lrpar)THEN
      num = nlakeregions
      values(1:num) = lregpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_mpar)THEN
      num = 12
      values(1:num) = monthpar(modparid(varindex)%parno,1:num)
    ENDIF

  END SUBROUTINE get_parametervalues

  !>Initiate accumulation variables for printout
  !>
  !>/b Consequences Module worldvar variables accdata, accdatahelp, accdata_classload, 
  !> accdata_basinload, taccdata, taccdatahelp, mapdata, mapdatahelp, maptime, tmap, 
  !> and accdataok are reset.
  !------------------------------------------------------------------
  SUBROUTINE initiate_output_routines(iout) 

    USE WORLDVAR, ONLY : accdata,      &  !OUT
                         accdatahelp,  &  !OUT
                         accdata_classload,  &  !OUT
                         accdata_basinload,  &  !OUT
                         taccdata,     &  !OUT
                         taccdatahelp, &  !OUT
                         mapdata,      &      !OUT!Accumulated data for mapping
                         maptime,      &      !OUT!Time for mapdata
                         mapdatahelp,  &  !OUT
                         tmap,         &  !OUT
                         writeload,    &
                         max_typeofoutput, &
                         accdataok  !OUT
    USE MODVAR, ONLY : numsubstances                     

    !Argument declarations
    INTEGER, INTENT(IN) :: iout(max_typeofoutput)      !<number of variables for printout

    accdata = 0.
    accdatahelp = 0.
    accdataok = .TRUE.
    IF(iout(2)>0)THEN
      mapdata = 0.
      mapdatahelp = 0.
      maptime = ''
      tmap = 1
    ENDIF
    IF(iout(3)>0)THEN
      taccdata = 0.
      taccdatahelp = 0.
    ENDIF
    IF(writeload)THEN
      IF(numsubstances>0)THEN
        accdata_classload = 0.
        accdata_basinload = 0.
      ENDIF
    ENDIF

  END SUBROUTINE initiate_output_routines

  !>Initiate output variables to missing value
  !>
  !>/b Consequences Module modvar variables firstoutstep and outvar is set.
  !-------------------------------------------------------------
  SUBROUTINE initiate_outvar(idt)

    USE MODVAR,   ONLY : outvar,        & !OUT
         missing_value, &
         firstoutstep
    USE WORLDVAR, ONLY : dtskip      

    INTEGER, INTENT(IN) :: idt  !<current timestep 

    firstoutstep = .FALSE. 
    IF(idt==dtskip+1) firstoutstep = .TRUE.      

    outvar = missing_value

  END SUBROUTINE initiate_outvar

  !>Set concentrations of waters with zero volume to missing_value
  !>
  !>/b Consequences Module modvar variables outvar maybe set.
  !------------------------------------------------------------------
  SUBROUTINE revise_outvar()

    USE MODVAR,   ONLY : outvar,        & !OUT
         outvarid,      &
         max_outvar,    &
         i_wmean,       &
         missing_value, &
         nsub

    !Local variables
    INTEGER i,isb

    DO i = 1, max_outvar
       IF(outvarid(i)%vartype==i_wmean)THEN
          DO isb = 1, nsub
             IF(outvar(isb,outvarid(i)%water)==0) outvar(isb,i) = missing_value
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE revise_outvar

  !>Write to file result from the nsim best MonteCarlo simulations
   !--------------------------------------------------------------------
  SUBROUTINE save_ensemble_simulations(dir,numpar,nperf,numcrit,&
       nsim,optcrit,performance,parameters)

    USE WORLDVAR, ONLY : writematlab,      &
         filename_best, &
         fileunit_temp

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    INTEGER, INTENT(IN) :: numpar             !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf              !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit            !<Number of optimised variables
    INTEGER, INTENT(IN) :: nsim               !<Number of simulations to be saved
    REAL, INTENT(IN) :: optcrit(nsim)         !<Optimation criterion for simulations
    REAL, INTENT(IN) :: performance(nsim,nperf,numcrit) !<Performance measures for simulations
    REAL, INTENT(IN) :: parameters(nsim,numpar) !<Parameter values for simulations
   
    !Local variables 
    INTEGER i

    !Open files
    OPEN(file=TRIM(dir)//filename_best,unit=fileunit_temp,status='unknown',form='formatted')
    !Write heading
    CALL write_ensemble_simulations_heading(fileunit_temp,numpar,nperf,numcrit,.FALSE.)
    !Write data
    DO i = 1,nsim
       CALL write_simulation_results(fileunit_temp,i,numpar,numpar,nperf,numcrit,&
            optcrit(i),performance(i,:,:),parameters(i,:),writematlab)
    ENDDO
    !End routine
    CLOSE(fileunit_temp)

  END SUBROUTINE save_ensemble_simulations

  !>Prepare a file to writes result from all MonteCarlo simulations
  !!
  !! @todo Why this procedure, why not let
  !! write_ensemble_simulations_heading open the file?
  !--------------------------------------------------------------------
  SUBROUTINE prepare_save_all_simulations(dir,filename,funit,numpar,nperf,numcrit)

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    CHARACTER(LEN=*), INTENT(IN) :: filename  !<Filename
    INTEGER, INTENT(IN) :: funit           !<Unit to be connected to file
    INTEGER, INTENT(IN) :: numpar          !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf           !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit         !<Number of optimised variables

    !Open file
    OPEN(file=TRIM(dir)//filename,unit=funit,status='unknown',form='formatted')
    !Write heading
    CALL write_ensemble_simulations_heading(funit,numpar,nperf,numcrit,.TRUE.)

  END SUBROUTINE prepare_save_all_simulations

  !>Write heading to file for MonteCarlo simulations
  !--------------------------------------------------------------------
  SUBROUTINE write_ensemble_simulations_heading(funit,numpar,nperf,numcrit,popflag)

    USE WORLDVAR, ONLY : performance_name,       &
                         writematlab,      &
                         optparid,      &
                         parindex, &
                         optim  
    USE MODVAR, ONLY :   modparid

    !Argument declarations
    INTEGER, INTENT(IN) :: funit    !<File unit
    INTEGER, INTENT(IN) :: numpar   !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf    !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit  !<Number of optimised variables
    LOGICAL, INTENT(IN) :: popflag  !<Flag for population columns included
    
    !Local variables
    INTEGER i,j,index
    INTEGER lt,lout
    CHARACTER (LEN=16)  t
    CHARACTER (LEN=800) outtxt     !heading
    CHARACTER (LEN=1)   sep        !separator

    !Write heading
    sep = ','
    IF(writematlab)THEN
       outtxt(1:4) = '%NO'//sep       !i
       outtxt(5:9) = 'CRIT'//sep     !optcrit
       lout = 9
    ELSE
       outtxt(1:3) = 'NO'//sep       !i
       outtxt(4:8) = 'CRIT'//sep     !optcrit
       lout = 8
    ENDIF
    DO i = 1,numcrit                !performance measures
       DO j = 1,nperf
          t  = performance_name(j)
          t  = ADJUSTL(t)
          lt = LEN_TRIM(t)
          outtxt(lout+1:lout+lt+1) = t(1:lt)//sep
          lout = lout+lt+1
       ENDDO
    ENDDO
    DO i = 1,numpar                 !parameters
       index = optparid(parindex(i,1))
       t = modparid(index)%shortname
       t  = ADJUSTL(t)
       lt = LEN_TRIM(t)
       outtxt(lout+1:lout+lt+1) = t(1:lt)//sep
       lout = lout+lt+1
    ENDDO
    !IF DREAM or DE-MC simulation, also print out the generation and population indeces
    IF(popflag .AND. (optim%task_DEMC .OR. optim%task_DREAM))THEN
       outtxt(lout+1:lout+4+1) = 'jpop'//sep
       lout = lout+4+1
       outtxt(lout+1:lout+4+1) = 'igen'//sep
       lout = lout+4+1
       outtxt(lout+1:lout+4+1) = 'iacc'//sep
       lout = lout+4+1
    ENDIF
    WRITE(funit,'(a)') outtxt(1:lout-1)

  END SUBROUTINE write_ensemble_simulations_heading

  !>Write the result from the last simulations to the file with all
  !>simulations' results.
  !--------------------------------------------------------------------
  SUBROUTINE write_simulation_results(funit,i,mpar,numpar,nperf,&
       numcrit,optcrit,performance,parameters,mlab,jpop,igen,iacc)

    INTEGER, INTENT(IN) :: funit         !<File unit
    INTEGER, INTENT(IN) :: i             !<Simulation number
    INTEGER, INTENT(IN) :: mpar          !<Dimension of parameters
    INTEGER, INTENT(IN) :: numpar        !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf         !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit       !<Number of optimised variables
    REAL, INTENT(IN) :: optcrit          !<Value of optimation criterion for the simulation
    REAL, INTENT(IN) :: performance(nperf,numcrit)  !<performance criteria for the simulation
    REAL, INTENT(IN) :: parameters(mpar)  !<Parameter values for the simulation
    LOGICAL, INTENT(IN) :: mlab           !<MATLAB format for print out
    INTEGER, OPTIONAL, INTENT(IN) :: jpop !<Population Index in DE-MC simulation
    INTEGER, OPTIONAL, INTENT(IN) :: igen !<Generation Index in DE-MC simulation
    INTEGER, OPTIONAL, INTENT(IN) :: iacc !<Acceptance Index in DE-MC simulation

    !Local variables
    REAL, ALLOCATABLE :: x(:)
    INTEGER n          !Number of columns with data to be written
    INTEGER ndec       !Number of decimals to be written (rounding to ndec decimals)

    !Write data
    ndec=8
    n = 1 + nperf*numcrit + numpar
    IF(PRESENT(iacc)) n = n+3
    IF(.NOT.ALLOCATED(x)) ALLOCATE(x(n))
    IF(PRESENT(iacc))THEN
       x = (/optcrit,performance(:,:),parameters(1:numpar),REAL(jpop),REAL(igen),REAL(iacc)/)
    ELSE
       x = (/optcrit,performance(:,:),parameters(1:numpar)/)
    ENDIF
    CALL write_dataline(funit,n,x,ndec,0,',',0,mlab,id=i)
    DEALLOCATE(x)

  END SUBROUTINE write_simulation_results

  !>Save the loads for the last year to files. Subroutine is called once per year.
  !>
  !>/b Consequences Module worldvar variables accdata_classload and 
  !>accdata_basinload are zeroed.
  !--------------------------------------------------------------------
  SUBROUTINE save_loadfiles(dir,year)

    USE WORLDVAR, ONLY : accdata_classload, &     !OUT
                         accdata_basinload, &     !OUT
                         fileunit_temp
    USE MODVAR, ONLY : basin,               &
                       numsubstances,       &
                       nclass,              &
                       i_in,i_on,i_sp,i_pp, &
                       nsub,                &
                       max_classoutvar,     &
                       max_basinoutvar,     &
                       loadheadings

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory
    INTEGER, INTENT(IN) :: year           !<Current year

    !Local variables
    INTEGER i       ! substance number       
    INTEGER j       ! class number
    INTEGER k       ! heading number in 'Headings'
    INTEGER l       ! heading number 'Prepheadings'
    INTEGER s       ! subbasin number
    INTEGER subnr   ! subbasin id
    INTEGER ffunit   ! file unit
    INTEGER nvalues   ! number of values in xload
    REAL xload(max_classoutvar*nclass+max_basinoutvar) !input to write_tab_sep
    CHARACTER (LEN=7) str !file name
    CHARACTER (LEN=10),DIMENSION(max_classoutvar*nclass+max_basinoutvar+1) :: headings ! Array of headings for writing to file

    !Preparation of heading
    headings(:)='0'
    WRITE(headings(1),'(a7,a)') loadheadings(1),CHAR(9)                     !subbasin number
    DO l=2,max_classoutvar+1
      DO j=1,nclass
        k=(l-2)*nclass+j+1 ! position in Headings
        WRITE(headings(k),'(a6,a1,I2.2,a)') loadheadings(l),'_',j,CHAR(9)     !class dependent
      ENDDO
    ENDDO
    DO l=max_classoutvar+2,max_classoutvar+max_basinoutvar+1 
      k=max_classoutvar*nclass+(l-max_classoutvar)
      WRITE(headings(k),'(a6,a)') loadheadings(l),CHAR(9)                   !subbasin dependent
    ENDDO

    !Save data to files
    ffunit=fileunit_temp
    DO i=1,numsubstances
      IF (i==i_in.OR.i==i_on.OR.i==i_sp.OR.i==i_pp)THEN
        IF (i==i_in) WRITE(str,'(I4,3a)') year,'_IN'
        IF (i==i_on) WRITE(str,'(I4,3a)') year,'_ON'
        IF (i==i_sp) WRITE(str,'(I4,3a)') year,'_SP'
        IF (i==i_pp) WRITE(str,'(I4,3a)') year,'_PP'
        OPEN(UNIT=ffunit,FILE=TRIM(dir)//TRIM(ADJUSTL(str))//'.txt',status='unknown',form='formatted')
        WRITE(ffunit,'(999a10)') headings     !Write headings to file
        DO s=1,nsub
          !Write data to a row array
          DO k=1,max_classoutvar
            xload(((k-1)*nclass+1):(k*nclass))=accdata_classload(1:nclass,k,i,s)
          ENDDO
          xload((max_classoutvar*nclass+1):(max_classoutvar*nclass)+max_basinoutvar)=accdata_basinload(i,:,s)
          nvalues=max_classoutvar*nclass+max_basinoutvar
          subnr = basin(s)%subid
          CALL write_dataline(ffunit,nvalues,xload,3,0,CHAR(9),0,.FALSE.,id=subnr) 
        ENDDO ! end subarea loop (s)
        CLOSE(ffunit)
      ENDIF
    ENDDO !end substance lop (i)

    !Reset accumulations to zero at end of year 
    accdata_classload(:,:,:,:) = 0.0  
    accdata_basinload(:,:,:)   = 0.0 

  END SUBROUTINE save_loadfiles

  !>Check input data from info, GeoData, BranchData, GeoClass and par files
  !--------------------------------------------------------------------
  SUBROUTINE checkindata_part1()

    USE WORLDVAR, ONLY : bdate,ndt,   &
                         sdate,       &
                         checkdata,   &
                         modeldir,    &
                         fileunit_temp, &
                         comment_str, &
                         vegtypereset
    USE MODVAR, ONLY : basin,         &
                       classbasin,    &
                       classdata,     &
                       cropdata,      &
                       classmodel,    &
                       slc_olake,     &
!                       vegtype,       &
                       pathsubid,     &
                       branchsubid,   &
                       soildepth,     &
                       streamdepth,   &
                       nlakeregions,  &
                       nregions,      &
                       nwqregions,    &
                       nclass,        &
                       nluse,nsoil,   &
                       conductN,conductP,conductC,  &
                       nsub_basemodel, &
                       modparid,  &
                       modeloption,  &
                       max_par,   &
                       m_gpar,m_bpar,  &
                       m_spar,m_lpar,  &
                       m_rpar,m_lrpar, &
                       m_ldpar,m_mpar, &
                       m_wqrpar,  &
                       p_growthstart, &
                       m_ilrpar, &
                       m_olrpar, & 
                       nilakeregions, & 
                       nolakeregions

    !Local variables
    REAL    basinarea
    REAL    classtest
    REAL    basintest
    INTEGER subidproc,classproc
    INTEGER i,k,bdim,ksource,kbranch
    INTEGER j
    INTEGER dim
    INTEGER dimcheck                !current parameter supposed size
    INTEGER ncomment                !number of lines with comments
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER varindex                !parameter index in array
    LOGICAL found
    CHARACTER (LEN=220) filename
    CHARACTER (LEN=10) varstr       !string with variable name
    CHARACTER(LEN=18000) line
    REAL,ALLOCATABLE :: values(:)               !values of parameter read from file

    WRITE(6,*)
    WRITE(6,*) 'Check indata part one'
    WRITE(6,*) '---------------------'
    
    !Check simulation period
    WRITE(6,*) 'Checking info.txt...'
    IF(ndt<=0) THEN
      WRITE(6,*) 'Number of timesteps less than one. Check bdate and edate in info.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(bdate%day<1.OR.bdate%day>31) WRITE(6,*) 'Check bdate day in info.txt'
    IF(bdate%month<1.OR.bdate%month>12) WRITE(6,*) 'Check bdate month in info.txt'
    IF(sdate%day<1.OR.sdate%day>31) WRITE(6,*) 'Check edate day in info.txt'
    IF(sdate%month<1.OR.sdate%month>12) WRITE(6,*) 'Check edate month in info.txt'
    
    !Check GeoData
    WRITE(6,*) 'Checking GeoData.txt...'
    basinarea = SUM(basin(:)%area)
    IF(basinarea<=0.)THEN 
      WRITE(6,*) 'Error: No area of catchment. Check GeoData.txt.'
      checkdata(1,2) = .TRUE.
    ENDIF
    basinarea = MINVAL(basin(:)%area)
    IF(basinarea<=0.)THEN
      DO i = 1,nsub_basemodel
        IF(basin(i)%area<=0.)THEN
          WRITE(6,*) 'Subbasin',basin(i)%subid,'has area',basin(i)%area,'.'
          IF(SUM(classbasin(i,:)%part)/=0.)THEN
            WRITE(6,*) 'and its sum of slc-fraction is not zero (Error).'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF  
      ENDDO
    ENDIF
    subidproc = MINVAL(basin(:)%subid)
    IF(subidproc==0)THEN
      WRITE(6,*) 'Error: Some subid are zero. Check GeoData.txt.'
      checkdata(1,2)=.TRUE.
    ENDIF
    subidproc = MAXVAL(basin(:)%subid)
    IF(subidproc>=100000000) WRITE(6,*) 'Subid is too large. Maximum 9999999 may be used. Check GeoData.txt.'
    basintest = MINVAL(basin(:)%slope)
    IF(basintest<0)THEN
      WRITE(6,*) 'Error: Some slope is less than zero. Check GeoData.txt.'
      checkdata(1,2)=.TRUE.
    ENDIF
    IF(conductN.OR.conductP.OR.conductC) THEN
      IF(nlakeregions<=0)THEN
        WRITE(6,*) 'Error: Max value of lake regions is zero. Lake regions are necessary for'
        WRITE(6,*) 'NPC simulation. Check GeoData.txt.'
        checkdata(1,2)=.TRUE.
      ELSE
        subidproc = MINVAL(basin(:)%lakeregion)
        IF(subidproc<=0)THEN
          WRITE(6,*) 'Error: Some lake regions are not set. Lake regions are necessary for NPC simulation. Check GeoData.txt.'
          checkdata(1,2)=.TRUE.
        ENDIF
      ENDIF
    ENDIF  
    IF(nregions<=0)THEN
      WRITE(6,*) 'Max value of parameter regions is zero. Check GeoData.txt.'
    ENDIF
    IF(conductN.OR.conductP) THEN
      IF(nwqregions<=0)THEN
        WRITE(6,*) 'Max value of water quality parameter regions is zero. Check GeoData.txt.'
      ENDIF
      subidproc = MINVAL(basin(:)%wqparregion)
      IF(subidproc<=0) WRITE(6,*) 'Some water quality parameter regions are not set. Check GeoData.txt.'
    ENDIF
    subidproc = MINVAL(basin(:)%parregion)
    IF(subidproc<=0) WRITE(6,*) 'Some parameter regions are not set. Check GeoData.txt.'
    subidproc = MINVAL(basin(:)%region)
    IF(subidproc<=0)THEN
      WRITE(6,*) 'Some regions are not set. Regions are necessary for NP simulation'
      WRITE(6,*) '(and irrigation). Check GeoData.txt.'
    ENDIF
    basinarea=SUM(basin(:)%rivlen)
    IF(basinarea<=0) WRITE(6,*) 'No main rivers? Check GeoData.txt (rivlen).'
    basinarea=MINVAL(basin(:)%rivlen)
    IF(basinarea<0) WRITE(6,*) 'Negative length of some rivers. That is not allowed. Check GeoData.txt'
    !Check slc-fractions
    DO i = 1,nsub_basemodel
      IF(SUM(classbasin(i,:)%part)/=1.) WRITE(6,*) 'Subbasin ',basin(i)%subid,' has not its sum of class-fractions equal to one.'
      IF(SUM(classbasin(i,:)%part)<0.99 .OR. SUM(classbasin(i,:)%part)>1.01) WRITE(6,*) 'Subbasin ',basin(i)%subid,' has sum of class-fractions deviating more then 1%.'
    ENDDO
    !Check olake depth
    IF(slc_olake>0)THEN
      DO i = 1,nsub_basemodel
        IF(classbasin(i,slc_olake)%part>0. .AND. basin(i)%lakedepth(2)<=0.)THEN
          WRITE(6,*) 'Subbasin ',basin(i)%subid,' has olake with depth zero. This is not allowed, model may crash.'
          WRITE(6,*) 'Check GeoData.txt, LakeData.txt and/or DamData.txt'
          checkdata(1,2) = .TRUE.
        ENDIF
      ENDDO
    ENDIF
    !Check linkage and subbasin order
    WRITE(6,*) 'Checking subbasin linkage...'
    DO i = 1,nsub_basemodel
      DO k = i-1,1,-1
        IF(basin(k)%subid==pathsubid(i)%main)THEN
          WRITE(6,*) 'Linking error: main branch (',pathsubid(i)%main,') of subbasin',basin(i)%subid,'comes before subbasin in GeoData.txt.'
          checkdata(1,2) = .TRUE.
        ENDIF
      ENDDO
      found = .FALSE.
      IF(pathsubid(i)%main<=0) found = .TRUE.
      IF(.NOT.found)THEN
        DO k = i+1,nsub_basemodel
          IF(pathsubid(i)%main==basin(k)%subid)THEN
            found=.TRUE.
            EXIT
          ENDIF  
        ENDDO
      ENDIF
      IF(.NOT.found) WRITE(6,*) 'Main branch (',pathsubid(i)%main,') of subbasin',basin(i)%subid,' not found. Check GeoData.txt'
    ENDDO  
    bdim=SIZE(branchsubid(:)%source)
    DO i = 1,bdim
      ksource = 0
      kbranch = 0
      DO k = 1,nsub_basemodel
        IF(branchsubid(i)%source==basin(k)%subid)  ksource = k
        IF(branchsubid(i)%branch==basin(k)%subid)  kbranch = k
      ENDDO
      IF(ksource==0) WRITE(6,*) 'Branch from subid',branchsubid(i)%source,' not in model. Check BranchData.txt.'
      IF(kbranch==0)THEN
        IF(branchsubid(i)%branch>0)  WRITE(6,*) 'Branch (',branchsubid(i)%branch,') of subbasin',branchsubid(i)%source,' not found. Check BranchData.txt.'
      ELSE
        IF(kbranch<ksource)THEN
          WRITE(6,*) 'Linking error: Branch (',branchsubid(i)%branch,') of subbasin'
          WRITE(6,*) branchsubid(i)%source,' lies upstream subbasin. Check BranchData.txt.'
          checkdata(1,2) = .TRUE.
        ENDIF
      ENDIF  
    ENDDO  
    IF(modeloption(p_growthstart)==1)THEN
      basinarea=SUM(basin(:)%latitude)
      IF(basinarea<=0) WRITE(6,*) 'Latitude is needed for using growth start model 1. No latitude given. Check GeoData.txt (latitude)?'
      basinarea=MINVAL(basin(:)%latitude)
      IF(basinarea<0) WRITE(6,*) 'Negative latitudes given. Check GeoData.txt'
    ENDIF
    
    !Check GeoClass
    WRITE(6,*) 'Checking GeoClass.txt...'
    IF(nclass<=0)THEN
      WRITE(6,*) 'Error: Zero classes. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    classproc=MINVAL(classdata(:)%luse)
    IF(classproc<=0)THEN
      WRITE(6,*) 'Error: No land use code given. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    classproc=MINVAL(classdata(:)%soil)
    IF(classproc<=0)THEN
      WRITE(6,*) 'Error: No soil type code given. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    DO i=1,nclass
      IF((classmodel(i)>=0.AND.classmodel(i)<=3).OR.classmodel(i)==11.OR.classmodel(i)==12)THEN
      ELSE
        WRITE(6,*) 'Unknown information in special class column for slc-class',i,'. Check GeoClass.txt.'
      ENDIF
    ENDDO
    classtest = MINVAL(soildepth(1,:))
    IF(classtest<=0.)THEN
      WRITE(6,*) 'Zero soil depths found for some class. Check GeoClass.txt'
    ENDIF
    classtest=MAXVAL(soildepth(:,:))
    IF(classtest<=0.)THEN
      WRITE(6,*) 'Error: No soil depths found. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    DO i = 1,nclass
      IF(vegtypereset)THEN
        WRITE(6,*) 'Error: vegetation type not given in GeoClass.txt for slc-class vegtype=1 (open) is used (when no indatacheck is off)'
        checkdata(1,2) = .TRUE.
      ENDIF
      IF(MAXVAL(soildepth(:,i))<streamdepth(i))THEN
        WRITE(6,*) 'Warning: Stream depth below deepest soillayer for class:',i
      ENDIF
    ENDDO

    !Check parameters
    WRITE(6,*) 'Checking par.txt...'
    filename=TRIM(modeldir)//'par.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    !Initiate some variables
    dim = MAX(1,nsub_basemodel,nsoil,nluse,nlakeregions,nregions,nwqregions,12)
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dim))
    !Read and check heading
    READ(fileunit_temp,'(a)',END=100) line
    READ(line,*,ERR=10) varstr
    varindex = 0
    DO j = 1,max_par
      IF(varstr==modparid(j)%shortname)THEN
        IF(modparid(j)%deptype/=m_ldpar)THEN 
          varindex = j
          EXIT
        ENDIF
      ENDIF
    ENDDO
    IF(varindex>0) WRITE(6,*) 'First line (heading) has parameter name. This parameter will not be used. Ok?'
10  CONTINUE
    !Read and check rest of parameters and their values    
    nline = 1
    ncomment = 0
    DO 
      nline = nline + 1
      READ(fileunit_temp,'(a)',END=100) line
      IF(line(1:2)==comment_str)THEN
        ncomment = ncomment + 1
        CYCLE
      ENDIF
      CALL read_parameterline(line,dim,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')THEN
        WRITE(6,*) 'No parameter name found on line',nline,'. End of par.txt file?'
        EXIT    !end of file
      ENDIF
      varindex = 0
      DO j = 1,max_par      !Find corresponding varindex
        IF(varstr==modparid(j)%shortname)THEN
          IF(modparid(j)%deptype/=m_ldpar)THEN 
            varindex = j
            EXIT
          ENDIF
        ENDIF
      ENDDO
      IF(varindex==0)THEN
        WRITE(6,*) 'Unknown variable (',TRIM(varstr),') on line',nline,'.'
      ELSE
        !Find and check variable dimension and set parameter values
        SELECT CASE(modparid(varindex)%deptype)
        CASE(m_gpar)
          dimcheck = 1  
        CASE(m_bpar)
          dimcheck = nsub_basemodel 
        CASE(m_spar)
          dimcheck = nsoil
        CASE(m_lpar)
          dimcheck = nluse
        CASE(m_lrpar)
          dimcheck = nlakeregions
        CASE(m_rpar)
          dimcheck = nregions
        CASE(m_wqrpar)
          dimcheck = nwqregions
        CASE(m_ilrpar)
          dimcheck = nilakeregions
        CASE(m_olrpar)
          dimcheck = nolakeregions
        CASE(m_mpar)
          dimcheck = 12
        END SELECT
        IF(dimcheck /= nvalues)THEN
          IF(modparid(varindex)%deptype==m_rpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is parameter region dependent, and extra values will be skipped.'
          ELSEIF(modparid(varindex)%deptype==m_lrpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is lake region dependent, and extra values will be skipped.'
          ELSEIF(modparid(varindex)%deptype==m_wqrpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is water quality region dependent, and extra values will be skipped.'
          ELSEIF(modparid(varindex)%deptype==m_ilrpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is ilake region dependent, and extra values will be skipped.'
          ELSEIF(modparid(varindex)%deptype==m_olrpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is olake region dependent, and extra values will be skipped.'
          ELSE
            WRITE(6,*) 'Error: wrong number of values in parameter file for parameter ',TRIM(varstr),'.'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
100 CONTINUE
    WRITE(6,*) 'Found',ncomment,' lines with comments in par.txt.'
    CLOSE(fileunit_temp)
    IF(ALLOCATED(values)) DEALLOCATE(values)

    !Check CropData.txt
    WRITE(6,*) 'Checking CropData.txt...'
    IF(modeloption(p_growthstart)==1)THEN
      filename=TRIM(modeldir)//'CropData.txt'
      INQUIRE(FILE = filename, EXIST= found)
      IF(.NOT.found)THEN
        WRITE(6,*) 'CropData.txt necessary for using growth start model 1. File not found.'
      ELSE
        basinarea=SUM(cropdata(:)%gddsow)
        IF(basinarea<=0) WRITE(6,*) 'Four columns needed for using growth start model 1, but gddsow not found. Check CropData.txt.'
        basinarea=MINVAL(cropdata(:)%gddsow)
        IF(basinarea<=0) WRITE(6,*) 'Negative or zero gddsow given. Check CropData.txt'
      ENDIF
    ENDIF

    WRITE(6,*) '---------------------'

  END SUBROUTINE checkindata_part1

  !>Check forcing data and other observations
  !------------------------------------------------------------
  SUBROUTINE checkindata_part2() 

    USE WORLDVAR, ONLY : simsequence,     &
                         tobsid,pobsid,   &
                         sfobsid,swobsid, &
                         uobsid,         &
                         rhobsid,         &
                         tminobsid,tmaxobsid, &
                         modeldir,        &
                         fileunit_temp,   &
                         filename_Qobs,   &
                         maxcharpath,     &
                         readsfobs,       &
                         readswobs, &
                         readwind, &
                         readhumid, &
                         readtminmaxobs,  &
                         bdate,sdate,     &
                         checkdata,       &
                         i_str,i_real,i_intg
    USE MODVAR, ONLY : forcing, &
                       basin,   &
                       nsub_basemodel  

    !Argument declarations
    
    !Local parameters
    INTEGER, PARAMETER :: maxcol = 3      !maximum number of column of ForcData.txt

    !Local variables
    CHARACTER(LEN=3)   seqnr
    CHARACTER(LEN=maxcharpath) filepath
    CHARACTER(LEN=10)  colstr(maxcol)
    TYPE(DateType) :: fbdate,fedate
    REAL    forcproc
    INTEGER i,j
    INTEGER mcols                         !number of columns in file
    INTEGER nrows                         !number of data rows in file
    INTEGER forcid
    INTEGER n_Result
    INTEGER itobsid                       !column for tobsid
    INTEGER itobselev                     !column for tobselev
    INTEGER :: code(maxcol)        !Code for column variable
    INTEGER :: rindex(maxcol)      !Index for column real variables
    INTEGER :: iindex(maxcol)      !Index for column integer variables
    INTEGER obsindex(nsub_basemodel)         !Index of observation stations
    INTEGER nobscol                !Number of columns with observations in file
    INTEGER numneg(nsub_basemodel)
    LOGICAL fileexist,found,notimefound
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file

    WRITE(6,*)
    WRITE(6,*) 'Check indata part two'
    WRITE(6,*) '---------------------'
    
    !Checking ForcKey and coupling between subbasin and observations
    WRITE(6,*) 'Checking ForcKey.txt...'
    IF(ALLOCATED(forcing))THEN
      forcproc = MINVAL(forcing(:)%tobselev)
      IF(forcproc<0.) WRITE(6,*) 'Elevation of forcing data (temperature) is below zero. Check ForcKey.txt.'
    ENDIF
    forcid = MINVAL(tobsid)
    IF(forcid<=0)THEN
      WRITE(6,*) 'Identification number of temperature forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
      checkdata(2,2) = .TRUE.
    ENDIF
    forcid = MINVAL(pobsid)
    IF(forcid<=0)THEN
      WRITE(6,*) 'Identification number of precipitation forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
      checkdata(2,2) = .TRUE.
    ENDIF
    IF(readsfobs)THEN
      forcid = MINVAL(sfobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of snowfall forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
    ENDIF    
    IF(readswobs)THEN
      forcid = MINVAL(swobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of shortwave radiation forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
    ENDIF
    IF(readwind)THEN
      forcid = MINVAL(uobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of wind speed forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
    ENDIF
    IF(readhumid)THEN
      forcid = MINVAL(rhobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of relative humidity forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
    ENDIF
    IF(readtminmaxobs)THEN
      forcid = MINVAL(tminobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of min temperature forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
      forcid = MINVAL(tmaxobsid)
      IF(forcid<=0)THEN
        WRITE(6,*) 'Identification number of max temperature forcing is zero for some subbasin. Check ForcKey.txt or GeoData.txt.'
        checkdata(2,2) = .TRUE.
      ENDIF
    ENDIF
    
    !Checking ForcData
    IF(simsequence>0)THEN                     !Set filename of files to be used
      WRITE(seqnr,'(I3.3)') simsequence
      filepath = TRIM(modeldir)//'ForcData_'//seqnr//'.txt'
    ELSE
      filepath = TRIM(modeldir)//'ForcData.txt'
    ENDIF
    INQUIRE(FILE=filepath,EXIST=fileexist)    !Check for file existance
    IF(fileexist)THEN
      WRITE(6,*) 'ForcData.txt exist, checking ForcData.txt...'
      CALL count_data_rows(fileunit_temp,filepath,1,nrows,n_Result)
      IF(n_Result/=0)THEN
        WRITE(6,*) 'Unknown error in ForcData.txt. Check file.'
        checkdata(2,2) = .TRUE.
      ELSE
        ALLOCATE(xi(nrows,maxcol))
        ALLOCATE(xr(nrows,maxcol))
        OPEN(UNIT=fileunit_temp,FILE=filepath,STATUS='old', ACTION='read',ERR=900)
        CALL read_column_headings(fileunit_temp,maxcol,colstr,mcols,n_Result)
        IF(n_Result/=0)THEN
          WRITE(6,*) 'Error in ForcData.txt (heading). Check file.'
          checkdata(2,2) = .TRUE.
        ELSE
          code = i_str
          itobsid = 0; itobselev = 0
          DO i = 1,mcols
            IF(colstr(i)(1:10)=='tobsid    ')THEN
              code(i) = i_intg
              itobsid = i
            ENDIF  
            IF(colstr(i)(1:10)=='tobselev  ')THEN
              code(i) = i_real
              itobselev = i
            ENDIF  
          ENDDO
          CALL read_basindata5(fileunit_temp,filepath,maxcol,nrows,mcols,code,rindex,iindex,xi,xr) 
          CLOSE(fileunit_temp)
          forcproc = MINVAL(xr(:,rindex(itobselev)))
          IF(forcproc<0.) WRITE(6,*) 'Elevation of forcing data (temperature) is below zero. Check ForcData.txt.'
          DO i = 1,nsub_basemodel
            found = .FALSE.
            DO j = 1,nrows
              IF(tobsid(i)==xi(j,iindex(itobsid)))THEN
                found = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF(.NOT.found)THEN
              WRITE(6,*) 'Temperature observation id (',tobsid(i),') for subbasin ',basin(i)%subid,' not found in ForcData. Check ForcData.txt'
              checkdata(2,2) = .TRUE.
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ELSE
      WRITE(6,*) 'ForcData.txt not found'
    ENDIF

    !Checking Tobs
    WRITE(6,*) 'Checking Tobs.txt...'
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      filepath = TRIM(modeldir)//'Tobs_'//seqnr//'.txt'
    ELSE
      filepath = TRIM(modeldir)//'Tobs.txt'
    ENDIF
    CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,tobsid,obsindex,nobscol,n_Result)
    IF(n_Result/=0)THEN
      WRITE(6,*) 'Error: Missmatch in observations between model and Tobs.txt. Check GeoData/ForcKey and Tobs.txt.'
      checkdata(2,2)=.TRUE.
    ENDIF 
    CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &   !(also checks timestep)
         sdate,fbdate,fedate,notimefound,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'Error: Temperature forcing data has problems with time period. Check Tobs.txt.'
      checkdata(2,2)=.TRUE.
    ENDIF

    !Checking Pobs
    WRITE(6,*) 'Checking Pobs.txt...'
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      filepath = TRIM(modeldir)//'Pobs_'//seqnr//'.txt'
    ELSE
      filepath = TRIM(modeldir)//'Pobs.txt'
    ENDIF
    CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,pobsid,obsindex,nobscol,n_Result)
    IF(n_Result/=0)THEN
      WRITE(6,*) 'Error: Missmatch in observations between model and Pobs.txt. Check GeoData/ForcKey and Pobs.txt.'
      checkdata(2,2)=.TRUE.
    ENDIF 
    CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
         sdate,fbdate,fedate,notimefound,n_Result)
    IF(n_Result.NE.0)THEN
      WRITE(6,*) 'Error: Precipitation forcing data has problems with time period. Check Pobs.txt.'
      checkdata(2,2)=.TRUE.
    ENDIF
    CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,n_Result)
    IF(n_Result==1)THEN
      WRITE(6,*) 'Negative precipitation data exist. Check Pobs.txt.'
      IF(SUM(numneg)>0)THEN 
        checkdata(2,2)=.TRUE.    !negative prec found for used stations
        WRITE(6,*) 'The following subbasins have negative precipitation:'
        DO i = 1,nsub_basemodel
          IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'pobsid:',pobsid(i)
        ENDDO
      ENDIF  
    ENDIF

    !Checking SFobs
    IF(readsfobs)THEN
      WRITE(6,*) 'Checking SFobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'SFobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'SFobs.txt'
      ENDIF
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath),EXIST=fileexist)
      IF(.NOT.fileexist)THEN
        WRITE(6,*) 'Error: Snowfall forcing data file is missing, while setting is set to read it.'
        WRITE(6,*) 'Error: ', TRIM(filepath)
        checkdata(2,2)=.TRUE.
      ELSE
        CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,sfobsid,obsindex,nobscol,n_Result)
        IF(n_Result/=0)THEN
          WRITE(6,*) 'Error: Missmatch in observations between model and SFobs.txt. Check GeoData/ForcKey and SFobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF 
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
             sdate,fbdate,fedate,notimefound,n_Result)
        IF(n_Result.NE.0)THEN
          WRITE(6,*) 'Error: Snowfall forcing data has problems with time period. Check SFobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,n_Result)
        IF(n_Result==1)THEN
          WRITE(6,*) 'Negative snowfall fraction data exist. Check SFobs.txt.'
          IF(SUM(numneg)>0)THEN 
            checkdata(2,2)=.TRUE.    !negative data found for used stations
            WRITE(6,*) 'The following subbasins have negative snowfall fraction:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'sfobsid:',sfobsid(i)
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    !Checking SWobs
    IF(readswobs)THEN
      WRITE(6,*) 'Checking SWobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'SWobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'SWobs.txt'
      ENDIF
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath),EXIST=fileexist)
      IF(.NOT.fileexist)THEN
        WRITE(6,*) 'Error: Shortwave radiation data file is missing, while setting is set to read it.'
        WRITE(6,*) 'Error: ', TRIM(filepath)
        checkdata(2,2)=.TRUE.
      ELSE
        CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,swobsid,obsindex,nobscol,n_Result)
        IF(n_Result/=0)THEN
          WRITE(6,*) 'Error: Missmatch in observations between model and SWobs.txt. Check GeoData/ForcKey and SWobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF 
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
             sdate,fbdate,fedate,notimefound,n_Result)
        IF(n_Result.NE.0)THEN
          WRITE(6,*) 'Error: Shortwave radiation data has problems with time period. Check SWobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,n_Result)
        IF(n_Result==1)THEN
          WRITE(6,*) 'Negative Shortwave radiation data exist. Check SWobs.txt.'
          IF(SUM(numneg)>0)THEN 
            checkdata(2,2)=.TRUE.    !negative data found for used stations
            WRITE(6,*) 'The following subbasins have negative shortwave radiation:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'swobsid:',swobsid(i)
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    !Checking Uobs
    IF(readwind)THEN
      WRITE(6,*) 'Checking Uobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'Uobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'Uobs.txt'
      ENDIF
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath),EXIST=fileexist)
      IF(.NOT.fileexist)THEN
        WRITE(6,*) 'Error: Wind speed data file is missing, while setting is set to read it.'
        WRITE(6,*) 'Error: ', TRIM(filepath)
        checkdata(2,2)=.TRUE.
      ELSE
        CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,uobsid,obsindex,nobscol,n_Result)
        IF(n_Result/=0)THEN
          WRITE(6,*) 'Error: Missmatch in observations between model and Uobs.txt. Check GeoData/ForcKey and Uobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF 
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
             sdate,fbdate,fedate,notimefound,n_Result)
        IF(n_Result.NE.0)THEN
          WRITE(6,*) 'Error: Wind speed data has problems with time period. Check Uobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,n_Result)
        IF(n_Result==1)THEN
          WRITE(6,*) 'Negative wind speed data exist. Check Uobs.txt.'
          IF(SUM(numneg)>0)THEN 
            checkdata(2,2)=.TRUE.    !negative data found for used stations
            WRITE(6,*) 'The following subbasins have negative wind speed:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'uobsid:',uobsid(i)
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    !Checking RHobs
    IF(readhumid)THEN
      WRITE(6,*) 'Checking RHobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'RHobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'RHobs.txt'
      ENDIF
      !Check that file exist
      INQUIRE(FILE=TRIM(filepath),EXIST=fileexist)
      IF(.NOT.fileexist)THEN
        WRITE(6,*) 'Error: Relative humidity data file is missing, while setting is set to read it.'
        WRITE(6,*) 'Error: ', TRIM(filepath)
        checkdata(2,2)=.TRUE.
      ELSE
        CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,rhobsid,obsindex,nobscol,n_Result)
        IF(n_Result/=0)THEN
          WRITE(6,*) 'Error: Missmatch in observations between model and RHobs.txt. Check GeoData/ForcKey and RHobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF 
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
             sdate,fbdate,fedate,notimefound,n_Result)
        IF(n_Result.NE.0)THEN
          WRITE(6,*) 'Error: Relative humidity data has problems with time period. Check RHobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,n_Result)
        IF(n_Result==1)THEN
          WRITE(6,*) 'Negative relative humidity data exist. Check Uobs.txt.'
          IF(SUM(numneg)>0)THEN 
            checkdata(2,2)=.TRUE.    !negative data found for used stations
            WRITE(6,*) 'The following subbasins have negative relative humidity:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'rhobsid:',rhobsid(i)
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    !Checking TMINobs and TMAXobs
    IF(readtminmaxobs)THEN
      WRITE(6,*) 'Checking TMINobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'TMINobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'TMINobs.txt'
      ENDIF
      CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,tminobsid,obsindex,nobscol,n_Result)
      IF(n_Result/=0)THEN
        WRITE(6,*) 'Error: Missmatch in observations between model and TMINobs.txt. Check GeoData/ForcKey and TMINobs.txt.'
        checkdata(2,2)=.TRUE.
      ENDIF 
      CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
           sdate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'Error: Temperature forcing data has problems with time period. Check TMINobs.txt.'
        checkdata(2,2)=.TRUE.
      ENDIF
      WRITE(6,*) 'Checking TMAXobs.txt...'
      IF(simsequence>0)THEN
        WRITE(seqnr,'(I3.3)') simsequence
        filepath = TRIM(modeldir)//'TMAXobs_'//seqnr//'.txt'
      ELSE
        filepath = TRIM(modeldir)//'TMAXobs.txt'
      ENDIF
      CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,tmaxobsid,obsindex,nobscol,n_Result)
      IF(n_Result/=0)THEN
        WRITE(6,*) 'Error: Missmatch in observations between model and TMAXobs.txt. Check GeoData/ForcKey and TMAXobs.txt.'
        checkdata(2,2)=.TRUE.
      ENDIF 
      CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
           sdate,fbdate,fedate,notimefound,n_Result)
      IF(n_Result.NE.0)THEN
        WRITE(6,*) 'Error: Temperature forcing data has problems with time period. Check TMAXobs.txt.'
        checkdata(2,2)=.TRUE.
      ENDIF
    ENDIF

    !Checking Qobs
    WRITE(6,*) 'Checking Qobs.txt...'
    filepath = TRIM(modeldir)//filename_Qobs
    INQUIRE(FILE=filepath,EXIST=fileexist)
    IF(fileexist)THEN
      CALL check_q_stn(fileunit_temp,filepath,nsub_basemodel,   &
                       basin(:)%subid,nobscol,obsindex,n_Result)
      IF(n_Result==0)THEN
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,   &
            bdate,sdate,fbdate,fedate,notimefound,n_Result)
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,1,numneg,n_Result)
        IF(n_Result==1)THEN
          WRITE(6,*) 'Negative discharge data exist. Check Qobs.txt.'
          IF(SUM(numneg)>0)THEN 
            WRITE(6,*) 'The following subbasins have negative discharge:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

    WRITE(6,*) '---------------------'
    RETURN

900   WRITE(6,*) 'Error open file: ', TRIM(filepath)
    checkdata(2,2) = .TRUE.
    WRITE(6,*) '---------------------'

  END SUBROUTINE checkindata_part2

  !>Check for negative data in file maybe also for missing values (Pobs.txt, Qobs.txt)
  !------------------------------------------------------------
  SUBROUTINE check_data_positive(filepath,notimefound,ncols,ns,oindex,miss,numneg,n_Result)

    USE WORLDVAR, ONLY : readmatlab,      &
                         fileunit_temp   
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: filepath  !<File (Pobs.txt)
    LOGICAL, INTENT(IN)  :: notimefound  !<Flag for time format: date or date and time
    INTEGER, INTENT(IN)  :: ncols        !<Number of data columns in file
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: oindex(ns)   !<Index for columns used in model set-up
    INTEGER, INTENT(IN)  :: miss         !<Flag for counting missing values, 0=count, 1=ignore
    INTEGER, INTENT(OUT) :: numneg(ns)   !<Number of negative values per subbasin
    INTEGER, INTENT(OUT) :: n_Result     !<Number of/code for found negative data
    
    !Local variables
    INTEGER timeform
    INTEGER i
    INTEGER :: neg(ncols)
    REAL :: y(ncols)                 !Data (one time step)
    CHARACTER(LEN=16)  d2,d3         !Date yyyy-mm-dd[ hh:mm]

    IF(notimefound)THEN
      timeform = 0
    ELSE
      timeform = 1
    ENDIF

    !Load precipitation forcing data
    OPEN(UNIT = fileunit_temp,FILE = filepath, STATUS = 'old', ACTION='read')
    READ(fileunit_temp,*)  !Skip heading
    neg = 0
    numneg = 0
    DO
      y = missing_value
      IF(readmatlab)THEN
        READ(fileunit_temp,*,END=900) d2,y
      ELSEIF(timeform==0)THEN
        READ(fileunit_temp,*,END=900) d2,y
      ELSEIF(timeform==1)THEN
        READ(fileunit_temp,*,END=900) d2,d3,y    
      ENDIF
      DO i = 1,ncols
        IF(miss==0)THEN
          IF(y(i)<0.) neg(i) = neg(i) + 1
        ELSEIF(miss==1)THEN
          IF(y(i)<0..AND.y(i)/=missing_value) neg(i) = neg(i) + 1
        ENDIF
      ENDDO  
    ENDDO
900 IF(SUM(neg)==0.)THEN
      n_Result = 0    !no negative data in file
    ELSE
      n_Result = 1
      DO i = 1,ns
        IF(oindex(i)>0)THEN
          numneg(i) = neg(oindex(i))
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE check_data_positive

  !>Check parameter file
  !------------------------------------------------------------
  SUBROUTINE checkindata_stop() 

    USE WORLDVAR, ONLY : checkdata

    IF(checkdata(1,2).OR.checkdata(2,2))THEN
      WRITE(6,*)
      WRITE(6,*) 'Simulation stopped'
      WRITE(6,*) 'Serious errors in indata found'
      STOP 1
    ENDIF

  END SUBROUTINE checkindata_stop

END MODULE
