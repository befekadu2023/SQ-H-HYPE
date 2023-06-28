!> \file main.f90
!> Contains main program of HYSS-HYPE

!> Main program for HYSS - Hydrological Simulation System
!> 
!> The system features hydrological simulation of soil and surface water system, 
!> calibration of parameters and criteria calculations, updating of flow and state 
!> to observed values, ensemble simulation, and more.
PROGRAM MAIN

!Copyright 2011-2015 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------

!Used modules
!      USE IFQWIN                               !for nobutton
      USE MODELMODULE, ONLY : model, &
                              model_version_information, &
                              initiate_output_variables, &
                              initiate_model,            &
                              initiate_model_state, &
                              initiate_model_parameters, &
                              load_modeldefined_input, &
                              set_modelconfig_from_parameters
      USE STATETYPE_MODULE
      USE WORLDVAR, ONLY :  writematlab,       &
                            writeload,         &
                            outvarinfo,        &
                            o_nout,            &
                            nacrit,            &
                            ndt,               &
                            fileunit_temp,     &
                            allocate_accumulation,     &
                            numqobsstn,        &
                            maxcharpath,       &
                            simsequence,       &
                            infodir,           &
                            modeldir,          &
                            resdir,            &
                            bdate,             &
                            sdate,             &
                            outstartdate,      &
                            dtskip,            &
                            yearid,            &
                            doopt,             &
                            numoptimpar,       &   
                            deallocate_worldvar,      &
                            deallocate_MCvariables,   &
                            optim,             &
                            bestMCoptcrit,     &
                            bestMCperformance, &
                            bestMCparameters,  &
                            maxperf,           &
                            maxsubass,         &
                            simsubmodel,       &
                            ibasemodel,     &
                            readsfobs,readswobs,  &
                            readwind, &
                            readhumid, &
                            readtminmaxobs, &
                            checkdata,  &
                            psdates,    &
                            lddates,    &
                            dddates,    &
                            optimStartTime, &      
                            optimFuncCall,  &
                            lineSearchCallCount
      USE MODVAR, ONLY : preci,tempi,qobsi,xobsi,  &
                         snowfraci,shortwavei,     &
                         windi,humidi, &
                         tmini,tmaxi,          &
                         nsub,ncrop,        &
                         nsub_basemodel,       &
                         nclass,numsubstances, &
                         naquifers,   &
                         maxsoillayers,   &
                         max_classoutvar,   &
                         max_basinoutvar,   &
                         doupdate,i_qar,i_war,  &
                         conductN,conductP,conductC,  &
                         i_t2, &
                         wetlandexist,glacierexist,  &
                         doirrigation,  &
                         modeloption,p_lakeriverice, &
                         p_growthstart,  &
                         allocate_outvar,   &
                         deallocate_modvar,  &
                         dimriverlag, &
                         timesteps_per_day
      USE COMPOUT, ONLY : compute_mapvar,           &
                          compute_outloads,         &
                          prepare_to_compute_crit,  &
                          calculate_criteria
      USE TIMEROUTINES, ONLY : calculate_time_for_model
      USE READWRITE_ROUTINES
      USE LIBDATE, ONLY : DateType, OPERATOR(.EQ.)
      USE DATAMODULE
      USE STATE_DATAMODULE, ONLY : initiate_state_for_submodel, &
                                   load_saved_state,  &
                                   finalize_outstate
								   
		!Modules added by BTW for state-space continuous simulation
		USE gear_GlobVARs
		USE GEAR_IMPLICIT
		USE uawp, ONLY : dist_max, norm_max, awp
							 
		USE t_dgls, ONLY : 	dgl						!!to be populated
		USE fgauss, ONLY : gauss
		!End of Modules added by BTW for state-space continuous simulation

      IMPLICIT NONE

!Parameter declarations
      INTEGER, PARAMETER :: maxoutstates = 10     !Max number of dates for saving state

!Variable declarations
      TYPE(DateType) d            !Current time (e.g. 2005-01-26 18:00)
      INTEGER idt                 !Current time step number
      INTEGER ivar                !Current output variable
      INTEGER iens                !Current ensemble being simulated
!      INTEGER nobutton            !No exit window
      LOGICAL pwrite              !Flag for periodend, time to write to file
      CHARACTER(LEN=maxcharpath+25) filename  !hyss filename
      CHARACTER(LEN=8)  :: logdate  !Date for log-file name
      CHARACTER(LEN=10) :: logtime  !Time for log-file name
      CHARACTER(LEN=3)  :: logseq   !Seqnr for log-file name
      INTEGER :: datim(8) 
      INTEGER :: oldyear            !yearid of last time step

      REAL, ALLOCATABLE :: par(:) !,parmin(:),parmax(:)
      REAL optcrit, condcrit, condthres
      REAL, ALLOCATABLE :: basincrit(:,:,:)   !R2, CC, RE, RSDE, QC, QR, STDC, STDR, MAE, RMSE, Bias, STDbias, KGE, KGEpartSTD, KGEpartMM, NRMSE per subbasin och kriterie
      REAL, ALLOCATABLE :: simperformance(:,:)   !rr2,sr2,wr2,rmae,sbias,rrve,wrve,rra,sra,meanRA,tau,medianr2,medianra,meanrs,meancc,mediankg,meanabsre
      INTEGER :: n_Result = 0
      INTEGER npar              !Number of parameters to be calibrated (couted in file)
      INTEGER :: nsubout        !Number of subbasins for basinoutput
      INTEGER :: nmapperiod     !Number of periods for map print out
      LOGICAL :: stateinput                               !Code for reading state
      INTEGER :: numoutstates                             !Number of dates for saving state
      TYPE(DateType) :: stateoutadate(maxoutstates)       !Dates for saving state
      TYPE(DateType) :: prestateoutdate(maxoutstates)     !Day before date for saving state (or 0)
      
!Variables for updating of Q and W
      LOGICAL :: quseobsallstations, quseobsnostations
      LOGICAL :: qarnostations, warnostations
      LOGICAL :: wendupdallstations, wendupdnostations
      CHARACTER(LEN=4) :: wobsvarname
      
!Model state variables and other saved variables declaration
      TYPE(SNOWICESTATETYPE) :: frozenstate
      TYPE(SOILSTATETYPE)    :: soilstate      
      TYPE(AQUIFERSTATETYPE) :: aquiferstate      
      TYPE(RIVERSTATETYPE)   :: riverstate      
      TYPE(LAKESTATETYPE)    :: lakestate
      TYPE(MISCSTATETYPE)    :: miscstate

!Program start

      CALL DATE_AND_TIME (logdate, logtime,values=datim)

!Current model domain
      CALL get_hyss_arguments(infodir,simsequence)   
      WRITE(logseq,'(I3.3)') simsequence
      WRITE(filename,'(a)') TRIM(infodir)//'hyss_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:4)//'.log'
      OPEN(UNIT=6,FILE=TRIM(filename))
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job start date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)
      WRITE(6,*) '---------------------------------------------------'
!      nobutton = SETEXITQQ(qwin$exitnopersist)    !nobutton version
      CALL model_version_information(6)   !Collect and print model version information
      CALL initiate_output_variables()    !Define model's output variables
      CALL initiate_model_parameters()    !Define model's parameters
!Information for this simulation
      CALL load_coded_info(infodir,n_Result,bdate,sdate,outstartdate,dtskip,ndt,numsubstances,  &
                           stateinput,maxoutstates,numoutstates,stateoutadate,prestateoutdate,  &
                           modeldir,resdir,nsubout, &
                           quseobsallstations,quseobsnostations,qarnostations,warnostations,    &
                           wendupdallstations,wendupdnostations,wobsvarname)
      IF(n_Result/=0)THEN
        STOP 1
      ENDIF

!Read input data
      CALL load_cropdata(modeldir,'CropData.txt',ncrop,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_basindata(modeldir,'GeoClass.txt','GeoData.txt',nsub_basemodel,n_Result,lakestate) !MH2017: pass lakestate for conditioning initialization
      IF(n_Result.NE.0) STOP 1
      CALL load_pointsourcedata(modeldir,'PointSourceData.txt',nsub_basemodel,n_Result) 
      IF(n_Result.NE.0) STOP 1
      CALL load_branchdata(modeldir,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_aquiferdata(modeldir,nsub_basemodel,naquifers,n_Result) 
      IF(n_Result.NE.0) STOP 1
      CALL load_glacierdata(modeldir,nsub_basemodel,n_Result) 
      IF(n_Result.NE.0) STOP 1
      IF(checkdata(1,1)) CALL checkindata_part1()
      IF(checkdata(2,1)) CALL checkindata_part2()
      CALL checkindata_stop()

      CALL load_submodel_info(infodir,simsubmodel,nsub,n_Result)   !allocation and initialisation of ibasemodel
      IF(n_Result/=0) STOP 1

      CALL load_observations(modeldir,nsub_basemodel,bdate,sdate,ndt,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_parameters(modeldir,nsub_basemodel,'par.txt')

      IF(simsubmodel)THEN
        CALL reform_inputdata_for_submodel(nsub_basemodel,nsub,ibasemodel)
      ENDIF
      CALL calculate_path(nsub)

!Read model defined input data and set model configuration
      CALL load_modeldefined_input(modeldir,nsub_basemodel,nsub,ibasemodel,bdate,sdate,n_Result)
      IF(n_Result.NE.0) STOP 1

      CALL set_model_configuration(nsub,infodir)

      CALL prepare_for_update(infodir,wobsvarname,quseobsallstations, &
              quseobsnostations,qarnostations,warnostations,wendupdallstations,wendupdnostations,nsub)

!Initialisations for memory allocation
      CALL allocate_model_states(nsub_basemodel,numsubstances,nclass,naquifers, &
              maxsoillayers,2,2,dimriverlag,timesteps_per_day,conductN,conductP,conductC, &
              i_t2>0,wetlandexist,glacierexist,modeloption(p_lakeriverice)>=1,doirrigation, &
              doupdate(i_qar).OR.doupdate(i_war),modeloption(p_growthstart)==1, &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      CALL allocate_outvar(nsub,nclass,numsubstances)
      CALL allocate_accumulation(outvarinfo(:,o_nout),nsub,nsubout,nclass,numsubstances,max_classoutvar,max_basinoutvar)

!Allocate local variables
      ALLOCATE(basincrit(nsub,maxsubass,nacrit))
      ALLOCATE(simperformance(maxperf,nacrit))

!Preparations for subbasin output
      CALL prepare_subbasin_output()

      CALL DATE_AND_TIME (values=datim)
      WRITE(6,*)
      WRITE(6,*) '---------------------------------------------------'
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
                 ' Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
                 '  time: ',datim(5),':',datim(6),':',datim(7)
      WRITE(6,*) '---------------------------------------------------'


!Optimization
      IF(doopt)THEN
        CALL load_optpar(infodir)       !Reads optpar.txt and set all numerical optimization variables accordingly

        IF(optim%task_MC)THEN
          WRITE(6,*) 'MonteCarlo simulation with',optim%nruns_MC,'runs'
          CALL MonteCarlo_simulation(resdir,optim%task_writeall,frozenstate,  &
               soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
          CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
          CALL save_respar(resdir,npar,nsub)                       !and to file
        ENDIF

        IF(optim%task_boundps)THEN
          WRITE(6,*) 'Reduce bounds of parameter space for continued MonteCarlo simulation'
          CALL bounded_MonteCarlo_simulation(optim%task_MC,frozenstate,soilstate, &
               aquiferstate,riverstate,lakestate,miscstate,npar)
          CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
          CALL save_respar(resdir,npar,nsub)                       !and to file
        ENDIF

        lineSearchCallCount = 0
        optimFuncCall = 0
        CALL cpu_time(optimStartTime)

        IF(optim%task_Scanning) CALL param_scanning(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

        IF(optim%task_stageMC)THEN
          CALL stage_MonteCarlo(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
          CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(1,:))    !save the best to modpar
          CALL save_respar(resdir,numoptimpar,nsub)                               !and file
        ENDIF
        
        IF(optim%task_BrentNew .OR. optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS)THEN
          ALLOCATE(par(numoptimpar))
          CALL linesearch_methods_calibration(numoptimpar,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,par)
          CALL set_optim_modpar(numoptimpar,numoptimpar,par)    !save the best to modpar
          CALL save_respar(resdir,numoptimpar,nsub)                            !and file
        ENDIF
        IF(optim%task_DEMC)THEN
          WRITE(6,*) 'DEMC Differential-Evolution Markov Chain, with',optim%DEMC_npop, &
          'populations, and ',optim%DEMC_ngen,'generations (',optim%DEMC_ngen-1,  &
          'evolution steps)',', in total',optim%DEMC_npop*(optim%DEMC_ngen),'simulations.'
          
          CALL DEMC_simulation(resdir,optim%task_writeall,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
          CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))         !save the MEDIAN parameters to modpar (1 row in bestMC)
          CALL save_respar(resdir,npar,nsub)        
        ENDIF
        IF(optim%task_DREAM)THEN
          WRITE(6,*) 'DREAM DiffeRential-Evolution Adaptive Metropolis simulation, with',optim%DREAM_npop,'chains, and ' &
          ,optim%DREAM_ngen,'generations (',optim%DREAM_ngen-1,'evolution steps)',', in total',optim%DREAM_npop*(optim%DREAM_ngen),'model iterations.'
          WRITE(6,*) ''
        !  WRITE(6,*) ' DREAM-OUTPUTS:'
        !  WRITE(6,*) '    1) the mean of the last 50% of the chains (after burnin) is reported in RESPAR.TXT'
        !  WRITE(6,*) '    2) a simulation with the RESPAR.TXT parameters is printed in ensemble output files no. 1'
        !  WRITE(6,*) '    3) simulations with the BEST chain members are printed in ensemble output files no. 2:nchains+1'
        !  WRITE(6,*) '    4) identification and replacement of outlayer chains is reported in DREAM_OUTLIER.txt'
        !  WRITE(6,*) '    5) Criterias and parameter values for all proposals are reported in ALLSIM.TXT (similar to DE-MC)'
        !  WRITE(6,*) '    6) Parameter chains (only accepted proposals) are reported in DREAM_PARCHAINS_###.TXT after each completed generation'
        !  WRITE(6,*) '    7) Gelman-Rubin convergence diagnostics for each parameter are reported in DREAM_GELMANRUBIN.TXT'
        !  WRITE(6,*) '    8) The tuned Crossover value probabilities are reported in DREAM_CRPROB.TXT'
        !
        !  WRITE(6,*) ''
          
          CALL DREAM_simulation(resdir,optim%task_writeall,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
          CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))         !save the MEAN parameters based on the last half of the DREAM chains to modpar (1 row in bestMC)
          CALL save_respar(resdir,npar,nsub)        
        ENDIF
        IF(optim%task_runens)THEN       !write best ensemble resultat to file
          CALL save_ensemble_simulations(resdir,numoptimpar,maxperf,nacrit, &
               optim%nruns_best,bestMCoptcrit,bestMCperformance,bestMCparameters)
        ENDIF
      ENDIF

      IF(.NOT.optim%task_writesim)THEN
!Ensemble loop, simulate all ensemble members
        ensemble_loop:    &
     &  DO iens = 1, optim%nruns_best

  !Simulation start
          CALL prepare_outputfiles(resdir,nsub,naquifers,iens,optim%task_runens,optim%task_writesim)

  !        IF(doopt .AND. (optim%task_runens .OR. optim%task_DEMC))THEN
          IF(doopt .AND. optim%task_runens)THEN
            CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(iens,:))   !set model parameters
          ENDIF
          CALL initiate_output_routines(outvarinfo(:,o_nout))     !All output accumulation variables zeroed
          optcrit=0.

          !Initial model calculations; initial states, parameters
          CALL set_modelconfig_from_parameters() !set some model config parameters, for instance basin%lakedepth from parameter inputs
          
          IF(simsubmodel)THEN
            CALL initiate_state_for_submodel(modeldir,ibasemodel,dimriverlag,stateinput,  &
                frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
          ELSE
            IF(stateinput) THEN
              CALL load_saved_state(modeldir,nsub,dimriverlag,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
              WRITE(6,*)
              WRITE(6,*) 'Loading saved state.'
            ELSE
              CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            ENDIF
            CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
          ENDIF
          oldyear = 0

  !Time Loop:
          time_loop:    &
     &    DO idt = 1,ndt

  !Get current input data
            CALL get_current_date_prec(idt,nsub,d,preci(1:nsub))      !TEST för bättre körtid/minnesanvändning
            CALL calculate_time_for_model(idt,d)
            CALL log_progress(oldyear,yearid)
            tempi(1:nsub) = get_current_temp(idt,nsub)
            IF(ALLOCATED(qobsi)) qobsi(1:nsub) = get_current_flow(idt,d,nsub,numqobsstn)    !m3/s
            IF(ALLOCATED(xobsi)) xobsi(:) = get_current_otherobs(idt,d)
            IF(readsfobs) snowfraci(1:nsub) = get_current_snowfallfrac(idt,nsub)
            IF(readswobs) shortwavei(1:nsub) = get_current_shortwave(idt,nsub)
            IF(readwind)  windi(1:nsub) = get_current_windspeed(idt,nsub)
            IF(readhumid) humidi(1:nsub) = get_current_humidity(idt,nsub)
            IF(readtminmaxobs) tmini(1:nsub) = get_current_tmin(idt,nsub)
            IF(readtminmaxobs) tmaxi(1:nsub) = get_current_tmax(idt,nsub)
            IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,n_Result) 
            IF(n_Result.NE.0) STOP 1
!            IF(ALLOCATED(lddates)) CALL get_current_lakedata(modeldir,'LakeData.txt',nsub,d,n_Result) 
            IF(n_Result.NE.0) STOP 1
!            IF(ALLOCATED(dddates)) CALL get_current_damdata(modeldir,'DamData.txt',nsub,d,n_Result) 
            IF(n_Result.NE.0) STOP 1

            CALL initiate_outvar(idt)
            CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,d,idt)
        
            DO ivar = 1,numoutstates    !Write state output
              IF (d.EQ.prestateoutdate(ivar)) THEN
                IF(simsubmodel)THEN
                  WRITE(6,*) 'State can not be saved when submodel is simulated'
                ELSE
                  CALL finalize_outstate(resdir,dimriverlag,stateoutadate(ivar),  &
                       frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
                ENDIF
              ENDIF
            ENDDO

            CALL revise_outvar()
            CALL prepare_to_compute_crit(d,idt,ndt)

  !Calculate and write time dependent output
            CALL write_subbasinfiles(idt,ndt,d)
            CALL compute_mapvar(d,idt,ndt,outvarinfo(2,o_nout),writematlab,nmapperiod)  !Save data for map output
            CALL write_timefiles(idt,ndt,iens,d)
            IF(writeload)THEN
              CALL compute_outloads(d,pwrite,idt,ndt)       !Write yearly load total for all subbasins
              IF(pwrite) CALL save_loadfiles(resdir,yearid)
            ENDIF

          ENDDO time_loop

  !Compute and write criteria
          IF(nacrit/=0) THEN
            CALL calculate_criteria(optcrit,basincrit,simperformance,condcrit,condthres)
            CALL write_simulation_assessment(resdir,iens,nacrit,optcrit,     &
                simperformance,optim%task_runens,condcrit,condthres)
            CALL write_subbasin_assessment(resdir,nsub,basincrit,iens,optim%task_runens)
          ENDIF

  !Save and close files or prepare them for next ensemble member simulation
          CALL close_outputfiles(nsub,naquifers,iens)
          IF(iens==optim%nruns_best)THEN
            CALL close_observations(modeldir)
          ELSE
            CALL reset_observations(modeldir,n_Result)
            IF(n_Result.NE.0) STOP 1
          ENDIF

  !Write results to files
          IF(outvarinfo(2,o_nout)>0) CALL save_mapfiles(resdir,nsub,nmapperiod,iens,optim%task_runens,optim%task_writesim)

        ENDDO ensemble_loop   !simulate ensembles
      ENDIF

!Deallocate variables
      CALL deallocate_worldvar()
      CALL deallocate_modvar(nsub)
      CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      IF(doopt.AND.optim%task_MC) CALL deallocate_MCvariables()
      DEALLOCATE(basincrit,simperformance)
      IF(ALLOCATED(par))    DEALLOCATE(par)

!Write stop time on log-file
      CALL DATE_AND_TIME (values=datim)
      WRITE(6,*)
      WRITE(6,*) '---------------------------------------------------'
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job finished date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)

      CLOSE(6)
      STOP 84
      END PROGRAM
