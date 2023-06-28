!> \file model_hype.f90
!> Contains module modelmodule.

!>Main module for hydrological model HYPE.
!>
!> Detailed information of HYPE is found in the <a href=../HYPEmodelDescription.odt>HYPE model description</a>
MODULE MODELMODULE
!The HYPE model (HYdrological Predictions for the Environment)

!Copyright 2011-2016 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!----------------------------------------------------------------------

!Used modules
  USE STATETYPE_MODULE  !model state variable types
  USE MODVAR            !model variables, general model
  USE HYPEVARIABLES     !model variables, HYPE model
  USE HYPE_WATERBALANCE !water balance variables and indices
!Subroutines also uses GENERAL_WATER_CONCENTRATION, GLACIER_SOILMODEL, SOILMODEL_DEFAULT, 
!SOIL_PROCESSES, NPC_SOIL_PROCESSES, SURFACEWATER_PROCESSES, NPC_SURFACEWATER_PROCESSES, 
!IRRIGATION_MODULE, REGIONAL_GROUNDWATER_MODULE

  IMPLICIT NONE
  PRIVATE

  ! Private subroutines
  !--------------------------------------------------------------------
  ! apply_quseobs 
  ! apply_wendupd
  ! apply_nutrientcorr
  ! apply_qarupd
  !--------------------------------------------------------------------
  PUBLIC :: model_version_information, &
            initiate_output_variables, &
            initiate_model_parameters, &
            initiate_model_state, &
            initiate_model,  &
            set_special_models, &
            set_special_parameters_from_model, &
            calculate_special_model_parameters,  &
            model, &
            load_modeldefined_input, &
            reload_modeldefined_observations, &
            open_modeldefined_outputfiles, &
            close_modeldefined_outputfiles, &
            set_modelconfig_from_parameters

  CONTAINS
  
  !>Information about the model version to print in log-file
  !-----------------------------------------------------------
  SUBROUTINE model_version_information(funit)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: funit !<fileunit for log-file
    
    WRITE(funit,600) '----------------------'
    WRITE(funit,600) ' HYPE version 4.10.7  '
    WRITE(funit,600) ' Branch CryoProcesses '
    WRITE(funit,600) '----------------------'
600 FORMAT(A22)
    
  END SUBROUTINE model_version_information
  
  !>Set variables holding output information from the HYPE model; 
  !!outvarid and loadheadings
  !-----------------------------------------------------------------
  SUBROUTINE initiate_output_variables()
  
    max_outvar_trunk = 323
    max_outvar = max_outvar_trunk + 46    !number of output variables
    
    ALLOCATE(outvarid(max_outvar))
    
    outvarid(1)  = idtype('crun','computed runoff'   ,i_sum , 0,'mm'   ,'mm'             ,'mapCRUN.txt' ,'timeCRUN.txt')
    outvarid(2)  = idtype('rrun','recorded runoff'   ,i_sum , 0,'mm'   ,'mm'             ,'mapRRUN.txt' ,'timeRRUN.txt')
    outvarid(3)  = idtype('prec','precipitation'     ,i_sum , 0,'mm'   ,'mm'             ,'mapPREC.txt' ,'timePREC.txt')
    outvarid(4)  = idtype('temp','air temperature'   ,i_mean , 0,'deg'  ,'degree Celsius','mapTEMP.txt' ,'timeTEMP.txt')
    outvarid(5)  = idtype('coT1','computed T1 runoff',i_wmean, 1,'?'    ,'?'             ,'mapCOT1.txt','timeCOT1.txt')
    outvarid(6)  = idtype('coT2','computed T2 runoff',i_wmean, 1,'deg' ,'degree Celsius' ,'mapCOT2.txt','timeCOT2.txt')
    outvarid(7)  = idtype('coT3','computed T3 runoff',i_wmean, 1,'?'    ,'?'             ,'mapCOT3.txt','timeCOT3.txt')
    outvarid(o_cprecT1) = idtype('cpT1','recorded T1 precip',i_wmean, 3,'?','?'          ,'mapCPT1.txt'  ,'timeCPT1.txt')
    outvarid(9)  = idtype('ceT1','computed T1 evap'  ,i_wmean,15,'?'    ,'?'             ,'mapCET1.txt'  ,'timeCET1.txt')
    outvarid(10) = idtype('soim','computed soil water',i_mean, 0,'mm'   ,'mm'            ,'mapSOIM.txt' ,'timeSOIM.txt')
    outvarid(11) = idtype('csT1','computed T1 soil'  ,i_wmean, 10,'?'   ,'?'             ,'mapCST1.txt'  ,'timeCST1.txt')
    outvarid(12) = idtype('csT2','computed T2 soil'  ,i_wmean, 10,'deg','degree Celsius' ,'mapCST2.txt'  ,'timeCST2.txt')
    outvarid(13) = idtype('csT3','computed T3 soil'  ,i_wmean, 10,'?'   ,'?'             ,'mapCST3.txt'  ,'timeCST3.txt')
    outvarid(14) = idtype('snow','snow pack'         ,i_mean , 0,'mm'   ,'mm water'      ,'mapSNOW.txt' ,'timeSNOW.txt')
    outvarid(15) = idtype('evap','subbasin evaporation',i_sum , 0,'mm'   ,'mm'            ,'mapEVAP.txt' ,'timeEVAP.txt')    !DG: this should be the total including iterception losses
    outvarid(16) = idtype('reT1','recorded T1 outflow olake',i_wmean, 54,'?'    ,'?'             ,'mapRET1.txt'  ,'timeRET1.txt')
    outvarid(17) = idtype('reT2','recorded T2 outflow olake',i_wmean, 54,'deg'  ,'degree Celsius','mapRET2.txt'  ,'timeRET2.txt') 
    outvarid(18) = idtype('reT3','recorded T3 outflow olake',i_wmean, 54,'?'    ,'?'             ,'mapRET3.txt'  ,'timeRET3.txt')
    outvarid(19) = idtype('gwat','computed grw level',i_mean , 0,'m'    ,'m'             ,'mapGWAT.txt'  ,'timeGWAT.txt')
    outvarid(20) = idtype('coIN','computed IN runoff',i_wmean, 1,'ug/L' ,'ug Inorg-N/L'  ,'mapCOIN.txt','timeCOIN.txt')
    outvarid(21) = idtype('coON','computed ON runoff',i_wmean, 1,'ug/L' ,'ug Org-N/L'    ,'mapCOON.txt','timeCOON.txt')
    outvarid(22) = idtype('coSP','computed SP runoff',i_wmean, 1,'ug/L' ,'ug SRP-P/L'    ,'mapCOSP.txt','timeCOSP.txt')
    outvarid(23) = idtype('coPP','computed PP runoff',i_wmean, 1,'ug/L' ,'ug PartP-P/L'  ,'mapCOPP.txt','timeCOPP.txt')
    outvarid(24) = idtype('reIN','recorded IN outflow olake',i_wmean, 54,'ug/L' ,'ug INORG-N/L'  ,'mapREIN.txt','timeREIN.txt')
    outvarid(25) = idtype('reON','recorded ON outflow olake',i_wmean, 54,'ug/L' ,'ug ORG-N/L'    ,'mapREON.txt','timeREON.txt')
    outvarid(26) = idtype('reSP','recorded SP outflow olake',i_wmean, 54,'ug/L' ,'ug SRP-P/L'    ,'mapRESP.txt','timeRESP.txt')
    outvarid(27) = idtype('rePP','recorded PP outflow olake',i_wmean, 54,'ug/L' ,'ug PartP-P/L'  ,'mapREPP.txt','timeREPP.txt')
    outvarid(28) = idtype('reTN','recorded TN outflow olake',i_wmean, 54,'ug/L' ,'ug Tot-N/L'    ,'mapRETN.txt','timeRETN.txt')
    outvarid(29) = idtype('reTP','recorded TP outflow olake',i_wmean, 54,'ug/L' ,'ug Tot-P/L'    ,'mapRETP.txt','timeRETP.txt')
    outvarid(30) = idtype('qerr','outflow error',i_mean , 0,'m3/s','m3/s' ,'mapQERR.txt','timeQERR.txt')
    outvarid(31) = idtype('cobc','computed outflow before correction',i_mean , 0,'m3/s','m3/s' ,'mapCOBC.txt','timeCOBC.txt')
    outvarid(32) = idtype('wtmp','computed water temp' ,i_mean , 0,'deg'  ,'degree Celsius','mapWTMP.txt' ,'timeWTMP.txt')
    outvarid(33) = idtype('werr','waterstage error',i_mean , 0,'m','m' ,'mapWERR.txt','timeWERR.txt')
    outvarid(34) = idtype('cwbc','computed wcom before correction',i_mean , 0,'m','m' ,'mapCWBC.txt','timeCWBC.txt')
    outvarid(35) = idtype('reOC','recorded TOC outflow',i_wmean, 54,'mg/L','mg org-C/L'    ,'mapREOC.txt','timeREOC.txt')
    outvarid(36) = idtype('csIN','computed IN soil'    ,i_wmean, 10,'ug/L','ug INORG-N/L'  ,'mapCSIN.txt','timeCSIN.txt')
    outvarid(37) = idtype('sfst','computed soil frost' ,i_mean , 0, 'cm'  ,'cm'            ,'mapSFST.txt','timeSFST.txt')
    outvarid(38) = idtype('stmp','computed soil temp'  ,i_mean , 0, 'deg' ,'degree Celcius','mapSTMP.txt','timeSTMP.txt')
    outvarid(39) = idtype('sdep','computed snow depth' ,i_mean , 0, 'cm'  ,'cm'            ,'mapSDEP.txt','timeSDEP.txt')
    outvarid(40) = idtype('epot','potential evap'      ,i_sum  , 0, 'mm'  ,'mm'            ,'mapEPOT.txt','timeEPOT.txt') !This should be improved to take into account intereption evap!
    outvarid(o_reepot) = idtype('repo','recorded pot. evap',i_sum ,0,'mm' ,'mm'   ,'mapREPO.txt','timeREPO.txt')
    outvarid(42) = idtype('eobs','recorded evaporation'       ,i_sum , 0,'mm'   ,'mm'    ,'mapEOBS.txt' ,'timeEOBS.txt')
    outvarid(43) = idtype('cprc','corr precipitation'     ,i_sum , 0,'mm'   ,'mm'            ,'mapCPRC.txt' ,'timeCPRC.txt')
    outvarid(44) = idtype('coOC','computed TOC runoff',i_wmean, 1,'mg/L'    ,'mg org-C/L'     ,'mapCOOC.txt','timeCOOC.txt')
    outvarid(45) = idtype('csOC','computed TOC soil'  ,i_wmean, 10,'mg/L'    ,'mg org-C/L'    ,'mapCSOC.txt','timeCSOC.txt')
    outvarid(46) = idtype('ccOC','computed TOC outflow',i_wmean, 53,'mg/L','mg org-C/L' ,'mapCCOC.txt','timeCCOC.txt')
    outvarid(47) = idtype('phC1','pool humusC soil1',i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC1.txt','timePHC1.txt')
    outvarid(48) = idtype('pfC1','pool fastC soil1',i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC1.txt','timePFC1.txt')
    outvarid(49) = idtype('phC2','pool humusC soil2',i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC2.txt','timePHC2.txt')
    outvarid(50) = idtype('pfC2','pool fastC soil2',i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC2.txt','timePFC2.txt')
    outvarid(51) = idtype('wcom','olake water stage',i_mean , 0,'m','meter' ,'mapWCOM.txt','timeWCOM.txt')
    outvarid(o_rewstr) = idtype('wstr','rec olake water st',i_mean , 0,'m','meter' ,'mapWSTR.txt','timeWSTR.txt')
    outvarid(53) = idtype('cout','comp outflow olake',i_mean , 0,'m3/s','m3/s' ,'mapCOUT.txt','timeCOUT.txt')
    outvarid(54) = idtype('rout','rec outflow olake',i_mean , 0,'m3/s','m3/s' ,'mapROUT.txt','timeROUT.txt')
    outvarid(55) = idtype('ccIN','comp conc IN olake',i_wmean, 53,'ug/L','ug InorgN-N/L' ,'mapCCIN.txt','timeCCIN.txt')
    outvarid(56) = idtype('ccON','comp conc ON olake',i_wmean, 53,'ug/L','ug OrgN-N/L' ,'mapCCON.txt','timeCCON.txt')
    outvarid(57) = idtype('ccSP','comp conc SP olake',i_wmean, 53,'ug/L','ug SRP-P/L' ,'mapCCSP.txt','timeCCSP.txt')
    outvarid(58) = idtype('ccPP','comp conc PP olake',i_wmean, 53,'ug/L','ug PartP-P/L' ,'mapCCPP.txt','timeCCPP.txt')
    outvarid(59) = idtype('rsnw','recorded snow depth',i_mean , 0,'cm','cm' ,'mapRSNW.txt','timeRSNW.txt')
    outvarid(60) = idtype('resf','recorded soil frost',i_mean , 0,'cm','cm' ,'mapRESF.txt','timeRESF.txt')
    outvarid(61) = idtype('regw','recorded grw level',i_mean , 0,'m','meter' ,'mapREGW.txt','timeREGW.txt')
    outvarid(63) = idtype('ccT1','comp conc T1 olake',i_wmean, 53,'?','?' ,'mapCCT1.txt','timeCCT1.txt')
    outvarid(64) = idtype('ccT2','comp conc T2 olake',i_wmean, 53,'deg'  ,'degree Celsius','mapCCT2.txt','timeCCT2.txt')
    outvarid(65) = idtype('ccT3','comp conc T3 olake',i_wmean, 53,'?','?' ,'mapCCT3.txt','timeCCT3.txt')
    outvarid(66) = idtype('coTN','computed TN runoff',i_wmean, 1,'ug/L' ,'ug Tot-N/L'    ,'mapCOTN.txt','timeCOTN.txt')
    outvarid(67) = idtype('coTP','computed TP runoff',i_wmean, 1,'ug/L' ,'ug Tot-P/L'    ,'mapCOTP.txt','timeCOTP.txt')
    outvarid(68) = idtype('pfN1','pool fastN soil1' ,i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPFN1.txt','timePFN1.txt')
    outvarid(69) = idtype('pfN2','pool fastN soil2' ,i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPFN2.txt','timePFN2.txt')
    outvarid(70) = idtype('pfN3','pool fastN soil3' ,i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPFN3.txt','timePFN3.txt')
    outvarid(71) = idtype('phN1','pool humusN soil1',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHN1.txt','timePHN1.txt')
    outvarid(72) = idtype('phN2','pool humusN soil2',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHN2.txt','timePHN2.txt')
    outvarid(73) = idtype('phN3','pool humusN soil3',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHN3.txt','timePHN3.txt')
    outvarid(74) = idtype('pIN1','pool InorgN soil1',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPIN1.txt','timePIN1.txt')
    outvarid(75) = idtype('pIN2','pool InorgN soil2',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPIN2.txt','timePIN2.txt')
    outvarid(76) = idtype('pIN3','pool InorgN soil3',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPIN3.txt','timePIN3.txt')
    outvarid(77) = idtype('ccTN','computed TN olake',i_wmean, 53,'ug/L' ,'ug Tot-N/L'    ,'mapCCTN.txt','timeCCTN.txt')
    outvarid(78) = idtype('ccTP','computed TP olake',i_wmean, 53,'ug/L' ,'ug Tot-P/L'    ,'mapCCTP.txt','timeCCTP.txt')
    outvarid(79) = idtype('cro1','computed runoff 1',i_sum  , 0,'mm'     ,'mm'           ,'mapCRO1.txt','timeCRO1.txt')
    outvarid(80) = idtype('cro2','computed runoff 2',i_sum  , 0,'mm'     ,'mm'           ,'mapCRO2.txt','timeCRO2.txt')
    outvarid(81) = idtype('cgwl','computed groundwater loss',i_mean , 0,'m3/s' ,'m3/s'  ,'mapCGWL.txt' ,'timeCGWL.txt')
    outvarid(82) = idtype('crod','computed rf drain',i_sum  , 0,'mm' ,'mm'           ,'mapCROD.txt' ,'timeCROD.txt')
    outvarid(83) = idtype('cros','computed surface rf'    ,i_sum , 0,'mm' ,'mm'      ,'mapCROS.txt' ,'timeCROS.txt')
    outvarid(84) = idtype('deni','denitrifikation'  ,i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapDENI.txt' ,'timeDENI.txt')
    outvarid(85) = idtype('crut','crop uptake'      ,i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapCRUT.txt' ,'timeCRUT.txt')
    outvarid(86) = idtype('faIN','fast to inorganic',i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapFAIN.txt' ,'timeFAIN.txt')
    outvarid(87) = idtype('atmd','atmospheric deposition N'      ,i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapATMD.txt' ,'timeATMD.txt')
    outvarid(88) = idtype('ppP1','pool partP soil1',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPPP1.txt','timePPP1.txt')   
    outvarid(89) = idtype('ppP2','pool partP soil2',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPPP2.txt','timePPP2.txt')
    outvarid(90) = idtype('ppP3','pool partP soil3',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPPP3.txt','timePPP3.txt')
    outvarid(91) = idtype('pSP1','pool SRP soil1',  i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPSP1.txt','timePSP1.txt')
    outvarid(92) = idtype('pSP2','pool SRP soil2',  i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPSP2.txt','timePSP2.txt')
    outvarid(93) = idtype('pSP3','pool SRP soil3',  i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPSP3.txt','timePSP3.txt')
    outvarid(94) = idtype('pfP1','pool fastP soil1',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPFP1.txt','timePFP1.txt')
    outvarid(95) = idtype('pfP2','pool fastP soil2',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPFP2.txt','timePFP2.txt')
    outvarid(96) = idtype('pfP3','pool fastP soil3',i_mean , 0,'kg/km2' ,'kg P/km2'    ,'mapPFP3.txt','timePFP3.txt')
    outvarid(o_cprecIN) = idtype('cpIN','recorded IN precip',i_wmean, 3,'ug/L','ug N/L','mapCPIN.txt'  ,'timeCPIN.txt')
    outvarid(o_cprecSP) = idtype('cpSP','recorded SP precip',i_wmean, 3,'ug/L','ug P/L','mapCPSP.txt'  ,'timeCPSP.txt')
    outvarid(o_reswe) = idtype('rswe','rec snow water eq.',i_mean , 0,'mm','mm' ,'mapRSWE.txt','timeRSWE.txt')
    outvarid(100) = idtype('acdf','accumulated volume error',i_sum , 0,'mm','mm' ,'mapACDF.txt','timeACDF.txt')
    outvarid(101) = idtype('cloc','comp local flow',   i_mean , 0, 'm3/s','m3/s' ,'mapCLOC.txt','timeCLOC.txt')
    outvarid(102) = idtype('clIN','comp local conc IN',i_wmean, 101,'ug/L','ug InorgN-N/L' ,'mapCLIN.txt','timeCLIN.txt')
    outvarid(103) = idtype('clON','comp local conc ON',i_wmean, 101,'ug/L','ug OrgN-N/L' ,'mapCLON.txt','timeCLON.txt')
    outvarid(104) = idtype('clSP','comp local conc SP',i_wmean, 101,'ug/L','ug SRP-P/L' ,'mapCLSP.txt','timeCLSP.txt')
    outvarid(105) = idtype('clPP','comp local conc PP',i_wmean, 101,'ug/L','ug PartP-P/L' ,'mapCLPP.txt','timeCLPP.txt')
    outvarid(106) = idtype('clTN','comp local conc TN',i_wmean, 101,'ug/L' ,'ug Tot-N/L'    ,'mapCLTN.txt','timeCLTN.txt')
    outvarid(107) = idtype('clTP','comp local conc TP',i_wmean, 101,'ug/L' ,'ug Tot-P/L'    ,'mapCLTP.txt','timeCLTP.txt')
    outvarid(108) = idtype('pfC3','pool fastC soil3',  i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC3.txt','timePFC3.txt')
    outvarid(109) = idtype('phC3','pool humusC soil3', i_mean , 0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC3.txt','timePHC3.txt')
    outvarid(110) = idtype('totN','comp load TN olake',i_sum , 0,'kg' ,'kg totN'    ,'mapTOTN.txt','timeTOTN.txt')
    outvarid(111) = idtype('totP','comp load TP olake',i_sum , 0,'kg' ,'kg totP'    ,'mapTOTP.txt','timeTOTP.txt')
    outvarid(112) = idtype('cinf','comp inflow to olake',i_mean , 0,'m3/s' ,'m3/s'    ,'mapCINF.txt','timeCINF.txt')
    outvarid(113) = idtype('rinf','rec inflow to olake',i_mean , 0,'m3/s' ,'m3/s'    ,'mapRINF.txt','timeRINF.txt')
    outvarid(114) = idtype('clrv','comp local river volume',i_mean , 0,'m3' ,'m3'    ,'mapCLRV.txt','timeCLRV.txt')
    outvarid(115) = idtype('cmrv','comp main river volume', i_mean , 0,'m3' ,'m3'    ,'mapCMRV.txt','timeCMRV.txt')
    outvarid(116) = idtype('atmp','atmospheric deposition TP',i_sum , 0,'kg/km2' ,'kg P/km2'  ,'mapATMP.txt' ,'timeATMP.txt')
    outvarid(117) = idtype('glcv','comp glacier ice volume', i_mean , 0,'km3' ,'km3'    ,'mapGLCV.txt','timeGLCV.txt')
    outvarid(118) = idtype('glca','comp glacier area', i_mean , 0,'km2' ,'km2'    ,'mapGLCA.txt','timeGLCA.txt')
    outvarid(119) = idtype('rtoN','rec load TN olake', i_sum , 0,'kg' ,'kg totN'    ,'mapRTON.txt','timeRTON.txt')
    outvarid(120) = idtype('rtoP','rec load TN olake', i_sum , 0,'kg' ,'kg totP'    ,'mapRTOP.txt','timeRTOP.txt')
    outvarid(121) = idtype('ctmp','corrected air temperature'   ,i_mean , 0,'deg'  ,'degree Celsius','mapCTMP.txt' ,'timeCTMP.txt') 
    outvarid(122) = idtype('irel','irrigation evap losses'   ,i_sum , 0,'m3'  ,'m3','mapIREL.txt' ,'timeIREL.txt') 
    outvarid(123) = idtype('coum','comp outflow main',i_mean , 0,'m3/s','m3/s' ,'mapCOUM.txt','timeCOUM.txt')
    outvarid(124) = idtype('irld','abstr local dam f irr'    ,i_sum , 0,'m3'  ,'m3','mapIRLD.txt' ,'timeIRLD.txt')
    outvarid(125) = idtype('irlr','abstr local river f irr'  ,i_sum , 0,'m3'  ,'m3','mapIRLR.txt' ,'timeIRLR.txt')
    outvarid(126) = idtype('irra','applied irrigation water' ,i_sum , 0,'m3'  ,'m3','mapIRRA.txt' ,'timeIRRA.txt')
    outvarid(127) = idtype('irrg','gwater abstracted f irr'  ,i_sum , 0,'m3'  ,'m3','mapIRRG.txt' ,'timeIRRG.txt') 
    outvarid(128) = idtype('irrs','irr abstr f other subbasins'   ,i_sum , 0,'m3'  ,'m3','mapIRRS.txt' ,'timeIRRS.txt')
    outvarid(129) = idtype('upsd','aver upstream soildf',i_mean, 0,'mm','mm','mapUPSD.txt','timeUPSD.txt')
    outvarid(130) = idtype('coub','comp outflow branch',i_mean , 0,'m3/s','m3/s' ,'mapCOUB.txt','timeCOUB.txt')
    outvarid(131) = idtype('smdf','soil moisture deficit'    ,i_mean , 0,'mm'   ,'mm'            ,'mapSMDF.txt' ,'timeSMDF.txt')
    outvarid(132) = idtype('phP1','pool humusP soil1',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHP1.txt','timePHP1.txt')
    outvarid(133) = idtype('phP2','pool humusP soil2',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHP2.txt','timePHP2.txt')
    outvarid(134) = idtype('phP3','pool humusP soil3',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPHP3.txt','timePHP3.txt')
    outvarid(135) = idtype('totC','comp load OC olake',i_sum , 0,'kg' ,'kg OC'    ,'mapTOTC.txt','timeTOTC.txt')
    outvarid(136) = idtype('pON1','pool orgN soil1',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPON1.txt','timePON1.txt')
    outvarid(137) = idtype('pON2','pool orgN soil2',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPON2.txt','timePON2.txt')
    outvarid(138) = idtype('pON3','pool orgN soil3',i_mean , 0,'kg/km2' ,'kg N/km2'    ,'mapPON3.txt','timePON3.txt')
    outvarid(139) = idtype('stm1','computed soillayer 1 temp'  ,i_mean , 0, 'deg' ,'degree Celcius','mapSTM1.txt','timeSTM1.txt')
    outvarid(140) = idtype('stm2','computed soillayer 2 temp'  ,i_mean , 0, 'deg' ,'degree Celcius','mapSTM2.txt','timeSTM2.txt')
    outvarid(141) = idtype('stm3','computed soillayer 3 temp'  ,i_mean , 0, 'deg' ,'degree Celcius','mapSTM3.txt','timeSTM3.txt')
    outvarid(142) = idtype('cpRF','precipitation as rain'     ,i_sum , 0,'mm'   ,'mm'             ,'mapCPRF.txt' ,'timeCPRF.txt')
    outvarid(143) = idtype('cpSF','precipitation as snow'     ,i_sum , 0,'mm'   ,'mm'             ,'mapCPSF.txt' ,'timeCPSF.txt')
    outvarid(144) = idtype('colv','volume in olake',i_mean , 0,'M(m3)' ,'M(m3)'    ,'mapCOLV.txt','timeCOLV.txt') !Computed Lake Volume  (basins summed to outlet if any) 
    outvarid(145) = idtype('cilv','volume in olake',i_mean , 0,'M(m3)' ,'M(m3)'    ,'mapCILV.txt','timeCILV.txt') !Computed ILake Volume 
    outvarid(146) = idtype('clbv','volume in olake',i_mean , 0,'M(m3)' ,'M(m3)'    ,'mapCLBV.txt','timeCLBV.txt') !Computed Olake Volume Computed (volumes for individual basins if any)
    outvarid(147) = idtype('cro3','computed runoff 3',i_sum  , 0,'mm'     ,'mm'           ,'mapCRO3.txt','timeCRO3.txt')
    !------------------------------------------------------------------------------------
    ! New Lake and River Ice, Snow, and Water Temperature model variables David(20130426)
    !------------------------------------------------------------------------------------
    !Lake ice and snow depth variables
    outvarid(148) = idtype('coli','comp olake ice depth'      ,i_mean , 0,'cm'  ,'cm','mapCOLI.txt' ,'timeCOLI.txt')
    outvarid(149) = idtype('cili','comp ilake ice depth'      ,i_mean , 0,'cm'  ,'cm','mapCILI.txt' ,'timeCILI.txt')
    outvarid(150) = idtype('colb','comp olake blackice depth' ,i_mean , 0,'cm'  ,'cm','mapCOLB.txt' ,'timeCOLB.txt')
    outvarid(151) = idtype('cilb','comp ilake blackice depth' ,i_mean , 0,'cm'  ,'cm','mapCILB.txt' ,'timeCILB.txt')
    outvarid(152) = idtype('cols','comp olake snow depth'     ,i_mean , 0,'cm'  ,'cm','mapCOLS.txt' ,'timeCOLS.txt')
    outvarid(153) = idtype('cils','comp ilake snow depth'     ,i_mean , 0,'cm'  ,'cm','mapCILS.txt' ,'timeCILS.txt')
    outvarid(154) = idtype('roli','rec. olake ice depth'      ,i_mean , 0,'cm'  ,'cm','mapROLI.txt' ,'timeROLI.txt')
    outvarid(155) = idtype('rili','rec. ilake ice depth'      ,i_mean , 0,'cm'  ,'cm','mapRILI.txt' ,'timeRILI.txt')
    outvarid(156) = idtype('rolb','rec. olake blackice depth' ,i_mean , 0,'cm'  ,'cm','mapROLB.txt' ,'timeROLB.txt')
    outvarid(157) = idtype('rilb','rec. ilake blackice depth' ,i_mean , 0,'cm'  ,'cm','mapRILB.txt' ,'timeRILB.txt')
    outvarid(158) = idtype('rols','rec. olake snow depth'     ,i_mean , 0,'cm'  ,'cm','mapROLS.txt' ,'timeROLS.txt')
    outvarid(159) = idtype('rils','rec. ilake snow depth'     ,i_mean , 0,'cm'  ,'cm','mapRILS.txt' ,'timeRILS.txt')
    
    !River ice and snow depth variables
    outvarid(160) = idtype('cmri','comp main river ice depth'      ,i_mean , 0,'cm'  ,'cm','mapCMRI.txt' ,'timeCMRI.txt')
    outvarid(161) = idtype('clri','comp local river ice depth'     ,i_mean , 0,'cm'  ,'cm','mapCLRI.txt' ,'timeCLRI.txt')
    outvarid(162) = idtype('cmrb','comp main river blackice depth' ,i_mean , 0,'cm'  ,'cm','mapCMRB.txt' ,'timeCMRB.txt')
    outvarid(163) = idtype('clrb','comp local river blackice depth',i_mean , 0,'cm'  ,'cm','mapCLRB.txt' ,'timeCLRB.txt')
    outvarid(164) = idtype('cmrs','comp main river snow depth'     ,i_mean , 0,'cm'  ,'cm','mapCMRS.txt' ,'timeCMRS.txt')
    outvarid(165) = idtype('clrs','comp local river snow depth'    ,i_mean , 0,'cm'  ,'cm','mapCLRS.txt' ,'timeCLRS.txt')
    outvarid(166) = idtype('rmri','rec. main river ice depth'      ,i_mean , 0,'cm'  ,'cm','mapRMRI.txt' ,'timeRMRI.txt')
    outvarid(167) = idtype('rlri','rec. local river ice depth'     ,i_mean , 0,'cm'  ,'cm','mapRLRI.txt' ,'timeRLRI.txt')
    outvarid(168) = idtype('rmrb','rec. main river blackice depth' ,i_mean , 0,'cm'  ,'cm','mapRMRB.txt' ,'timeRMRB.txt')
    outvarid(169) = idtype('rlrb','rec. local river blackice depth',i_mean , 0,'cm'  ,'cm','mapRLRB.txt' ,'timeRLRB.txt')
    outvarid(170) = idtype('rmrs','rec. main river snow depth'     ,i_mean , 0,'cm'  ,'cm','mapRMRS.txt' ,'timeRMRS.txt')
    outvarid(171) = idtype('rlrs','rec. local river snow depth'    ,i_mean , 0,'cm'  ,'cm','mapRLRS.txt' ,'timeRLRS.txt')
    
    ! Lake and river temperature variables (surface, upper layer, lower layer, mean water temp)
    outvarid(172) = idtype('olst','comp olake surface temp'         ,i_mean , 0,'deg'  ,'degree Celsius','mapOLST.txt' ,'timeOLST.txt')
    outvarid(173) = idtype('olut','comp olake upper temp'           ,i_mean , 0,'deg'  ,'degree Celsius','mapOLUT.txt' ,'timeOLUT.txt')
    outvarid(174) = idtype('ollt','comp olake lower temp'           ,i_mean , 0,'deg'  ,'degree Celsius','mapOLLT.txt' ,'timeOLLT.txt')
    outvarid(175) = idtype('olwt','comp olake mean  temp'           ,i_mean , 0,'deg'  ,'degree Celsius','mapOLWT.txt' ,'timeOLWT.txt')
    outvarid(176) = idtype('ilst','comp ilake surface temp'         ,i_mean , 0,'deg'  ,'degree Celsius','mapILST.txt' ,'timeILST.txt')
    outvarid(177) = idtype('ilwt','comp ilake mean  temp'           ,i_mean , 0,'deg'  ,'degree Celsius','mapILWT.txt' ,'timeILWT.txt')
    outvarid(178) = idtype('lrst','comp local river surface temp'   ,i_mean , 0,'deg'  ,'degree Celsius','mapLRST.txt' ,'timeLRST.txt')
    outvarid(179) = idtype('lrwt','comp local river mean  temp'     ,i_mean , 0,'deg'  ,'degree Celsius','mapLRWT.txt' ,'timeLRWT.txt')
    outvarid(180) = idtype('mrst','comp main  river surface temp'   ,i_mean , 0,'deg'  ,'degree Celsius','mapMRST.txt' ,'timeMRST.txt')
    outvarid(181) = idtype('mrwt','comp main  river mean  temp'     ,i_mean , 0,'deg'  ,'degree Celsius','mapMRWT.txt' ,'timeMRWT.txt')
    
    !MODIS Water surface temperatures
    outvarid(182) = idtype('rolt','rec. olake surface temp'         ,i_mean , 0,'deg'  ,'degree Celsius','mapROLT.txt' ,'timeROLT.txt')
    outvarid(183) = idtype('rilt','rec. ilake surface temp'         ,i_mean , 0,'deg'  ,'degree Celsius','mapRILT.txt' ,'timeRILT.txt')
    outvarid(184) = idtype('rmrt','rec. main river surface temp'    ,i_mean , 0,'deg'  ,'degree Celsius','mapRMRT.txt' ,'timeRMRT.txt')
    
    !Additional extra outputs from old water temperature model
    outvarid(185) = idtype('mrto','OLD comp main river mean temp'   ,i_mean , 0,'deg'  ,'degree Celsius','mapMRTO.txt' ,'timeMRTO.txt')
    outvarid(186) = idtype('lrto','OLD comp local river mean temp'  ,i_mean , 0,'deg'  ,'degree Celsius','mapLRTO.txt' ,'timeLRTO.txt')
    outvarid(187) = idtype('ilto','OLD comp ilake mean temp'        ,i_mean , 0,'deg'  ,'degree Celsius','mapILTO.txt' ,'timeILTO.txt')
    outvarid(188) = idtype('olto','OLD comp olake mean temp'        ,i_mean , 0,'deg'  ,'degree Celsius','mapOLTO.txt' ,'timeOLTO.txt')
    
    !Snowcover fraction, computed and recorded (for land fractions only)
    outvarid(189) = idtype('cfsc','comp frac snowcover area'       ,i_mean , 0,'-'  ,'fraction','mapCFSC.txt' ,'timeCFSC.txt')
    outvarid(190) = idtype('rfsc','rec frac snowcover area'        ,i_mean , 0,'-'  ,'fraction','mapRFSC.txt' ,'timeRFSC.txt')
    outvarid(191) = idtype('smax','comp snowmax in winter'         ,i_mean , 0,'mm' ,'mm','mapSMAX.txt' ,'timeSMAX.txt')
    outvarid(192) = idtype('rfse','rec frac snowcover area error'  ,i_mean , 0,'-'  ,'fraction','mapRFSE.txt' ,'timeRFSE.txt')
    outvarid(193) = idtype('rfsm','rec frac snowcover multi'       ,i_mean , 0,'-'  ,'fraction','mapRFSM.txt' ,'timeRFSM.txt')
    outvarid(194) = idtype('rfme','rec frac snowcover multi error' ,i_mean , 0,'-'  ,'fraction','mapRFME.txt' ,'timeRFME.txt')
    outvarid(195) = idtype('levp','land evaporation'               ,i_sum  , 0,'mm' ,'mm'      ,'mapLEVP.txt' ,'timeLEVP.txt')
    outvarid(196) = idtype('som2','soil water upper 2 l',i_mean , 0,'mm' ,'mm'            ,'mapSOM2.txt' ,'timeSOM2.txt')
    outvarid(197) = idtype('sml1','soil moisture lay 1' ,i_mean , 0,'mm' ,'mm'            ,'mapSML1.txt' ,'timeSML1.txt')
    outvarid(198) = idtype('sml2','soil moisture lay 2' ,i_mean , 0,'mm' ,'mm'            ,'mapSML2.txt' ,'timeSML2.txt')
    outvarid(199) = idtype('sml3','soil moisture lay 3' ,i_mean , 0,'mm' ,'mm'            ,'mapSML3.txt' ,'timeSML3.txt')
    outvarid(200) = idtype('stsw','standing soil water' ,i_mean , 0,'mm' ,'mm'            ,'mapSTSW.txt' ,'timeSTSW.txt')
    outvarid(201) = idtype('smrz','soil moisture root z',i_mean , 0,'mm' ,'mm'            ,'mapSMRZ.txt' ,'timeSMRZ.txt')
    outvarid(202) = idtype('sm13','soil moisture sl 1-3',i_mean , 0,'mm' ,'mm'            ,'mapSM13.txt' ,'timeSM13.txt')
    outvarid(203) = idtype('clCO','comp local conc OC',i_wmean, 101,'mg/L' ,'mg org-C/L'  ,'mapCLCO.txt','timeCLCO.txt')    !clOC occupied (cloc)
    outvarid(204) = idtype('lrdp','local river depth'    ,i_mean , 0,'m'  ,'m','mapLRDP.txt' ,'timeLRDP.txt')
    outvarid(205) = idtype('mrdp','main river depth'     ,i_mean , 0,'m'  ,'m','mapMRDP.txt' ,'timeMRDP.txt')
    
    outvarid(206) = idtype('icpe','interception evap'    ,i_sum , 0,'mm'  ,'mm','mapICPE.txt' ,'timeICPE.txt') !better call it evaporation

    !Load recorded concentrations and computed flow
    outvarid(207) = idtype('rlIN','load rec IN com flow',i_sum, 0,'kg' ,'kg Inorg-N'  ,'mapRLIN.txt','timeRLIN.txt')
    outvarid(208) = idtype('rlON','load rec ON com flow',i_sum, 0,'kg' ,'kg Org-N'    ,'mapRLON.txt','timeRLON.txt')
    outvarid(209) = idtype('rlSP','load rec SP com flow',i_sum, 0,'kg' ,'kg SRP-P'    ,'mapRLSP.txt','timeRLSP.txt')
    outvarid(210) = idtype('rlPP','load rec PP com flow',i_sum, 0,'kg' ,'kg PartP-P'  ,'mapRLPP.txt','timeRLPP.txt')
    outvarid(211) = idtype('rlTN','load rec TN com flow',i_sum, 0,'kg' ,'kg Tot-N'    ,'mapRLTN.txt','timeRLTN.txt')
    outvarid(212) = idtype('rlTP','load rec TP com flow',i_sum, 0,'kg' ,'kg Tot-P'    ,'mapRLTP.txt','timeRLTP.txt')
    outvarid(213) = idtype('rlOC','load rec OC com flow',i_sum, 0,'kg' ,'kg org-C'    ,'mapRLOC.txt','timeRLOC.txt')

    outvarid(214) = idtype('srff','sm frac of field cap'     ,i_mean , 0,'-'  ,'-','mapSRFF.txt' ,'timeSRFF.txt')
    outvarid(215) = idtype('smfd','sm frac of depth'         ,i_mean , 0,'-'  ,'-','mapSMFD.txt' ,'timeSMFD.txt')
    outvarid(216) = idtype('srfd','sm root frac of depth'    ,i_mean , 0,'-'  ,'-','mapSRFD.txt' ,'timeSRFD.txt')
    outvarid(217) = idtype('smfp','sm frac of pore v'        ,i_mean , 0,'-'  ,'-','mapSMFP.txt' ,'timeSMFP.txt')
    outvarid(218) = idtype('srfp','sm root frac of pore v'   ,i_mean , 0,'-'  ,'-','mapSRFP.txt' ,'timeSRFP.txt')
    outvarid(219) = idtype('upfp','aver upstream smfp',i_mean, 0,'mm','mm','mapUPFP.txt','timeUPFP.txt')
    outvarid(o_xobsm)   = idtype('xom0','recorded mean var 0',i_mean , 0,'?','?' ,'mapXOM0.txt','timeXOM0.txt')
    outvarid(o_xobsm+1) = idtype('xom1','recorded mean var 1',i_mean , 0,'?','?' ,'mapXOM1.txt','timeXOM1.txt')
    outvarid(o_xobsm+2) = idtype('xom2','recorded mean var 2',i_mean , 0,'?','?' ,'mapXOM2.txt','timeXOM2.txt')
    outvarid(o_xobsm+3) = idtype('xom3','recorded mean var 3',i_mean , 0,'?','?' ,'mapXOM3.txt','timeXOM3.txt')
    outvarid(o_xobsm+4) = idtype('xom4','recorded mean var 4',i_mean , 0,'?','?' ,'mapXOM4.txt','timeXOM4.txt')
    outvarid(o_xobsm+5) = idtype('xom5','recorded mean var 5',i_mean , 0,'?','?' ,'mapXOM5.txt','timeXOM5.txt')
    outvarid(o_xobsm+6) = idtype('xom6','recorded mean var 6',i_mean , 0,'?','?' ,'mapXOM6.txt','timeXOM6.txt')
    outvarid(o_xobsm+7) = idtype('xom7','recorded mean var 7',i_mean , 0,'?','?' ,'mapXOM7.txt','timeXOM7.txt')
    outvarid(o_xobsm+8) = idtype('xom8','recorded mean var 8',i_mean , 0,'?','?' ,'mapXOM8.txt','timeXOM8.txt')
    outvarid(o_xobsm+9) = idtype('xom9','recorded mean var 9',i_mean , 0,'?','?' ,'mapXOM9.txt','timeXOM9.txt')
    outvarid(o_xobss)   = idtype('xos0','recorded sum var 0',i_sum  , 0,'?','?' ,'mapXOS0.txt','timeXOS0.txt')
    outvarid(o_xobss+1) = idtype('xos1','recorded sum var 1',i_sum  , 0,'?','?' ,'mapXOS1.txt','timeXOS1.txt')
    outvarid(o_xobss+2) = idtype('xos2','recorded sum var 2',i_sum  , 0,'?','?' ,'mapXOS2.txt','timeXOS2.txt')
    outvarid(o_xobss+3) = idtype('xos3','recorded sum var 3',i_sum  , 0,'?','?' ,'mapXOS3.txt','timeXOS3.txt')
    outvarid(o_xobss+4) = idtype('xos4','recorded sum var 4',i_sum  , 0,'?','?' ,'mapXOS4.txt','timeXOS4.txt')
    outvarid(o_xobss+5) = idtype('xos5','recorded sum var 5',i_sum  , 0,'?','?' ,'mapXOS5.txt','timeXOS5.txt')
    outvarid(o_xobss+6) = idtype('xos6','recorded sum var 6',i_sum  , 0,'?','?' ,'mapXOS6.txt','timeXOS6.txt')
    outvarid(o_xobss+7) = idtype('xos7','recorded sum var 7',i_sum  , 0,'?','?' ,'mapXOS7.txt','timeXOS7.txt')
    outvarid(o_xobss+8) = idtype('xos8','recorded sum var 8',i_sum  , 0,'?','?' ,'mapXOS8.txt','timeXOS8.txt')
    outvarid(o_xobss+9) = idtype('xos9','recorded sum var 9',i_sum  , 0,'?','?' ,'mapXOS9.txt','timeXOS9.txt')

    outvarid(240) = idtype('aqin','aquifer recharge'   ,i_sum , 0,'m3'  ,'m3','mapAQIN.txt' ,'timeAQIN.txt')
    outvarid(241) = idtype('aqut','aquifer return flow',i_sum , 0,'m3'  ,'m3','mapAQUT.txt' ,'timeAQUT.txt')
    outvarid(242) = idtype('aqwl','aquifer water depth',i_mean, 0,'m'   ,'m','mapAQWL.txt' ,'timeAQWL.txt')

    outvarid(243) = idtype('uppr','aver upstream prec',i_sum, 0,'mm','mm','mapUPPR.txt','timeUPPR.txt')
    outvarid(244) = idtype('upev','aver upstream evap',i_sum, 0,'mm','mm','mapUPEV.txt','timeUPEV.txt')  !Total upstream incl. interception losses
    outvarid(245) = idtype('upsn','aver upstream snoww',i_mean, 0,'mm','mm','mapUPSN.txt','timeUPSN.txt')
    outvarid(246) = idtype('upso','aver upstream soilw',i_mean, 0,'mm','mm','mapUPSO.txt','timeUPSO.txt')
    outvarid(247) = idtype('upro','specific discharge',i_sum, 0,'mm','mm','mapUPRO.txt','timeUPRO.txt')

    outvarid(248) = idtype('coic','comp olake ice cover',i_mean, 0,'-','-','mapCOIC.txt','timeCOIC.txt')
    outvarid(249) = idtype('ciic','comp ilake ice cover',i_mean, 0,'-','-','mapCIIC.txt','timeCIIC.txt')
    outvarid(250) = idtype('cmic','comp m river ice cov',i_mean, 0,'-','-','mapCMIC.txt','timeCMIC.txt')
    outvarid(251) = idtype('clic','comp l river ice cov',i_mean, 0,'-','-','mapCLIC.txt','timeCLIC.txt')

    !Glacier outputs and corresponding WGMS observations, using the following 4-letter code structure:
    !  xGnn, where x  = r or c (recorded or computed) and nn = 2-letter code identifying the glacier variable
    outvarid(252) = idtype('cgmb','comp. glacier mass balance',i_mean, 0,'mm','mm','mapCGMB.txt','timeCGMB.txt')
    outvarid(253) = idtype('rgmb','rec. glacier mass balance',i_mean, 0,'mm','mm','mapRGMB.txt','timeRGMB.txt')
    outvarid(254) = idtype('cgma','area used in com. mass balance',i_mean, 0,'km2','km2','mapCGMA.txt','timeCGMA.txt')
    outvarid(255) = idtype('rgma','area used in rec. mass balance',i_mean, 0,'km2','km2','mapRGMA.txt','timeRGMA.txt')
    outvarid(256) = idtype('rgmp','rec. mass balance period',i_mean, 0,'days','days','mapRGMP.txt','timeRGMP.txt')
    
    !Glacier Xobs variables not yet included in the code:
    ! ---------------------------------------------------
    !xGEL Equilibrium Line Altitude  ELA                	(m)
    !xGAA Accumulation Area Ratio of total area         	(%)
    !xGGA Area Survey year (GTC, geodetic thicness changes)	(km2)
    !xGAC Area Change (GTC)                                 (1000 m2)
    !xGTC Thickness Change (GTC)                            (mm)
    !xGVC Volume Change (GTC)                               (1000 m3)
    !xGSP GTC Survey period, days since Reference survey    (days)

    !Sovjet Snow Coarse observations in forest and open areas, and corresponding model variables:
    outvarid(257) = idtype('S105','fsusnow fsc surr. open',i_mean, 0,'0-10','0-10','mapS105.txt','timeS105.txt')
    outvarid(258) = idtype('S106','fsusnow fsc course open',i_mean, 0,'0-10','0-10','mapS106.txt','timeS106.txt')
    outvarid(259) = idtype('S108','fsusnow mean depth open',i_mean, 0,'cm','cm','mapS108.txt','timeS108.txt')
    outvarid(260) = idtype('S111','fsusnow mean density open',i_mean, 0,'g/cm3','g/cm3','mapS111.txt','timeS111.txt')
    outvarid(261) = idtype('S114','fsusnow snow water eq. open',i_mean, 0,'mm','mm','mapS114.txt','timeS114.txt')
    outvarid(262) = idtype('S205','fsusnow fsc surr. forest',i_mean, 0,'0-10','0-10','mapS205.txt','timeS205.txt')
    outvarid(263) = idtype('S206','fsusnow fsc course forest',i_mean, 0,'0-10','0-10','mapS206.txt','timeS206.txt')
    outvarid(264) = idtype('S208','fsusnow mean depth forest',i_mean, 0,'cm','cm','mapS208.txt','timeS208.txt')
    outvarid(265) = idtype('S211','fsusnow mean density forest',i_mean, 0,'g/cm3','g/cm3','mapS211.txt','timeS211.txt')
    outvarid(266) = idtype('S214','fsusnow snow water eq. forest ',i_mean, 0,'mm','mm','mapS214.txt','timeS214.txt')
    outvarid(267) = idtype('C106','comp. fsc course open',i_mean, 0,'0-10','0-10','mapC106.txt','timeC106.txt')
    outvarid(268) = idtype('C108','comp. mean depth open',i_mean, 0,'cm','cm','mapC108.txt','timeC108.txt')
    outvarid(269) = idtype('C111','comp. mean density open',i_mean, 0,'g/cm3','g/cm3','mapC111.txt','timeC111.txt')
    outvarid(270) = idtype('C114','comp. snow water eq. open',i_mean, 0,'mm','mm','mapC114.txt','timeC114.txt')
    outvarid(271) = idtype('C206','comp. fsc course forest',i_mean, 0,'0-10','0-10','mapC206.txt','timeC206.txt')
    outvarid(272) = idtype('C208','comp. mean depth forest',i_mean, 0,'cm','cm','mapC208.txt','timeC208.txt')
    outvarid(273) = idtype('C211','comp. mean density forest',i_mean, 0,'g/cm3','g/cm3','mapC211.txt','timeC211.txt')
    outvarid(274) = idtype('C214','comp. snow water eq. forest ',i_mean, 0,'mm','mm','mapC214.txt','timeC214.txt')

    !Additional upstream water balance components (snowfall and rainfall)
    outvarid(275) = idtype('upsf','aver upstream snowfall',i_sum, 0,'mm','mm','mapUPSF.txt','timeUPSF.txt')
    outvarid(276) = idtype('uprf','aver upstream rainfall',i_sum, 0,'mm','mm','mapUPRF.txt','timeUPRF.txt')
    outvarid(277) = idtype('uppe','aver upstream pot evap',i_sum, 0,'mm','mm','mapUPPE.txt','timeUPPE.txt')
  
    outvarid(278) = idtype('evsn','snow evaporation',i_sum , 0,'mm' ,'mm','mapEVSN.txt' ,'timeEVSN.txt')   !This is now excluding glacier evap!
    outvarid(279) = idtype('evpt','total evaporation',i_sum , 0,'mm' ,'mm','mapEVPT.txt' ,'timeEVPT.txt')  !This is now redundant, EVAP is already=Total
    outvarid(280) = idtype('clwc','cleaned wcom',i_mean , 0,'m','meter' ,'mapCLWC.txt','timeCLWC.txt')
    outvarid(281) = idtype('clws','cleaned wstr',i_mean , 0,'m','meter' ,'mapCLWS.txt','timeCLWS.txt')
    outvarid(282) = idtype('wcoa','waterstage adjusted', i_mean , 0,'m','meter' ,'mapWCOA.txt','timeWCOA.txt')
    outvarid(283) = idtype('wtm0','comp water temp >0' , i_mean , 0,'deg'  ,'degree Celsius','mapWTM0.txt' ,'timeWTM0.txt')
    outvarid(284) = idtype('psim','sim. precipitation', i_sum , 0,'mm','mm','mapPSIM.txt' ,'timePSIM.txt')

    !Soil load output variables; 1-12 soillayer 1-2, 13-24 soillayer 3, odd=gross load, even=net load, 25-36 soillayer3+tiledrain
    outvarid(285) = idtype('sl07','soill 1-2 grs ld SP',i_sum , 0,'kg' ,'kg','mapSL07.txt' ,'timeSL07.txt')
    outvarid(286) = idtype('sl08','soill 1-2 net ld SP',i_sum , 0,'kg' ,'kg','mapSL08.txt' ,'timeSL08.txt')
    outvarid(287) = idtype('sl09','soill 1-2 grs ld PP',i_sum , 0,'kg' ,'kg','mapSL09.txt' ,'timeSL09.txt')
    outvarid(288) = idtype('sl10','soill 1-2 net ld PP',i_sum , 0,'kg' ,'kg','mapSL10.txt' ,'timeSL10.txt')
    outvarid(289) = idtype('sl11','soill 1-2 grs ld TP',i_sum , 0,'kg' ,'kg','mapSL11.txt' ,'timeSL11.txt')
    outvarid(290) = idtype('sl12','soill 1-2 net ld TP',i_sum , 0,'kg' ,'kg','mapSL12.txt' ,'timeSL12.txt')
    outvarid(291) = idtype('sl13','soill 3 grs ld IN',i_sum , 0,'kg' ,'kg','mapSL13.txt' ,'timeSL13.txt')
    outvarid(292) = idtype('sl14','soill 3 net ld IN',i_sum , 0,'kg' ,'kg','mapSL14.txt' ,'timeSL14.txt')
    outvarid(293) = idtype('sl15','soill 3 grs ld ON',i_sum , 0,'kg' ,'kg','mapSL15.txt' ,'timeSL15.txt')
    outvarid(294) = idtype('sl16','soill 3 net ld ON',i_sum , 0,'kg' ,'kg','mapSL16.txt' ,'timeSL16.txt')
    outvarid(295) = idtype('sl17','soill 3 grs ld TN',i_sum , 0,'kg' ,'kg','mapSL17.txt' ,'timeSL17.txt')
    outvarid(296) = idtype('sl18','soill 3 net ld TN',i_sum , 0,'kg' ,'kg','mapSL18.txt' ,'timeSL18.txt')
    outvarid(297) = idtype('sl19','soill 3 grs ld SP',i_sum , 0,'kg' ,'kg','mapSL19.txt' ,'timeSL19.txt')
    outvarid(298) = idtype('sl20','soill 3 net ld SP',i_sum , 0,'kg' ,'kg','mapSL20.txt' ,'timeSL20.txt')
    outvarid(299) = idtype('sl21','soill 3 grs ld PP',i_sum , 0,'kg' ,'kg','mapSL21.txt' ,'timeSL21.txt')
    outvarid(300) = idtype('sl22','soill 3 net ld PP',i_sum , 0,'kg' ,'kg','mapSL22.txt' ,'timeSL22.txt')
    outvarid(301) = idtype('sl23','soill 3 grs ld TP',i_sum , 0,'kg' ,'kg','mapSL23.txt' ,'timeSL23.txt')
    outvarid(302) = idtype('sl24','soill 3 net ld TP',i_sum , 0,'kg' ,'kg','mapSL24.txt' ,'timeSL24.txt')
    outvarid(303) = idtype('den3','denitrification sl3'  ,i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapDEN3.txt' ,'timeDEN3.txt')
    outvarid(304) = idtype('denz','denitrification sl12'  ,i_sum , 0,'kg/km2' ,'kg N/km2'  ,'mapDENZ.txt' ,'timeDENZ.txt')
    outvarid(305) = idtype('sl01','soill 1-2 grs ld IN',i_sum , 0,'kg' ,'kg','mapSL01.txt' ,'timeSL01.txt')
    outvarid(306) = idtype('sl02','soill 1-2 net ld IN',i_sum , 0,'kg' ,'kg','mapSL02.txt' ,'timeSL02.txt')
    outvarid(307) = idtype('sl03','soill 1-2 grs ld ON',i_sum , 0,'kg' ,'kg','mapSL03.txt' ,'timeSL03.txt')
    outvarid(308) = idtype('sl04','soill 1-2 net ld ON',i_sum , 0,'kg' ,'kg','mapSL04.txt' ,'timeSL04.txt')
    outvarid(309) = idtype('sl05','soill 1-2 grs ld TN',i_sum , 0,'kg' ,'kg','mapSL05.txt' ,'timeSL05.txt')
    outvarid(310) = idtype('sl06','soill 1-2 net ld TN',i_sum , 0,'kg' ,'kg','mapSL06.txt' ,'timeSL06.txt')
    outvarid(311) = idtype('sl25','soill 3+t grs ld IN',i_sum , 0,'kg' ,'kg','mapSL25.txt' ,'timeSL25.txt')
    outvarid(312) = idtype('sl26','soill 3+t net ld IN',i_sum , 0,'kg' ,'kg','mapSL26.txt' ,'timeSL26.txt')
    outvarid(313) = idtype('sl27','soill 3+t grs ld ON',i_sum , 0,'kg' ,'kg','mapSL27.txt' ,'timeSL27.txt')
    outvarid(314) = idtype('sl28','soill 3+t net ld ON',i_sum , 0,'kg' ,'kg','mapSL28.txt' ,'timeSL28.txt')
    outvarid(315) = idtype('sl29','soill 3+t grs ld TN',i_sum , 0,'kg' ,'kg','mapSL29.txt' ,'timeSL29.txt')
    outvarid(316) = idtype('sl30','soill 3+t net ld TN',i_sum , 0,'kg' ,'kg','mapSL30.txt' ,'timeSL30.txt')
    outvarid(317) = idtype('sl31','soill 3+t grs ld SP',i_sum , 0,'kg' ,'kg','mapSL31.txt' ,'timeSL31.txt')
    outvarid(318) = idtype('sl32','soill 3+t net ld SP',i_sum , 0,'kg' ,'kg','mapSL32.txt' ,'timeSL32.txt')
    outvarid(319) = idtype('sl33','soill 3+t grs ld PP',i_sum , 0,'kg' ,'kg','mapSL33.txt' ,'timeSL33.txt')
    outvarid(320) = idtype('sl34','soill 3+t net ld PP',i_sum , 0,'kg' ,'kg','mapSL34.txt' ,'timeSL34.txt')
    outvarid(321) = idtype('sl35','soill 3+t grs ld TP',i_sum , 0,'kg' ,'kg','mapSL35.txt' ,'timeSL35.txt')
    outvarid(322) = idtype('sl36','soill 3+t net ld TP',i_sum , 0,'kg' ,'kg','mapSL36.txt' ,'timeSL36.txt')
    outvarid(323) = idtype('clwa','cleaned wcoa',i_mean , 0,'m','meter' ,'mapCLWC.txt','timeCLWC.txt')

!Additional for CryoProcesses
    
    !Prec.correction (to separate interception losses from preciptation)     
    outvarid(max_outvar_trunk+1) = idtype('pcor','precipitation correction',i_sum , 0,'mm'  ,'mm','mapPCOR.txt' ,'timePCOR.txt') !323
    
    !Evapotraspiration components
    outvarid(max_outvar_trunk+2) = idtype('evgl','glacier evaporation',i_sum , 0,'mm' ,'mm','mapEVGL.txt' ,'timeEVGL.txt')
    outvarid(max_outvar_trunk+3) = idtype('evsv','soil+veg evaporation',i_sum , 0,'mm' ,'mm','mapEVSV.txt' ,'timeEVSV.txt')
 
    !Water balance components, normalized per unit area, local and upstream, storages(mm) and flows (mm/ts)
    !------------------------------------------------------------------------------------------------------
    !STORAGES
    !
    !--SNOW   
    !---local:     snow,  snowi,                                 outvar(14),         mm (per unit basin landarea)
    !              snwi,  snowi*landarea(i)/area(land+localw),   outvar(???),        mm (per unit basin area, including LOCAL rivers and lakes)
    !              snwo,  snowi*landarea(i)/basin(i)%area,       outvar(???),        mm (per unit basin area, including ALL lakes and rivers)
    !
    !---upstream:  upsn,  upstream_snow,                         outvar(245),        mm (per unit upstream area, including ALL lakes and rivers)
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+4) = idtype('snwi','snowpack pr.basin(ilake)'         ,i_mean , 0,'mm'   ,'mm water'      ,'mapSNWI.txt' ,'timeSNWI.txt')
    outvarid(max_outvar_trunk+5) = idtype('snwo','snowpack pr.basin(olake)'         ,i_mean , 0,'mm'   ,'mm water'      ,'mapSNWO.txt' ,'timeSNWO.txt')
    !------------------------------------------------------------------------------------------------------
    !--SOIL
    !---local:     soim, soili,                                  outvar(10),         mm (per unit basin landarea)
    !              soii, soili*landarea(i)/area(land+localw),    outvar(???),        mm (per unit basin area, inluding LOCAL lakes and rivers)
    !              soio, soili*landarea(i)/basin(i)%area,        outvar(???),        mm (per unit basin area, inluding ALL lakes and rivers)
    !
    !---upstream:  upso, upstream_soil,                          outvar(246),        mm (per unit upstream area, inluding lakes and rivers
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+6) = idtype('soii','soilwater pr.basin(ilake)'       ,i_mean , 0,'mm' ,'mm','mapSOII.txt' ,'timeSOII.txt')
    outvarid(max_outvar_trunk+7) = idtype('soio','soilwater pr.basin(olake)'       ,i_mean , 0,'mm' ,'mm','mapSOIO.txt' ,'timeSOIO.txt')
    !------------------------------------------------------------------------------------------------------
    !--STREAMS (LOCAL and ALL)
    !---local:     clrv,   riverstate%water...                   outvar(114),       m3 local river volume
    !              cmrv,                                         outvar(115),       m3 main river volume
    !              clrw,   clrv*1000/(landarea+(locriv+ilake))   outvar(???),       mm (per unit land and local surface water area in basin)
    !              criw,   (cmrv+clrv)*1000/basin(i)%area        outvar(???),       mm (per unit land and total surface water area in basin)
    !
    !---upstream:  upri,   upstream_river,                       outvar(???),       mm (per unit upstream area, including lakes and rivers
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+8) = idtype('clrw','local riverwater pr.basin(ilake)',i_mean, 0, 'mm', 'mm', 'mapCLRW.txt', 'timeCLRW.txt')
    outvarid(max_outvar_trunk+9) = idtype('criw','river water pr.basin(olake)'     ,i_mean, 0, 'mm', 'mm', 'mapCRIW.txt', 'timeCRIW.txt')
    outvarid(max_outvar_trunk+10) = idtype('upri','aver upstream riverw'            ,i_mean, 0, 'mm', 'mm', 'mapUPRI.txt', 'timeUPRI.txt') !332
    !------------------------------------------------------------------------------------------------------
    !--LAKES   (ilake, olake, total)
    !---local:     cilv,   lakebasinvol(1)/1000000.,             outvar(145),       Mm3 ilake volume 
    !              clbv,   lakebasinvol(2)/1000000.,             outvar(146),       Mm3 olake volume
    !              cimm,   lakebasinvol(1)*1000/(land+lriv+il),  outvar(???),       mm  ilake volume (per land and local surface water area) 
    !              clmm,   lakebasinvol(1+2)*1000/basin%area,    outvar(???),       mm  lake volume  (per subbasin area) 
    !
    !---upstream:  upla,   upstream_lake,                        outvar(???),       mm  lake volume (per upstream basin area, including lakes and streams)
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+11) = idtype('cilw','ilakew pr.basin(ilake)' ,i_mean , 0,'mm' ,'mm'            ,'mapCILW.txt' ,'timeCILW.txt')
    outvarid(max_outvar_trunk+12) = idtype('claw','olakew pr.basin(olake)' ,i_mean , 0,'mm' ,'mm'            ,'mapCLAW.txt' ,'timeCLAW.txt')
    outvarid(max_outvar_trunk+13) = idtype('upla','aver upstream lakew'    ,i_mean , 0,'mm','mm','mapUPLA.txt','timeUPLA.txt')
    !------------------------------------------------------------------------------------------------------
    !--GLACIERS
    !---local:     glcv,   frozenstate%glacvol(i)*1E-9,             outvar(117),       km3 glacier volume
    !              glwi,   frozenstate%glacvol(i)*1000/(la+lr+il),  outvar(???),       mm glacier volume (per land+local water area)
    !              glwo,   frozenstate%glacvol(i)*1000/basin%area,  outvar(???),       mm glacier volume (per basin area)
    !
    !---upstream:  upgl,   upstream_glacier,                        outvar(???),       mm glacier water (per upstream basin area, including lake and rivers)
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+14) = idtype('glwi','glacierw pr.basin(ilake)',i_mean, 0,'mm' ,'mm'            ,'mapGLWI.txt' ,'timeGLWI.txt')
    outvarid(max_outvar_trunk+15) = idtype('glwo','glacierw pr.basin(olake)',i_mean, 0,'mm' ,'mm'            ,'mapGLWO.txt' ,'timeGLWO.txt')
    outvarid(max_outvar_trunk+16) = idtype('upgl','aver upstream glacierw'  ,i_mean, 0,'mm' ,'mm'            ,'mapUPGL.txt' ,'timeUPGL.txt')
    !------------------------------------------------------------------------------------------------------
    !-FLOWS
    !------------------------------------------------------------------------------------------------------
    !--Runoff
    !---specific runoff(land), specific runoff(ilake), specific runoff(olake, upstream area)
    !
    !     upro,   specific runoff (olake, upstream area),        outvar(247)
    !     uprr,   rec. specific runoff (olake, upstream area),   outvar(???)
    !     crun,   specific runoff (local, land),                 outvar(1)
    !     crui,   specific runoff (local, ilake),                outvar(???)
    !     uprn,   specific runoff(land,upstream),                outvar
    !------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+17) = idtype('uprr','rec. specific discharge(upstream)' ,i_sum, 0,'mm','mm','mapUPRR.txt','timeUPRR.txt')
    outvarid(max_outvar_trunk+18) = idtype('crui','comp. runoff (below ilake)'        ,i_sum, 0,'mm','mm','mapCRUI.txt','timeCRUI.txt')
    outvarid(max_outvar_trunk+19) = idtype('uprn','aver upstr.spec.runoff (land)'     ,i_sum, 0,'mm','mm','mapUPRN.txt','timeUPRN.txt')
    !------------------------------------------------------------------------------------------------------------------
    !--Snowmelt, glaciermelt, snowevaporation, glacierevaporation, soil+vegevaporation   
    !------------------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+20) = idtype('upsm','aver upstream snow melt'      ,i_sum, 0,'mm','mm','mapUPSM.txt','timeUPSM.txt') !342
    outvarid(max_outvar_trunk+21) = idtype('upgm','aver upstream glacier melt'   ,i_sum, 0,'mm','mm','mapUPGM.txt','timeUPGM.txt')
    outvarid(max_outvar_trunk+22) = idtype('upse','aver upstream snow evap'      ,i_sum, 0,'mm','mm','mapUPSE.txt','timeUPSE.txt')
    outvarid(max_outvar_trunk+23) = idtype('upge','aver upstream glacier evap'   ,i_sum, 0,'mm','mm','mapUPGE.txt','timeUPGE.txt')
    outvarid(max_outvar_trunk+24) = idtype('upve','aver upstream soil+veg evap'  ,i_sum, 0,'mm','mm','mapUPVE.txt','timeUPVE.txt')
    outvarid(max_outvar_trunk+25) = idtype('uple','aver upstream lake evap'      ,i_sum, 0,'mm','mm','mapUPLE.txt','timeUPLE.txt')
    outvarid(max_outvar_trunk+26) = idtype('upre','aver upstream river evap'     ,i_sum, 0,'mm','mm','mapUPRE.txt','timeUPRE.txt')
    outvarid(max_outvar_trunk+27) = idtype('upfs','aver upstream snow cover'     ,i_mean, 0,'-','-','mapUPFS.txt','timeUPFS.txt')
    
    !------------------------------------------------------------------------------------------------------------------
    !--Precipitation
    !---local, local(ilake+lriver), upstream
    !      uppr, upstream_prec       outvarid(243) = idtype('uppr','aver upstream prec',i_sum, 0,'mm','mm','mapUPPR.txt','timeUPPR.txt')
    !      upsf, upstream snowfall   outvarid(275) = idtype('upsf','aver upstream snowfall',i_sum, 0,'mm','mm','mapUPSF.txt','timeUPSF.txt')
    !      uprf, upstream rainfall   outvarid(276) = idtype('uprf','aver upstream rainfall',i_sum, 0,'mm','mm','mapUPRF.txt','timeUPRF.txt')
    !------------------------------------------------------------------------------------------------------------------
    !--Evapotranspiration
    !---local, local(ilake+lriver), upstream
    !    outvarid(244) = idtype('upev','aver upstream evap',i_sum, 0,'mm','mm','mapUPEV.txt','timeUPEV.txt')
    !    outvarid(277) = idtype('uppe','aver upstream pot evap',i_sum, 0,'mm','mm','mapUPPE.txt','timeUPPE.txt')
    !------------------------------------------------------------------------------------------------------------------
    !--Remote sensning snow observations
    !---snow water equivalent:
    !              rswe,  -,                                     outvar(o_reswe),    mm (per unit basin landarea)
    !              ucsn,  upstream_compsnow,                     outvar(???),        mm (per unit upstream landarea with observations)
    !              ursn,  upstream_recsnow,                      outvar(???),        mm (per unit upstream landarea with observations)
    !                     upstream_recsnowarea
    !                     upstream_recsnowmissing
    !---fractional snow cover area:
    !                     upstream_compfsc
    !                     upstream_recfsc
    !                     upstream_recfscarea
    !                     upstream_recfscerror
    !                     upstream_recfscmissing
    !------------------------------------------------------------------------------------------------------------------
    outvarid(max_outvar_trunk+28) = idtype('ucsn','aver upstream comp snow'              ,i_mean , 0,'mm'   ,'mm water' ,'mapUCSN.txt' ,'timeUCSN.txt')
    outvarid(max_outvar_trunk+29) = idtype('ursn','aver upstream rec snow'               ,i_mean , 0,'mm'   ,'mm water' ,'mapURSN.txt' ,'timeURSN.txt')
!    outvarid(296) = idtype('uasn','frac nonmissing upstream rec snow'    ,i_mean, 0,'-'     ,'-'        ,'mapUASN.txt' ,'timeUASN.txt')
!    outvarid(297) = idtype('umsn','frac missing upstream rec. snow'      ,i_mean, 0,'-'     ,'-'        ,'mapUMSN.txt' ,'timeUMSN.txt')
!    outvarid(298) = idtype('ucsc','aver upstream comp snowfrac'          ,i_mean, 0,'-'     ,'-'        ,'mapUCSC.txt' ,'timeUCSC.txt')
!    outvarid(299) = idtype('ursc','aver upstream rec snowfrac'           ,i_mean, 0,'-'     ,'-'        ,'mapURSC.txt' ,'timeURSC.txt')
!    outvarid(300) = idtype('uasc','frac nonerror upstream rec. snowfrac' ,i_mean, 0,'-'     ,'-'        ,'mapUASC.txt' ,'timeUASC.txt')
!    outvarid(301) = idtype('urfe','frac error upstream rec. snowfrac'    ,i_mean, 0,'-'     ,'-'        ,'mapURFE.txt' ,'timeURFE.txt')
!    outvarid(302) = idtype('umsc','frac missing upstream rec. snowfrac'  ,i_mean, 0,'-'     ,'-'        ,'mapUMSC.txt' ,'timeUMSC.txt')

    !upstream average temperature
    outvarid(max_outvar_trunk+30) = idtype('upte','average upstream temp'                ,i_mean, 0,'deg'   ,'deg'      ,'mapUPTE.txt' ,'timeUPTE.txt') !352
    
    !water source concentrations (fraction of basin outflow)
    outvarid(max_outvar_trunk+31) = idtype('wfsm','basin.outflow.snowmelt.frac'         ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFSN.txt' ,'timeWFSN.txt')
    outvarid(max_outvar_trunk+32) = idtype('wfgm','basin.outflow.glaciermelt.frac'      ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFGL.txt' ,'timeWFGL.txt')
    outvarid(max_outvar_trunk+33) = idtype('wfrn','basin.outflow.rainfall.frac'         ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFRN.txt' ,'timeWFRN.txt')
    outvarid(max_outvar_trunk+34) = idtype('wfli','basin.outflow.lakeini.frac'          ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFLI.txt' ,'timeWFLI.txt')
    outvarid(max_outvar_trunk+35) = idtype('wfri','basin.outflow.riverini.frac'         ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFRI.txt' ,'timeWFRI.txt')
    outvarid(max_outvar_trunk+36) = idtype('wfsi','basin.outflow.soilini.frac'          ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFSI.txt' ,'timeWFSI.txt')
    outvarid(max_outvar_trunk+37) = idtype('wfsr','basin.outflow.surfrunoff.frac'       ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFSR.txt' ,'timeWFSR.txt')
    outvarid(max_outvar_trunk+38) = idtype('wfdr','basin.outflow.drainage.frac'         ,i_mean, 0,'m3/s'   ,'m3/s'      ,'mapWFDR.txt' ,'timeWFDR.txt')

    !Simulated and Observed Solar and Net Radiation
    outvarid(max_outvar_trunk+39) = idtype('ntrd','net radiation'                 ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapNTRD.txt' ,'timeNTRD.txt')
    outvarid(max_outvar_trunk+40) = idtype('swrd','sw radiation'                  ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapSWRD.txt' ,'timeSWRD.txt') !362
    outvarid(max_outvar_trunk+41) = idtype('upnr','upstream netradiation'         ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapUPNR.txt' ,'timeUPNR.txt')
    outvarid(max_outvar_trunk+42) = idtype('upsw','upstream solarradiation'       ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapUPSW.txt' ,'timeUPSW.txt')
    outvarid(max_outvar_trunk+43) = idtype('rntr','rec. net radiation'            ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapRNTR.txt' ,'timeRNTR.txt')
    outvarid(max_outvar_trunk+44) = idtype('rswr','rec. sw radiation'             ,i_mean, 0,'MJm-2d-1'   ,'MJm-2d-1'      ,'mapRSWR.txt' ,'timeRSWR.txt')

    !Additional temperature diagnostics
    outvarid(max_outvar_trunk+45) = idtype('cdpt','comp deep temp'                ,i_mean, 0,'deg'   ,'deg'      ,'mapCDPT.txt' ,'timeCDPT.txt')
    outvarid(max_outvar_trunk+46) = idtype('cstu','comp soil temp upper 2 layer'  ,i_mean, 0,'deg'   ,'deg'      ,'mapCSTU.txt' ,'timeCSTU.txt') !368


    !Max length: 1st==4, 2nd=20, 3rd=Int, 4th=Int, 5th=5, 6th=20, 7th=20, 8th=20
    !i_mean=medel per dag, i_sum=summa över vald period (eller år för hela sim perioden)
    !i_wmean=volymviktat med vattnet i kolumnen efteråt
    
    !Set variable exchange for criterion calculation (outvarid-index)
    ALLOCATE(changecritvar(3,2))
    changecritvar = RESHAPE(SOURCE=(/51,o_rewstr,282,280,281,323/),SHAPE=(/3,2/))
      
    !Preparation of heading for load files
    loadheadings=['subid ','WetAtm','DryAtm','Fertil','PDecay','RuralA','GrwSIn','IrrSrc','Runoff','RuralB','Urban1','Urban2','Urban3','Rgrwmr','Rgrwol','A     ','B     ','C     ','D     ','E     ','F     ','G     ','H     ','I     ','J     ','K     ','L     ','M     ','N     ']
    
  END SUBROUTINE initiate_output_variables
  
  !>Set variable holding model parameter information with parameters of HYPE model; modparid
  !------------------------------------------------------------------------------------------
  SUBROUTINE initiate_model_parameters()
  !>Note that only small letters allowed in parameter names
  !>modparid form: name, dependence type, index (defined in hypevar)
  
    max_par = 300          !maximum number of parameters
    ALLOCATE(modparid(max_par))
    modparid = modparidtype('',0,0) !default initialisation
 
    modparid(n_lp)= modparidtype('lp        ', m_gpar , m_lp)
    modparid(n_cevpa)   = modparidtype('cevpam    ', m_gpar , m_cevpam)
    modparid(n_cevpp)   = modparidtype('cevpph    ', m_gpar , m_cevpph)
    modparid(4)   = modparidtype('deadl     ', m_gpar , m_deadl)
    modparid(5)   = modparidtype('fastn0    ', m_gpar , m_fastN0)
    modparid(6)   = modparidtype('fastp0    ', m_gpar , m_fastP0)
    modparid(7)   = modparidtype('init1     ', m_gpar , m_iniT1)
    modparid(8)   = modparidtype('init2     ', m_gpar , m_iniT2)
    modparid(9)   = modparidtype('init3     ', m_gpar , m_iniT3)
    modparid(n_rivv)  = modparidtype('rivvel    ', m_gpar , m_rivvel)
    modparid(n_damp)  = modparidtype('damp      ', m_gpar , m_damp)
    modparid(n_tcalt)  = modparidtype('tcalt     ', m_gpar , m_tcalt)
    modparid(13)  = modparidtype('gratk     ', m_gpar , m_grat1)
    modparid(14)  = modparidtype('gratp     ', m_gpar , m_grat2)
    modparid(15)  = modparidtype('denitrlu  ', m_lpar , m_denitrlu)
    modparid(16)  = modparidtype('denitwrm  ', m_gpar , m_denitwr)
    modparid(17)  = modparidtype('maxwidth  ', m_gpar , m_maxwidth)
    modparid(18)  = modparidtype('dissolfn  ', m_lpar , m_dissolfN)
    modparid(19)  = modparidtype('dissolfp  ', m_lpar , m_dissolfP)
    modparid(20)  = modparidtype('litterdays', m_gpar , m_littdays)
    modparid(n_tcea)  = modparidtype('tcelevadd ', m_gpar , m_tcelevadd)
    modparid(n_pcem)  = modparidtype('pcelevmax ', m_gpar , m_pcelevmax)
    modparid(n_tcorr)  = modparidtype('tempcorr  ', m_rpar , m_tcadd)   
    modparid(n_pcea)  = modparidtype('pcelevadd ', m_gpar , m_pcelevadd)
    modparid(25)  = modparidtype('fastlake  ', m_gpar , m_fastlake) 
    modparid(26)  = modparidtype('epotdist  ', m_gpar , m_epotdist)
    modparid(27)  = modparidtype('qmean     ', m_gpar , m_qmean)
    modparid(28)  = modparidtype('pcaddg    ', m_gpar , m_pcaddg)
    modparid(n_pcet)  = modparidtype('pcelevth  ', m_gpar , m_pcelevth)
    modparid(30)  = modparidtype('tpmean    ', m_lrpar, m_tpmean)
    modparid(31)  = modparidtype('ttmp      ', m_lpar , m_ttmp)
    modparid(32)  = modparidtype('cmlt      ', m_lpar , m_cmlt)
    modparid(33)  = modparidtype('cevp      ', m_lpar , m_cevp)
    modparid(34)  = modparidtype('frost     ', m_lpar , m_cfrost)
    modparid(35)  = modparidtype('srrcs     ', m_lpar , m_srrcs)
    modparid(36)  = modparidtype('ttpd      ', m_gpar , m_ttpd)    
    modparid(37)  = modparidtype('ttpi      ', m_gpar , m_ttpi)    
    modparid(38)  = modparidtype('dissolhn  ', m_lpar , m_dissolhN)
    modparid(39)  = modparidtype('dissolhp  ', m_lpar , m_dissolhP)
    modparid(45)  = modparidtype('humusn0   ', m_lpar , m_humusN0) 
    modparid(46)  = modparidtype('partp0    ', m_lpar , m_partP0)  
    modparid(47)  = modparidtype('wcfc      ', m_spar , m_wcfc)
    modparid(48)  = modparidtype('wcwp      ', m_spar , m_wcwp)
    modparid(49)  = modparidtype('wcep      ', m_spar , m_wcep)
    modparid(50)  = modparidtype('wcfc1     ', m_spar , m_wcfc1)
    modparid(51)  = modparidtype('wcwp1     ', m_spar , m_wcwp1)
    modparid(52)  = modparidtype('wcep1     ', m_spar , m_wcep1)
    modparid(53)  = modparidtype('wcfc2     ', m_spar , m_wcfc2)
    modparid(54)  = modparidtype('wcwp2     ', m_spar , m_wcwp2)
    modparid(55)  = modparidtype('wcep2     ', m_spar , m_wcep2)
    modparid(56)  = modparidtype('wcfc3     ', m_spar , m_wcfc3)
    modparid(57)  = modparidtype('wcwp3     ', m_spar , m_wcwp3)
    modparid(58)  = modparidtype('wcep3     ', m_spar , m_wcep3)
    modparid(59)  = modparidtype('gldepi    ', m_gpar , m_gldepi)
    modparid(60)  = modparidtype('trrcs     ', m_spar , m_trrcs)
    modparid(61)  = modparidtype('mperc1    ', m_spar , m_perc1)
    modparid(62)  = modparidtype('mperc2    ', m_spar , m_perc2)
    modparid(65)  = modparidtype('sswcorr   ', m_gpar , m_sswcorr)   
    modparid(66)  = modparidtype('depthrel  ', m_lpar , m_depthrel)
    modparid(67)  = modparidtype('regirr    ', m_gpar , m_regirr)
    modparid(68)  = modparidtype('immdepth  ', m_gpar , m_immdep)
    modparid(69)  = modparidtype('iwdfrac   ', m_gpar , m_iwdfrac)
    modparid(70)  = modparidtype('irrdemand ', m_gpar , m_wdpar)  
    modparid(71)  = modparidtype('minerfn   ', m_lpar , m_minerfN)
    modparid(72)  = modparidtype('minerfp   ', m_lpar , m_minerfP)
    modparid(73)  = modparidtype('degradhn  ', m_lpar , m_degradhN)
    modparid(74)  = modparidtype('deepmem   ', m_gpar , m_deepmem)
    modparid(75)  = modparidtype('surfmem   ', m_lpar , m_surfmem)
    modparid(76)  = modparidtype('freuc     ', m_spar , m_freuc)  
    modparid(77)  = modparidtype('freuexp   ', m_spar , m_freuexp)
    modparid(78)  = modparidtype('freurate  ', m_spar , m_freurate)
    modparid(80)  = modparidtype('wprodn    ', m_gpar , m_wprodn) 
    modparid(81)  = modparidtype('sedon     ', m_gpar , m_sedon)
    modparid(82)  = modparidtype('sedpp     ', m_gpar , m_sedpp)
    modparid(83)  = modparidtype('sedexp    ', m_gpar , m_sedexp)
    modparid(84)  = modparidtype('rcgrw     ', m_gpar , m_rcgrw)
    modparid(85)  = modparidtype('rrcs1     ', m_spar , m_rrcs1)
    modparid(86)  = modparidtype('rrcs2     ', m_spar , m_rrcs2)
    modparid(n_rrcs3)  = modparidtype('rrcs3     ', m_gpar , m_rrcs3)
    modparid(88)  = modparidtype('hnhalf    ', m_lpar , m_hNhalf)
    modparid(89)  = modparidtype('pphalf    ', m_lpar , m_pPhalf)
    modparid(90)  = modparidtype('locsoil   ', m_gpar , m_locsoil)
    modparid(91)  = modparidtype('bufffilt  ', m_lpar , m_filtPbuf)
    modparid(92)  = modparidtype('innerfilt ', m_lpar , m_filtPinner)
    modparid(93)  = modparidtype('otherfilt ', m_lpar , m_filtPother)
    modparid(94)  = modparidtype('drydeppp  ', m_lpar , m_drypp)
    modparid(95)  = modparidtype('wetdepsp  ', m_gpar , m_wetsp)
    modparid(96)  = modparidtype('rivvel1   ', m_lrpar, m_velpar1)
    modparid(97)  = modparidtype('rivvel2   ', m_lrpar, m_velpar2)
    modparid(98)  = modparidtype('rivvel3   ', m_lrpar, m_velpar3)
    modparid(99)  = modparidtype('rivwidth1 ', m_lrpar, m_widpar1)
    modparid(100) = modparidtype('rivwidth2 ', m_lrpar, m_widpar2)
    modparid(101) = modparidtype('rivwidth3 ', m_lrpar, m_widpar3)
    modparid(102) = modparidtype('srrate    ', m_spar , m_srrate)
    modparid(103) = modparidtype('macrate   ', m_spar , m_macrate) 
    modparid(104) = modparidtype('mactrinf  ', m_spar , m_mactrinf)
    modparid(105) = modparidtype('mactrsm   ', m_spar , m_mactrsm) 
    modparid(106) = modparidtype('soilcoh   ', m_spar , m_soilcoh) 
    modparid(107) = modparidtype('soilerod  ', m_spar , m_soilerod)
    modparid(108) = modparidtype('wprodp    ', m_gpar , m_wprodp)
    modparid(109) = modparidtype('sreroexp  ', m_gpar , m_sreroexp)
    modparid(110) = modparidtype('snowdensdt', m_gpar , m_dsndens)
    modparid(111) = modparidtype('sfrost    ', m_spar , m_sfrost)
    modparid(112) = modparidtype('macrofilt ', m_spar , m_macfilt)
    modparid(113) = modparidtype('fertdays  ', m_gpar , m_fertdays)
    modparid(114) = modparidtype('humusp0   ', m_lpar , m_humusP0)   
    modparid(115) = modparidtype('hphalf    ', m_lpar , m_hPhalf)   
    modparid(116) = modparidtype('degradhp  ', m_lpar , m_degradhP)
    modparid(n_cevpc) = modparidtype('cevpcorr  ', m_rpar , m_cevpcorr)
    modparid(118) = modparidtype('minc      ', m_gpar , m_minc)
    modparid(119) = modparidtype('humusc1   ', m_lpar , m_humusC1)
    modparid(120) = modparidtype('fastc1    ', m_lpar , m_fastC1)
    modparid(121) = modparidtype('klh       ', m_gpar , m_crate1)
    modparid(122) = modparidtype('klo       ', m_gpar , m_crate2)
    modparid(123) = modparidtype('kho       ', m_gpar , m_crate3)
    modparid(124) = modparidtype('tcobselev ', m_gpar , m_tcobselev)
    modparid(125) = modparidtype('ripz      ', m_lpar , m_ripz)
    modparid(126) = modparidtype('ripe      ', m_gpar , m_ripe)
    modparid(127) = modparidtype('sedoc     ', m_gpar , m_sedoc)
    modparid(128) = modparidtype('koc       ', m_gpar , m_crate5)
    modparid(129) = modparidtype('kcgwreg   ', m_gpar , m_crate6)
    modparid(130) = modparidtype('rips      ', m_gpar , m_rips)
    modparid(131) = modparidtype('snowdens0 ', m_gpar , m_sndens0)
    modparid(132) = modparidtype('pprelmax  ', m_gpar , m_pprelmax)
    modparid(133) = modparidtype('pprelexp  ', m_gpar , m_pprelexp)
    modparid(134) = modparidtype('tnmean    ', m_lrpar, m_tnmean)
    modparid(135) = modparidtype('tocmean   ', m_lrpar, m_tocmean)
    modparid(136) = modparidtype('humusc2   ', m_lpar , m_humusC2)
    modparid(137) = modparidtype('humusc3   ', m_lpar , m_humusC3)
    modparid(138) = modparidtype('fastc2    ', m_lpar , m_fastC2)
    modparid(139) = modparidtype('fastc3    ', m_lpar , m_fastC3)
    modparid(n_rrcsc) = modparidtype('rrcscorr  ', m_rpar , m_rrcscorr)
    modparid(141) = modparidtype('kof       ', m_gpar , m_crate9)
    modparid(142) = modparidtype('koflim    ', m_gpar , m_crate10)
    modparid(143) = modparidtype('grata     ', m_gpar , m_grat3)
    modparid(144) = modparidtype('limqprod  ', m_gpar , m_limprod)
    modparid(145) = modparidtype('partp1    ', m_lpar , m_partP1)
    modparid(146) = modparidtype('partp2    ', m_lpar , m_partP2)
    modparid(147) = modparidtype('partp3    ', m_lpar , m_partP3)
    modparid(148) = modparidtype('wprodc    ', m_gpar , m_wprodc)
    modparid(149) = modparidtype('denitwrl  ', m_gpar,  m_denitwrl) !denitw local river
    modparid(150) = modparidtype('qmean     ', m_ldpar, m_ldqmean)
    modparid(151) = modparidtype('tpmean    ', m_ldpar, m_ldtpmean)
    modparid(152) = modparidtype('tnmean    ', m_ldpar, m_ldtnmean)
    modparid(153) = modparidtype('tocmean   ', m_ldpar, m_ldtocmean)
    modparid(154) = modparidtype('wprodn    ', m_ldpar, m_ldwprodn) 
    modparid(155) = modparidtype('sedon     ', m_ldpar, m_ldsedon)
    modparid(156) = modparidtype('sedoc     ', m_ldpar, m_ldsedoc)
    modparid(157) = modparidtype('sedpp     ', m_ldpar, m_ldsedpp)
    modparid(158) = modparidtype('wprodp    ', m_ldpar, m_ldwprodp) 
    modparid(159) = modparidtype('limqprod  ', m_ldpar, m_ldlimprod)
    modparid(160) = modparidtype('deeplake  ', m_gpar , m_deeplake)  
    modparid(161) = modparidtype('denitwl   ', m_gpar , m_denitwl)   
    modparid(162) = modparidtype('denitwl   ', m_ldpar, m_lddenitwl) 
    modparid(163) = modparidtype('wprodc    ', m_ldpar, m_ldwprodc)  
    modparid(164) = modparidtype('prodpp    ', m_ldpar, m_ldprodpp)  
    modparid(165) = modparidtype('prodsp    ', m_ldpar, m_ldprodsp)  
    modparid(166) = modparidtype('deeplake  ', m_ldpar, m_lddeeplake)
    modparid(167) = modparidtype('fastlake  ', m_ldpar, m_ldfastlake)
    modparid(168) = modparidtype('laketemp  ', m_gpar , m_laketemp)  
    modparid(169) = modparidtype('deadm     ', m_gpar , m_deadm)    
    modparid(170) = modparidtype('incorr    ', m_wqrpar, m_incorr) 
    modparid(171) = modparidtype('oncorr    ', m_wqrpar, m_oncorr) 
    modparid(172) = modparidtype('phoscorr  ', m_wqrpar, m_phoscorr)  
    modparid(173) = modparidtype('ratcorr   ', m_rpar , m_ratcorr)   
    modparid(174) = modparidtype('ponatm    ', m_lpar , m_ponatm)    
    modparid(175) = modparidtype('preccorr  ', m_rpar , m_preccorr)
    modparid(176) = modparidtype('cirrsink  ', m_rpar , m_cirrsink)  
    modparid(177) = modparidtype('irrcomp   ', m_gpar , m_irrcomp)   
    modparid(178) = modparidtype('ocsoimsat ', m_lpar , m_ocsoim)     
    modparid(179) = modparidtype('ocsoimslp ', m_lpar , m_ocsmslp)     
    modparid(180) = modparidtype('pirrs     ', m_rpar , m_pirrs)     
    modparid(181) = modparidtype('pirrg     ', m_rpar , m_pirrg)     
    modparid(182) = modparidtype('onpercred ', m_lpar , m_onpercred) 
    modparid(183) = modparidtype('pppercred ', m_lpar , m_pppercred) 
    modparid(184) = modparidtype('onconc0   ', m_lpar , m_onconc0)   
    modparid(185) = modparidtype('ppconc0   ', m_lpar , m_ppconc0)   
    modparid(186) = modparidtype('occonc0   ', m_lpar , m_occonc0)  
    modparid(187) = modparidtype('snalbmin  ', m_lpar , m_snalbmin) 
    modparid(188) = modparidtype('snalbmax  ', m_lpar , m_snalbmax) 
    modparid(189) = modparidtype('snalbkexp ', m_lpar , m_snalbkexp)
    modparid(190) = modparidtype('cmrad     ', m_lpar , m_cmrad)  
    modparid(191) = modparidtype('pcluse    ', m_lpar , m_pcluse)
    modparid(192) = modparidtype('aloadconst', m_gpar , m_atmload)
    
    !Water T2 temperature parameters
    modparid(193) = modparidtype('t2trriver ', m_gpar , m_t2trriver)  !temp flow from air to river
    modparid(194) = modparidtype('t2trlake  ', m_gpar , m_t2trlake)   !temp flow from air to lake
    modparid(195) = modparidtype('krelflood ', m_gpar , m_krelflood)  !some flood control dam parameter
    modparid(196) = modparidtype('upper2deep', m_gpar , m_upper2deep) !temp flow from upper to lower lake layer

    !Lake ice parameters
    modparid(197) = modparidtype('licewme   ', m_gpar , m_licewme)    !lake ice, water melt efficiency
    modparid(198) = modparidtype('licetf    ', m_gpar , m_licetf)     !lake ice, freezing temperature     
    modparid(199) = modparidtype('licesndens', m_gpar , m_licesndens) !lake ice, snow compaction parameter
    modparid(200) = modparidtype('licekika  ', m_gpar , m_licekika)   !lake ice, ki/ka     
    modparid(201) = modparidtype('licekexp  ', m_gpar , m_licekexp)   !lake ice, ks = ki*(dsnow(dice)^kexp, 1.88 
    modparid(202) = modparidtype('licetmelt ', m_gpar , m_licetmelt)  !lake ice, degreeday factor for ice melt
    modparid(203) = modparidtype('licewcorr ', m_gpar , m_licewcorr)  !lake ice, snow fall reduction for winddrift 

    !River ice parameters
    modparid(204) = modparidtype('ricewme   ', m_gpar , m_ricewme)    !river ice, water melt efficiency
    modparid(205) = modparidtype('ricetf    ', m_gpar , m_ricetf)     !river ice, freezing temperature     
    modparid(206) = modparidtype('ricesndens', m_gpar , m_ricesndens) !river ice, snow compaction parameter
    modparid(207) = modparidtype('ricekika  ', m_gpar , m_ricekika)   !river ice, ki/ka     
    modparid(208) = modparidtype('ricekexp  ', m_gpar , m_ricekexp)   !river ice, ki/ks
    modparid(209) = modparidtype('ricetmelt ', m_gpar , m_ricetmelt)  !river ice, degreeday factor for ice melt
    
    !snow cover area parameters
    modparid(210) = modparidtype('fscmax    ', m_gpar, m_fscmax)      !maximum snow cover area (0.95 [-])
    modparid(211) = modparidtype('fscmin    ', m_gpar, m_fscmin)      !minimum fsc             (0.001 [-])
    modparid(212) = modparidtype('fsclim    ', m_gpar, m_fsclim)      !fsc limit for onset of snowmax (0.001 [-])
    modparid(213) = modparidtype('fscdistmax', m_lpar, m_fscdistmax)  !maximum snow distribution factor (0.8 [mm^-1])
    modparid(214) = modparidtype('fscdist0  ', m_lpar, m_fscdist0)    !minimum snow distribution factor (0.6 [mm^-1])
    modparid(215) = modparidtype('fscdist1  ', m_lpar, m_fscdist1)    !elev_std coefficient for snow distr factor (0.001 [m^-1])
    modparid(216) = modparidtype('fsck1     ', m_gpar, m_fsck1)       !time constant in decrease of snowmax during melt (0.2 [-])
    modparid(217) = modparidtype('fsckexp   ', m_gpar, m_fsckexp)     !exponential coefficient in decrease of snowmax during melt (1e-6 [s^-1])
    
    !parameters for radiation and optional potential evaporation calculations
    modparid(218) = modparidtype('krs       ', m_gpar, m_krs)         !Hargreaves adjustment factor in estimation of shortwave radiation from tmin and tmax 0.16 (inland) 0.19 (coastal)
    modparid(219) = modparidtype('jhtadd    ', m_gpar, m_jhtadd)      !PET parameter in Jensen-Haise/McGuinness following Oudin et al (2005), recommended value 5
    modparid(220) = modparidtype('jhtscale  ', m_gpar, m_jhtscale)    !PET parameter in Jensen-Haise/McGuinness following Oudin et al (2005), recommended value 100
    modparid(221) = modparidtype('alfapt    ', m_gpar, m_alfapt)      !Priestly-taylor coefficient, recommended value 1.26
    modparid(222) = modparidtype('mwind     ', m_gpar, m_mwind)       !Mean windspeed, to replace missing data in FAO Penman-Monteith pet-equation, suggested value 2 m/s
    modparid(223) = modparidtype('kc        ', m_lpar, m_kc)          !Landuse dependent crop coefficient, used to scale the optional reference PET estimates (Hargreaves, Jensen, Priestly-Taylor,Penman-Monteith....)
    modparid(224) = modparidtype('alb       ', m_lpar, m_alb)         !Landuse dependent albedo, used for net radiation calculation (Priestly-Taylor,Penman-Monteith), suggested value 0.23

    !Glacier parameters (previous non-optional constants in the code)
    modparid(225) = modparidtype('glacvcoef ', m_gpar, m_glacvcoef)   !Coefficient glacier volume-area relationship (default 0.205)
    modparid(226) = modparidtype('glacvexp  ', m_gpar, m_glacvexp)    !Exponent glacier volume-area relationship (default 1.375)
    modparid(227) = modparidtype('glacdens  ', m_gpar, m_glacdens)    !Glacier density (default 0.85 m3 water/m3 ice)
    modparid(228) = modparidtype('glacvcoef1', m_gpar, m_glacvcoef1)  !Coefficient glacier volume-area relationship, type 1 glacier (default 1.701)
    modparid(229) = modparidtype('glacvexp1 ', m_gpar, m_glacvexp1)   !Exponent glacier volume-area relationship, type 1 glacier (default 1.25)
    modparid(230) = modparidtype('glac2arlim', m_gpar, m_glac2arlim)  !Area limit to separate glacer type 1 (Glaciers or Small) from type 2 (Ice Caps or Large)

    modparid(231) = modparidtype('rcgrwst   ', m_spar, m_rcgrwst)     !Deep percolation to aquifer 
    modparid(232) = modparidtype('aqretcor  ', m_rpar, m_aqretcorr)   !Aquifer return flow adjustment parameter
    modparid(233) = modparidtype('aqdelcor  ', m_rpar, m_aqdelcorr)   !Aquifer percolation delay adjustment parameter
    modparid(234) = modparidtype('aqpercor  ', m_rpar, m_aqpercorr)   !Aquifer percolation adjustment parameter
    modparid(235) = modparidtype('denitaq   ', m_gpar, m_denitaq)     !Aquifer denitrification 
    
    !Additional Lake and River water temperature parameters, to be replacing the t2trlake and t2trriver parameters
    modparid(236) = modparidtype('tcfriver  ', m_gpar, m_tcfriver)    !air-riverwater heat flow, temperature difference coefficient
    modparid(237) = modparidtype('scfriver  ', m_gpar, m_scfriver)    !air-riverwater heat flow, solar radiation coefficient
    modparid(238) = modparidtype('ccfriver  ', m_gpar, m_ccfriver)    !air-riverwater heat flow, constant coefficient
    modparid(239) = modparidtype('lcfriver  ', m_gpar, m_lcfriver)    !air-riverwater heat flow, linear coefficient
    modparid(240) = modparidtype('tcflake   ', m_gpar, m_tcflake)      !air-lakewater heat flow, temperature difference coefficient
    modparid(241) = modparidtype('scflake   ', m_gpar, m_scflake)      !air-lakewater heat flow, solar radiation coefficient
    modparid(242) = modparidtype('ccflake   ', m_gpar, m_ccflake)      !air-lakewater heat flow, constant coefficient
    modparid(243) = modparidtype('lcflake   ', m_gpar, m_lcflake)      !air-lakewater heat flow, linear coefficient
    modparid(244) = modparidtype('stbcorr1  ', m_gpar, m_stbcorr1)    !parameter for stability correction
    modparid(245) = modparidtype('stbcorr2  ', m_gpar, m_stbcorr2)    !parameter for stability correction
    modparid(246) = modparidtype('stbcorr3  ', m_gpar, m_stbcorr3)    !parameter for stability correction
    modparid(247) = modparidtype('zwind     ', m_gpar, m_zwind)       !wind observation level
    modparid(248) = modparidtype('zwish     ', m_gpar, m_zwish)       !wanted wind observation level for PET
    modparid(249) = modparidtype('zpdh      ', m_gpar, m_zpdh)        !zero plane displacement height for wind
    modparid(250) = modparidtype('roughness ', m_gpar, m_roughness)   !surface roughness for wind
    
    !Additional precipitation correction; orographic (pcelevstd) and phase (pcusnow, pcurain) impact on undercatch
    modparid(251) = modparidtype('pcelevstd ', m_gpar, m_pcelevstd)    !fractional increase in precipitation per 100m elev.std
    modparid(n_pcur) = modparidtype('pcurain   ', m_gpar, m_pcurain)      !undercatch correction factor for rain
    modparid(n_pcus) = modparidtype('pcusnow   ', m_gpar, m_pcusnow)      !undercatch correction factor for snow

    !Snow evaporation, refreeze, liquid water content holding capacity
    modparid(254) = modparidtype('fepotsnow ', m_lpar, m_fepotsnow)   !fraction of potential evaporation used for snow evaporation
    modparid(255) = modparidtype('cmrefr    ', m_gpar, m_cmrefr)      !snow refreeze efficiency (fraction of degree day factor cmlt)
    modparid(256) = modparidtype('fsceff    ', m_gpar, m_fsceff)      !efficiency of fractional snow cover to reduce melt and evap
    
    !Separate glacier melt parameters and sublimation/evaporation parameters
    modparid(257) = modparidtype('glacalb   ', m_gpar, m_glacalb)           !glacier ice albedo (0.35)
    modparid(258) = modparidtype('glacttmp  ', m_gpar, m_glacttmp)
    modparid(259) = modparidtype('glaccmlt  ', m_gpar, m_glaccmlt)
    modparid(260) = modparidtype('glaccmrad ', m_gpar, m_glaccmrad)
    modparid(261) = modparidtype('glaccmrefr', m_gpar, m_glaccmrefr)
    modparid(262) = modparidtype('fepotglac ', m_gpar, m_fepotglac)
    
    modparid(263) = modparidtype('kthrflood ', m_gpar, m_kthrflood)      !Threshold inflow over which a flood control dam save water (fraction of max inflow)
    modparid(264) = modparidtype('klowflood ', m_gpar, m_klowflood)      !Threshold level for extra flood control releases (fraction of regvol, typical 1/3)
    modparid(265) = modparidtype('monthlapse', m_mpar, m_mlapse)
    modparid(266) = modparidtype('limsedon  ', m_gpar, m_limsedon)
    modparid(267) = modparidtype('limsedpp  ', m_gpar, m_limsedpp)
    
   !Soil temperature control on soil evapotranspiration
    modparid(268) = modparidtype('ttrig      ', m_lpar , m_ttrig)
    modparid(269) = modparidtype('treda      ', m_lpar , m_treda)
    modparid(270) = modparidtype('tredb      ', m_lpar , m_tredb)
    
    !Interception loss - fixed fraction of rainfall and snowfall
    modparid(271) = modparidtype('sncluse    ', m_lpar , m_sncluse)
    modparid(272) = modparidtype('rncluse    ', m_lpar , m_rncluse)
 
    !Frozen soils infiltration
    modparid(273) = modparidtype('bfroznsoil', m_spar, m_bfroznsoil) ! empirical coefficient = 2.10 & 1.14 for prairie & forest soils, respectively
    
    !ilake and olake regional parameters
    modparid(274) = modparidtype('ilrratk   ', m_ilrpar , m_ilrrat1)
    modparid(275) = modparidtype('ilrratp   ', m_ilrpar , m_ilrrat2)
    modparid(276) = modparidtype('ilrldep   ', m_ilrpar , m_ilrldep)
    modparid(277) = modparidtype('ilricatch ', m_ilrpar , m_ilricatch)
    modparid(278) = modparidtype('olrratk   ', m_olrpar , m_olrrat1)
    modparid(279) = modparidtype('olrratp   ', m_olrpar , m_olrrat2)
    modparid(280) = modparidtype('olrldep   ', m_olrpar , m_olrldep)
    
    !additional general olake/ilake parameters (used only if not set by geodata/lakedata/damdata or by olakeregion parameters
    modparid(281) = modparidtype('gldepo    ', m_gpar   , m_gldepo)
    modparid(282) = modparidtype('gicatch   ', m_gpar   , m_gicatch)

    !glacier initial volume scaling coefficicent (kg/m2/year)
    modparid(283) = modparidtype('glacinimb ', m_gpar   , m_glacinimb)
    
    !MM2016: correction factors from A-HYPE parameters to Hudson
    modparid(284) = modparidtype('kc_corr   ', m_gpar , m_kccorr)
    modparid(285) = modparidtype('fc_corr   ', m_gpar , m_fccorr)
    modparid(286) = modparidtype('wp_corr   ', m_gpar , m_wpcorr)
    modparid(287) = modparidtype('fpsno_corr', m_gpar , m_fpsnocorr)
    modparid(288) = modparidtype('srrc_corr ', m_gpar , m_srrccorr)
    modparid(289) = modparidtype('rrc_corr  ', m_gpar , m_rrccorr)
    modparid(290) = modparidtype('ep_corr   ', m_gpar , m_epcorr)
    modparid(291) = modparidtype('cmlt_corr ', m_gpar , m_cmltcorr)
    modparid(292) = modparidtype('surfm_corr   ', m_gpar , m_surfmcorr)
    modparid(293) = modparidtype('deprl_corr ', m_gpar , m_deprlcorr)


    !MM2016: Groundwater runoff for frozen soils
    modparid(294) = modparidtype('logsatmp', m_spar, m_logsatmp)
    modparid(295) = modparidtype('bcosby', m_spar, m_bcosby)
    
    !MH2017: river variables
    modparid(296) = modparidtype('nc_lim1', m_gpar, m_nclim1)
    modparid(297) = modparidtype('nc_lim2', m_gpar, m_nclim2)
    modparid(298) = modparidtype('nc_fact1', m_gpar, m_ncfact1)
    modparid(299) = modparidtype('nc_fact2', m_gpar, m_ncfact2)
    modparid(300) = modparidtype('nc_fact3', m_gpar, m_ncfact3)
    

  END SUBROUTINE initiate_model_parameters

!>Set flags for special HYPE models that will be used
!------------------------------------------------------------
  SUBROUTINE set_special_models(nc,arr,glacexist,irrexist)

  USE IRRIGATION_MODULE, ONLY : irrtype,  &
                                check_for_irrigated_classes
  !Argument declarations
  INTEGER,INTENT(IN) :: nc         !<dimension, number of classes
  INTEGER,INTENT(IN) :: arr(nc)    !<soilmodels
  LOGICAL, INTENT(OUT) :: glacexist   !<status of glacier model
  LOGICAL, INTENT(OUT) :: irrexist    !<status of irrigation model

  INTEGER i

    !Check for glacier soilmodel
    glacexist = .FALSE.
    DO i=1,nc
      IF(arr(i)==glacier_model) glacexist=.TRUE.
    ENDDO
  
    !Check for irrigation included in model set-up
    irrexist = .TRUE.
    IF(.NOT.ALLOCATED(irrigationsystem))THEN
      irrexist = .FALSE.      !Irrigation information not provided in MgmtData
    ELSEIF(.NOT.ALLOCATED(cropirrdata))THEN
      irrexist = .FALSE.      !Irrigation information not provided in CropData
    ELSE
      IF(.NOT.ALLOCATED(irrtype))THEN
        ALLOCATE(irrtype(nc))
        CALL check_for_irrigated_classes(nc)
      ENDIF
      IF(SUM(irrtype)==0)THEN
        irrexist = .FALSE.
        DEALLOCATE(irrtype)
      ENDIF
    ENDIF

  END SUBROUTINE set_special_models
  
!>Set HYPE parameters needed by HYSS
!------------------------------------------------------------
  SUBROUTINE set_special_parameters_from_model(velindex,dampindex)

  !Argument declarations
  INTEGER, INTENT(OUT) :: velindex  !<index of rivvel in modparid
  INTEGER, INTENT(OUT) :: dampindex !<index of damp in modparid

    velindex = n_rivv
    dampindex = n_damp

  END SUBROUTINE set_special_parameters_from_model
  
!>Calculate special HYPE parameters needed by HYSS
!------------------------------------------------------------
  SUBROUTINE calculate_special_model_parameters(nsubb,optrivvel,optdamp,dimriver)

  USE SURFACEWATER_PROCESSES, ONLY : calculate_landarea_riverlength
  USE HYPE_INDATA, ONLY : clust_group,catdes, &
                          ndes_reg, idx_catdes, wght_catdes, indx_par
  
  !Argument declarations
  INTEGER, INTENT(IN)  :: nsubb     !<dimension, number of subbasins (submodel)
  REAL, INTENT(IN)     :: optrivvel !<lower river velocity boundary (optpar)
  REAL, INTENT(IN)     :: optdamp   !<lower damp boundary (optpar)
  INTEGER, INTENT(OUT) :: dimriver  !<maximum size river lag needed

  INTEGER i,k,itype
  REAL maxtrans
  REAL landarea(nsubb)
  REAL rivlength(2,nsubb)   !river length (m)
  REAL totaltime,transtime  !total time in river and time in river train (translation) (in timesteps)
  REAL localrivvel,localdamp

    !Calculate local and main river length
    CALL calculate_landarea_riverlength(nsubb,landarea,rivlength)
    
    !Calculate river translation time, and the maximum, for every subbasin
    maxtrans = 0.
    localrivvel = genpar(m_rivvel)
    localdamp = genpar(m_damp)
    IF(.NOT.conductregest)THEN
      localrivvel = MIN(localrivvel,optrivvel)    !This will give the "worst" combination
      localdamp = MIN(localdamp,optdamp)
      IF(localrivvel==0)THEN !Test for missing rivvel in model set-up
        dimriver = 1
        WRITE(6,*) 'Warning: Parameter rivvel not set. River flow will take one time step.'
        RETURN
      ENDIF
      totaltime = MAXVAL(rivlength) / localrivvel / seconds_per_timestep
      transtime = (1. - localdamp) * totaltime
      dimriver = INT(transtime) + 1
    ELSE
      !If regional parameter estimations is used calculate rivvel and damp
      DO i = 1,nsub
        IF(indx_par(n_rivv).gt.0)THEN
          localrivvel=0.0
          DO k = 1,ndes_reg(clust_group(i),indx_par(n_rivv))
            localrivvel = localrivvel + wght_catdes(clust_group(i),indx_par(n_rivv),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rivv),k))
          ENDDO
        ENDIF
        IF(indx_par(n_damp).gt.0) THEN
          localdamp=0.0
          DO k = 1,ndes_reg(clust_group(i),indx_par(n_damp))
            localdamp = localdamp + wght_catdes(clust_group(i),indx_par(n_damp),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_damp),k))
          ENDDO
        ENDIF   
        localrivvel = MIN(localrivvel,optrivvel)    !This will give the "worst" combination
        localdamp = MIN(localdamp,optdamp)
        IF(localrivvel==0)THEN !Test for missing rivvel in model set-up
          dimriver = 1
          WRITE(6,*) 'Warning: Parameter rivvel not set. River flow will take one time step.'
          RETURN
        ENDIF
        DO itype = 1,2          
          totaltime = rivlength(itype,i) / localrivvel / seconds_per_timestep
          transtime = (1. - localdamp) * totaltime
          maxtrans = MAX(maxtrans,transtime)
        ENDDO
      ENDDO
      dimriver = INT(maxtrans) + 1
    ENDIF

  END SUBROUTINE calculate_special_model_parameters
  
  !>Set some model configuration data from input parameters
  !--------------------------------------------------------
  SUBROUTINE set_modelconfig_from_parameters()
    !Variable declarations
    INTEGER i     !loop-variable

    !Set basin ilake and olake lakedepth/icatch from parameter inputs (if not set already by geodata/lakedata/damdata input)
    IF(slc_ilake>0)THEN
      DO i = 1,nsub
        !use ilake region depth parameter if available, otherwise general gldepi parameter
        IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))THEN
          basin(i)%lakedepth(1) = ilregpar(m_ilrldep,basin(i)%ilakeregion)
        ELSE
          basin(i)%lakedepth(1) = genpar(m_gldepi)
        ENDIF
        !use general or ilake region ICATCH parameter, if NOT set by geodata
        IF(basin(i)%ilakecatch.LE.0)THEN
          IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))THEN
            basin(i)%ilakecatch = ilregpar(m_ilricatch,basin(i)%ilakeregion)
          ELSE
            basin(i)%ilakecatch = genpar(m_gicatch)
          ENDIF
          IF(basin(i)%ilakecatch.LE.0)basin(i)%ilakecatch = 1.
        ENDIF
      ENDDO
    ENDIF
    IF(slc_olake>0)THEN
      DO i = 1,nsub
        !Check if lakedepth(2) is not yet set (=0) by GEODATA/LAKEDATA/DAMDATA=0
        IF(basin(i)%lakedepth(2).LE.0.)THEN
          !Use olake region parameter if available, otherwise general olake depth parameter gldepo 
          IF(basin(i)%olakeregion.GT.0 .AND. ALLOCATED(olregpar))THEN
            basin(i)%lakedepth(2) = olregpar(m_olrldep,basin(i)%olakeregion)
          ELSE
            basin(i)%lakedepth(2) = genpar(m_gldepo)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE set_modelconfig_from_parameters
  
  !>Initiates the model state for a simulation with no saved states. 
  !-----------------------------------------
  SUBROUTINE initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE GLACIER_SOILMODEL, ONLY :   initiate_glacier_state
    USE SOIL_PROCESSES, ONLY :      initiate_soil_water_state
    USE NPC_SOIL_PROCESSES, ONLY :  initiate_soil_npc_state
    USE SURFACEWATER_PROCESSES, ONLY : sum_upstream_area, &
                                       calculate_landarea_riverlength
    USE NPC_SURFACEWATER_PROCESSES, ONLY : initiate_river_npc_state,     &
                                           initiate_lake_npc_state
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : initiate_aquifer_state

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Variable declarations
    INTEGER i     !loop-variable

    !Parameter declarations
    REAL, PARAMETER :: seconds_per_day = 86400.  

    !Existence of soillayers
    IF(MAXVAL(soilthick(3,1:nclass))>0.)THEN
      nummaxlayers = 3
    ELSEIF(MAXVAL(soilthick(2,1:nclass))>0.)THEN
      nummaxlayers = 2
    ELSE
      nummaxlayers = 1
    ENDIF
    
    !Initiate states to zero
    CALL initiate_state_zero(numsubstances,naquifers,i_t2>0,wetlandexist, &
                             glacierexist,modeloption(p_lakeriverice)>=1, &
                             doirrigation,doupdate(i_qar).OR.doupdate(i_war), &
                             modeloption(p_growthstart)==1, &
                             frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Initiate soil water state variables and soil water processes help parameters.
    CALL initiate_soil_water_state(soilstate)
    
    !Initialize glacier volume
    IF(glacierexist) CALL initiate_glacier_state(nclass,classmodel,frozenstate)
    
    !Initialize lake state for pure water model simulation to total lake volume
    DO i = 1,nsub
!      IF(slc_ilake>0) lakestate%water(1,i) = genpar(m_gldepi)*1000.  !ilake water stage (mm)
!      IF(slc_ilake>0)THEN
!        lakestate%water(1,i) = genpar(m_gldepi)*1000.  !ilake water stage (mm)
!        IF(basin(i)%ilakeregion.GT.0 .AND. ALLOCATED(ilregpar))lakestate%water(1,i) = ilregpar(m_ilrldep,basin(i)%ilakeregion)*1000.
!      ENDIF
      IF(slc_ilake>0) lakestate%water(1,i) = basin(i)%lakedepth(1)*1000.  !ilake water stage (mm)
      IF(slc_olake>0)THEN
        !Check if lakedepth(2)=0 and olakeregion parameter olrldep>0
!        if(basin(i)%lakedepth(2).LE.0. .AND. basin(i)%olakeregion.GT.0 .AND. ALLOCATED(olregpar)) basin(i)%lakedepth(2) = olregpar(m_olrldep,basin(i)%olakeregion)

        lakestate%water(2,i) = basin(i)%lakedepth(2) * 1000.            !ordinary olake water stage (mm)
        IF(ALLOCATED(damindex))THEN
          IF(damindex(i)>0)THEN
            IF(dam(damindex(i))%purpose==3) THEN
              lakestate%water(2,i) = 0                               !dam for flood control start empty
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
    !Initialize river states
    IF(.NOT. ALLOCATED(landarea))     ALLOCATE(landarea(nsub))
    IF(.NOT. ALLOCATED(riverlength))  ALLOCATE(riverlength(2,nsub))
    IF(.NOT. ALLOCATED(deadriver))    ALLOCATE(deadriver(2,nsub))
    IF(.NOT. ALLOCATED(upstreamarea)) ALLOCATE(upstreamarea(nsub))
    CALL sum_upstream_area(nsub,upstreamarea)             !calculates the area upstream of the outlet point of catchment
    CALL calculate_landarea_riverlength(nsub,landarea,riverlength)
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadriver(1,i) = genpar(m_deadl) * (basin(i)%area/1.0E6) * riverlength(1,i)
      deadriver(2,i) = genpar(m_deadm) * (upstreamarea(i)/1.0E6) * riverlength(2,i)
      riverstate%water(:,i) = deadriver(:,i)    !initialize river damping box to deadvolume (m3)
    ENDDO
    
    !Initiate average discharge variable, Qmeani = 365 day mean Q before damping  
    DO i = 1,nsub
      IF(ALLOCATED(lakedatapar))THEN
        riverstate%Qmean(1,i) = lakedatapar(lakedataparindex(i,1),m_ldqmean)* 0.001 * basin(i)%area / (365. * seconds_per_day)
        riverstate%Qmean(2,i) = lakedatapar(lakedataparindex(i,2),m_ldqmean)* 0.001 * upstreamarea(i) / (365. * seconds_per_day)
      ELSE
        riverstate%Qmean(1,i) = genpar(m_Qmean)* 0.001 * basin(i)%area / (365. * seconds_per_day)
        riverstate%Qmean(2,i) = genpar(m_Qmean)* 0.001 * upstreamarea(i) / (365. * seconds_per_day)
      ENDIF
    ENDDO
   
    !Initialize soil nutrient variables
    CALL initiate_soil_npc_state(conductN,conductP,conductC,nummaxlayers,soilstate)
    
    !Initialize soil for tracer concentration
    IF(i_t1>0)THEN
      soilstate%conc(i_t1,1,:,:)  = genpar(m_iniT1)
      IF(nummaxlayers>1) soilstate%conc(i_t1,2,:,:) = genpar(m_iniT1)
      IF(nummaxlayers==3) soilstate%conc(i_t1,3,:,:) = genpar(m_iniT1)
    ENDIF
    IF(i_t2>0)THEN
      soilstate%conc(i_t2,1,:,:)  = genpar(m_iniT2)
      soilstate%temp(1,:,:) = genpar(m_iniT2)
      soilstate%deeptemp(:,:)    = genpar(m_iniT2)
      IF(nummaxlayers>1)THEN
        soilstate%conc(i_t2,2,:,:) = genpar(m_iniT2)
        soilstate%temp(2,:,:) = genpar(m_iniT2) 
      ENDIF
      IF(nummaxlayers==3)THEN
        soilstate%conc(i_t2,3,:,:) = genpar(m_iniT2)
        soilstate%temp(3,:,:) = genpar(m_iniT2)
      ENDIF
    ENDIF
    IF(i_t3>0)THEN
      soilstate%conc(i_t3,1,:,:)  = genpar(m_iniT3)
      IF(nummaxlayers>1) soilstate%conc(i_t3,2,:,:) = genpar(m_iniT3)
      IF(nummaxlayers==3) soilstate%conc(i_t3,3,:,:) = genpar(m_iniT3)
    ENDIF
    
    !tracers for keeping track of water origin
    IF(i_si>0)THEN 
      soilstate%conc(i_si,1,:,:)  = 1.
      IF(nummaxlayers>1) soilstate%conc(i_si,2,:,:) = 1.
      IF(nummaxlayers==3) soilstate%conc(i_si,3,:,:) = 1.
    ENDIF
    
    !Initiate river variables      
    CALL initiate_river_npc_state(conductN,conductP,conductC,conductT,conductWS,riverstate)
    
    !Initialize lake concentrations
    CALL initiate_lake_npc_state(conductN,conductP,conductC,conductT,conductWS,slc_ilake>0,slc_olake>0,lakestate)
    
    !Initiate wetland variable
    IF(wetlandexist)THEN
      !IF(conductT2) riverstate%cwetland(i_t2,:,:)=genpar(m_iniT2)
      IF(conductT) riverstate%cwetland(i_t2,:,:)=genpar(m_iniT2)   !Bug? conductT2 is never set, DG 20151125
      IF(conductWS) riverstate%cwetland(i_si,:,:)=1.
    ENDIF

    !Initiate aquifer state variables
    IF(modeloption(p_deepgroundwater)==2) CALL initiate_aquifer_state(aquiferstate)
    
  END SUBROUTINE initiate_model_state

  !>\brief Initiates model variables and parameters for a simulation. 
  !>This include 
  !>calculating HYPE variables for holding model set-up information and more...
  !--------------------------------------------------------------------------
  SUBROUTINE initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE GLACIER_SOILMODEL, ONLY :   initiate_glacier
    USE SOIL_PROCESSES, ONLY : initiate_soil_water
    USE SURFACEWATER_PROCESSES, ONLY : sum_upstream_area, &
                                       set_general_rating_k,  &
                                       calculate_landarea_riverlength
    USE NPC_SURFACEWATER_PROCESSES, ONLY : set_lake_slowwater_maxvolume
    USE IRRIGATION_MODULE, ONLY :   initiate_irrigation
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : initiate_regional_groundwater_flow, &
                                            initiate_aquifer_model, &
                                            calculate_delayed_water
    USE ATMOSPHERIC_PROCESSES, ONLY : calculate_class_wind_transformation_factor, &
                                      set_atmospheric_parameters_corrections
    USE HYPE_INDATA, ONLY : clust_group,catdes, & 
                            ndes_reg, idx_catdes, wght_catdes, indx_par, &
                            deallocate_regest_input_variables
    USE HYPEVARIABLES, ONLY : m_surfmcorr, m_deprlcorr
    
    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Variable declarations
    INTEGER i,j,k     !loop-variable
    INTEGER itype     !loop-variable, river type
    REAL help
    REAL vel,part,totaltime,kt  !parameters and help variables for river transport
    REAL delayedwater(naquifers)        !water in delay to reach aquifer (m3)
    REAL localrivvol(nsub)  !help variable for calculating local river initial volume
    REAL mainrivvol(nsub)   !help variable for calculating main river initial volume

    !Parameter declarations
    REAL, PARAMETER :: seconds_per_day = 86400.  

    
    !Flag for radiation and humidity calculations
    IF(modeloption(p_snowmelt).GE.2 .OR. modeloption(p_petmodel).GT.1 .OR. modeloption(p_lakeriverice).GE.1)THEN
      calcSWRAD = .TRUE.
    ELSE
      calcSWRAD = .FALSE.
    ENDIF
    IF(modeloption(p_petmodel).GT.1)THEN
      calcVAPOUR = .TRUE.
    ELSE
      calcVAPOUR = .FALSE.
    ENDIF
    IF(modeloption(p_petmodel).EQ.5)THEN
      calcWIND = .TRUE.
      IF(ALLOCATED(windi)) CALL calculate_class_wind_transformation_factor(windtrans)
    ELSE
      calcWIND = .FALSE.
    ENDIF

    !Set precipitation parameter corrections
    CALL set_atmospheric_parameters_corrections()
    
    !Existence of soillayers
    IF(MAXVAL(soilthick(3,1:nclass))>0)THEN
      nummaxlayers = 3
    ELSEIF(MAXVAL(soilthick(2,1:nclass))>0)THEN
      nummaxlayers = 2
    ELSE
      nummaxlayers = 1
    ENDIF
    
    !Set potential evaporation correction (not glacier?)
    IF(.NOT.ALLOCATED(basincevpcorr)) ALLOCATE(basincevpcorr(nsub))
    basincevpcorr = 1.
    DO i = 1,nsub
      IF(basin(i)%parregion>0)THEN
        basincevpcorr(i) = 1. + regpar(m_cevpcorr,basin(i)%parregion)
      ENDIF
      !Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        IF(indx_par(n_cevpc).gt.0) THEN
          basincevpcorr(i) = 1.
          DO k=1,ndes_reg(clust_group(i),indx_par(n_cevpc))
            basincevpcorr(i) = basincevpcorr(i) + wght_catdes(clust_group(i),indx_par(n_cevpc),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_cevpc),k))
          ENDDO
        ENDIF
      ENDIF
    ENDDO
       
    !Initiate soil water processes help parameters.
    CALL initiate_soil_water()
    
    !Initialize glacier parameters and type
    IF(glacierexist) CALL initiate_glacier(nclass,classmodel,frozenstate)

    !Initiate soil temperature parameters and variables
    avertemp =(/5.,10.,20.,30./)*timesteps_per_day    !Number of timestep over which meantemp is calculated
    IF(.NOT.ALLOCATED(soilmem)) ALLOCATE(soilmem(maxsoillayers,nclass))
    soilmem = 0.
    DO j= 1,nclass
      DO k = 1,maxsoillayers
        IF(k>1)THEN
          soilmem(k,j) = timesteps_per_day*genpar(m_surfmcorr)*landpar(m_surfmem,classdata(j)%luse)*EXP(genpar(m_deprlcorr)*landpar(m_depthrel,classdata(j)%luse)*(soildepth(k-1,j)+(soilthick(k,j) / 2.)))
        ELSE  
          soilmem(k,j) = timesteps_per_day*genpar(m_surfmcorr)*landpar(m_surfmem,classdata(j)%luse)*EXP(genpar(m_deprlcorr)*landpar(m_depthrel,classdata(j)%luse)*(soilthick(k,j) / 2.))
        ENDIF
      ENDDO
    ENDDO
    
    !Initialize lake and river
    IF(.NOT. ALLOCATED(landarea))     ALLOCATE(landarea(nsub))
    IF(.NOT. ALLOCATED(riverlength))  ALLOCATE(riverlength(2,nsub))
    IF(.NOT. ALLOCATED(deadriver))    ALLOCATE(deadriver(2,nsub))
    IF(.NOT. ALLOCATED(upstreamarea)) ALLOCATE(upstreamarea(nsub))
    CALL sum_upstream_area(nsub,upstreamarea)             !calculates the area upstream of the outlet point of catchment
    CALL calculate_landarea_riverlength(nsub,landarea,riverlength)
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadriver(1,i) = genpar(m_deadl) * (basin(i)%area/1.0E6) * riverlength(1,i)
      deadriver(2,i) = genpar(m_deadm) * (upstreamarea(i)/1.0E6) * riverlength(2,i)
    ENDDO
    IF(.NOT. ALLOCATED(deadwidth))    ALLOCATE(deadwidth(2,nsub))
    IF(.NOT. ALLOCATED(ratingk))      ALLOCATE(ratingk(2,nsub))
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadwidth(1,i) = SQRT(genpar(m_deadl) * (basin(i)%area/1.0E6)*0.1)   !(m)  !width=10*depth, width = sqrt(area/10)
      deadwidth(2,i) = SQRT(genpar(m_deadm) * (upstreamarea(i)/1.0E6)*0.1)   !(m)  !width=10*depth, width = sqrt(area/10)
    ENDDO
    CALL set_general_rating_k(2,nsub,landarea,upstreamarea,ratingk)
    
    !Calculate river transport time parameters
    IF(.NOT. ALLOCATED(transtime)) ALLOCATE(transtime(2,nsub))
    IF(.NOT. ALLOCATED(ttstep))    ALLOCATE(ttstep(2,nsub))
    IF(.NOT. ALLOCATED(ttpart))    ALLOCATE(ttpart(2,nsub))
    IF(.NOT. ALLOCATED(riverrc))   ALLOCATE(riverrc(2,nsub))
    vel = genpar(m_rivvel)                          !peak velocity of river (m/s)
    IF((.NOT.conductregest).AND.vel==0)THEN
      transtime = 0.
      ttstep    = 0
      ttpart    = 0.
      riverrc   = 1.
    ELSE
      part = genpar(m_damp)                           !part of delay from damping
      DO i = 1,nsub
        !>Replace parameter values with regional parameter estimates
        IF(conductregest)THEN
          IF(indx_par(n_rivv).gt.0) THEN
            vel=0.0
            DO k=1,ndes_reg(clust_group(i),indx_par(n_rivv))
              vel = vel + wght_catdes(clust_group(i),indx_par(n_rivv),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rivv),k))
            ENDDO
          ENDIF
          IF(indx_par(n_damp).gt.0) THEN
            part=0.0
            DO k=1,ndes_reg(clust_group(i),indx_par(n_damp))
              part = part + wght_catdes(clust_group(i),indx_par(n_damp),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_damp),k))
            ENDDO
          ENDIF
        ENDIF
        
        DO itype = 1,2          
          totaltime = riverlength(itype,i) / vel / seconds_per_timestep  !t=s/v, river length from GeoData (time step)
          transtime(itype,i) = (1. - part) * totaltime
          ttstep(itype,i) = INT(transtime(itype,i))
          ttpart(itype,i) = transtime(itype,i) - REAL(ttstep(itype,i))
          kt = part * totaltime                    !damptime, k in equation q=(1/k)*S
          IF(kt>0)THEN
            riverrc(itype,i) = 1. - kt + kt * exp(-1./kt)  !recession coefficient in equation Q=r*S
          ELSE
            riverrc(itype,i) = 1.   !Safe for all lake subbasin
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    
    !Check of dimension river queue
    help = SIZE(riverstate%qqueue,DIM=1) !number of elements (0:ml)
!    WRITE(6,*) 'Test: dim,ttstep',help,MAXVAL(ttstep)
    DO i = 1,nsub
      DO itype=1,2
        IF(help-1<ttstep(itype,i))THEN
          WRITE(6,*) 'Error: dimension queue',itype,i,ttstep(itype,i),help
        ENDIF
      ENDDO
    ENDDO
    
    !Set river bankful variables
    IF(conductN.OR.conductP) CALL set_Qvariables_for_bankfulflow(nsub,riverstate%Q365)

    !Initialize lake parameters
!    IF(conductN.OR.conductP.OR.conductC.OR.conductT)  &
!      CALL set_lake_slowwater_maxvolume(nsub,genpar(m_gldepi),basin(:)%lakedepth(2),lakedataparindex,lakedatapar(:,m_lddeeplake),slc_ilake>0,slc_olake>0)
    IF(conductN.OR.conductP.OR.conductC.OR.conductT)  &
      CALL set_lake_slowwater_maxvolume(nsub,basin(:)%lakedepth(1),basin(:)%lakedepth(2),lakedataparindex,lakedatapar(:,m_lddeeplake),slc_ilake>0,slc_olake>0)
 
    !Initiate regional groundwater flow/aquifer
    IF(modeloption(p_deepgroundwater)==1)THEN
      CALL initiate_regional_groundwater_flow(nsub,numsubstances,slc_ilake,slc_olake,slc_lriver,slc_mriver)
    ELSEIF(modeloption(p_deepgroundwater)==2)THEN
      CALL initiate_aquifer_model(nsub,numsubstances,naquifers)
    ENDIF
        
    !Allocate and initiate output variable
    IF(.NOT. ALLOCATED(accdiff)) ALLOCATE(accdiff(nsub))
    accdiff = 0.
    
    !Allocate nutrient (NP) load variables
    IF(conductload)THEN
       IF(.NOT.ALLOCATED(Latmdep))  ALLOCATE(Latmdep(nclass,2,numsubstances))   !Substances in order IN, ON, SP, PP
       IF(.NOT.ALLOCATED(Lcultiv))  ALLOCATE(Lcultiv(nclass,2,numsubstances))   !1=fertiliser, 2=plantdecay
       IF(.NOT.ALLOCATED(Lirrsoil)) ALLOCATE(Lirrsoil(nclass,numsubstances))    !irrigation on soil
       IF(.NOT.ALLOCATED(Lrurala))  ALLOCATE(Lrurala(nclass,numsubstances))     !rural a
       IF(.NOT.ALLOCATED(Lstream))  ALLOCATE(Lstream(nclass,numsubstances))     !runoff from soil to stream
       IF(.NOT.ALLOCATED(Lpathway)) ALLOCATE(Lpathway(numsubstances,13))        ! class independent, 13 is points A to M along flow pathway.
       IF(.NOT.ALLOCATED(Lbranch))  ALLOCATE(Lbranch(numsubstances))            !part of outflow
       IF(.NOT.ALLOCATED(Lgrwmr))   ALLOCATE(Lgrwmr(numsubstances))             !reg.grw to main river
       IF(.NOT.ALLOCATED(Lgrwol))   ALLOCATE(Lgrwol(numsubstances))             !reg.grw to olake
    ENDIF 
    IF(numsubstances>0)THEN
      IF(.NOT.ALLOCATED(Lruralb))  ALLOCATE(Lruralb(numsubstances))            !rural b
      IF(.NOT.ALLOCATED(Lpoints))  ALLOCATE(Lpoints(numsubstances,max_pstype)) !point sources 1=ps1, 2=ps2, 3=ps3, ..
      IF(.NOT.ALLOCATED(Lgrwsoil)) ALLOCATE(Lgrwsoil(nclass,numsubstances))       !regional groundwaterflow to soil
      IF(.NOT.ALLOCATED(Lgrwclass)) ALLOCATE(Lgrwclass(nclass,numsubstances,nsub))       !regional groundwater outflow from soil
    ENDIF
    
    !Initiate irrigation
    IF(doirrigation) CALL initiate_irrigation(nsub,nclass)

    !Initate calculations for water balance output
    IF(conductwb)THEN
      CALL initiate_waterbalance_output(nsub)
      !set initial wb_stores
      wbstores = 0.
      DO i = 1,nsub
        DO j = 1,nclass
          wbstores(w_snow,i) = wbstores(w_snow,i) + frozenstate%snow(j,i)*classbasin(i,j)%part  !only land?
          IF(classmodel(j)==glacier_model) wbstores(w_glacier,i) = frozenstate%glacvol(i)*genpar(m_glacdens)
          IF(classmodel(j)==0.OR.classmodel(j)==glacier_model)THEN
            wbstores(w_soil1,i) = wbstores(w_soil1,i) + soilstate%water(1,j,i)*classbasin(i,j)%part
            wbstores(w_soil2,i) = wbstores(w_soil2,i) + soilstate%water(2,j,i)*classbasin(i,j)%part
            wbstores(w_soil3,i) = wbstores(w_soil3,i) + soilstate%water(3,j,i)*classbasin(i,j)%part
          ENDIF
          IF(doirrigation) wbstores(w_irrcanal,i) = wbstores(w_irrcanal,i) + miscstate%nextirrigation(j,i)*classbasin(i,j)%part
        ENDDO
        localrivvol(i) = riverstate%water(1,i) + (SUM(riverstate%qqueue(1:ttstep(1,i),1,i)) + riverstate%qqueue(ttstep(1,i)+1,1,i) * ttpart(1,i))
        mainrivvol(i)  = riverstate%water(2,i) + (SUM(riverstate%qqueue(1:ttstep(2,i),2,i)) + riverstate%qqueue(ttstep(2,i)+1,2,i) * ttpart(2,i))
      ENDDO
      wbstores(w_snow,:)   = wbstores(w_snow,:)  * basin(:)%area *1.E-3  !m3
      wbstores(w_soil1,:)  = wbstores(w_soil1,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_soil2,:)  = wbstores(w_soil2,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_soil3,:)  = wbstores(w_soil3,:) * basin(:)%area *1.E-3  !m3
      IF(doirrigation) wbstores(w_irrcanal,:) = wbstores(w_irrcanal,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_iriver,:) = localrivvol(:) !m3
      wbstores(w_mriver,:) = mainrivvol(:)  !m3
      IF(slc_ilake>0) wbstores(w_ilake,:) = lakestate%water(1,:)  !mm
      IF(slc_olake>0)wbstores(w_olake,:)  = lakestate%water(2,:)  !mm
      IF(numsubstances>0)THEN
        IF(slc_ilake>0) wbstores(w_ilake,:) = wbstores(w_ilake,:) + lakestate%slowwater(1,:)  !mm
        IF(slc_olake>0) wbstores(w_olake,:) = wbstores(w_olake,:) + lakestate%slowwater(2,:)  !mm
      ENDIF
      IF(slc_ilake>0) wbstores(w_ilake,:) = wbstores(w_ilake,:) * basin(:)%area *classbasin(:,slc_ilake)%part *1.E-3  !m3
      IF(slc_olake>0) wbstores(w_olake,:) = wbstores(w_olake,:) * basin(:)%area *classbasin(:,slc_olake)%part *1.E-3  !m3
      CALL calculate_delayed_water(aquiferstate,naquifers,delayedwater)
      IF(naquifers>0) wbstores(w_aquifer,1:naquifers) = aquiferstate%water + aquiferstate%nextoutflow + delayedwater
      CALL print_initial_waterbalance_stores(nsub,naquifers)
    ENDIF
    
    IF(conductregest) CALL deallocate_regest_input_variables()

  END SUBROUTINE initiate_model
  
  !>Model subroutine for HYPE.
  !!Calculates what happen during one timestep.
  !---------------------------------------------------
  SUBROUTINE model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,d,idt)
    
    USE STATETYPE_MODULE
    USE NPC_SURFACEWATER_PROCESSES, ONLY : add_dry_deposition_to_river,   &
                                           add_dry_deposition_to_lake,   &
                                           np_processes_in_river,  &
                                           np_processes_in_lake,   &
                                           oc_processes_in_river,  &
                                           oc_processes_in_lake,   &
                                           add_diffuse_source_to_local_river,  &
                                           add_point_sources_to_main_river,    &
                                           calculate_river_wetland
    USE SURFACEWATER_PROCESSES, ONLY : add_precipitation_to_river, &
                                       calculate_river_evaporation, &
                                       calculate_actual_lake_evaporation, &
                                       calculate_water_temperature,     &
                                       set_water_temperature, &
                                       calculate_river_characteristics, &
                                       translation_in_river,            &
                                       point_abstraction_from_main_river, &
                                       point_abstraction_from_outlet_lake, &
                                       calculate_outflow_from_lake,     &
                                       check_outflow_from_lake,         &
                                       remove_outflow_from_lake,        &
                                       calculate_flow_within_lake,      &
                                       calculate_olake_waterstage,      &
                                       calculate_regamp_adjusted_waterstage, &
!                                       calculate_olake_waterlevel_for_outflow,      &
!                                       calculate_cleaned_olake_waterstage,      &
                                       calculate_branched_flow,&
                                       calculate_lake_volume, &
                                       T2_processes_in_river, &
                                       T2_processes_in_lake, &
                                       ice_processes_in_river, &
                                       ice_processes_in_lake, &
                                       add_T2_concentration_in_precipitation_on_water,  &
                                       get_rivertempvol
    USE NPC_SOIL_PROCESSES, ONLY : set_class_precipitation_concentration_and_load,  &
                                   croprotation_soilpoolaverage
    USE IRRIGATION_MODULE, ONLY : initiate_timestep_irrigation,  &
                                  calculate_irrigation
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : calculate_regional_groundwater_flow,  &
                                            calculate_current_flow_from_aquifer, &
                                            add_regional_groundwater_flow_to_olake, &
                                            calculate_river_groundwaterflow_removal, &
                                            calculate_aquifer,  &
                                            add_aquifer_flow_to_river, &
                                            calculate_delayed_water, &
                                            calculate_aquifer_waterlevel
    USE SOILMODEL_DEFAULT, ONLY : soilmodel_0
    USE GLACIER_SOILMODEL, ONLY : soilmodel_3, calculate_glacier_massbalance
    USE ATMOSPHERIC_PROCESSES, ONLY : calculate_class_atmospheric_forcing,  &
                                calculate_subbasin_temperature, &
                                calculate_subbasin_precipitation, &
                                calculate_rain_snow_from_precipitation, & 
                                calculate_extraterrestrial_radiation, &
                                set_precipitation_concentration,  &
                                calculate_daylength
    USE SOIL_PROCESSES, ONLY : calculate_potential_evaporation
    USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water,         &
                                            error_remove_water,   &
                                            add_water
    USE HYPE_INDATA, ONLY : get_current_xoms, &
                            num_xoms, &
                            xom, xos
    USE LIBDATE, ONLY : datetype !allows assigning date to special dams
    USE MODVAR, ONLY : dam, damindex !MH2017: allows assigning date, yesterday's flow, to special dams
    
    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    TYPE(datetype),INTENT(IN)            :: d           !<MH2017: date, used for special dams
    INTEGER,INTENT(IN)                   :: idt

    !Variable declarations
    INTEGER i,j,k    !loop-variables: subbasin, class, div.
    INTEGER itype    !loop-variable: ilake=1,olake=2
    INTEGER numrotations   !Number of crop rotation class groups
    INTEGER counter !loop-variable: updates dam(damindex(i))%qin14day
    
    REAL a        !area fraction of class
    REAL asum, asumseclayer, asumthdlayer    !sum of land area parts
    REAL divasum, divasum2, divasum3         !inverse of area sums
    REAL divasum4, divasum5                  !subbasin area - upstream area relation
    REAL incorr, oncorr    !correction of inorganic and organic nitrogen level
    REAL phoscorr     !correction of phosphorus level
    
    !Help variables for class summation to subbasin output
    REAL gwati,eacti,epoti,melti,runoffi,snowi,soili,landeacti,evapsnowi
    REAL runoff1i,runoff2i,runoff3i, runoffdi,runoffsri
    REAL plantuptakei, nitrifi, denitrifi, denitrif12i, denitrif3i
    REAL atmdepi, atmdep2i                !atmospheric deposition of TN and TP on soil (kg/km2)
    REAL crunoffi(numsubstances),csoili(numsubstances)         !(mg/L)
    REAL soilfrosti,soiltempi,snowdepthi
    REAL soilltmpi(maxsoillayers)
    REAL cevapi,icpevapi  
    REAL csoili1(numsubstances),csoili2(numsubstances),csoili3(numsubstances)         !(mg/L)
    REAL soili1,soili2,soili3,soili12,standsoili
    REAL srffi,smfdi,smfpi,srfdi,srfpi
    REAL fastNp1i,fastNp2i,fastNp3i,humusNp1i,humusNp2i,humusNp3i             !pools of nitrogen
    REAL partPp1i, partPp2i, partPp3i, fastPp1i, fastPp2i, fastPp3i, humusPp1i, humusPp2i, humusPp3i
    REAL humusCp1i,humusCp2i,humusCp3i,soilorgCp1i,soilorgCp2i,soilorgCp3i    !pools of organic carbon
    REAL fastCp1i,fastCp2i,fastCp3i
    REAL corrpreci,corrtempi, qcinfli 
    REAL snowfalli, rainfalli  !Precipiation as snow, rain (mm)
    REAL runoffrzi,crunoffrzi(numsubstances),runoffnorzi,crunoffnorzi(numsubstances)
    REAL soildefi
    REAL snowcovi,snowmaxi
!    REAL pcorricep(nsub)  DG20151125: stop treating interception loss as precipitation correction -> it's evaporation
    REAL soillgrossload(3,numsubstances), soillnetload(3,numsubstances)   !1=sl 1-2, 2=sl 3, 3=sl3+tile
!Additional in CryoProcesses branch/DG 20151125:
    REAL deeptempi,soil2tempi
    REAL landevapsnowi,landevapglaci,evapglaci,evapsoili,evaplakei,evapriveri
    REAL glacmelti
    REAL swradi,netradi
    REAL pcorrbasin(nsub)   !precipitation corrections to subbasin scale 
    REAL snowfallthroughi, rainfallthroughi  !Precipiation as snow, rain (mm)
    REAL pcorr,pcorri       !precipitation corrections to classes and summation over subasin
        
    !Help variables for summation of upstream output variables
    REAL upstream_prec(nsub),upstream_evap(nsub),upstream_epot(nsub)
    REAL upstream_snow(nsub),upstream_soil(nsub),upstream_soildef(nsub),upstream_smfp(nsub)
    REAL upstream_snowfall(nsub),upstream_rainfall(nsub)
!Additional in CryoProcesses branch/DG 20151125:
    REAL upstream_snowcov(nsub)
    REAL upstream_river(nsub),upstream_lake(nsub),upstream_glacier(nsub)
    REAL upstream_recsnow(nsub),upstream_compsnow(nsub)
    !REAL upstream_recsnowarea(nsub),upstream_recsnowmissing(nsub)
    !REAL upstream_compfsc(nsub),upstream_recfsc(nsub)
    !REAL upstream_recfscarea(nsub),upstream_recfscerror(nsub),upstream_recfscmissing(nsub)
    REAL upstream_temp(nsub)
    REAL upstream_evapsnow(nsub),upstream_evapsoil(nsub),upstream_evapglac(nsub),upstream_evapriver(nsub),upstream_evaplake(nsub)
    REAL upstream_snowmelt(nsub), upstream_glacmelt(nsub)
    REAL upstream_runoff(nsub)
    REAL upstream_swrad(nsub), upstream_netrad(nsub)
    
    !Current forcing variables for estimating radiation for snow and evaporation calculations
    REAL radexti(nsub)           !Extraterrestrial solar radiation [MJ/m2/day] for current time step
    REAL tmin,tmax,rhmin,netrad,actvap,satvap,wind
    REAL sffrac    !sffrac = snow fall fraction of precipitation [-]
    REAL swrad     !shortwave radiation downward [MJ/m2/day]
    REAL daylength  !day length (hours)
    
    !Variables for class values
    REAL classarea            !class area (km2)
    REAL classheight          !class height (möh)
    REAL lakeareakm2
    REAL gwat,melt,epot,evap,evapsnow,epotsnow   !calculated class values
    REAL surfaceflow,icpevap
    REAL nitrif, denitrif(maxsoillayers)
    REAL cropuptake          !calculate class values, uptake of IN
    REAL runoffd,crunoffd (numsubstances)         !calculated class values
    REAL crunoff1(numsubstances),crunoff2(numsubstances), crunoff3(numsubstances)   !<calculated class values (mg/L)
    REAL csrunoff(numsubstances)
    REAL cevap(numsubstances),cevapsnow(numsubstances)   !calculated class values (mg/L)
    REAL cprec(numsubstances),cprecj(numsubstances)   !concentration of precipitation, subbasin and class
    REAL evapl,cevapl(numsubstances)  !evaporation lake
    REAL concarr(numsubstances)   !help variable for concentration
    REAL temp                 !class elevation corrected temperature
    REAL prec                 !class elevation corrected precipitation
    REAL smdef    !soil moisture deficit
    REAL frostdepth,snowdepth                    !soil frost variables 
    REAL glac_part        !glacier fraction of glacier class (-)
    REAL snowfall, rainfall  !Precipiation as snow, rain (mm)
    REAL snowplain_area_fraction, old_snowfield_part
    REAL atmdepload1(numsubstances),atmdepload2(numsubstances)  !wet and drydep
    REAL irrload(numsubstances)
    REAL grwloadlake(numsubstances),grwloadmr(numsubstances)
    REAL rgrwload(numsubstances)
    REAL ruralaload(numsubstances)
    REAL cultivload(2,numsubstances)
    REAL infiltrationflows(7) !fraction from snow melt,infiltration,surface runoff,macroporeflow to soil layer 1,2,and 3,fraction from glacier melt
    REAL glacierflows(2)      !flow snow to glacier, precipitation on glacier
    REAL evapflows(4)         !flow evapotranspiration (land classes), 1=sl1, 2=sl2,3=snow,4=glacier
    REAL runofflows(7)        !different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3,7=saturated surface runoff
    REAL crunofflows(numsubstances,6) !concentration of different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3
    REAL verticalflows(6)     !vertical flows:1-2=percolation,3-4=upwelling due to rural,5-6=upwelling due to reg. grw flows
    REAL cverticalflows(2,numsubstances) !<concentration of vertical flows:1-2=percolation
    REAL horizontalflows(3)   !horizontal flows:1-3=recieved rural load flow
    REAL horizontalflows2(3,nsub)   !horizontal flows:1-3=division of regional groundwater flows to grwdown
    REAL cruralflow(numsubstances)   !concentration of rural load flow
    REAL nutrientloadflows(4) !1=rural load flow to local stream,2=point source to main river,3=abstraction from main river,4=abstraction from outlet lake
    REAL irrigationflows(7)   !withdrawal from 1=deep groundwater(unlimited),2=internal lake,3=outlet lake,4=main river volume and inflow,5=local evaporation losses (local sources),6=unlimited external source (alternative),7=akvifer(modelled)
    REAL regionalirrflows(4,nsub) !withdrawal from 1=outlet lake,2=main river volume and inflow,3=evaporation losses olake as regional source,4=evaporation losses main river as regional source
    REAL regionalirrevap(nsub)    !local evaporation losses (regional source)
!Additional for CryoProcesses branch /DG20151125
    REAL glacmelt, glacevap                    !glacier melt and sublimation
!    REAL csnowfall(numsubstances),crainfall(numsubstances) 
!    REAL cicpevap(numsubstances) !concentration in interception evaporation
    REAL ctemp_sm,ctemp_gm,ctemp_rn,ctemp_li,ctemp_ri,ctemp_si
    REAL lakeareai(2)                          !lake area (m2) separately for ilake and olake
    REAL snowfallthrough,rainfallthrough       !throughfall (snow and rain)
    REAL csnowfallthrough(numsubstances),crainfallthrough(numsubstances) !concentration in throughfall (snowfall and rainfall), subbasin and class
    REAL cicpevap(numsubstances) !interception evaporation, water and it's concentration of substances
    
    !Variables for routing
    REAL accinflow(nsub)                  !accumulated upstream inflow to subbasin (m3/s)
    REAL acccinflow(numsubstances,nsub)   !concentration of upstream inflow (mg/L)
    REAL mainflow,branchflow              !divided discharge for output (m3/s)
    REAL qinmm,qin,cin(numsubstances)
    REAL inflowpart !part of runoff to ilake
    REAL qupstream, cupstream(numsubstances),qunitfactor,outflowm3s,outflowmm
    REAL transq,transc(numsubstances),dampq,dampc(numsubstances)
    REAL lakeoutflow(2)                   !Outflow of lakes (m3/s)
    REAL,ALLOCATABLE :: clakeoutflow(:,:)    !Concentration of outflow of lakes (mg/L)
    REAL outflowtemp,wstlaketemp                        !variables for lakeoutflow/wstlake before update is done
    REAL olakewstold                         !wst olake last time step (mm)
    REAL lakewst(2)   !lake water stage for ilake and olake (mm)
    REAL wstlake,wstlakeadj      !temporary variables for olake water stage output
    REAL w0ref        !olake reference level for water stage
    REAL lakebasinvol(2)   !volume of ilake & olake (m3)
    REAL lakevol           !volume for single olake and lakes composed of basins(sum)  (m3)
    REAL lakevolsum(nbasinlakes)    !accumulation of lake basin volumes to outlet (m3)
    REAL grwout(4,nsub),grwout2             !Outflow from groundwater for print out (all, and resp soillayer)(main river)
    REAL grwtomr,grwtool                    !Flow of regional groundwater to main river and olake
    REAL riverQbank,riverarea(2),riverdepth,riverareakm2
    REAL evapr,cevapr(numsubstances)
    REAL lakearea             !lake area (m2) for lake
    REAL lakeareatemp,lakeareatemp2         !lake area (m2) for outlet lake or whole lake for last lakebasin
    REAL flow1000m3ts         !olake outflow in 1000 m3/ts
    
    !Variables for irrigation 
    REAL pwneedi                      !irrigation water demand for all irrigated classes (m3)
    REAL gwremi                       !groudwater removed for irrigation (mm)
    REAL ldremi,lrremi,rsremi         !abstraction surface water for irrigation (m3)
    REAL irrevap(nsub)                !accumulation variable for irrigation field and network losses (m3)
    REAL irrappl(2),irrappli          !irrappl = applied irrigation (mm), for summation basin output, irrapli is the sum of these
    REAL irrsinkmass(numsubstances)   !irrsinkmass = mass of substances in the irrigation sink(kg for mg/L concentrations)
    
    !Variables for the lake and river ice model
    REAL lakesnowdepth(2),riversnowdepth(2),lakesurftemp(2), riversurftemp(2)
    REAL ctemp_T2
    REAL meanrivertemp,totrivervol    !mean temperature and total river volume over all water course elements
    REAL freezeuparea !freezeuparea for river and lakes
    INTEGER lakefreezeupday(2), lakebreakupday(2), riverfreezeupday(2), riverbreakupday(2)
    
    REAL delayedwater(naquifers)        !water in delay to reach aquifer (m3) 
    REAL aquiferoutflow(naquifers)
    REAL aqremflow(nsub)            !water removed from soillayers for aquifer
    REAL aqoutside(naquifers)       !water from aquifer to outside model system
    REAL aqirrloss(naquifers)       !water for irrigation from aquifer
    
    !Variables for error handling
    INTEGER status
    CHARACTER(LEN=80) :: errstring(1)  !error message for location of remove_water call
    PARAMETER (errstring = (/'outflow from river damping box'/))
    
    !Variables for comparing forest and open snow cover variables with the FSUS data (Former Soviet Union Snow coarse data)
    REAL snowdepth_forest,snowdepth_open,snowcover_forest,snowcover_open,snowdensity_forest,snowdensity_open,swe_forest,swe_open
    REAL area_forest, area_open
    
    !MM2016 temp
    REAL outflowshit, adjnoncontr(nsub)
   

    !Start of subroutine calculations
    !> \b Algorithm \n
    IF(conductwb)  wbflows=0.  !zeroing the timestep for all flows
    
    !>Initial calculation of regional groundwater
    IF(numsubstances>0) Lgrwclass=0.
    grwout = 0.
    IF(modeloption(p_deepgroundwater)==1 .OR. &
       modeloption(p_deepgroundwater)==2)THEN
      CALL calculate_regional_groundwater_flow(soilstate,grwout,Lgrwclass)    !regional groundwater leaving the soil
      outvar(:,240) = grwout(1,:)*1.E-3   !will be replaced later for aquifermodel
    ENDIF
    IF(conductwb)THEN
      wbflows(w_rgrwof1,:) = grwout(2,:)
      wbflows(w_rgrwof2,:) = grwout(3,:)
      wbflows(w_rgrwof3,:) = grwout(4,:)
    ENDIF
    IF(modeloption(p_deepgroundwater)==2)THEN     !prepare for adding regional groundwater to model
      CALL calculate_current_flow_from_aquifer(naquifers,numsubstances,aquiferstate,aqoutside)    
      IF(conductwb) wbflows(w_rgrwtoos,1:naquifers) = aqoutside
    ENDIF
        
    !>Get current observations for HYPE specific indata
    IF(conductxoms) CALL get_current_xoms(currentdate,nsub)
    
    !>Calculate corrected atmospheric forcing for all subbasins
    outvar(:,4) = tempi       !original input data
    outvar(:,3) = preci       !original input data
    CALL calculate_subbasin_temperature(nsub,month,tempi)

!DG20151125: stop treating interception losses as precipitation correction!
!    CALL calculate_subbasin_precipitation(nsub,outvar(:,4),preci,pcorricep)
!    outvar(:,206) = pcorricep
    CALL calculate_subbasin_precipitation(nsub,outvar(:,4),preci,pcorrbasin)
    outvar(:,max_outvar_trunk+1) = pcorrbasin
    
    IF(ALLOCATED(tmini))THEN
      CALL calculate_subbasin_temperature(nsub,month,tmini)
      CALL calculate_subbasin_temperature(nsub,month,tmaxi)
    ENDIF
    
    !Calculate extra terrestrial radiation for all subbasins, if needed
    IF(calcSWRAD)THEN
      CALL calculate_extraterrestrial_radiation(nsub,dayno,radexti)
    ENDIF
    
    !Initiations for main subbasin calculation loop
    accinflow = 0.            
    acccinflow = 0.
    irrevap = 0. 
    aqirrloss = 0.
    lakevolsum = 0.
    IF(conductwb) wbirrflows = 0.
    regionalirrflows=0.
    horizontalflows2=0.
    CALL initiate_timestep_irrigation(numsubstances,miscstate)
    IF(.NOT.ALLOCATED(clakeoutflow)) ALLOCATE(clakeoutflow(numsubstances,2))
    upstream_prec = 0.; upstream_evap = 0.; upstream_snow = 0.; upstream_soil = 0.
    upstream_soildef = 0.; upstream_smfp = 0.
    upstream_snowfall = 0.; upstream_rainfall = 0.; upstream_epot = 0.
!CryoProcesses branch /DG20151125
    upstream_river = 0.; upstream_lake = 0.; upstream_glacier = 0.
    upstream_recsnow = 0. ; upstream_compsnow = 0.
   ! upstream_recsnowarea = 0. ;upstream_recsnowmissing = 0.
   ! upstream_compfsc = 0. ;upstream_recfsc = 0.
   ! upstream_recfscarea = 0. ;upstream_recfscerror = 0. ;upstream_recfscmissing = 0.
    upstream_temp = 0.
    upstream_evapsnow = 0. ; upstream_evapsoil = 0.; upstream_evapriver = 0.; upstream_evaplake = 0.; upstream_evapglac = 0.
    upstream_runoff=0.;upstream_snowcov=0.;
    upstream_swrad=0 ; upstream_netrad=0.
    
    !MM2016: non-contributing area
    adjnoncontr(:) = basin(:)%noncontr

    !>Main subbasin-loop, subbasins calculated in flow order
    subbasinloop:  &
    DO i = 1,nsub
       
      !Initiate variables for calculation of subbasin averages
      gwati=0.; eacti=0.; cevapi=0.; epoti=0.; melti=0.; snowi=0.; soili=0.
      runoffi=0.; runoffdi=0.; runoffsri=0.
      runoff1i=0.; runoff2i=0.; runoff3i=0.
      runoffrzi=0.;crunoffrzi=0.;runoffnorzi=0.;crunoffnorzi=0.
      denitrifi=0.; nitrifi=0.; plantuptakei=0.; atmdepi=0.; atmdep2i=0.
      denitrif12i=0.; denitrif3i=0.
      soilfrosti=0.; soiltempi=0.; snowdepthi=0.; soilltmpi=0.
      csoili=0.; crunoffi=0. 
      soili1=0.; soili2=0.; soili3=0.; soili12=0.; standsoili=0.
      csoili1=0.; csoili2=0.; csoili3=0.
      soildefi=0.; icpevapi=0.
      corrpreci=0.; corrtempi=0.
      snowfalli=0.;rainfalli=0.
      fastNp1i=0.;fastNp2i=0.;fastNp3i=0.
      humusNp1i=0.;humusNp2i=0.;humusNp3i=0.
      partPp1i=0.;partPp2i=0.;partPp3i=0. 
      fastPp1i=0.;fastPp2i=0.;fastPp3i=0.
      humusPp1i=0.;humusPp2i=0.;humusPp3i=0.
      humusCp1i=0.;humusCp2i=0.;humusCp3i=0.
      soilorgCp1i=0.;soilorgCp2i=0.;soilorgCp3i=0.
      fastCp1i=0.;fastCp2i=0.;fastCp3i=0.
      qcinfli=0.
      irrappli=0.;pwneedi=0.
      ldremi=0.;lrremi=0.;gwremi=0.;rsremi=0.
      snowcovi=0.;snowmaxi=0.
      srffi=0.;smfdi=0.;smfpi=0.;srfdi=0.;srfpi=0.
      asum=0.;asumseclayer=0.;asumthdlayer=0.
      glacierflows=0.;nutrientloadflows=0.;irrigationflows=0.
      evapsnowi=0.
      !snow variables to compare with FSUS data on forest and open land
      snowdepth_forest = 0.;snowdepth_open = 0.;snowcover_forest = 0.;snowcover_open = 0.
      snowdensity_forest = 0.;snowdensity_open = 0.;swe_forest = 0.;swe_open = 0.
      area_forest = 0.;area_open = 0.
      soillgrossload =0.;soillnetload=0.
      !additional for CryoProcesses
      deeptempi=0.; soil2tempi=0.
      glacmelti=0.;evapglaci=0.;evapsoili=0.;evaplakei=0.;evapriveri=0.
      netradi = 0. ; swradi = 0.
      snowfallthroughi=0.;rainfallthroughi=0.
      pcorri=0.
      
      !Initiate other subbasin variables
      IF(conductload)THEN
        Latmdep=0.
        Lcultiv=0.
        Lirrsoil=0.
        Lrurala=0.
        Lgrwsoil=0.
        Lgrwol=0.
        Lgrwmr=0.
        IF(numsubstances>0) Lstream=0.
      ENDIF
       
      !Short notation for parameters not dependent on class
      IF(basin(i)%wqparregion>0)THEN
        incorr = (1. + wqregpar(m_incorr,basin(i)%wqparregion))     !Correction of inorgnic nitrogen
        oncorr = (1. + wqregpar(m_oncorr,basin(i)%wqparregion))     !Correction of organic nitrogen
        phoscorr = 1. + wqregpar(m_phoscorr,basin(i)%wqparregion)   !Correction of phosphorus level
      ELSE
        incorr    = 1.
        oncorr    = 1.
        phoscorr  = 1.
      ENDIF
       
      !Calculation of mean air temperature
      miscstate%temp5(i) = miscstate%temp5(i) + (tempi(i) - miscstate%temp5(i)) / avertemp(1)
      miscstate%temp10(i) = miscstate%temp10(i) + (tempi(i) - miscstate%temp10(i)) / avertemp(2)
      miscstate%temp20(i) = miscstate%temp20(i) + (tempi(i) - miscstate%temp20(i)) / avertemp(3)
      miscstate%temp30(i) = miscstate%temp30(i) + (tempi(i) - miscstate%temp30(i)) / avertemp(4)

      !Set subbasin precipitation concentrations
      CALL set_precipitation_concentration(i,numsubstances,cprec)
       
      !Calculate daylength for temperature dependent growing season
      daylength = 0.
      IF(modeloption(p_growthstart)==1)  &
        CALL calculate_daylength(dayno,basin(i)%latitude,daylength)

      !>Main class-loop for calculation of soil water and substances (land classes)
      DO j=1,nclass
        a=classbasin(i,j)%part
        classarea = a * basin(i)%area * 1.0E-6    !km2
        classheight = basin(i)%elev+classbasin(i,j)%deltah  !möh
        IF(a>0)THEN
          irrappl = 0.
          cultivload = 0.
          irrload = 0.
          ruralaload = 0.
          atmdepload1 = 0.
          atmdepload2 = 0.
          rgrwload = 0.
          old_snowfield_part = 1.
          !> \li Calculate class forcing data
!DG20151125: stop treating interception evaporation as precipitation correction!
!          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
!                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,netrad,wind,sffrac, &
                  snowfall,rainfall,pcorr)
               
          CALL set_class_precipitation_concentration_and_load(numsubstances,classarea,outvar(i,3),temp,prec,cprec,cprecj,atmdepload1)
  
          !> \li Calculate soil processes            
          IF(classmodel(j)==0)THEN
            CALL soilmodel_0(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,classheight,prec,cprecj,temp,  & 
                 daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,  &
                 frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,    &
                 cropuptake,nitrif,denitrif,melt,epot,gwat,frostdepth,snowdepth,irrappl,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                 pwneedi,snowfall,rainfall,cultivload,irrload,ruralaload,rgrwload,  &
                 atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                 infiltrationflows,evapflows,runofflows,crunofflows,verticalflows,  &
                 cverticalflows,horizontalflows,horizontalflows2,evapsnow,cevapsnow,cruralflow, &
                 snowfallthrough,rainfallthrough,csnowfallthrough,crainfallthrough,icpevap,cicpevap, &
                 tmin,tmax)
            snowplain_area_fraction = 1.
          ELSEIF(j==slc_lriver)THEN
            CYCLE
          ELSEIF(j==slc_mriver)THEN
            CYCLE
          ELSEIF(classmodel(j)==glacier_model)THEN
            !only one glacier class is allowed
            CALL soilmodel_3(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,prec,cprecj,temp, & 
                  daylength,basincevpcorr(i),basinrrcscorr(i),phoscorr, &
                  frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,  &
                  cropuptake,nitrif,denitrif,melt,epot,gwat,frostdepth,snowdepth,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                  glac_part,old_snowfield_part,snowfall,rainfall,cultivload,ruralaload, &
                  rgrwload,atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                  infiltrationflows,glacierflows,evapflows,runofflows,crunofflows,  &
                  verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cevapsnow,cruralflow, &
                  snowfallthrough,rainfallthrough,csnowfallthrough,crainfallthrough,icpevap,cicpevap,glacmelt)
            snowplain_area_fraction = 1. - glac_part
            outvar(i,117) = frozenstate%glacvol(i)*1.E-9    !glacier ice volume (km3)
            outvar(i,118) = glac_part * classarea           !glacier area (km2)
            IF(xobsindex(256,i).GT.0) outvar(i,256) = xobsi(xobsindex(256,i)) !Glacier mass balance period
            IF(outvar(i,256) /= missing_value)THEN   !If there is a Mass balance period to evaluate:
              IF(xobsindex(253,i).GT.0) outvar(i,253) = xobsi(xobsindex(253,i)) !recorded mass balance
              IF(xobsindex(255,i).GT.0) outvar(i,255) = xobsi(xobsindex(255,i)) !recorded mass balance area
              CALL calculate_glacier_massbalance(i,outvar(i,256),outvar(i,252),outvar(i,254))   !calculate simulated glacier mass balance and area
            ENDIF
          ELSEIF(j==slc_olake)THEN
            CYCLE
          ELSEIF(j==slc_ilake)THEN
            CYCLE
          ELSE
            CALL soilmodel_0(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,classheight,prec,cprecj,temp,  & 
                 daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,  &
                 frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,    &
                 cropuptake,nitrif,denitrif,melt,epot,gwat,frostdepth,snowdepth,irrappl,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                 pwneedi,snowfall,rainfall,cultivload,irrload,ruralaload,rgrwload,  &
                 atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                 infiltrationflows,evapflows,runofflows,crunofflows,verticalflows,  &
                 cverticalflows,horizontalflows,horizontalflows2,evapsnow,cevapsnow,cruralflow, & 
                 snowfallthrough,rainfallthrough,csnowfallthrough,crainfallthrough,icpevap,cicpevap, &
                 tmin,tmax) 
            snowplain_area_fraction = 1.
          ENDIF
             
          !Set load variables for printout
          IF(conductload)THEN 
            Lcultiv(j,:,:) = cultivload
            Lirrsoil(j,:)  = irrload
            Lrurala(j,:)   = ruralaload
            Latmdep(j,1,:) = atmdepload1
            Latmdep(j,2,:) = atmdepload2
            Lgrwsoil(j,:)  = rgrwload
          ENDIF
             
          !> \li Accumulate variables for mean over subbasin and layers
          ! TODO: gör till fler arrayer?
          asum = asum + a                                         !landarea fraction
          corrpreci = corrpreci + prec*a
          corrtempi = corrtempi + temp*a
          rainfalli = rainfalli + rainfall * a 
          snowfalli = snowfalli + snowfall * a
          !new variables for throughfall (precipitation-interception)
          rainfallthroughi = rainfallthroughi + rainfallthrough * a 
          snowfallthroughi = snowfallthroughi + snowfallthrough * a
          !and some further for CryoProcesses
          swradi = swradi + swrad*a
          netradi= netradi+ netrad*a
          pcorri = pcorri + pcorr*a
          !---------------------------------------------------------
          snowi = snowi + frozenstate%snow(j,i)*a*snowplain_area_fraction  !zero snow on glacier area
          snowdepthi = snowdepthi + snowdepth*a*snowplain_area_fraction    !-"-
          melti = melti+melt*a                                             !(melt rescaled in glacier_soilmodel)
          snowmaxi = snowmaxi+frozenstate%snowmax(j,i)*a*snowplain_area_fraction   !zero snow on glacier area
          soili1 = soili1+soilstate%water(1,j,i)*a
          csoili1(:) = csoili1(:) + soilstate%conc(:,1,j,i)*soilstate%water(1,j,i)*a     !accumulation of amount
          runoff1i = runoff1i + runofflows(1)*a
          runoffdi = runoffdi + runoffd*a
          runoffsri = runoffsri + surfaceflow*a
          runoffnorzi = runoffnorzi + runoffd*a + surfaceflow*a
          crunoffnorzi(:) = crunoffnorzi(:) + crunoffd(:)*runoffd*a + csrunoff(:)*surfaceflow*a
          plantuptakei = plantuptakei + cropuptake*a
          nitrifi = nitrifi + nitrif*a        
          denitrifi = denitrifi + SUM(denitrif(:))*a
          denitrif12i = denitrif12i + denitrif(1)*a + denitrif(2)*a
          denitrif3i = denitrif3i + denitrif(3)*a  
          IF(i_in>0) atmdepi = atmdepi + (atmdepload1(i_in)+atmdepload2(i_in))*a/classarea
          IF(i_sp>0) atmdep2i = atmdep2i + (atmdepload1(i_sp)+atmdepload2(i_sp))*a/classarea
          epoti = epoti+epot*a                                   !total potential evap (dry evap,  excluding interception, n.b.!)
          eacti = eacti+evap*a                                   !total actual evapotranspiration, now including evp from interception, n.b.!
          !evapsnowi = evapsnowi +evapsnow*a                     !sublimation from snow and glaciers
          
          !Alternative evaporation components for CryoProcesses
          evapsoili = evapsoili + (evap-evapsnow-icpevap)*a     !evaporation from soillayer 1+2
          IF(classmodel(j)==glacier_model)THEN
            glacmelti = glacmelti+glacmelt*a                                 !glaciermelt
            glacevap  = evapflows(4)*(1.-snowplain_area_fraction)            !glacier sublimation
            evapglaci = evapglaci + glacevap*a                               !sum glacier sublimation
            evapsnowi = evapsnowi +(evapsnow-glacevap)*a                     !sublimation from snow (remove glacier subl)
          ELSE
            evapsnowi = evapsnowi +evapsnow*a                                !sublimation from snow
          ENDIF
               
          IF(i_t1>0) cevapi = cevapi+cevap(i_t1)*evap*a   !accumulation of amount
          icpevapi = icpevapi + icpevap*a
          gwati = gwati+gwat*a
          soilfrosti = soilfrosti + frostdepth*a
          soiltempi = soiltempi + (soilstate%temp(1,j,i)*soilthick(1,j)+soilstate%temp(2,j,i)*soilthick(2,j)+soilstate%temp(3,j,i)*soilthick(3,j))/soildepth(3,j)*a  !CP121114
          !Additional for CryoProcess
          soil2tempi = soil2tempi + (soilstate%temp(1,j,i)*soilthick(1,j)+soilstate%temp(2,j,i)*soilthick(2,j))/(soilthick(1,j)+soilthick(2,j))*a  !DG150901
          deeptempi = deeptempi + soilstate%deeptemp(j,i)*a  !DG150901
          !------------------------------------------------------------
          soilltmpi(1) = soilltmpi(1) + soilstate%temp(1,j,i)*a
          irrappli = irrappli + (irrappl(1)+irrappl(2))*a
          soildefi = soildefi + smdef*a
          IF(soilstate%water(1,j,i)>pwmm(1,j))THEN
            standsoili = standsoili + (soilstate%water(1,j,i)-pwmm(1,j))*a
            srffi = srffi + (pwmm(1,j)+soilstate%water(2,j,i)-wpmm(1,j)-wpmm(2,j))/(fcmm(1,j)+fcmm(2,j))*a
            smfdi = smfdi + (pwmm(1,j)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(1000.*soildepth(3,j))*a
            smfpi = smfpi + (pwmm(1,j)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(pwmm(1,j)+pwmm(2,j)+pwmm(3,j))*a
            srfdi = srfdi + (pwmm(1,j)+soilstate%water(2,j,i))/(1000.*soildepth(2,j))*a
            srfpi = srfpi + (pwmm(1,j)+soilstate%water(2,j,i))/(pwmm(1,j)+pwmm(2,j))*a
          ELSE
            srffi = srffi + (soilstate%water(1,j,i)+soilstate%water(2,j,i)-wpmm(1,j)-wpmm(2,j))/(fcmm(1,j)+fcmm(2,j))*a
            smfdi = smfdi + (soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(1000.*soildepth(3,j))*a
            smfpi = smfpi + (soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(pwmm(1,j)+pwmm(2,j)+pwmm(3,j))*a
            srfdi = srfdi + (soilstate%water(1,j,i)+soilstate%water(2,j,i))/(1000.*soildepth(2,j))*a
            srfpi = srfpi + (soilstate%water(1,j,i)+soilstate%water(2,j,i))/(pwmm(1,j)+pwmm(2,j))*a
          ENDIF
          IF(conductN)THEN
            fastNp1i = fastNp1i + soilstate%fastN(1,j,i)*a
            humusNp1i = humusNp1i + soilstate%humusN(1,j,i)*a
          ENDIF
          IF(conductP)THEN
            humusPp1i = humusPp1i + soilstate%humusP(1,j,i)*a
            partPp1i = partPp1i + soilstate%partP(1,j,i)*a
            fastPp1i = fastPp1i + soilstate%fastP(1,j,i)*a
          ENDIF
          IF(conductC)THEN
            humusCp1i = humusCp1i + soilstate%humusC(1,j,i)*a
            fastCp1i  = fastCp1i  + soilstate%fastC(1,j,i)*a
          ENDIF

          !Soil load output, mm*mg/L
          IF(numsubstances>0)THEN
            soillgrossload(1,:)=soillgrossload(1,:)+a*cruralflow(:)*horizontalflows(1)+a*cruralflow(:)*horizontalflows(2)+  &
              (atmdepload1+atmdepload2+cultivload(1,:)+cultivload(2,:))/basin(i)%area*1.E6
            soillnetload(1,:)=soillnetload(1,:)+a*cverticalflows(2,:)*verticalflows(2)+a*crunofflows(:,1)*runofflows(1)+  &
              a*crunofflows(:,2)*runofflows(2)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)+   &
              a*csrunoff(:)*surfaceflow
            soillgrossload(2,:)=soillgrossload(2,:)+a*cverticalflows(2,:)*verticalflows(2)+a*cruralflow(:)*horizontalflows(3)
            soillnetload(2,:)=soillnetload(2,:)+a*crunofflows(:,3)*runofflows(3)+a*crunofflows(:,6)*runofflows(6)
            soillgrossload(3,:)=soillgrossload(3,:)+a*cverticalflows(2,:)*verticalflows(2)+a*cruralflow(:)*horizontalflows(3)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)
            soillnetload(3,:)=soillnetload(3,:)+a*crunofflows(:,3)*runofflows(3)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)+a*crunofflows(:,6)*runofflows(6)
          ENDIF

          IF(soilthick(2,j)>0)THEN
            asumseclayer = asumseclayer + a
            soili2 = soili2+soilstate%water(2,j,i)*a
            csoili2(:) = csoili2(:) + soilstate%conc(:,2,j,i)*soilstate%water(2,j,i)*a   !accumulation of amount
            runoff2i = runoff2i + runofflows(2)*a
            soilltmpi(2) = soilltmpi(2) + soilstate%temp(2,j,i)*a
            IF(conductN)THEN
              fastNp2i = fastNp2i + soilstate%fastN(2,j,i)*a
              humusNp2i = humusNp2i + soilstate%humusN(2,j,i)*a
            ENDIF
            IF(conductP)THEN
              humusPp2i = humusPp2i + soilstate%humusP(2,j,i)*a
              partPp2i = partPp2i + soilstate%partP(2,j,i)*a
              fastPp2i = fastPp2i + soilstate%fastP(2,j,i)*a
            ENDIF
            IF(conductC)THEN
              humusCp2i = humusCp2i + soilstate%humusC(2,j,i)*a
              fastCp2i  = fastCp2i  + soilstate%fastC(2,j,i)*a
            ENDIF
            IF(soilthick(3,j)>0)THEN
              asumthdlayer = asumthdlayer + a
              soili3 = soili3+soilstate%water(3,j,i)*a
              csoili3(:) = csoili3(:) + soilstate%conc(:,3,j,i)*soilstate%water(3,j,i)*a   !accumulation of amount
              runoff3i = runoff3i + runofflows(3)*a
              soilltmpi(3) = soilltmpi(3) + soilstate%temp(3,j,i)*a
              IF(conductN)THEN
                fastNp3i = fastNp3i + soilstate%fastN(3,j,i)*a
                humusNp3i = humusNp3i + soilstate%humusN(3,j,i)*a
              ENDIF
              IF(conductP)THEN
                humusPp3i = humusPp3i + soilstate%humusP(3,j,i)*a
                partPp3i = partPp3i + soilstate%partP(3,j,i)*a
                fastPp3i = fastPp3i + soilstate%fastP(3,j,i)*a
              ENDIF
              IF(conductC)THEN
                humusCp3i = humusCp3i + soilstate%humusC(3,j,i)*a
                fastCp3i  = fastCp3i  + soilstate%fastC(3,j,i)*a
              ENDIF
              soili = soili+soilstate%water(1,j,i)*a+soilstate%water(2,j,i)*a+soilstate%water(3,j,i)*a
              csoili(:) = csoili(:) + soilstate%conc(:,1,j,i)*soilstate%water(1,j,i)*a + soilstate%conc(:,2,j,i)*soilstate%water(2,j,i)*a + soilstate%conc(:,3,j,i)*soilstate%water(3,j,i)*a   !accumulation of amount
              runoffi=runoffi+runofflows(1)*a+runofflows(2)*a+runofflows(3)*a+surfaceflow*a+runoffd*a
              crunoffi(:) = crunoffi(:) + crunoff1(:)*runofflows(1)*a + crunoff2(:)*runofflows(2)*a + crunoff3(:)*runofflows(3)*a + csrunoff(:)*surfaceflow*a + crunoffd(:)*runoffd*a    !accumulation of amount
              runoffrzi=runoffrzi+runofflows(1)*a+runofflows(2)*a+runofflows(3)*a
              crunoffrzi(:) = crunoffrzi(:) + crunoff1(:)*runofflows(1)*a + crunoff2(:)*runofflows(2)*a + crunoff3(:)*runofflows(3)*a    !accumulation of amount
            ELSE
              soili = soili+soilstate%water(1,j,i)*a+soilstate%water(2,j,i)*a
              csoili(:) = csoili(:) + soilstate%conc(:,1,j,i)*soilstate%water(1,j,i)*a + soilstate%conc(:,2,j,i)*soilstate%water(2,j,i)*a   !accumulation of amount
              runoffi = runoffi+runofflows(1)*a+runofflows(2)*a+surfaceflow*a+runoffd*a
              crunoffi(:) = crunoffi(:) + crunoff1(:)*runofflows(1)*a + crunoff2(:)*runofflows(2)*a + csrunoff(:)*surfaceflow*a + crunoffd(:)*runoffd*a    !accumulation of amount
              runoffrzi = runoffrzi+runofflows(1)*a+runofflows(2)*a
              crunoffrzi(:) = crunoffrzi(:) + crunoff1(:)*runofflows(1)*a + crunoff2(:)*runofflows(2)*a    !accumulation of amount
            ENDIF
          ELSE
            soili = soili+soilstate%water(1,j,i)*a
            csoili(:) = csoili(:) + soilstate%conc(:,1,j,i)*soilstate%water(1,j,i)*a
            runoffi=runoffi+runofflows(1)*a+surfaceflow*a+runoffd*a
            crunoffi(:) = crunoffi(:) + crunoff1(:)*runofflows(1)*a + csrunoff(:)*surfaceflow*a + crunoffd(:)*runoffd*a    !accumulation of amount
            runoffrzi = runoffrzi+runofflows(1)*a
            crunoffrzi(:) = crunoffrzi(:) + crunoff1(:)*runofflows(1)*a     !accumulation of amount
          ENDIF
          !Loads from soil: soil runoff, surface runoff, drain pipe runoff
          IF(conductload)THEN
            Lstream(j,:) = (runofflows(1)*crunoff1(:)+runofflows(2)*crunoff2(:)+  &
                            runofflows(3)*crunoff3(:)+surfaceflow*csrunoff(:)+    &
                            crunoffd(:)*runoffd)*classarea    !Load in runoff (kg/timestep) (incl. surface runoff and tile runoff)
          ENDIF
          !Snow cover FSC (must be rescaled with the fraction of ilake+olake before output)
          snowcovi = snowcovi+frozenstate%snowcov(j,i)*a*snowplain_area_fraction   !+0% snow on glacier
          !Forest and Open field snow variables defined by vegtype
          IF(vegtype(j).EQ.1)THEN
            !OPEN FIELD LANDUSE
            snowdepth_open   = snowdepth_open + snowdepth*a*snowplain_area_fraction
            IF(snowdepth.GT.0)THEN
              snowdensity_open = snowdensity_open + 0.1 * frozenstate%snow(j,i)/snowdepth * a*snowplain_area_fraction
            ENDIF
            snowcover_open   = snowcover_open + 10. * frozenstate%snowcov(j,i)*a*snowplain_area_fraction
            swe_open         = swe_open + frozenstate%snow(j,i)*a*snowplain_area_fraction
            area_open        = area_open + a*snowplain_area_fraction
          ELSEIF(vegtype(j).EQ.2)THEN
            !FOREST LANDUSE
            snowdepth_forest = snowdepth_forest + snowdepth*a*snowplain_area_fraction
            IF(snowdepth.GT.0)THEN
              snowdensity_forest = snowdensity_forest + 0.1 * frozenstate%snow(j,i)/snowdepth * a*snowplain_area_fraction
            ENDIF
            snowcover_forest = snowcover_forest + 10. * frozenstate%snowcov(j,i)*a*snowplain_area_fraction
            swe_forest       = swe_forest + frozenstate%snow(j,i)*a*snowplain_area_fraction
            area_forest      = area_forest + a*snowplain_area_fraction
          ENDIF
          !DG20151126: I did not try to update these waterbalance flows to the new suggested treatment of interception and throughfall!!
          IF(conductwb)THEN
            wbflows(w_sfallsnow,i)= wbflows(w_sfallsnow,i) + snowfall * a * old_snowfield_part
            wbflows(w_smeltsl1,i) = wbflows(w_smeltsl1,i) + infiltrationflows(1) * infiltrationflows(2) * a
            wbflows(w_smeltsr,i)  = wbflows(w_smeltsr,i)  + infiltrationflows(1) * infiltrationflows(3) * a
            wbflows(w_smeltmp1,i) = wbflows(w_smeltmp1,i) + infiltrationflows(1) * infiltrationflows(4) * a
            wbflows(w_smeltmp2,i) = wbflows(w_smeltmp2,i) + infiltrationflows(1) * infiltrationflows(5) * a
            wbflows(w_smeltmp3,i) = wbflows(w_smeltmp3,i) + infiltrationflows(1) * infiltrationflows(6) * a
            wbflows(w_infrain,i)  = wbflows(w_infrain,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(2) * a
            wbflows(w_rainsr,i)   = wbflows(w_rainsr,i)  + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(3) * a
            wbflows(w_rainmp1,i)  = wbflows(w_rainmp1,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(4) * a
            wbflows(w_rainmp2,i)  = wbflows(w_rainmp2,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(5) * a
            wbflows(w_rainmp3,i)  = wbflows(w_rainmp3,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(6) * a
            wbflows(w_evap1,i)    = wbflows(w_evap1,i) + evapflows(1) * a
            wbflows(w_evap2,i)    = wbflows(w_evap2,i) + evapflows(2) * a
            wbflows(w_evap3,i)    = wbflows(w_evap3,i) + evapflows(3) * a * old_snowfield_part
            wbflows(w_gwrunf1,i)  = wbflows(w_gwrunf1,i) + runofflows(1) * a
            wbflows(w_gwrunf2,i)  = wbflows(w_gwrunf2,i) + runofflows(2) * a
            wbflows(w_gwrunf3,i)  = wbflows(w_gwrunf3,i) + runofflows(3) * a
            wbflows(w_tile1,i)    = wbflows(w_tile1,i) + runofflows(4) * a
            wbflows(w_tile2,i)    = wbflows(w_tile2,i) + runofflows(5) * a
            wbflows(w_tile3,i)    = wbflows(w_tile3,i) + runofflows(6) * a
            wbflows(w_surfrf,i)   = wbflows(w_surfrf,i) + runofflows(7) * a
            wbflows(w_perc1,i)    = wbflows(w_perc1,i) + verticalflows(1) * a
            wbflows(w_perc2,i)    = wbflows(w_perc2,i) + verticalflows(2) * a
            wbflows(w_upwell1,i)  = wbflows(w_upwell1,i) + (verticalflows(3)+verticalflows(5)) * a
            wbflows(w_upwell2,i)  = wbflows(w_upwell2,i) + (verticalflows(4)+verticalflows(6)) * a
            wbflows(w_rural1,i)   = wbflows(w_rural1,i) + horizontalflows(1) * a
            wbflows(w_rural2,i)   = wbflows(w_rural2,i) + horizontalflows(2) * a
            wbflows(w_rural3,i)   = wbflows(w_rural3,i) + horizontalflows(3) * a
            wbirrflows(w_apply1,i)  = wbirrflows(w_apply1,i) + irrappl(1) * a
            wbirrflows(w_apply2,i)  = wbirrflows(w_apply2,i) + irrappl(2) * a
            IF(classmodel(j)==glacier_model)THEN
              wbstores(w_glacier,i) = frozenstate%glacvol(i)*genpar(m_glacdens)
              wbflows(w_stoice,i)   = glacierflows(1)
              wbflows(w_precglac,i) = glacierflows(2)
              wbflows(w_gmeltsl1,i) = infiltrationflows(7) * infiltrationflows(2) * a
              wbflows(w_gmeltsr,i)  = infiltrationflows(7) * infiltrationflows(3) * a
              wbflows(w_gmeltmp1,i) = infiltrationflows(7) * infiltrationflows(4) * a
              wbflows(w_gmeltmp2,i) = infiltrationflows(7) * infiltrationflows(5) * a
              wbflows(w_gmeltmp3,i) = infiltrationflows(7) * infiltrationflows(6) * a
              wbflows(w_evap4,i)    = wbflows(w_evap4,i) + evapflows(4) * a * (1. - old_snowfield_part)
            ENDIF
          ENDIF
        ENDIF !a>0
      ENDDO   !j
        
      !Calculate average of soil pools for classes with crops in rotation
      IF(dayno==360)THEN
        numrotations = MAXVAL(rotation(:))    !LÄGG i modvar!
        IF(numrotations>0)THEN
          IF(i_sp>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusP)
          IF(i_sp>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%partP)
          IF(i_in>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusN)
          IF(i_oc>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusC)
        ENDIF
      ENDIF
        
      !Calculate subbasin mean over land classes for concentrations that was accumulated in the class-loop
      IF(soili>0.0)     csoili(:) = csoili(:) / soili
      IF(soili1>0.0)    csoili1(:) = csoili1(:) / soili1
      IF(nummaxlayers>1)THEN
        IF(soili2>0.0)   csoili2(:) = csoili2(:) / soili2
        IF(soili3>0.0)   csoili3(:) = csoili3(:) / soili3
      ENDIF
      IF(runoffi>0.0) crunoffi(:) = crunoffi(:) / runoffi
      !Land area volumes for soil variables
      IF(asum>0)THEN
        divasum = 1. / asum
      ELSE
        divasum = 0.
      ENDIF
      IF(asumseclayer>0)THEN
        divasum2 = 1. / asumseclayer
      ELSE
        divasum2 = 0.
      ENDIF
      IF(asumthdlayer>0)THEN
        divasum3 = 1. / asumthdlayer
      ELSE
        divasum3 = 0.
      ENDIF
      IF(conductwb)THEN
        wbstores(w_snow,i) = snowi * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil1,i) = soili1 * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil2,i) = soili2 * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil3,i) = soili3 * basin(i)%area * 1.E-3  !m3
        !After this point soilstate%water and frozenstate%snow may not be changed, or the water balance will be off.
      ENDIF
      !Calculate subbasin mean over land area (or second or third soillayer area)
      snowi = snowi * divasum
      soili = soili * divasum
      soili12 = (soili1 + soili2) * divasum       !Note! calc before scaling soili1 and soili2
      soili1 = soili1 * divasum
      standsoili = standsoili * divasum
      runoffi = runoffi * divasum
      runoff1i = runoff1i * divasum
      runoffdi = runoffdi * divasum
      runoffsri = runoffsri * divasum 
      plantuptakei= plantuptakei * divasum
      nitrifi = nitrifi * divasum    
      denitrifi = denitrifi * divasum
      atmdepi =atmdepi * divasum
      atmdep2i =atmdep2i * divasum
      melti = melti * divasum              !OBS only snow melt on land
      glacmelti = glacmelti * divasum     !OBS only glaciermelt on land, DG20151126 for CryoProcess
      !Epoti and eacti is calculated for whole basin area (incl lake), as is corrpreci and corrtempi
      landeacti = eacti * divasum    !eacti for land area
      gwati = gwati * divasum
      soilfrosti = soilfrosti * divasum
      soiltempi = soiltempi * divasum
      deeptempi = deeptempi * divasum   !DG20151126 for CryoProcess
      soil2tempi = soil2tempi * divasum  !DG20151126 for CryoProcess
      soilltmpi(1) = soilltmpi(1) * divasum
      snowdepthi = snowdepthi * divasum
      irrappli = irrappli * divasum 
      soildefi = soildefi * divasum 
      srffi = srffi * divasum
      smfdi = smfdi * divasum
      smfpi = smfpi * divasum
      srfdi = srfdi * divasum
      srfpi = srfpi * divasum
      fastNp1i = fastNp1i * divasum
      humusNp1i = humusNp1i * divasum
      humusPp1i = humusPp1i * divasum
      partPp1i = partPp1i * divasum
      fastPp1i = fastPp1i * divasum
      humusCp1i = humusCp1i * divasum
      fastCp1i = fastCp1i * divasum
      snowcovi = snowcovi * divasum
      snowmaxi = snowmaxi * divasum
      IF(nummaxlayers>1)THEN
        soili2    = soili2 * divasum2
        runoff2i  = runoff2i * divasum2
        soilltmpi(2) = soilltmpi(2) * divasum2
        fastNp2i  = fastNp2i * divasum2
        humusNp2i = humusNp2i * divasum2
        humusPp2i = humusPp2i * divasum2
        partPp2i  = partPp2i * divasum2
        fastPp2i  = fastPp2i * divasum2
        humusCp2i = humusCp2i * divasum2
        fastCp2i  = fastCp2i * divasum2
      ENDIF
      IF(nummaxlayers==3)THEN
        soili3    = soili3 * divasum3
        runoff3i  = runoff3i * divasum3
        soilltmpi(3) = soilltmpi(3) * divasum3
        fastNp3i  = fastNp3i * divasum3
        humusNp3i = humusNp3i * divasum3
        humusPp3i = humusPp3i * divasum3
        partPp3i  = partPp3i * divasum3
        fastPp3i  = fastPp3i * divasum3
        humusCp3i = humusCp3i * divasum3
        fastCp3i  = fastCp3i * divasum3
      ENDIF
      IF(conductwb)THEN
        wbflows(w_sfallsnow,i)= wbflows(w_sfallsnow,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltsl1,i) = wbflows(w_smeltsl1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltsr,i)  = wbflows(w_smeltsr,i)  * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp1,i) = wbflows(w_smeltmp1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp2,i) = wbflows(w_smeltmp2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp3,i) = wbflows(w_smeltmp3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_infrain,i)  = wbflows(w_infrain,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainsr,i)   = wbflows(w_rainsr,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp1,i)  = wbflows(w_rainmp1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp2,i)  = wbflows(w_rainmp2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp3,i)  = wbflows(w_rainmp3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap1,i)    = wbflows(w_evap1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap2,i)    = wbflows(w_evap2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap3,i)    = wbflows(w_evap3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf1,i)  = wbflows(w_gwrunf1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf2,i)  = wbflows(w_gwrunf2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf3,i)  = wbflows(w_gwrunf3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile1,i)    = wbflows(w_tile1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile2,i)    = wbflows(w_tile2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile3,i)    = wbflows(w_tile3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_surfrf,i)   = wbflows(w_surfrf,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_perc1,i)    = wbflows(w_perc1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_perc2,i)    = wbflows(w_perc2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_upwell1,i)  = wbflows(w_upwell1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_upwell2,i)  = wbflows(w_upwell2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural1,i)   = wbflows(w_rural1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural2,i)   = wbflows(w_rural2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural3,i)   = wbflows(w_rural3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbirrflows(w_apply1,i)  = wbirrflows(w_apply1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbirrflows(w_apply2,i)  = wbirrflows(w_apply2,i) * basin(i)%area * 1.E-3   !m3/ts
        IF(glacierexist)THEN
          wbflows(w_gmeltsl1,i) = wbflows(w_gmeltsl1,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltsr,i)  = wbflows(w_gmeltsr,i)  * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp1,i) = wbflows(w_gmeltmp1,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp2,i) = wbflows(w_gmeltmp2,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp3,i) = wbflows(w_gmeltmp3,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_evap4,i)    = wbflows(w_evap4,i) * basin(i)%area * 1.E-3   !m3/ts
        ENDIF
      ENDIF
      !Forest and open land snow variables
      IF(area_open.GT.0.)THEN
        !OPEN LAND
        area_open=1./area_open
        snowdepth_open   = snowdepth_open * area_open
        snowdensity_open = snowdensity_open * area_open
        snowcover_open   = snowcover_open * area_open
        swe_open         = swe_open * area_open
      ELSE
        snowdepth_open   = -9999.
        snowdensity_open = -9999.
        snowcover_open   = -9999.
        swe_open         = -9999.
      ENDIF
      IF(area_forest.GT.0.)THEN
        !FOREST LAND
        area_forest=1./area_forest
        snowdepth_forest   = snowdepth_forest * area_forest
        snowdensity_forest = snowdensity_forest * area_forest
        snowcover_forest   = snowcover_forest * area_forest
        swe_forest         = swe_forest * area_forest
      ELSE
        snowdepth_forest   = -9999.
        snowdensity_forest = -9999.
        snowcover_forest   = -9999.
        swe_forest         = -9999.
      ENDIF  

      !>Routing in rivers and lakes
      !This means that the primary unit is m3/s below. The waterlevel of the lakes are measured in mm though.

      !Remember olake waterstage for last time step, to be used in updating
      olakewstold=lakestate%water(2,i)
!      IF(conductN.OR.conductP.OR.conductC.OR.conductT) olakewstold=olakewstold+lakestate%slowwater(2,i)
      IF(conductN.OR.conductP.OR.conductC.OR.conductT.OR.conductWS) olakewstold=olakewstold+lakestate%slowwater(2,i) !DG20151126 CryoProcesses
        
      !Temperature in river and lakes
      IF(modeloption(p_swtemperature)==0)THEN
        CALL calculate_water_temperature(i,tempi(i),riverstate,lakestate)
      ENDIF
        
      !Inflow from upstream basins has been accumulated. Calculate the concentration.
      qupstream = accinflow(i)
      IF(qupstream>0)THEN
        cupstream = acccinflow(:,i) / qupstream
      ELSE
        cupstream = 0.
      ENDIF

      !Initiations
      lakefreezeupday=0
      lakebreakupday=0
      riverfreezeupday=0
      riverbreakupday=0
      lakeareai=0        !DG20151126 CryoProcesses
        
      !Local or main river and lakes, local (=1) and coupled to upstream/main (=2)
      routingloop:  &
      DO itype = 1,2
           
        !Initiation of river calculations
        riverarea(itype) = 0.
        IF(itype==1 .AND. slc_lriver>0 .OR. itype==2 .AND. slc_mriver>0)THEN
          IF(itype==1) j=slc_lriver
          IF(itype==2) j=slc_mriver
          a = classbasin(i,j)%part
        ELSE
          a = 0.
          j = 0
        ENDIF
              
        !Calculate inflow to river
        IF(itype==1)THEN
          IF( (runoffrzi + runoffnorzi)<genpar(m_nclim1) )THEN
            adjnoncontr(i) = genpar(m_ncfact1)*basin(i)%noncontr
            IF(adjnoncontr(i)>0.99) adjnoncontr(i) = 0.99
          ELSEIF( (runoffrzi + runoffnorzi)<genpar(m_nclim2) )THEN
            adjnoncontr(i) = genpar(m_ncfact2)*basin(i)%noncontr
          ELSE
            adjnoncontr(i) = genpar(m_ncfact3)*basin(i)%noncontr
          ENDIF
          !!adjnoncontr(i) = 0.0 ! !! basin(i)%noncontr   !!! TEMP, REMOVE: test using only published NCA or zero NCA
          qin   = (runoffrzi + runoffnorzi) * basin(i)%area * (1. - adjnoncontr(i)) / (seconds_per_timestep * 1000.)      !m3/s
          IF(runoffrzi+runoffnorzi>0.)THEN
            cin   = (crunoffrzi(:) + crunoffnorzi(:)) / (runoffrzi + runoffnorzi)              !mg/L (crunoffrzi/crunoffnorzi is amount)
          ELSE
            cin = 0.
          ENDIF
          IF(conductload) Lpathway(:,1) = cin(:) * qin * seconds_per_timestep / 1000.   !Total load from land, A (kg/time step)
        ELSEIF(itype==2)THEN
          qin   = lakeoutflow(1) + qupstream                       !m3/s
          IF(qin>0)THEN
            cin = (clakeoutflow(:,1)*lakeoutflow(1) + cupstream(:)*qupstream) / qin
          ELSE
            cin = 0.
          ENDIF
        ENDIF
           
        !Keep T2 concentration in memory, in case water is being extracted/added by irrigation and/or sources/wetlands
        IF(i_t2>0)ctemp_T2 = cin(i_t2)
        IF(i_sm>0) ctemp_sm = cin(i_sm) !DG20151126 CryoProcesses
        IF(i_gm>0) ctemp_gm = cin(i_gm) 
        IF(i_rn>0) ctemp_rn = cin(i_rn) 
        IF(i_li>0) ctemp_li = cin(i_sm)
        IF(i_ri>0) ctemp_ri = cin(i_gm) 
        IF(i_si>0) ctemp_si = cin(i_rn) !DG20151126 CryoProcesses 
        
         
        !Irrigation water removal from sources and calculate irrigation to be applied the next day
        IF(doirrigation .AND. itype==2)THEN
          CALL calculate_irrigation(i,naquifers,genpar(m_regirr),regpar(m_pirrs,basin(i)%parregion),      &
                 regpar(m_pirrg,basin(i)%parregion),qin,lakestate%water(:,i),riverstate%water(itype,i),aquiferstate%water,pwneedi,irrevap,gwremi,     &
                 ldremi,lrremi,rsremi,cin,regpar(m_cirrsink,basin(i)%parregion),irrsinkmass,  &
                 genpar(m_irrcomp),lakestate%conc(:,1,i),lakestate%conc(:,2,i),riverstate%conc(:,itype,i),aquiferstate%conc, &
                 soilstate,miscstate%nextirrigation,miscstate%cnextirrigation,irrigationflows,regionalirrflows,regionalirrevap, &
                 aqirrloss)
          IF(conductwb)THEN
            wbirrflows(w_wdfromdg,i) = irrigationflows(1)
            wbirrflows(w_wdfromil,i) = irrigationflows(2)
            wbirrflows(w_wdfromol,i) = irrigationflows(3)
            wbirrflows(w_wdfrommr,i) = irrigationflows(4)
            wbirrflows(w_wdoutside,i) = irrigationflows(6)
            wbirrflows(w_rgrwtoir,i) = irrigationflows(7)
            wbirrflows(w_evapirrc,i) = wbirrflows(w_evapirrc,i) + irrigationflows(5)  !local sources
            wbirrflows(w_evapirrc,:) = wbirrflows(w_evapirrc,:) + regionalirrevap(:)  !regional sources, local network
          ENDIF
        ENDIF
           
        !Add point sources and local diffuse source to river (in)flow
        IF(itype==1)THEN
          CALL add_diffuse_source_to_local_river(i,qin,cin,Lruralb,nutrientloadflows(1))
          IF(conductload) Lpathway(:,2) = cin(:) * qin * seconds_per_timestep * 1.E-3  !Load after adding rural b (kg/timestep), point B
        ENDIF
        IF(itype==2)THEN
          CALL add_point_sources_to_main_river(i,qin,cin,Lpoints,nutrientloadflows(2))
        ENDIF
          
        !Add aquifer inflow to main river
        IF(modeloption(p_deepgroundwater)==2.AND.itype==2)THEN
          CALL add_aquifer_flow_to_river(i,numsubstances,qin,cin,grwtomr,grwloadmr)
          IF(conductwb) wbflows(w_rgrwtomr,i) = grwtomr
          outvar(i,241) = grwtomr
          IF(conductload) Lgrwmr = grwloadmr
        ENDIF

        IF(itype==2.AND.conductload) Lpathway(:,9) = cin(:) * qin * seconds_per_timestep * 1.E-3  !Load after adding point sources and aquifer inflow (kg/timestep), point I

        !Calculate river wetland          
        IF(wetlandexist)  CALL calculate_river_wetland(i,itype,numsubstances,miscstate%temp5(i),miscstate%temp30(i),qin,cin,riverstate%cwetland(:,itype,i))
        IF(conductload)THEN
          IF(itype==1) Lpathway(:,3) = qin * cin * seconds_per_timestep * 1.E-3   !Total load after local river wetland (kg/timestep), point C
          IF(itype==2) Lpathway(:,10) = qin * cin * seconds_per_timestep * 1.E-3  !Total load after main river wetland (kg/timestep), point J
        ENDIF
           
        !Restore T2 concentration from temp, in case water was extracted/added by irrigation and/or sources/wetlands
        !(later, we should calculate ice and T2 in river wetlands explicitly)
        IF(i_t2>0) cin(i_t2) = ctemp_T2
        IF(i_sm>0) cin(i_sm) = ctemp_sm
        IF(i_gm>0) cin(i_gm) = ctemp_gm
        IF(i_rn>0) cin(i_rn) = ctemp_rn
        IF(i_li>0) cin(i_sm) = ctemp_li
        IF(i_ri>0) cin(i_gm) = ctemp_ri
        IF(i_si>0) cin(i_rn) = ctemp_si
        
        !Translation (delay) in river (in)flow
        CALL translation_in_river(i,itype,qin,cin,transq,transc,riverstate)
           
        !Calculate river processes in river water volume
        !The river water volume is a storage for damping river flow, mixing of concentrations 
        !with deadvolume and applying processes for nutrient and other substances
         
        !Add river inflow to riverboxi
        CALL add_water(numsubstances,riverstate%water(itype,i),riverstate%conc(:,itype,i),transq * seconds_per_timestep,transc)
         
        !Calculate river dimensions, velocity and mean flow for use in substance processes calculation
        CALL calculate_river_characteristics(i,itype,transq,conductN.OR.conductP,riverstate,riverdepth,riverarea(itype),riverQbank)
           
        !Calculate precipitation, atmospheric deposition and evaporation of river with area
        IF(a>0)THEN      !river has area
          riverarea(itype) = a * basin(i)%area       !m2
          riverareakm2 = riverarea(itype)*1.E-6      !km2
          atmdepload1 = 0.
          atmdepload2 = 0.
            
          !Forcing data and atmospheric deposition on river 
!          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
!                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,netrad,wind,sffrac, &
                  snowfall,rainfall,pcorr)
          IF(i_rn>0)THEN
            cprecj(i_rn)=1.-sffrac
            cprecj(i_sm)=sffrac
          ENDIF
          CALL set_class_precipitation_concentration_and_load(numsubstances,riverareakm2,outvar(i,3),temp,prec,cprec,cprecj,atmdepload1)
          CALL add_dry_deposition_to_river(i,riverarea(itype),itype,atmdepload2,    &
               load(i)%indrydep(vegtype(j)),landpar(m_drypp,classdata(j)%luse), &
               riverstate)
          IF(conductload)THEN
            Latmdep(j,1,:) = atmdepload1
            Latmdep(j,2,:) = atmdepload2
          ENDIF
               
          !Snowfall/rainfall on river
          !CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac)
            
          !"Temperature" T2 concentration in precipitation
          IF(i_t2>0)THEN
            CALL get_rivertempvol(i,itype,riverstate,meanrivertemp,totrivervol)  !Get total river water volume and mean T2 temperature
            CALL add_T2_concentration_in_precipitation_on_water(prec,temp,  &  !Set correct T2 temperature concentration in precipitation
                  snowfall,rainfall,meanrivertemp,cprecj(i_t2),frozenstate%rivericecov(itype,i))
          ENDIF

          !Add precipitation with substances to river water
          IF(prec>0) CALL add_precipitation_to_river(i,itype,riverarea(itype),prec,cprecj,riverstate)
              
          !Evaporation of river, taking partial ice-cover into account if icemodel is enabled
          CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
          epot = epot * basincevpcorr(i)
          IF(modeloption(p_lakeriverice)>0)THEN
            epot = epot * (1. - frozenstate%rivericecov(itype,i))
          ENDIF
          IF(epot>0.)THEN
            CALL calculate_river_evaporation(i,j,itype,numsubstances,riverarea(itype),temp,epot,evapr,cevapr,riverstate)
          ELSE
            evapr = 0.
          ENDIF
            
          !Accumulate for mean over subbasin (accumulation for cevapi)
          corrpreci = corrpreci + prec*a
          corrtempi = corrtempi + temp*a
          rainfalli = rainfalli + rainfall * a
          snowfalli = snowfalli + snowfall * a
          epoti = epoti+epot*a
          eacti = eacti+evapr*a
          evapriveri = evapriveri + evapr*a
          netradi = netradi + netrad*a
          swradi = swradi + swrad*a
          IF(i_t1>0) cevapi = cevapi+cevapr(i_t1)*evapr*a
!          icpevapi = icpevapi + icpevap*a
          pcorri = pcorri + pcorr*a

          !Set water balance output
          IF(conductwb)THEN
            IF(itype==1)THEN
              wbflows(w_piriver,i) = prec * riverarea(itype)*1.E-3     !m3
              wbflows(w_eiriver,i) = evapr * riverarea(itype)*1.E-3
            ENDIF
            IF(itype==2)THEN
              wbflows(w_pmriver,i) = prec * riverarea(itype)*1.E-3
              wbflows(w_emriver,i) = evapr * riverarea(itype)*1.E-3
            ENDIF
          ENDIF
            
        ENDIF !(a>0)    
           
        !Calculate river ice and T2 processes
        IF(modeloption(p_lakeriverice)>0 .OR. i_t2>0)THEN
            
          !If no river class: Assume same precipitation and temperature conditions as for lake class
          IF(a==0)THEN
            IF(itype==1) j=slc_ilake
            IF(itype==2) j=slc_olake
!            CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
!                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
            CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,netrad,wind,sffrac, &
                  snowfall,rainfall,pcorr)
!            IF(modeloption(p_lakeriverice)>0)THEN
!              CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac)
!            ENDIF
          ENDIF
              
          !Calculate river T2 processes only
          IF(i_t2>0) CALL T2_processes_in_river(i,itype,temp,swrad,riversurftemp,riverarea(itype),frozenstate,riverstate,riverfreezeupday,freezeuparea)
          
          !River ice calculations
          IF(modeloption(p_lakeriverice)>0) CALL ice_processes_in_river(i,itype,classdata(j)%luse,snowfall,temp, & 
                                        riversurftemp,riversnowdepth,riverarea(itype),swrad,         &
                                        frozenstate,riverstate,riverfreezeupday,riverbreakupday,freezeuparea)
          IF(modeloption(p_swtemperature)==1)THEN
            CALL set_water_temperature(itype*2-1,i,riverstate,lakestate)
          ENDIF
        ENDIF

        !Regional groundwater flow from main river
        IF(itype==2 .AND. modeloption(p_deepgroundwater)==2)THEN
          CALL calculate_river_groundwaterflow_removal(i,j,basin(i)%subid,numsubstances,riverarea(itype),riverstate,grwout2)
          IF(conductwb) wbflows(w_rgrwofmr,i) = grwout2
        ENDIF
        
        !Abstraction of water from main river           
        IF(itype==2) CALL point_abstraction_from_main_river(i,itype,riverstate,nutrientloadflows(3))

        !Calculate substances processes in river
        CALL np_processes_in_river(i,itype,riverarea(itype),riverdepth,transq,riverQbank,       &
                (2.-incorr)*genpar(m_denitwr),   &
                (2.-incorr)*genpar(m_denitwrl),  &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodn),    &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodp),    &
                genpar(m_sedexp),genpar(m_limsedpp),riverstate)
        CALL oc_processes_in_river(i,itype,riverarea(itype),riverdepth,   &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodc),genpar(m_limsedpp),riverstate)
           
        !Calculate and remove outflow from riverboxi
        dampq = (riverrc(itype,i)*(riverstate%water(itype,i) - deadriver(itype,i)))
        IF(dampq<0.) dampq=0.     !safe for irrigation
        dampc = riverstate%conc(:,itype,i)
        CALL remove_water(riverstate%water(itype,i),numsubstances,riverstate%conc(:,itype,i),dampq,dampc,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,itype)
        IF(conductwb)THEN
          IF(itype==1) wbflows(w_irtomr,i) = dampq  !in case of no ilake
          IF(itype==2) wbflows(w_mrtool,i) = dampq
        ENDIF
        
        !Load after river transformations
        IF(conductload)THEN
          IF(itype==1) Lpathway(:,4) = dampc * dampq * 1.E-3   !Load at point C (downstream of local river), point D
          IF(itype==2) Lpathway(:,11) = dampc * dampq * 1.E-3  !Load at point I (downstream of main river), point K
        ENDIF
        dampq = dampq / seconds_per_timestep
           
        !Initiation of lake calculations
        a = 0.
        IF(itype==1 .AND. slc_ilake>0 .OR. itype==2 .AND. slc_olake>0)THEN
          IF(itype==1) j=slc_ilake
          IF(itype==2) j=slc_olake
          a = classbasin(i,j)%part
        ENDIF
           
        !Lake calculations
        IF(a>0)THEN      !lake exist
          lakearea = a * basin(i)%area       !m2
          lakeareai(itype) = lakearea !DG CryoProcesses
          lakeareakm2 = lakearea*1.E-6       !km2
          qunitfactor = seconds_per_timestep * 1000. / lakearea     !m3/s->mm/timestep
          atmdepload1 = 0.
          atmdepload2 = 0.
            
          !Forcing data and atmospheric deposition on lake 
!          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
!                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,netrad,wind,sffrac, &
                  snowfall,rainfall,pcorr)
          CALL set_class_precipitation_concentration_and_load(numsubstances,lakeareakm2,outvar(i,3),temp,prec,cprec,cprecj,atmdepload1)
          !DG cryoprocesses
          IF(i_rn>0)THEN
            cprecj(i_rn)=1.-sffrac
            cprecj(i_sm)=sffrac
          ENDIF
          CALL add_dry_deposition_to_lake(i,lakeareakm2,itype,atmdepload2,   &         !Add atmospheric deposition
               load(i)%indrydep(vegtype(j)),landpar(m_drypp,classdata(j)%luse), &
               lakestate)
          IF(conductload)THEN
            Latmdep(j,1,:) = atmdepload1
            Latmdep(j,2,:) = atmdepload2
          ENDIF
 
          !!Snowfall/rainfall on lake
          !CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac) 
            
          !"Temperature" T2 deposition in precipitation
          IF(i_t2>0) CALL add_T2_concentration_in_precipitation_on_water(prec,temp, &
                           snowfall,rainfall,lakestate%uppertemp(itype,i),cprecj(i_t2), &
                           frozenstate%lakeicecov(itype,i))

          !Add precipitation with substances to lake water
          IF(prec>0) CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),prec,cprecj)
            
          !Evaporation of lake, taking partial lakeice cover into account, if lakeice model is enabled
          CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
          epot = epot * basincevpcorr(i)
          IF(modeloption(p_lakeriverice)>0)THEN
            epot = epot * (1. - frozenstate%lakeicecov(itype,i))
          ENDIF
          IF(epot.GT.0.)THEN
            CALL calculate_actual_lake_evaporation(i,j,itype,numsubstances,temp,epot,evapl,cevapl,lakestate)
          ELSE
            evapl = 0.
          ENDIF
          !Accumulate for mean over subbasin (accumulation for cevapi)
          corrpreci = corrpreci + prec*a
          corrtempi = corrtempi + temp*a
          rainfalli = rainfalli + rainfall * a 
          snowfalli = snowfalli + snowfall * a 
          epoti = epoti+epot*a
          eacti = eacti+evapl*a
          evaplakei = evaplakei+evapl*a
          netradi = netradi + netrad*a
          swradi = swradi + swrad*a
          IF(i_t1>0) cevapi = cevapi+cevapl(i_t1)*evapl*a
!          icpevapi = icpevapi + icpevap*a
          IF(itype==2) qcinfli = qcinfli + (prec-evapl)*lakearea/seconds_per_timestep/1000.    !m3/s
          pcorri = pcorri + pcorr * a
          
          !Set water balance output
          IF(conductwb)THEN
            IF(itype==1)THEN
              wbflows(w_pilake,i) = prec * lakearea*1.E-3     !m3
              wbflows(w_eilake,i) = evapl * lakearea*1.E-3
            ENDIF
            IF(itype==2)THEN
              wbflows(w_polake,i) = prec * lakearea*1.E-3
              wbflows(w_eolake,i) = evapl * lakearea*1.E-3
            ENDIF
          ENDIF

          !Calculate and add inflow to lake
          inflowpart = 1.
          IF(itype==1)  inflowpart = basin(i)%ilakecatch * (1. - adjnoncontr(i))
          qin = dampq * inflowpart                              !m3/s
          IF(itype==2) qcinfli = qcinfli + qin
          qinmm = qin * qunitfactor                             !mm/d
          cin = dampc
          IF(conductload)THEN
            IF(itype==1)THEN
              Lpathway(:,6) = cin * qin * seconds_per_timestep * 1.E-3                 !Load at point F, local flow to ilake (kg/timestep)
              Lpathway(:,5) = Lpathway(:,4) - Lpathway(:,6) !Load at point E, bypassed local flow
            ENDIF
            IF(itype==2) Lpathway(:,12) = cin * qin * seconds_per_timestep * 1.E-3     !Load at point L inflow to olake (kg/timestep)
          ENDIF
          IF(conductwb .AND. itype==1 .AND. qin>0)THEN
            wbflows(w_irtoil,i) = qin * seconds_per_timestep
            wbflows(w_irtomr,i) = wbflows(w_irtomr,i) - wbflows(w_irtoil,i)
          ENDIF
          
          !MH2017: skips addition of inflow for special dam
          !intent: still calc npc flow, but overwrite inflow
          lakewst(itype)=lakestate%water(itype,i)
          CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),qinmm,cin)
          skipadd:IF(damindex(i)>0.AND.itype==2)THEN
            IF(dam(damindex(i))%purpose>=5)THEN
              lakestate%water(itype,i)=lakewst(itype)
            ENDIF
          ENDIF skipadd
          
            
          !Add regional groundwater inflow to olake
          IF(modeloption(p_deepgroundwater)==1)THEN
            CALL add_regional_groundwater_flow_to_olake(i,itype,numsubstances,  &
               qunitfactor,lakestate%water(itype,i),lakestate%conc(:,itype,i),  &
               qcinfli,grwloadlake,grwtool)
            IF(conductload)THEN
              Lgrwol(:) = grwloadlake
              Lpathway(:,12) = Lpathway(:,12) + grwloadlake
            ENDIF
            IF(conductwb) wbflows(w_rgrwtool,i) = grwtool
          ENDIF

          !Calculate T2 Lake Temperature Processes
          IF(i_t2>0) CALL T2_processes_in_lake(i,itype,temp,swrad,lakesurftemp,lakearea,frozenstate,lakestate,lakefreezeupday,freezeuparea)

          !Calculate Lake Ice Processes
          IF(modeloption(p_lakeriverice)>0) CALL ice_processes_in_lake(i,itype,classdata(j)%luse,snowfall,temp, & 
                                        lakesurftemp,lakesnowdepth,swrad,frozenstate,lakestate, & 
                                        lakefreezeupday,lakebreakupday,lakearea,freezeuparea) 
          IF(modeloption(p_swtemperature)==1)THEN
            CALL set_water_temperature(itype*2,i,riverstate,lakestate)
          ENDIF
          !Abstraction of water from outlet lake
          IF(itype==2) CALL point_abstraction_from_outlet_lake(i,itype,qunitfactor,lakestate,nutrientloadflows(4))

          !Calculate substances processes in lakes
          CALL np_processes_in_lake(i,itype,lakearea,    &
               (2.-incorr)*lakedatapar(lakedataparindex(i,itype),m_lddenitwl),  &
               lakedatapar(lakedataparindex(i,itype),m_ldwprodn),   &
               lakedatapar(lakedataparindex(i,itype),m_ldwprodp),   &
               (2.-oncorr)*lakedatapar(lakedataparindex(i,itype),m_ldsedon),    &
               lakedatapar(lakedataparindex(i,itype),m_ldsedpp),genpar(m_limsedon),genpar(m_limsedpp),lakestate)
          CALL oc_processes_in_lake(i,itype,lakearea,lakedatapar(lakedataparindex(i,itype),m_ldwprodc),    &
               genpar(m_limsedpp),lakedatapar(lakedataparindex(i,itype),m_ldsedoc),lakestate)
          
          !Calculate and remove outflow from lake
          lakewst(itype)=lakestate%water(itype,i)
          
          !MH2017: if special dam, store date, iteration, updates qin14day
          storespecbefore:IF(itype==2.AND.damindex(i)>0)THEN
            IF(dam(damindex(i))%purpose>=5)THEN
              dam(damindex(i))%curmonth = d%Month
              dam(damindex(i))%dayofmonth = d%Day
              dam(damindex(i))%idt=idt
              DO counter = 13,1,-1
                dam(damindex(i))%qin14day(counter+1) = dam(damindex(i))%qin14day(counter)
              ENDDO
              dam(damindex(i))%qin14day(1) = qin
            ENDIF
            IF(dam(damindex(i))%condswitch==1)THEN !conditioning required
              DO counter = 1,15
                dam(damindex(i))%wsl14day(counter) = lakestate%wsl14day(dam(damindex(i))%condindex,counter) !may need to TRANSFORM
                dam(damindex(i))%qo14day(counter) = lakestate%qo14day(dam(damindex(i))%condindex,counter)
              ENDDO
            ENDIF
          ENDIF storespecbefore
          
          
          IF(conductN.OR.conductP.OR.conductC.OR.conductT.OR.conductWS) lakewst(itype)=lakewst(itype)+lakestate%slowwater(itype,i)   !For TNPC-simulation with divisioned lake
          CALL calculate_outflow_from_lake(i,itype,qin,lakearea,lakewst(itype),qunitfactor,outflowm3s,outflowmm,lakestate,branchflow)
          CALL check_outflow_from_lake(i,itype,lakearea,lakewst(itype),qunitfactor,outflowm3s,outflowmm,lakestate,branchflow,qin)
          lakeoutflow(itype) = outflowm3s
          outflowshit = outflowmm+branchflow*qunitfactor !MM2016: include dam branchflow for removal from lake
          !MH2017: may be irrelevent, currently only utilized in NPC processes
          CALL remove_outflow_from_lake(i,itype,numsubstances,outflowshit,basin(i)%subid,concarr,lakestate) !no dam
      
          
          clakeoutflow(:,itype) = concarr
          
          !MH2017: if special dam, store inflow/outflow for next timestep, and WSL/Qo for conditioning
          !CHECK: how this interacts with branches, qoutprev may be redundant, so might qinprev (qin14day above)
          storespecafter:IF(itype==2.AND.damindex(i)>0)THEN
            IF(dam(damindex(i))%purpose>=5)THEN
              dam(damindex(i))%qinprev = qin
              dam(damindex(i))%qoutprev = outflowm3s + branchflow
            ENDIF            
            IF(dam(damindex(i))%purpose==6)THEN !store inflow,reservoir outflow
              !timing of qin7day() is different from qin14day(), DO NOT interchange them
              DO counter = 6,1,-1
                dam(damindex(i))%resqo7day(counter+1) = dam(damindex(i))%resqo7day(counter)
                dam(damindex(i))%qin7day(counter+1) = dam(damindex(i))%qin7day(counter)
              ENDDO
              IF(dam(damindex(i))%branchswt==0)THEN !reservoir = main
                dam(damindex(i))%resqo7day(1) = outflowm3s
              ELSEIF(dam(damindex(i))%branchswt==1)THEN !reservoir = branch
                dam(damindex(i))%resqo7day(1) = branchflow
              ENDIF
              dam(damindex(i))%qin7day(1) = qin
            ENDIF
          ENDIF storespecafter
          
          
          
          !Calculate load from ilake
          IF(conductload.AND.itype==1) Lpathway(:,7) = outflowm3s * clakeoutflow(:,itype) * seconds_per_timestep * 1.E-3 !Load at point G, outflow ilake (kg/timestep)
          IF(conductwb .AND. itype==1) wbflows(w_iltomr,i) = outflowm3s * seconds_per_timestep
            
          !Lake water stage for print out 
          lakewst(itype) = lakestate%water(itype,i)
          IF(ALLOCATED(lakestate%slowwater)) lakewst(itype) = lakewst(itype) + lakestate%slowwater(itype,i)
          
          !Add flow from river not going through local lake
          lakeoutflow(itype) = lakeoutflow(itype) + dampq*(1.-inflowpart)
          IF(inflowpart<1 .AND. lakeoutflow(itype)>0.)THEN
            clakeoutflow(:,itype) = (clakeoutflow(:,itype)*outflowm3s + dampq*(1.-inflowpart)*dampc)/(lakeoutflow(itype))
          ENDIF
          
          !Flow between lake parts for divided lake
          CALL calculate_flow_within_lake(i,itype,basin(i)%subid,lakestate)
            
        ELSE    !No lake, outflow = inflow
          lakeoutflow(itype) = dampq       !m3/s
          IF(itype==2) qcinfli = qcinfli + dampq
          clakeoutflow(:,itype) = dampc
          IF(conductload)THEN
            IF(itype==1)THEN
              Lpathway(:,5)  = Lpathway(:,4)           !Load at point E, bypassed local flow
              Lpathway(:,6)  = 0.                      !Load at point F
              Lpathway(:,7)  = 0.                      !Load at point G
            ENDIF
            IF(itype==2) Lpathway(:,12) = lakeoutflow(itype) * clakeoutflow(:,itype) * seconds_per_timestep * 1.E-3  !Load at point L flow in main river (kg/timestep)
          ENDIF
          lakewst(itype) = missing_value   
          wstlaketemp    = missing_value
          lakesurftemp(itype) = missing_value
          lakesnowdepth(itype) = missing_value
        ENDIF   !a>0
           
        IF(conductload)THEN
          IF(itype==1) Lpathway(:,8)  = lakeoutflow(itype) * clakeoutflow(:,itype) * seconds_per_timestep * 1.E-3  !Load at point H, outflow ilake + bypassed water (kg/timestep)
          IF(itype==2) Lpathway(:,13) = lakeoutflow(itype) * clakeoutflow(:,itype) * seconds_per_timestep * 1.E-3  !Load at point M, outflow of subbasin (kg/timestep)
        ENDIF

        storecond:IF(lakestate%useswitch(i)==1.AND.itype==2)THEN !used for another dam's conditioning
          DO counter = 14,1,-1
            lakestate%wsl14day(i,counter+1) = lakestate%wsl14day(i,counter)
            lakestate%qo14day(i,counter+1) = lakestate%qo14day(i,counter)
          ENDDO
          lakestate%qo14day(i,1) = lakeoutflow(itype) + branchflow !MH2017: does branchflow get initialized to zero?
          IF(damindex(i)>0)THEN !get w0ref from dam()
            lakestate%wsl14day(i,1) = (lakestate%slowwater(itype,i) + lakestate%water(itype,i)) / 1000. + dam(damindex(i))%w0ref - basin(i)%lakedepth(2)
          ELSEIF(lakeindex(i)>0)THEN !get w0ref from lake()
            lakestate%wsl14day(i,1) = (lakestate%slowwater(itype,i) + lakestate%water(itype,i)) / 1000. + lake(lakeindex(i))%w0ref - basin(i)%lakedepth(2)
          ENDIF !note: what if conditioning is taken from something that is not a dam OR a lake? Then it won't have a WSL anyway
        ENDIF storecond !maybe make a dam()%wsl variable


        !>Update modelled subbasin outflow with observations or otherwise
        IF(itype==2) THEN
          outflowtemp = lakeoutflow(itype)            
          wstlaketemp = lakewst(itype)            
          IF(doupdate(i_war))THEN
            lakeareatemp = classbasin(i,slc_olake)%part*basin(i)%area
            CALL apply_warupd(i,lakeareatemp,olakewstold,lakewst(itype),lakeoutflow(itype),miscstate%updatestationsarcorr(i),lakestate,branchflow)
          ENDIF
          IF(doupdate(i_quseobs)) CALL apply_quseobs(i,lakeoutflow(itype))
          IF(doupdate(i_qar))     CALL apply_qarupd(i,outflowtemp,lakeoutflow(itype),miscstate%updatestationsarcorr(i))
        ENDIF
        IF(itype==2)THEN
          IF(doupdate(i_tpcorr))  CALL apply_nutrientcorr(updatetpcorr(i),clakeoutflow(i_sp,itype),clakeoutflow(i_pp,itype))          
          IF(doupdate(i_tncorr))  CALL apply_nutrientcorr(updatetncorr(i),clakeoutflow(i_in,itype),clakeoutflow(i_on,itype))          
        ELSEIF(itype==1)THEN  
          IF(doupdate(i_tploccorr))  CALL apply_nutrientcorr(updatetploccorr(i),clakeoutflow(i_sp,itype),clakeoutflow(i_pp,itype))          
          IF(doupdate(i_tnloccorr))  CALL apply_nutrientcorr(updatetnloccorr(i),clakeoutflow(i_in,itype),clakeoutflow(i_on,itype))          
        ENDIF
        IF(itype==2.AND.doupdate(i_wendupd))THEN
          lakeareatemp = classbasin(i,slc_olake)%part*basin(i)%area
          IF(ALLOCATED(lakestate%slowwater).AND.lakedatapar(lakedataparindex(i,itype),m_lddeeplake)>0.)THEN
            CALL apply_wendupd(i,itype,lakeareatemp,lakewst(itype),lakestate,slowlakeini(itype,i))
          ELSEIF(ALLOCATED(lakestate%slowwater))THEN
            CALL apply_wendupd(i,itype,lakeareatemp,lakewst(itype),lakestate)
          ELSE
            CALL apply_wendupd(i,itype,lakeareatemp,lakewst(itype),lakestate)
          ENDIF
        ENDIF
         
      ENDDO routingloop
        
      !Calculate division of discharge in branches
      CALL calculate_branched_flow(i,lakeoutflow(2),mainflow,branchflow)
      IF(conductload) Lbranch = branchflow * clakeoutflow(:,2) * seconds_per_timestep * 1.E-3  !Load in branch (kg/timestep)
      IF(conductwb)THEN
        wbflows(w_oltomb,i) = mainflow * seconds_per_timestep   !m3/timestep
        wbflows(w_oltob,i)  = branchflow * seconds_per_timestep
        wbflows(w_rural4,i) = nutrientloadflows(1)
        wbflows(w_pstomr,i) = nutrientloadflows(2)
        wbflows(w_abstmr,i) = nutrientloadflows(3)
        wbflows(w_abstol,i) = nutrientloadflows(4)
      ENDIF
        
      !Accumulate flow to downstream subbasins
      IF(lakeoutflow(2)>0)THEN
        IF(mainflow>0 .AND. path(i)%main>0)THEN
          accinflow(path(i)%main) = accinflow(path(i)%main) + mainflow
          acccinflow(:,path(i)%main) = acccinflow(:,path(i)%main) + clakeoutflow(:,2) * mainflow
        ENDIF
        IF(branchflow>0)THEN
          IF(branchdata(branchindex(i))%branch>0)THEN
            accinflow(branchdata(branchindex(i))%branch) = accinflow(branchdata(branchindex(i))%branch) + branchflow
            acccinflow(:,branchdata(branchindex(i))%branch) = acccinflow(:,branchdata(branchindex(i))%branch) + clakeoutflow(:,2) * branchflow
          ENDIF
        ENDIF
      ENDIF
   
      !Calculate subbasin mean over all classes
      IF(eacti>0.0 .AND. i_t1>0)   cevapi = cevapi / eacti
      
      !Calculate lake volume for print out
      lakebasinvol=0.
      lakevol=0.
      IF(slc_ilake>0)THEN
        a = classbasin(i,slc_ilake)%part
        IF(a>0)THEN
          lakearea = a * basin(i)%area
          IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
            CALL calculate_lake_volume(1,i,lakearea,lakestate%water(1,i)+lakestate%slowwater(1,i),lakebasinvol,lakevol,lakevolsum)
          ELSE
            CALL calculate_lake_volume(1,i,lakearea,lakestate%water(1,i),lakebasinvol,lakevol,lakevolsum)
          ENDIF
        ENDIF
      ENDIF
      IF(slc_olake>0)THEN
        a = classbasin(i,slc_olake)%part
        IF(a>0)THEN
          lakearea = a * basin(i)%area
          IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
            CALL calculate_lake_volume(2,i,lakearea,lakestate%water(2,i)+lakestate%slowwater(2,i),lakebasinvol,lakevol,lakevolsum)
          ELSE
            CALL calculate_lake_volume(2,i,lakearea,lakestate%water(2,i),lakebasinvol,lakevol,lakevolsum)
          ENDIF
        ENDIF
      ENDIF
        
      !>Set subbasin output variables (outvar) for possible print out
      flow1000m3ts=lakeoutflow(2)*seconds_per_timestep*1.E-3  !m3/s -> 1000m3/ts
      outvar(i,1)=runoffi
      IF(xobsindex(2,i)>0) outvar(i,2)=xobsi(xobsindex(2,i))   !local recorded runoff, mm/d 
      IF(i_t1>0) outvar(i,5)=crunoffi(i_t1)
      IF(i_t2>0) outvar(i,6)=crunoffi(i_t2)
      IF(i_t3>0) outvar(i,7)=crunoffi(i_t3)
      IF(i_t1>0.AND.xobsindex(o_cprecT1,i)>0) outvar(i,o_cprecT1)=xobsi(xobsindex(o_cprecT1,i))
      IF(i_t1>0) outvar(i,9)=cevapi
      outvar(i,10)=soili
      IF(i_t1>0) outvar(i,11)=csoili(i_t1)
      IF(i_t2>0) outvar(i,12)=csoili(i_t2)
      IF(i_t3>0) outvar(i,13)=csoili(i_t3)
      outvar(i,14)=snowi
      outvar(i,15)=eacti      
      outvar(i,195)=landeacti
      
!      outvar(i,206)=outvar(i,206) + icpevapi   !interception losses [mm] (preccorr and pcluse)
      outvar(i,206) = icpevapi                  !interception evaporation [mm] now calculated as pcluse (or sncluse and rncluse) and separated from preccorr(outvarid 323)

      outvar(i,max_outvar_trunk+1) = outvar(i,max_outvar_trunk+1) + pcorri    !precipitation corrections [mm]
      
      outvar(i,278)=evapsnowi                   !evaporatio from snow, separated from glacier evporation [m]
      
!      outvar(i,279)=outvar(i,206) + eacti  !total evaporation (incl interception losses) [mm]
      outvar(i,279)= eacti  !total evaporation (incl interception losses) [mm]  - this is now same as EVAP (outvarid 15)

      IF(xobsindex(16,i)>0) outvar(i,16)=xobsi(xobsindex(16,i))   !T1 in runoff
      IF(xobsindex(17,i)>0) outvar(i,17)=xobsi(xobsindex(17,i))   !T2 in runoff
      IF(xobsindex(18,i)>0) outvar(i,18)=xobsi(xobsindex(18,i))   !T3 in runoff
      outvar(i,19)=gwati                            !groundwater level [m]
      IF(i_in>0) outvar(i,20)=crunoffi(i_in)*1000.  !IN conc in runoff [ug/L]
      IF(i_on>0) outvar(i,21)=crunoffi(i_on)*1000.  !ON conc in runoff [ug/L]
      IF(i_sp>0) outvar(i,22)=crunoffi(i_sp)*1000.  !SRP conc in runoff [ug/L]
      IF(i_pp>0) outvar(i,23)=crunoffi(i_pp)*1000.  !PP conc in runoff [ug/L]
      IF(xobsindex(24,i)>0)THEN
        outvar(i,24)=xobsi(xobsindex(24,i))   !rec IN in outflow, ug/L
        IF(xobsi(xobsindex(24,i)).NE.missing_value) outvar(i,207)=xobsi(xobsindex(24,i))*flow1000m3ts*1.E-3  !rec IN load in outflow, kg/ts
      ENDIF  
      IF(xobsindex(25,i)>0)THEN
        outvar(i,25)=xobsi(xobsindex(25,i))   !rec ON in outflow, ug/L
        IF(xobsi(xobsindex(25,i)).NE.missing_value) outvar(i,208)=xobsi(xobsindex(25,i))*flow1000m3ts*1.E-3  !rec ON load in outflow, kg/ts
      ENDIF
      IF(xobsindex(26,i)>0)THEN
        outvar(i,26)=xobsi(xobsindex(26,i))   !rec SP in outflow, ug/L
        IF(xobsi(xobsindex(26,i)).NE.missing_value) outvar(i,209)=xobsi(xobsindex(26,i))*flow1000m3ts*1.E-3  !rec SP load in outflow, kg/ts
      ENDIF
      IF(xobsindex(27,i)>0)THEN
        outvar(i,27)=xobsi(xobsindex(27,i))   !rec PP in outflow, ug/L
        IF(xobsi(xobsindex(27,i)).NE.missing_value) outvar(i,210)=xobsi(xobsindex(27,i))*flow1000m3ts*1.E-3  !rec PP load in outflow, kg/ts
      ENDIF
      IF(xobsindex(28,i)>0)THEN
        outvar(i,28)=xobsi(xobsindex(28,i))   !rec TN in outflow, ug/L
        IF(xobsi(xobsindex(28,i)).NE.missing_value) outvar(i,211)=xobsi(xobsindex(28,i))*flow1000m3ts*1.E-3  !rec TN load in outflow, kg/TS
      ENDIF
      IF(xobsindex(29,i)>0)THEN
        outvar(i,29)=xobsi(xobsindex(29,i))   !rec TP in outflow, ug/L 
        IF(xobsi(xobsindex(29,i)).NE.missing_value) outvar(i,212)=xobsi(xobsindex(29,i))*flow1000m3ts*1.E-3  !rec TP load in outflow, kg/ts
      ENDIF

      IF(ALLOCATED(qobsi)) THEN
        IF(qobsi(i)/=-9999) THEN
          outvar(i,30) = outflowtemp-qobsi(i)                 !Daily error in Q
        ENDIF
      ENDIF
      outvar(i,31)=outflowtemp
      outvar(i,32) = riverstate%temp(2,i)                     !Water temperature (out of subbasin)
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0) outvar(i,32) = lakestate%temp(2,i)
      ENDIF
      outvar(i,283) = outvar(i,32)                            !Water temperature limited to zero.
      IF(outvar(i,283)<0.) outvar(i,283) = 0.
 
      !Lake and River Ice, Snow, and T2 Water Temperature Model<<<
      IF(modeloption(p_lakeriverice)>0)THEN
        !Lake ice variables
        IF(slc_olake>0)THEN
          outvar(i,148) = frozenstate%lakeice(2,i)               !'coli' 'comp olake ice depth'
          outvar(i,150) = frozenstate%lakebice(2,i)              !'colb','comp olake blackice depth'
          outvar(i,152) = lakesnowdepth(2)            !'cols','comp olake snow depth'
          IF(xobsindex(154,i)>0)THEN
            outvar(i,154) = xobsi(xobsindex(154,i))  !'roli','rec. olake ice depth'
          ENDIF
          IF(xobsindex(156,i)>0)THEN
            outvar(i,156) = xobsi(xobsindex(156,i))  !'rolb','rec. olake blackice depth'
          ENDIF
          IF(xobsindex(158,i)>0)THEN
            outvar(i,158) = xobsi(xobsindex(158,i))  !'rols','rec. olake snow depth'
          ENDIF
        ENDIF
        IF(slc_ilake>0)THEN
          outvar(i,149) = frozenstate%lakeice(1,i)               !'cili' 'comp ilake ice depth'
          outvar(i,151) = frozenstate%lakebice(1,i)              !'cilb','comp ilake blackice depth'
          outvar(i,153) = lakesnowdepth(1)            !'cils','comp ilake snow depth'
          IF(xobsindex(155,i)>0)THEN
            outvar(i,155) = xobsi(xobsindex(155,i))  !'rili','rec. ilake ice depth'
          ENDIF
          IF(xobsindex(157,i)>0)THEN
            outvar(i,157) = xobsi(xobsindex(157,i))  !'rilb','rec. ilake blackice depth'
          ENDIF
          IF(xobsindex(159,i)>0)THEN
            outvar(i,159) = xobsi(xobsindex(159,i))  !'rils','rec. ilake snow depth'
          ENDIF
        ENDIF
        !River ice and snow depth variables
        outvar(i,160) = frozenstate%riverice(2,i) !'cmri','comp main river ice depth'
        outvar(i,161) = frozenstate%riverice(1,i) !'clri','comp local river ice depth' 
        outvar(i,162) = frozenstate%riverbice(2,i) !'cmrb','comp main river blackice depth'
        outvar(i,163) = frozenstate%riverbice(1,i) !'clrb','comp local river blackice depth'
        outvar(i,164) = riversnowdepth(2) !'cmrs','comp main river snow depth' 
        outvar(i,165) = riversnowdepth(1) !'clrs','comp local river snow depth'
        outvar(i,248) = frozenstate%lakeicecov(2,i) !'coic','comp olake ice cover'
        outvar(i,249) = frozenstate%lakeicecov(1,i) !'ciic','comp ilake ice cover' 
        outvar(i,250) = frozenstate%rivericecov(2,i) !'cmic','comp main river ice cover'
        outvar(i,251) = frozenstate%rivericecov(1,i) !'clic','comp local river ice cover'
        IF(xobsindex(166,i)>0)THEN
          outvar(i,166) = xobsi(xobsindex(166,i))  !!'rmri','rec. main river ice depth'
        ENDIF
        IF(xobsindex(167,i)>0)THEN
          outvar(i,167) = xobsi(xobsindex(167,i))  !'rlri','rec. local river ice depth'
        ENDIF
        IF(xobsindex(168,i)>0)THEN
          outvar(i,168) = xobsi(xobsindex(168,i))  !'rmrb','rec. main river blackice depth'
        ENDIF
        IF(xobsindex(169,i)>0)THEN
          outvar(i,169) = xobsi(xobsindex(169,i))  !'rlrb','rec. local river blackice depth'
        ENDIF
        IF(xobsindex(170,i)>0)THEN
          outvar(i,170) = xobsi(xobsindex(170,i))  !'rmrs','rec. main river snow depth'
        ENDIF
        IF(xobsindex(171,i)>0)THEN
          outvar(i,171) = xobsi(xobsindex(171,i))  !'rlrs','rec. local river snow depth'
        ENDIF
      ENDIF
      IF(i_t2>0)THEN
      ! Lake temperature variables (surface, upper layer, lower layer, mean water temp)
        IF(slc_olake>0)THEN
          outvar(i,172) = lakesurftemp(2)     !'olst','comp olake surface temp'  
          outvar(i,173) = lakestate%uppertemp(2,i)  !'olut','comp olake upper temp'  
          outvar(i,174) = lakestate%lowertemp(2,i)   !'ollt','comp olake lower temp' 
          outvar(i,175) = lakestate%concslow(i_t2,2,i) !'olwt','comp olake mean  temp'
        ENDIF
        IF(slc_ilake>0)THEN
          outvar(i,176) = lakesurftemp(1)     !'ilst','comp olake surface temp'   
          outvar(i,177) = lakestate%concslow(i_t2,1,i) !'ilwt','comp olake mean  temp'  
        ENDIF
      ! River temperature variables (surface, mean, local and main)
        outvar(i,178) = riversurftemp(1)       !lrst','comp local river surface temp'
        CALL get_rivertempvol(i,1,riverstate,meanrivertemp,totrivervol)
        outvar(i,179) = meanrivertemp  !'lrwt','comp local river mean  temp'
        outvar(i,180) = riversurftemp(2)       !'mrst','comp main  river surface temp'
        CALL get_rivertempvol(i,2,riverstate,meanrivertemp,totrivervol)
        outvar(i,181) = meanrivertemp
      ENDIF
      !Additional variables from old water temp model
      outvar(i,185)= riverstate%temp(2,i) !main river temperature
      outvar(i,186)= riverstate%temp(1,i) !local river temp
      IF(slc_olake>0)THEN
        outvar(i,187)= lakestate%temp(2,i) !olake temp
      ENDIF
      IF(slc_ilake>0)THEN
        outvar(i,188)= lakestate%temp(1,i) !ilake temp
      ENDIF
      !Recorded Olake water surface temp
      IF(xobsindex(182,i)>0) outvar(i,182) = xobsi(xobsindex(182,i)) !'rolt','rec. olake surface temp'
      IF(xobsindex(183,i)>0) outvar(i,183) = xobsi(xobsindex(183,i)) !'rilt','rec. ilake surface temp' 
      IF(xobsindex(184,i)>0) outvar(i,184) = xobsi(xobsindex(184,i)) !'rmrt','rec. main river surface temp'
      !END of Ice and T2 temperature block<<<
      !Block with Remote Sensing Snow variables<<<
      outvar(i,189) = snowcovi  !Fractional snow cover area(-)
      IF(xobsindex(190,i)>0) outvar(i,190) = xobsi(xobsindex(190,i)) !Recorded Fractional snow cover area(-)
      outvar(i,191) = snowmaxi  !snowmax during winter(-)
      IF(xobsindex(192,i)>0) outvar(i,192) = xobsi(xobsindex(192,i)) !Recorded Fractional snow cover area error (-)
      IF(xobsindex(193,i)>0) outvar(i,193) = xobsi(xobsindex(193,i)) !Recorded Fractional snow cover multi (-)
      IF(xobsindex(194,i)>0) outvar(i,194) = xobsi(xobsindex(194,i)) !Recorded Fractional snow cover multi error (-)
      !END snow remote sensing block<<<
      
      !Output variables for olake water stage and lake volume
      IF(slc_olake>0)THEN
        lakeareatemp = classbasin(i,slc_olake)%part*basin(i)%area
        IF(lakeareatemp>0.)THEN
          CALL calculate_olake_waterstage(i,wstlaketemp,lakeareatemp,lakeareatemp2,wstlake,lakestate,w0ref)
          outvar(i,34) = wstlake + w0ref                               !olake water stage before updating (m)
          CALL calculate_olake_waterstage(i,lakewst(2),lakeareatemp,lakeareatemp2,wstlake,lakestate,w0ref)
          outvar(i,280) = wstlake                         !olake water stage (not in w-reference-system)
          outvar(i,51) = wstlake + w0ref    !olake water stage (in w-reference-system)
          CALL calculate_regamp_adjusted_waterstage(i,lakeareatemp,wstlake,wstlakeadj)
          IF(wstlakeadj/=missing_value) outvar(i,323) = wstlakeadj    !regamp adjusted olake water stage (not in w-reference-system)
          IF(wstlakeadj/=missing_value) outvar(i,282) = wstlakeadj + w0ref    !regamp adjusted olake water stage (in w-reference-system)
          IF(xobsindex(o_rewstr,i)>0)THEN
            IF(xobsi(xobsindex(o_rewstr,i))/=missing_value)THEN
              outvar(i,281) = xobsi(xobsindex(o_rewstr,i)) - w0ref       !recorded olake waterstage (cleaned from w0ref)
              outvar(i,33) = outvar(i,34) - xobsi(xobsindex(o_rewstr,i)) !error in W
            ENDIF
          ENDIF
          outvar(i,146) = lakebasinvol(2)/1000000. !OLake Volume Computed (10^6 m3)
          IF(lakevol .NE. missing_value) THEN !LakeBasin Volume Computed (10^6 m3)
            outvar(i,144) = lakevol/1000000. 
          ELSE
            outvar(i,144) = lakevol
          ENDIF
        ENDIF
      ENDIF
      
      IF(xobsindex(35,i)>0)THEN
        outvar(i,35)=xobsi(xobsindex(35,i))   !DOC in outflow, mg/L 
        IF(xobsi(xobsindex(35,i)).NE.missing_value) outvar(i,213)=xobsi(xobsindex(35,i))*flow1000m3ts  !rec DOC load outflow, kg/ts
      ENDIF
      IF(i_in>0) outvar(i,36)=csoili(i_in)*1000.                  !IN soil concentration, ug/L
      outvar(i,37)=soilfrosti                                     !soil frost depth
      outvar(i,38)=soiltempi                                      !soil temperature
      outvar(i,39)=snowdepthi                                     !snow depth depth
      outvar(i,40)=epoti                                          !potential evaporation used
      IF(xobsindex(o_reepot,i)>0) outvar(i,o_reepot)=xobsi(xobsindex(o_reepot,i))   !Observed potential evaporation
      IF(xobsindex(42,i)>0) outvar(i,42)=xobsi(xobsindex(42,i))   !Observed evaporation
      outvar(i,43)=corrpreci                                      !Corrected precipitation (mm)
      
!      outvar(i,284)=outvar(i,43) + outvar(i,206)                   !Simulated precipitation (mm)
      outvar(i,284)=outvar(i,43)                                !Simulated precipitation (mm) - this is now same as corrpreci since interception losses are not included in corrpreci!!

      IF(i_oc>0) outvar(i,44) = crunoffi(i_oc)                  !DOC conc in runoff, mg/L
      IF(i_oc>0) outvar(i,45) = csoili(i_oc)                    !DOC conc in soil, mg/L
      IF(i_oc>0) outvar(i,46) = clakeoutflow(i_oc,2)            !DOC conc outflow olake, mg/L
      IF(i_oc>0) outvar(i,47) = humusCp1i                       !pool humusC soil1,kg/km2
      IF(i_oc>0) outvar(i,48) = fastCp1i                        !pool fastC soil1,kg/km2
      IF(i_oc>0) outvar(i,49) = humusCp2i                       !pool humusC soil2,kg/km2
      IF(i_oc>0) outvar(i,50) = fastCp2i                        !pool fastC soil2,kg/km2
      IF(xobsindex(o_rewstr,i)>0) outvar(i,o_rewstr)=xobsi(xobsindex(o_rewstr,i))      !recorded olake waterstage
      outvar(i,53) = lakeoutflow(2)                          !computed outflow olake (m3/s)
      IF(ALLOCATED(qobsi)) outvar(i,54)=qobsi(i)             !recorded outflow olake (m3/s)
      IF(i_in>0) outvar(i,55) = clakeoutflow(i_in,2)*1000.   !comp conc IN outflow olake (ug/l)
      IF(i_on>0) outvar(i,56) = clakeoutflow(i_on,2)*1000.   !comp conc ON outflow olake (ug/l)
      IF(i_sp>0) outvar(i,57) = clakeoutflow(i_sp,2)*1000.   !comp conc SRP outflow olake (ug/l)
      IF(i_pp>0) outvar(i,58) = clakeoutflow(i_pp,2)*1000.   !comp conc PP outflow olake (ug/l)
      IF(xobsindex(59,i)>0) outvar(i,59) = xobsi(xobsindex(59,i))              !rec snow depth
      IF(xobsindex(60,i)>0) outvar(i,60) = xobsi(xobsindex(60,i))              !rec soil frost
      IF(xobsindex(61,i)>0) outvar(i,61) = xobsi(xobsindex(61,i))              !rec grw level
      IF(xobsindex(62,i)>0) outvar(i,62) = xobsi(xobsindex(62,i))              !rec T1 conc in outflow
      IF(i_t1>0) outvar(i,63) = clakeoutflow(i_t1,2)                        !comp conc T1 conc in outflow
      IF(i_t2>0) outvar(i,64) = clakeoutflow(i_t2,2)                        !comp conc T2 conc in outflow
      IF(i_t3>0) outvar(i,65) = clakeoutflow(i_t3,2)                        !comp conc T3 conc in outflow
      IF(i_in>0.AND.i_on>0) outvar(i,66)=(crunoffi(i_in) + crunoffi(i_on))*1000.  !TN conc in runoff, ug/L
      IF(i_sp>0.AND.i_pp>0) outvar(i,67)=(crunoffi(i_sp) + crunoffi(i_pp))*1000.  !TP conc in runoff, ug/L
      IF(i_in>0.AND.i_on>0) outvar(i,68)=fastNp1i                         !fastN pool soil1
      IF(i_in>0.AND.i_on>0) outvar(i,69)=fastNp2i                         !fastN pool soil2
      IF(i_in>0.AND.i_on>0) outvar(i,70)=fastNp3i                         !fastN pool soil3
      IF(i_in>0.AND.i_on>0) outvar(i,71)=humusNp1i                        !humusN pool soil1
      IF(i_in>0.AND.i_on>0) outvar(i,72)=humusNp2i                        !humusN pool soil2
      IF(i_in>0.AND.i_on>0) outvar(i,73)=humusNp3i                        !humusN pool soil3
      IF(i_in>0) outvar(i,74)=soili1*csoili1(i_in)             !inorgN pool soil1
      IF(i_in>0) outvar(i,75)=soili2*csoili2(i_in)             !inorgN pool soil2
      IF(i_in>0) outvar(i,76)=soili3*csoili3(i_in)             !inorgN pool soil3
      IF(i_in>0.AND.i_on>0) outvar(i,77)=(clakeoutflow(i_in,2) + clakeoutflow(i_on,2))*1000.  !TN conc in outflow olake, ug/L
      IF(i_sp>0.AND.i_pp>0) outvar(i,78)=(clakeoutflow(i_sp,2) + clakeoutflow(i_pp,2))*1000.  !TP conc in outflow olake, ug/L
      outvar(i,79)=runoff1i
      IF(nummaxlayers>1) outvar(i,80)=runoff2i
      IF(modeloption(p_deepgroundwater)==1)THEN
        IF(path(i)%grw1==0)THEN
          outvar(i,81)=grwout(1,i)/seconds_per_timestep*0.001          !loss of water via groundwater flow from the model system, m3/s
        ELSE
          outvar(i,81)=0.
        ENDIF
      ENDIF
      outvar(i,82)=runoffdi 
      outvar(i,83)=runoffsri
      outvar(i,84)=denitrifi
      outvar(i,85)=plantuptakei
      outvar(i,86)=nitrifi  
      outvar(i,87)=atmdepi                                     !atmospheric deposition of N on soil
      IF(i_sp>0.AND.i_pp>0) outvar(i,88)=partPp1i              !partP pool soil1
      IF(i_sp>0.AND.i_pp>0) outvar(i,89)=partPp2i              !partP pool soil2
      IF(i_sp>0.AND.i_pp>0) outvar(i,90)=partPp3i              !partP pool soil3
      IF(i_sp>0) outvar(i,91)=soili1*csoili1(i_sp)             !SRP pool soil1
      IF(i_sp>0) outvar(i,92)=soili2*csoili2(i_sp)             !SRP pool soil2
      IF(i_sp>0) outvar(i,93)=soili3*csoili3(i_sp)             !SRP pool soil3
      IF(i_sp>0.AND.i_pp>0) outvar(i,94)=fastPp1i              !fastP pool soil1
      IF(i_sp>0.AND.i_pp>0) outvar(i,95)=fastPp2i              !fastP pool soil2
      IF(i_sp>0.AND.i_pp>0) outvar(i,96)=fastPp3i              !fastP pool soil3
      IF(i_in>0.AND.xobsindex(o_cprecIN,i)>0) outvar(i,o_cprecIN)=xobsi(xobsindex(o_cprecIN,i))               !recorded concentration of IN in precipitation
      IF(i_sp>0.AND.xobsindex(o_cprecSP,i)>0) outvar(i,o_cprecSP)=xobsi(xobsindex(o_cprecSP,i))               !recorded concentration of SP in precipitation
      IF(xobsindex(o_reswe,i)>0) outvar(i,o_reswe)=xobsi(xobsindex(o_reswe,i))               !recorded snow water equivalent
      IF(firstoutstep) accdiff(i) = 0.  
      IF(ALLOCATED(qobsi))THEN
        IF(qobsi(i)>=0.)THEN
          accdiff(i) = accdiff(i) + (lakeoutflow(2)-qobsi(i))/upstreamarea(i)*seconds_per_timestep*1.E3              !accumulated volume error (mm)
          outvar(i,100) = accdiff(i)
        ENDIF
      ENDIF
      outvar(i,101) = lakeoutflow(1)
      IF(i_in>0) outvar(i,102) = clakeoutflow(i_in,1)*1000.
      IF(i_on>0) outvar(i,103) = clakeoutflow(i_on,1)*1000.
      IF(i_sp>0) outvar(i,104) = clakeoutflow(i_sp,1)*1000.
      IF(i_pp>0) outvar(i,105) = clakeoutflow(i_pp,1)*1000.
      IF(i_in>0.AND.i_on>0) outvar(i,106) = (clakeoutflow(i_in,1) + clakeoutflow(i_on,1))*1000.
      IF(i_sp>0.AND.i_pp>0) outvar(i,107) = (clakeoutflow(i_sp,1) + clakeoutflow(i_pp,1))*1000.
      IF(i_oc>0) outvar(i,108) = fastCp3i                         !fastC pool soil3
      IF(i_oc>0) outvar(i,109) = humusCp3i                        !humusC pool soil3
      IF(i_in>0.AND.i_on>0) outvar(i,110) = (clakeoutflow(i_in,2) + clakeoutflow(i_on,2))*flow1000m3ts  !TN load (kg/timestep)
      IF(i_sp>0.AND.i_pp>0) outvar(i,111) = (clakeoutflow(i_sp,2) + clakeoutflow(i_pp,2))*flow1000m3ts  !TP load (kg/timestep)
      outvar(i,112) = qcinfli                                          !total inflow (minus evaporation) to olake (m3/s)
      IF(xobsindex(113,i)>0) outvar(i,113) = xobsi(xobsindex(113,i))  !recorded inflow (-"-) (m3/s)
      outvar(i,114) = riverstate%water(1,i) + (SUM(riverstate%qqueue(1:ttstep(1,i),1,i)) + riverstate%qqueue(ttstep(1,i)+1,1,i) * ttpart(1,i))
      outvar(i,115) = riverstate%water(2,i) + (SUM(riverstate%qqueue(1:ttstep(2,i),2,i)) + riverstate%qqueue(ttstep(2,i)+1,2,i) * ttpart(2,i))
      outvar(i,116) = atmdep2i      !atmospheric deposition of TP on soil
      IF((i_in>0.OR.i_sp>0).AND.ALLOCATED(qobsi)) THEN
        IF(qobsi(i)/=missing_value) THEN
          IF(i_in>0.AND.xobsindex(28,i)>0)THEN
            IF(xobsi(xobsindex(28,i))/=missing_value) outvar(i,119) = xobsi(xobsindex(28,i))*qobsi(i)*seconds_per_timestep*1.E-6  !rec TN load (kg/timestep)
          ENDIF
          IF(i_sp>0.AND.xobsindex(29,i)>0)THEN
            IF(xobsi(xobsindex(29,i))/=missing_value) outvar(i,120) = xobsi(xobsindex(29,i))*qobsi(i)*seconds_per_timestep*1.E-6  !rec TP load (kg/timestep)
          ENDIF
        ENDIF
      ENDIF
      outvar(i,121) = corrtempi      !Corrected temperature (mm)
      IF(i==nsub) outvar(:,122) = irrevap(:)     !Irrigation losses (m3) Note! can change for later subbasins
      outvar(i,123) = mainflow       !Discharge in main channel (m3/s)
      outvar(i,124) = ldremi         !Abstraction local dam for irrigation (m3)
      outvar(i,125) = lrremi         !Abstraction local river for irrigation (m3)
      outvar(i,126) = irrappli*landarea(i)*1.E-3       !Applied irrigation (m3)
      outvar(i,127) = gwremi         !Ground water removed for irrigation (m3)
      outvar(i,128) = rsremi         !Abstraction regional surface water for irrigation (m3)
      outvar(i,130) = branchflow     !Discharge in branched river (m3/s)
      outvar(i,131) = soildefi       !Soil moisture deficit (mm)
      IF(i_sp>0.AND.i_pp>0) outvar(i,132)=humusPp1i              !humusP pool soil1
      IF(i_sp>0.AND.i_pp>0) outvar(i,133)=humusPp2i              !humusP pool soil2
      IF(i_sp>0.AND.i_pp>0) outvar(i,134)=humusPp3i              !humusP pool soil3
      IF(i_oc>0) outvar(i,135) = clakeoutflow(i_oc,2)*flow1000m3ts  !TOC load (kg/timestep)
      IF(i_on>0) outvar(i,136)=soili1*csoili1(i_on)             !orgN pool soil1
      IF(i_on>0) outvar(i,137)=soili2*csoili2(i_on)             !orgN pool soil2
      IF(i_on>0) outvar(i,138)=soili3*csoili3(i_on)             !orgN pool soil3
      outvar(i,139)=soilltmpi(1)                                !soil layer 1 temperature
      IF(nummaxlayers>1) outvar(i,140)=soilltmpi(2)             !soil layer 2 temperature
      IF(nummaxlayers>2) outvar(i,141)=soilltmpi(3)             !soil layer 3 temperature
      outvar(i,142)=rainfalli             !Part of precipitation (corrected) that falls as rain
      outvar(i,143)=snowfalli             !Part of precipitation (corrected) that falls as snow
      IF(slc_ilake>0)THEN
        outvar(i,145) = lakebasinvol(1)/1000000. !ILake Volume Computed (10^6 m3)
        ELSE
        outvar(i,145) = missing_value
      ENDIF
      IF(nummaxlayers>2) outvar(i,147) = runoff3i
      outvar(i,196) = soili12
      outvar(i,197) = soili1 - standsoili
      outvar(i,198) = soili2
      outvar(i,199) = soili3
      outvar(i,200) = standsoili
      outvar(i,201) = soili12 - standsoili
      outvar(i,202) = soili - standsoili
      IF(i_oc>0) outvar(i,203) = clakeoutflow(i_oc,1)   !orgC local flow
      IF(riverarea(1)>0) outvar(i,204) = outvar(i,114)/riverarea(1)   !local river depth [m]
      IF(riverarea(2)>0) outvar(i,205) = outvar(i,115)/riverarea(2)   !main river depth [m]
      outvar(i,214) = srffi
      outvar(i,215) = smfdi
      outvar(i,216) = srfdi
      outvar(i,217) = smfpi
      outvar(i,218) = srfpi
      DO k=0,9
        IF(xobsindex(o_xobsm+k,i)>0)THEN
          outvar(i,o_xobsm+k)=xobsi(xobsindex(o_xobsm+k,i))      !recorded meantype variable
        ELSEIF(num_xoms(1)>=k+1)THEN
          outvar(i,o_xobsm+k)=xom(k+1,i)
        ENDIF
        IF(xobsindex(o_xobss+k,i)>0)THEN
          outvar(i,o_xobss+k)=xobsi(xobsindex(o_xobss+k,i))      !recorded sumtype variable
        ELSEIF(num_xoms(2)>=k+1)THEN
          outvar(i,o_xobss+k)=xos(k+1,i)
        ENDIF
      ENDDO
      
      !Accumulation of upstream output variables
      upstream_prec(i) = upstream_prec(i) + corrpreci * basin(i)%area / upstreamarea(i)
      upstream_evap(i) = upstream_evap(i) + eacti * basin(i)%area / upstreamarea(i)
      upstream_epot(i) = upstream_epot(i) + epoti * basin(i)%area / upstreamarea(i)
      upstream_snow(i) = upstream_snow(i) + snowi * landarea(i) / upstreamarea(i)
      upstream_soil(i) = upstream_soil(i) + soili * landarea(i) / upstreamarea(i)
      upstream_soildef(i) = upstream_soildef(i) + soildefi * landarea(i) / upstreamarea(i)
      upstream_smfp(i) = upstream_smfp(i) + smfpi * landarea(i) / upstreamarea(i)
      upstream_snowfall(i) = upstream_snowfall(i) + snowfalli * basin(i)%area / upstreamarea(i)
      upstream_rainfall(i) = upstream_rainfall(i) + rainfalli * basin(i)%area / upstreamarea(i)
      !For CryoProcesses
      upstream_netrad(i) = upstream_netrad(i) + netradi * basin(i)%area / upstreamarea(i)
      upstream_swrad(i) = upstream_swrad(i) + swradi * basin(i)%area / upstreamarea(i)
      upstream_temp(i) = upstream_temp(i) + corrtempi * basin(i)%area / upstreamarea(i)
      upstream_river(i) = upstream_river(i) + 1000. * (outvar(i,114)+outvar(i,115)) / upstreamarea(i)
      upstream_lake(i) = upstream_lake(i) + 1000. * (lakebasinvol(1)+lakebasinvol(2)) / upstreamarea(i)
      upstream_glacier(i) = upstream_glacier(i) + 1000. * frozenstate%glacvol(i) * genpar(m_glacdens) / upstreamarea(i)
      upstream_snowcov(i)  = upstream_snowcov(i)   + snowcovi  * landarea(i) / upstreamarea(i)
      upstream_snowmelt(i) = upstream_snowmelt(i)  + melti     * landarea(i) / upstreamarea(i)
      upstream_glacmelt(i) = upstream_glacmelt(i)  + glacmelti * landarea(i) / upstreamarea(i)
      upstream_evapsoil(i) = upstream_evapsoil(i)  + evapsoili * basin(i)%area / upstreamarea(i)
      upstream_evapsnow(i) = upstream_evapsnow(i)  + evapsnowi * basin(i)%area / upstreamarea(i)
      upstream_evapglac(i) = upstream_evapglac(i)  + evapglaci * basin(i)%area / upstreamarea(i)
      upstream_evaplake(i) = upstream_evaplake(i)  + evaplakei * basin(i)%area / upstreamarea(i)
      upstream_evapriver(i)= upstream_evapriver(i) + evapriveri* basin(i)%area / upstreamarea(i)
      upstream_runoff(i)   = upstream_runoff(i)    + runoffi   * landarea(i) / upstreamarea(i)
      !upstream rempote sensing observed and corresponding modelled snow and snowcover for evaluation:
      !  (accumulate only if observation is available, but keep track of the fraction of observed area)
      ! a. snow water equivalent
      IF(xobsindex(o_reswe,i).GT.0)THEN
        IF(xobsi(xobsindex(o_reswe,i))/=missing_value)THEN
          upstream_compsnow(i)    = upstream_compsnow(i)    + snowi * landarea(i) / upstreamarea(i)
          upstream_recsnow(i)     = upstream_recsnow(i)     + xobsi(xobsindex(o_reswe,i)) * landarea(i) / upstreamarea(i)
!          upstream_recsnowarea(i) = upstream_recsnowarea(i) + landarea(i) / upstreamarea(i)
!        ELSE
!          upstream_recsnowmissing(i) = upstream_recsnowmissing(i) + landarea(i) / upstreamarea(i)
        ENDIF
!      ELSE
!        upstream_recsnowmissing(i) = upstream_recsnowmissing(i) + landarea(i) / upstreamarea(i)
      ENDIF
      !b. fractional snow cover area
!      IF(xobsindex(190,i).GT.0)THEN
!        IF(xobsi(xobsindex(190,i))/=missing_value)THEN
!          upstream_compfsc(i)     = upstream_compfsc(i)     + snowcovi * landarea(i) / upstreamarea(i)
!          upstream_recfsc(i)      = upstream_recfsc(i)      + xobsi(xobsindex(190,i)) * landarea(i) / upstreamarea(i)
!          upstream_recfscarea(i)  = upstream_recfscarea(i)  + landarea(i) / upstreamarea(i)
!          upstream_recfscerror(i) = upstream_recfscerror(i) + xobsi(xobsindex(192,i)) * landarea(i) / upstreamarea(i)
!        ELSE
!          upstream_recfscmissing(i)  = upstream_recfscmissing(i)  + landarea(i) / upstreamarea(i)
!        ENDIF
!      ELSE
!        upstream_recfscmissing(i)  = upstream_recfscmissing(i)  + landarea(i) / upstreamarea(i)
!      ENDIF
      
      !accumulate upstream variables on downstream MAIN
      IF(path(i)%main>0)THEN
        divasum4 = upstreamarea(i) / upstreamarea(path(i)%main)
        IF(ALLOCATED(branchdata))THEN
          IF(branchindex(i)>0)THEN
            IF(branchdata(branchindex(i))%mainpart<1.) divasum4 = branchdata(branchindex(i))%mainpart * divasum4  !OBS: No consideration of varying branchflow fraction
          ENDIF
        ENDIF
        upstream_prec(path(i)%main) = upstream_prec(path(i)%main) + upstream_prec(i) * divasum4
        upstream_evap(path(i)%main) = upstream_evap(path(i)%main) + upstream_evap(i) * divasum4
        upstream_epot(path(i)%main) = upstream_epot(path(i)%main) + upstream_epot(i) * divasum4
        upstream_snow(path(i)%main) = upstream_snow(path(i)%main) + upstream_snow(i) * divasum4
        upstream_soil(path(i)%main) = upstream_soil(path(i)%main) + upstream_soil(i) * divasum4
        upstream_soildef(path(i)%main) = upstream_soildef(path(i)%main) + upstream_soildef(i) * divasum4
        upstream_smfp(path(i)%main) = upstream_smfp(path(i)%main) + upstream_smfp(i) * divasum4
        upstream_snowfall(path(i)%main) = upstream_snowfall(path(i)%main) + upstream_snowfall(i) * divasum4 ! DG 2014113
        upstream_rainfall(path(i)%main) = upstream_rainfall(path(i)%main) + upstream_rainfall(i) * divasum4 ! DG 2014113
        !For CryoProcesses
        upstream_netrad(path(i)%main) = upstream_netrad(path(i)%main) + upstream_netrad(i) * divasum4
        upstream_swrad(path(i)%main) = upstream_swrad( path(i)%main) + upstream_swrad(i) * divasum4
        upstream_temp(path(i)%main) = upstream_temp(path(i)%main) + upstream_temp(i) * divasum4
        upstream_river(path(i)%main) = upstream_river(path(i)%main) + upstream_river(i) * divasum4 
        upstream_lake(path(i)%main) = upstream_lake(path(i)%main) + upstream_lake(i) * divasum4 
        upstream_glacier(path(i)%main) = upstream_glacier(path(i)%main) + upstream_glacier(i) * divasum4 
        upstream_snowcov(path(i)%main)  = upstream_snowcov(path(i)%main)  + upstream_snowcov(i)  * divasum4 
        upstream_snowmelt(path(i)%main) = upstream_snowmelt(path(i)%main) + upstream_snowmelt(i) * divasum4 
        upstream_glacmelt(path(i)%main) = upstream_glacmelt(path(i)%main) + upstream_glacmelt(i) * divasum4 
        upstream_evapsoil(path(i)%main) = upstream_evapsoil(path(i)%main) + upstream_evapsoil(i) * divasum4
        upstream_evapsnow(path(i)%main) = upstream_evapsnow(path(i)%main) + upstream_evapsnow(i) * divasum4
        upstream_evapglac(path(i)%main) = upstream_evapglac(path(i)%main) + upstream_evapglac(i) * divasum4
        upstream_evaplake(path(i)%main) = upstream_evaplake(path(i)%main) + upstream_evaplake(i) * divasum4
        upstream_evapriver(path(i)%main)= upstream_evapriver(path(i)%main)+ upstream_evapriver(i)* divasum4
        upstream_runoff(path(i)%main)   = upstream_runoff(path(i)%main)   + upstream_runoff(i)   * divasum4
        IF(xobsindex(o_reswe,i).GT.0)THEN
          IF(xobsi(xobsindex(o_reswe,i))/=missing_value)THEN
            upstream_compsnow(path(i)%main)    = upstream_compsnow(path(i)%main)    + upstream_compsnow(i) * divasum4
            upstream_recsnow(path(i)%main)     = upstream_recsnow(path(i)%main)     + upstream_recsnow(i) * divasum4
!            upstream_recsnowarea(path(i)%main) = upstream_recsnowarea(path(i)%main) + upstream_recsnowarea(i) * divasum4
!          ELSE
!            upstream_recsnowmissing(path(i)%main) = upstream_recsnowmissing(path(i)%main) + upstream_recsnowmissing(i) * divasum4
          ENDIF
!        ELSE
!          upstream_recsnowmissing(path(i)%main) = upstream_recsnowmissing(path(i)%main) + upstream_recsnowmissing(i) * divasum4
        ENDIF
!        IF(xobsindex(190,i).GT.0)THEN
!          IF(xobsi(xobsindex(190,i))/=missing_value)THEN
!            upstream_compfsc(path(i)%main)     = upstream_compfsc(path(i)%main)     + upstream_compfsc(i) * divasum4
!            upstream_recfsc(path(i)%main)      = upstream_recfsc(path(i)%main)      + upstream_recfsc(i) * divasum4
!            upstream_recfscarea(path(i)%main)  = upstream_recfscarea(path(i)%main)  + upstream_recfscarea(i) * divasum4
!            upstream_recfscerror(path(i)%main) = upstream_recfscerror(path(i)%main) + upstream_recfscerror(i) * divasum4
!          ELSE
!            upstream_recfscmissing(path(i)%main)  = upstream_recfscmissing(path(i)%main)  + upstream_recfscmissing(i) * divasum4
!          ENDIF
!        ELSE
!          upstream_recfscmissing(path(i)%main)  = upstream_recfscmissing(path(i)%main)  + upstream_recfscmissing(i) * divasum4
!        ENDIF
      ENDIF
      !accumulate upstream variables on downstream BRANCH
      IF(ALLOCATED(branchdata))THEN
        IF(branchindex(i)>0)THEN
          IF(branchdata(branchindex(i))%branch>0)THEN
            divasum5 = (1.- MIN(1.,branchdata(branchindex(i))%mainpart)) * upstreamarea(i) / upstreamarea(branchdata(branchindex(i))%branch)  !OBS: No consideration of varying branchflow fraction
            IF(divasum5>0.)THEN
              upstream_prec(branchdata(branchindex(i))%branch) = upstream_prec(branchdata(branchindex(i))%branch) + upstream_prec(i) * divasum5
              upstream_evap(branchdata(branchindex(i))%branch) = upstream_evap(branchdata(branchindex(i))%branch) + upstream_evap(i) * divasum5
              upstream_epot(branchdata(branchindex(i))%branch) = upstream_epot(branchdata(branchindex(i))%branch) + upstream_epot(i) * divasum5
              upstream_snow(branchdata(branchindex(i))%branch) = upstream_snow(branchdata(branchindex(i))%branch) + upstream_snow(i) * divasum5
              upstream_soil(branchdata(branchindex(i))%branch) = upstream_soil(branchdata(branchindex(i))%branch) + upstream_soil(i) * divasum5
              upstream_soildef(branchdata(branchindex(i))%branch) = upstream_soildef(branchdata(branchindex(i))%branch) + upstream_soildef(i) * divasum5
              upstream_smfp(branchdata(branchindex(i))%branch) = upstream_smfp(branchdata(branchindex(i))%branch) + upstream_smfp(i) * divasum5
              upstream_snowfall(branchdata(branchindex(i))%branch) = upstream_snowfall(branchdata(branchindex(i))%branch) + upstream_snowfall(i) * divasum5
              upstream_rainfall(branchdata(branchindex(i))%branch) = upstream_rainfall(branchdata(branchindex(i))%branch) + upstream_rainfall(i) * divasum5
              !For CryoProcesses
              upstream_netrad(branchdata(branchindex(i))%branch) = upstream_netrad(branchdata(branchindex(i))%branch) + upstream_netrad(i) * divasum5
              upstream_swrad(branchdata(branchindex(i))%branch) = upstream_swrad(branchdata(branchindex(i))%branch) + upstream_swrad(i) * divasum5
              upstream_temp(branchdata(branchindex(i))%branch) = upstream_temp(branchdata(branchindex(i))%branch) + upstream_temp(i) * divasum5
              upstream_river(branchdata(branchindex(i))%branch) = upstream_river(branchdata(branchindex(i))%branch) + upstream_river(i) * divasum5
              upstream_lake(branchdata(branchindex(i))%branch) = upstream_lake(branchdata(branchindex(i))%branch) + upstream_lake(i) * divasum5
              upstream_glacier(branchdata(branchindex(i))%branch) = upstream_glacier(branchdata(branchindex(i))%branch) + upstream_glacier(i) * divasum5

              upstream_snowcov(branchdata(branchindex(i))%branch)  = upstream_snowcov(branchdata(branchindex(i))%branch)  + upstream_snowcov(i)  * divasum5
              upstream_snowmelt(branchdata(branchindex(i))%branch) = upstream_snowmelt(branchdata(branchindex(i))%branch) + upstream_snowmelt(i) * divasum5
              upstream_glacmelt(branchdata(branchindex(i))%branch) = upstream_glacmelt(branchdata(branchindex(i))%branch) + upstream_glacmelt(i) * divasum5
              upstream_evapsoil(branchdata(branchindex(i))%branch) = upstream_evapsoil(branchdata(branchindex(i))%branch) + upstream_evapsoil(i) * divasum5
              upstream_evapsnow(branchdata(branchindex(i))%branch) = upstream_evapsnow(branchdata(branchindex(i))%branch) + upstream_evapsnow(i) * divasum5
              upstream_evapglac(branchdata(branchindex(i))%branch) = upstream_evapglac(branchdata(branchindex(i))%branch) + upstream_evapglac(i) * divasum5
              upstream_evaplake(branchdata(branchindex(i))%branch) = upstream_evaplake(branchdata(branchindex(i))%branch) + upstream_evaplake(i) * divasum5
              upstream_evapriver(branchdata(branchindex(i))%branch)= upstream_evapriver(branchdata(branchindex(i))%branch)+ upstream_evapriver(i)* divasum5
              upstream_runoff(branchdata(branchindex(i))%branch)   = upstream_runoff(branchdata(branchindex(i))%branch)   + upstream_runoff(i)   * divasum5

              IF(xobsindex(o_reswe,i).GT.0)THEN
                IF(xobsi(xobsindex(o_reswe,i))/=missing_value)THEN
                  upstream_compsnow(branchdata(branchindex(i))%branch)    = upstream_compsnow(branchdata(branchindex(i))%branch)    + upstream_compsnow(i) * divasum5
                  upstream_recsnow(branchdata(branchindex(i))%branch)     = upstream_recsnow(branchdata(branchindex(i))%branch)     + upstream_recsnow(i) * divasum5
!                  upstream_recsnowarea(branchdata(branchindex(i))%branch) = upstream_recsnowarea(branchdata(branchindex(i))%branch) + upstream_recsnowarea(i) * divasum5
!                ELSE
!                  upstream_recsnowmissing(branchdata(branchindex(i))%branch) = upstream_recsnowmissing(branchdata(branchindex(i))%branch) + upstream_recsnowmissing(i) * divasum5
                ENDIF
!              ELSE
!                upstream_recsnowmissing(branchdata(branchindex(i))%branch) = upstream_recsnowmissing(branchdata(branchindex(i))%branch) + upstream_recsnowmissing(i) * divasum5
              ENDIF
!              IF(xobsindex(190,i).GT.0)THEN
!                IF(xobsi(xobsindex(190,i))/=missing_value)THEN
!                  upstream_compfsc(branchdata(branchindex(i))%branch)     = upstream_compfsc(branchdata(branchindex(i))%branch)     + upstream_compfsc(i) * divasum5
!                  upstream_recfsc(branchdata(branchindex(i))%branch)      = upstream_recfsc(branchdata(branchindex(i))%branch)      + upstream_recfsc(i) * divasum5
!                  upstream_recfscarea(branchdata(branchindex(i))%branch)  = upstream_recfscarea(branchdata(branchindex(i))%branch)  + upstream_recfscarea(i) * divasum5
!                  upstream_recfscerror(branchdata(branchindex(i))%branch) = upstream_recfscerror(branchdata(branchindex(i))%branch) + upstream_recfscerror(i) * divasum5
!                ELSE
!                  upstream_recfscmissing(branchdata(branchindex(i))%branch)  = upstream_recfscmissing(branchdata(branchindex(i))%branch)  + upstream_recfscmissing(i) * divasum5
!                ENDIF
!              ELSE
!                upstream_recfscmissing(branchdata(branchindex(i))%branch)  = upstream_recfscmissing(branchdata(branchindex(i))%branch)  + upstream_recfscmissing(i) * divasum5
!              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      outvar(i,243) = upstream_prec(i)
      outvar(i,244) = upstream_evap(i)
      outvar(i,277) = upstream_epot(i)
      outvar(i,245) = upstream_snow(i)
      outvar(i,246) = upstream_soil(i)
      outvar(i,129) = upstream_soildef(i)
      outvar(i,219) = upstream_smfp(i)
      outvar(i,247) = 1000.*(lakeoutflow(2)*seconds_per_timestep) / upstreamarea(i) !Specific discharge (upstream runoff)
      outvar(i,275) = upstream_snowfall(i)
      outvar(i,276) = upstream_rainfall(i)
          
      !Set subbasin load output variables for possible printout for Source Apportionment program
      IF(conductload) THEN
        outvar_classload(1:nclass,1:2,:,i)  = Latmdep(1:nclass,1:2,:)
        outvar_classload(1:nclass,3:4,:,i)  = Lcultiv(1:nclass,1:2,:)
        outvar_classload(1:nclass,5,:,i)    = Lrurala(1:nclass,:)
        outvar_classload(1:nclass,6,:,i)    = Lgrwsoil(1:nclass,:)
        outvar_classload(1:nclass,7,:,i)    = Lirrsoil(1:nclass,:)
        outvar_classload(1:nclass,8,:,i)    = Lstream(1:nclass,:)
        outvar_basinload(:,1,i)             = Lruralb(:)
        outvar_basinload(:,2:4,i)           = Lpoints(:,1:3)     !Note: Load-files and SourceApp need to change if max_pstype is increased
        outvar_basinload(:,5,i)             = Lgrwmr(:)
        outvar_basinload(:,6,i)             = Lgrwol(:)
        outvar_basinload(:,7:19,i)          = Lpathway(:,1:13)
        outvar_basinload(:,20,i)            = Lbranch(:)
      ENDIF
      
      !Set water balance flows and stores for print out (current subbasin)
      IF(conductwb)THEN
        wbstores(w_ilake,i)  = lakebasinvol(1) !m3
        wbstores(w_olake,i)  = lakebasinvol(2) !m3
      ENDIF
        
      !FSUS snow outputs, forest open land
      IF(xobsindex(257,i)>0) outvar(i,257) = xobsi(xobsindex(257,i)) !'S105','fsusnow fsc surr. open'
      IF(xobsindex(258,i)>0) outvar(i,258) = xobsi(xobsindex(258,i)) !'S106','fsusnow fsc course open'
      IF(xobsindex(259,i)>0) outvar(i,259) = xobsi(xobsindex(259,i)) !'S108','fsusnow mean depth open'
      IF(xobsindex(260,i)>0) outvar(i,260) = xobsi(xobsindex(260,i)) !'S111','fsusnow mean density open'
      IF(xobsindex(261,i)>0) outvar(i,261) = xobsi(xobsindex(261,i)) !'S114','fsusnow snow water eq. open'
      IF(xobsindex(262,i)>0) outvar(i,262) = xobsi(xobsindex(262,i)) !'S205','fsusnow fsc surr. forest'
      IF(xobsindex(263,i)>0) outvar(i,263) = xobsi(xobsindex(263,i)) !'S206','fsusnow fsc course forest'
      IF(xobsindex(264,i)>0) outvar(i,264) = xobsi(xobsindex(264,i)) !'S208','fsusnow mean depth forest'
      IF(xobsindex(265,i)>0) outvar(i,265) = xobsi(xobsindex(265,i)) !'S211','fsusnow mean density forest'
      IF(xobsindex(266,i)>0) outvar(i,266) = xobsi(xobsindex(266,i)) !'S214','fsusnow snow water eq. forest'
      IF(area_open.GT.0)THEN
         outvar(i,267) =snowcover_open     !'C106','comp. fsc course open'
         outvar(i,268) =snowdepth_open     !'C108','comp. mean depth open'
         outvar(i,269) =snowdensity_open   !'C111','comp. mean density open'
         outvar(i,270) =swe_open           !'C114','comp. snow water eq. open'
      ENDIF
      IF(area_forest.GT.0)THEN
         outvar(i,271) =snowcover_forest   !'C206','comp. fsc course forest'
         outvar(i,272) =snowdepth_forest   !'C208','comp. mean depth forest'
         outvar(i,273) =snowdensity_forest !'C211','comp. mean density forest'
         outvar(i,274) =swe_forest         !'C214','comp. snow water eq. forest'
      ENDIF
      
      !Soil layer load
      IF(numsubstances>0.)THEN
        IF(i_in>0)THEN
          outvar(i,305) = soillgrossload(1,i_in)*basin(i)%area*1.E-6    !mm*mg/L*m2/1000000=kg
          outvar(i,306) = soillnetload(1,i_in)*basin(i)%area*1.E-6
          outvar(i,307) = soillgrossload(1,i_on)*basin(i)%area*1.E-6
          outvar(i,308) = soillnetload(1,i_on)*basin(i)%area*1.E-6
          outvar(i,309) = outvar(i,305) + outvar(i,307)
          outvar(i,310) = outvar(i,306) + outvar(i,308)
          outvar(i,291) = soillgrossload(2,i_in)*basin(i)%area*1.E-6
          outvar(i,292) = soillnetload(2,i_in)*basin(i)%area*1.E-6
          outvar(i,293) = soillgrossload(2,i_on)*basin(i)%area*1.E-6
          outvar(i,294) = soillnetload(2,i_on)*basin(i)%area*1.E-6
          outvar(i,295) = outvar(i,291) + outvar(i,293)
          outvar(i,296) = outvar(i,292) + outvar(i,294)
          outvar(i,303) = denitrif3i*basin(i)%area*1.E-6    !denitrification soil layer 3 (kg)
          outvar(i,304) = denitrif12i*basin(i)%area*1.E-6   !denitrification soil layer 1 and 2 (kg)
          outvar(i,311) = soillgrossload(3,i_in)*basin(i)%area*1.E-6
          outvar(i,312) = soillnetload(3,i_in)*basin(i)%area*1.E-6
          outvar(i,313) = soillgrossload(3,i_on)*basin(i)%area*1.E-6
          outvar(i,314) = soillnetload(3,i_on)*basin(i)%area*1.E-6
          outvar(i,315) = outvar(i,311) + outvar(i,313)
          outvar(i,316) = outvar(i,312) + outvar(i,314)
        ENDIF
        IF(i_sp>0)THEN
          outvar(i,285) = soillgrossload(1,i_sp)*basin(i)%area*1.E-6    
          outvar(i,286) = soillnetload(1,i_sp)*basin(i)%area*1.E-6    
          outvar(i,287) = soillgrossload(1,i_pp)*basin(i)%area*1.E-6    
          outvar(i,288) = soillnetload(1,i_pp)*basin(i)%area*1.E-6    
          outvar(i,289) = outvar(i,285) + outvar(i,287)
          outvar(i,290) = outvar(i,286) + outvar(i,288)
          outvar(i,297) = soillgrossload(2,i_sp)*basin(i)%area*1.E-6
          outvar(i,298) = soillnetload(2,i_sp)*basin(i)%area*1.E-6
          outvar(i,299) = soillgrossload(2,i_pp)*basin(i)%area*1.E-6
          outvar(i,300) = soillnetload(2,i_pp)*basin(i)%area*1.E-6
          outvar(i,301) = outvar(i,297) + outvar(i,299)
          outvar(i,302) = outvar(i,298) + outvar(i,300)
          outvar(i,317) = soillgrossload(3,i_sp)*basin(i)%area*1.E-6
          outvar(i,318) = soillnetload(3,i_sp)*basin(i)%area*1.E-6
          outvar(i,319) = soillgrossload(3,i_pp)*basin(i)%area*1.E-6
          outvar(i,320) = soillnetload(3,i_pp)*basin(i)%area*1.E-6
          outvar(i,321) = outvar(i,317) + outvar(i,319)
          outvar(i,322) = outvar(i,318) + outvar(i,320)
        ENDIF
      ENDIF
      
      !Additional outvars for CryoProcesses
      outvar(i,max_outvar_trunk+2) = evapglaci
      outvar(i,max_outvar_trunk+3) = evapsoili
      
      !some storages and flows per (land+ilake+localriver) area
      IF((landarea(i)+riverarea(1)+lakeareai(1)).GT.0.)THEN
        outvar(i,max_outvar_trunk+4) = snowi * landarea(i) / (landarea(i)+riverarea(1)+lakeareai(1))        !snow storage
        outvar(i,max_outvar_trunk+6) = soili * landarea(i)/(landarea(i)+riverarea(1)+lakeareai(1))          !soil storage
        outvar(i,max_outvar_trunk+8) = outvar(i,114) * 1000./(landarea(i)+riverarea(1)+lakeareai(1))        !local river storage
        outvar(i,max_outvar_trunk+11) = lakebasinvol(1) * 1000./(landarea(i)+riverarea(1)+lakeareai(1))      !internal lake storage
        outvar(i,max_outvar_trunk+14) = frozenstate%glacvol(i)*1000./(landarea(i)+riverarea(1)+lakeareai(1)) !glacier storage
        outvar(i,max_outvar_trunk+16) = lakeoutflow(1) * 1000. * seconds_per_timestep / (landarea(i)+riverarea(1)+lakeareai(1))     !internal lake outflow
      ELSE
        outvar(i,max_outvar_trunk+4)=0.
        outvar(i,max_outvar_trunk+6)=0.
        outvar(i,max_outvar_trunk+8)=0.
        outvar(i,max_outvar_trunk+11)=0.
        outvar(i,max_outvar_trunk+14)=0.
        outvar(i,max_outvar_trunk+16)=0.
      ENDIF
      !some storages per basin area
      outvar(i,max_outvar_trunk+5) = snowi * landarea(i) / basin(i)%area
      outvar(i,max_outvar_trunk+7) = soili * landarea(i) / basin(i)%area
      outvar(i,max_outvar_trunk+9) = (outvar(i,114)+outvar(i,115)) * 1000. / basin(i)%area
      outvar(i,max_outvar_trunk+12) = (lakebasinvol(1)+lakebasinvol(2)) * 1000. / basin(i)%area
      outvar(i,max_outvar_trunk+5) = frozenstate%glacvol(i) * 1000. / basin(i)%area

      outvar(i,max_outvar_trunk+10) = upstream_river(i)
      outvar(i,max_outvar_trunk+13) = upstream_lake(i)
      outvar(i,max_outvar_trunk+16) = upstream_glacier(i)
      !recorded, specific runoff per upstream area
      IF(ALLOCATED(qobsi))THEN
        IF(qobsi(i)/=missing_value)THEN
          outvar(i,max_outvar_trunk+17)=qobsi(i) * 1000. * seconds_per_timestep / upstreamarea(i)
        ELSE
          outvar(i,max_outvar_trunk+17)=missing_value
        ENDIF
      ENDIF

    !------------------------------------------------------------------------------------------------------------------
    !--Snowmelt, glaciermelt, snowevaporation, glacierevaporation, soil+vegevaporation   
    !------------------------------------------------------------------------------------------------------------------
    outvar(i,max_outvar_trunk+19) = upstream_runoff(i)
    outvar(i,max_outvar_trunk+20) = upstream_snowmelt(i)
    outvar(i,max_outvar_trunk+21) = upstream_glacmelt(i)
    outvar(i,max_outvar_trunk+22) = upstream_evapsnow(i)
    outvar(i,max_outvar_trunk+23) = upstream_evapglac(i)
    outvar(i,max_outvar_trunk+24) = upstream_evapsoil(i)
    outvar(i,max_outvar_trunk+25) = upstream_evaplake(i)
    outvar(i,max_outvar_trunk+26) = upstream_evapriver(i)
    outvar(i,max_outvar_trunk+27) = upstream_snowcov(i)
      
      !recorded and computed upstream area average snow (per land area with observations)
      outvar(i,max_outvar_trunk+28) = upstream_compsnow(i)
      outvar(i,max_outvar_trunk+29) = upstream_recsnow(i)
!      outvar(i,296) = upstream_recsnowarea(i)
!      outvar(i,297) = upstream_recsnowmissing(i)
!      outvar(i,298) = upstream_compfsc(i)
!      outvar(i,299) = upstream_recfsc(i)
!      outvar(i,300) = upstream_recfscarea(i)
!      outvar(i,301) = upstream_recfscerror(i)
!      outvar(i,302) = upstream_recfscmissing(i)
      
    !Upstream air temperature
    outvar(i,max_outvar_trunk+30) = upstream_temp(i)
    
    !Water source fractions
    IF(i_sm>0)outvar(i,max_outvar_trunk+31) = lakeoutflow(2) * clakeoutflow(i_sm,2) ! snow melt
    IF(i_gm>0)outvar(i,max_outvar_trunk+32) = lakeoutflow(2) * clakeoutflow(i_gm,2) ! glacier melt
    IF(i_rn>0)outvar(i,max_outvar_trunk+33) = lakeoutflow(2) * clakeoutflow(i_rn,2) ! rainfall
    IF(i_si>0)outvar(i,max_outvar_trunk+34) = lakeoutflow(2) * clakeoutflow(i_li,2) ! initial lake
    IF(i_li>0)outvar(i,max_outvar_trunk+35) = lakeoutflow(2) * clakeoutflow(i_ri,2) ! initial river
    IF(i_ri>0)outvar(i,max_outvar_trunk+36) = lakeoutflow(2) * clakeoutflow(i_si,2) ! initial soil
    IF(i_sr>0)outvar(i,max_outvar_trunk+37) = lakeoutflow(2) * clakeoutflow(i_sr,2) ! surface runoff
    IF(i_dr>0)outvar(i,max_outvar_trunk+38) = lakeoutflow(2) * clakeoutflow(i_dr,2) ! drainage

    !Radiation
    outvar(i,max_outvar_trunk+39) = netradi
    outvar(i,max_outvar_trunk+40) = swradi
    outvar(i,max_outvar_trunk+41) = upstream_netrad(i)
    outvar(i,max_outvar_trunk+42) = upstream_swrad(i)
    IF(xobsindex(max_outvar_trunk+43,i)>0) outvar(i,max_outvar_trunk+43) = xobsi(xobsindex(max_outvar_trunk+43,i)) ! recorded net radiation
    IF(xobsindex(max_outvar_trunk+44,i)>0) outvar(i,max_outvar_trunk+44) = xobsi(xobsindex(max_outvar_trunk+44,i)) ! recorded shortwave radiation

    !Additional soil temperatures
    outvar(i,max_outvar_trunk+45)=deeptempi                                      !deep soil temperature
    outvar(i,max_outvar_trunk+46)=soil2tempi                                     !soil temperature (mean layer 1+2, used for frost depth calculation)
   
    !>End main subbasin-loop   
    ENDDO subbasinloop
     
    !Calculate aquifer delay and outflow
    IF(modeloption(p_deepgroundwater)==2)THEN
      CALL calculate_aquifer(aquiferstate,aquiferoutflow,aqremflow,aqirrloss)
      outvar(:,240) = aqremflow
    ENDIF
    IF(naquifers>0)THEN
      outvar(:,242) = calculate_aquifer_waterlevel(nsub,naquifers,aquiferstate)
    ENDIF
     
    !Set water balance flows and stores for print out (all subbasins)
    IF(conductwb)THEN
      wbflows(w_rgrwto1,:)  = horizontalflows2(1,:)     !regional groundwater flow from this subbasin's groundwater reservoir
      wbflows(w_rgrwto2,:)  = horizontalflows2(2,:)
      wbflows(w_rgrwto3,:)  = horizontalflows2(3,:)
      IF(modeloption(p_deepgroundwater)==1) wbflows(w_rgrwtoos,:) = outvar(:,81) * seconds_per_timestep
      wbirrflows(w_wdregol,:)  = regionalirrflows(1,:)  !regional sources to this subbasin
      wbirrflows(w_wdregmr,:)  = regionalirrflows(2,:)  !regional sources to this subbasin
      wbirrflows(w_evapregol,:) = regionalirrflows(3,:)
      wbirrflows(w_evapregmr,:) = regionalirrflows(4,:)
      wbstores(w_iriver,:) = outvar(:,114)
      wbstores(w_mriver,:) = outvar(:,115)
      IF(doirrigation)THEN
        wbstores(w_irrcanal,:) = 0.
        DO i = 1,nsub
          DO j = 1,nclass
            wbstores(w_irrcanal,i) = wbstores(w_irrcanal,i) + miscstate%nextirrigation(j,i)*classbasin(i,j)%part
          ENDDO
        ENDDO
        wbstores(w_irrcanal,:) = wbstores(w_irrcanal,:) * basin(:)%area *1.E-3  !m3
      ENDIF  
      IF(naquifers>0)THEN
        CALL calculate_delayed_water(aquiferstate,naquifers,delayedwater)
        wbstores(w_aquifer,1:naquifers) = aquiferstate%water + aquiferstate%nextoutflow + delayedwater
      ENDIF
    ENDIF
     
    !Print out water balance flow 
    IF(conductwb) CALL print_waterbalance_timestep(nsub,naquifers,currentdate)
     
  END SUBROUTINE model

  !Private subroutines, may be moved, for updating   
  !----------------------------------------------------------------
  
  !>Update outflow of subbasin to observed value
  !-------------------------------------------------
  SUBROUTINE apply_quseobs(i,simflow)

    USE MODVAR, ONLY : qobsi,   &
                       missing_value,  &
                       updatestations
     
    !Argument declaration
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    REAL, INTENT(INOUT) :: simflow   !<simulated outflow of subbasin !!(lakeoutflow)
     
    IF(ALLOCATED(updatestations))THEN 
      IF(updatestations(i))THEN
        IF(ALLOCATED(qobsi))THEN
          IF(qobsi(i)/=missing_value)  simflow = qobsi(i)          
        ENDIF
      ENDIF
    ENDIF
     
  END SUBROUTINE apply_quseobs
   
  !>Update outflow of subbasin from observed value with AR method
  !----------------------------------------------------------------
  SUBROUTINE apply_qarupd(i,simflow,corrFlow,arcorr)

    USE MODVAR, ONLY : qobsi,  &
                       missing_value,       &
                       updatestationsqar,   &
                       updatestationsarfact
     
    !Argument declaration
    INTEGER, INTENT(IN)            :: i           !<index of current subbasin
    REAL, INTENT(IN)               :: simflow     !<simulated outflow of subbasin
    REAL, INTENT(INOUT)            :: corrFlow    !<updated outflow (lakeoutflow)
    REAL, INTENT(INOUT)            :: arcorr      !<current AR-error (state-variable)
     
    IF(ALLOCATED(updatestationsqar)) THEN 
      IF(updatestationsqar(i)) THEN
        IF(ALLOCATED(qobsi)) THEN
          IF(qobsi(i)/=missing_value) THEN  !calculates the error
            arcorr =  simflow - qobsi(i)
          ELSE  !no observation, using AR
            arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
            corrFlow = simflow - arcorr
            IF(corrFlow<0.) corrFlow=0.
          ENDIF
        ELSE  !no observation, using AR
          arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
          corrFlow = simflow - arcorr
          IF(corrFlow<0.) corrFlow=0.
        ENDIF
      ENDIF
    ENDIF
  
  END SUBROUTINE apply_qarupd
   
  !>Update outflow of subbasin from observed waterstage value with AR method
  !-------------------------------------------------------------------------
  SUBROUTINE apply_warupd(i,lakeareain,wstold,corrWst,corrFlow,arcorr,lakestate,&
                                        branchflow)

    USE MODVAR, ONLY : xobsi,  &
                       xobsindex, &
                       missing_value,       &
                       wobsvar, &
                       updatestationswar,   &
                       updatestationsarfact
    USE SURFACEWATER_PROCESSES, ONLY : calculate_olake_waterstage, &
                                       calculate_flow_from_lake_waterstage
  
    !Argument declaration
    INTEGER, INTENT(IN)            :: i           !<index of current subbasin
    REAL, INTENT(IN)               :: lakeareain  !<olake area of subbasin (m2)
    REAL, INTENT(IN)               :: wstold      !<simulated lake water end of last time step (mm)
    REAL, INTENT(INOUT)            :: corrWst     !<IN: simulated lake water, OUT: updated lake water (lakewst - for print out only) (mm)
    REAL, INTENT(INOUT)            :: corrFlow    !<updated outflow (lakeoutflow, m3/s)
    REAL, INTENT(INOUT)            :: arcorr      !<current AR-error (state-variable, mm)
    TYPE(lakestatetype),INTENT(IN) :: lakestate   !<Lake state
    REAL,INTENT(INOUT)             :: branchflow  !<outflow from dam branch (used subsequently in calculate_branch_flow)

    !Local variables     
    REAL corrWstold     !updated water stage last time step (mm)
    REAL errorw   !current error in wst in mm
    REAL lakearea !lake area of outlet lake (m2) (whole lake of last lakebasin)
    REAL qoutold  !outflow at old waterstage
    REAL qoutnew  !outflow at new waterstage
    REAL wstobs   !observed lake waterstage wstr from Xobs in w-ref system (m)
    REAL wstm     !lake water stage in local system (m)
    REAL w0ref    !waterstage reference level (m)
    
    !>\b Algorithm \n
    lakearea = lakeareain

    IF(ALLOCATED(updatestationswar)) THEN 
      IF(updatestationswar(i)) THEN
        !>Calculate updated waterstage last timestep
        corrWstold = wstold - arcorr
        !>Calculate error and new arcorr-factor
        IF(ALLOCATED(xobsi))THEN
          IF(xobsindex(wobsvar,i)>0)THEN
            wstobs = xobsi(xobsindex(wobsvar,i))  !Get current wst observation, m in wref-system
          ELSE
            wstobs = missing_value
          ENDIF
          IF(wstobs/=missing_value)THEN
            !Calculate current error of water stage
            CALL calculate_olake_waterstage(i,corrWst,lakeareain,lakearea,wstm,lakestate,w0ref)
            errorw = (wstobs-w0ref-wstm)*1000.         !mm
            IF(ALLOCATED(lakebasinindex))THEN   !For last lake basin, transform change to local deltaw
              IF(lakebasinindex(i)>0)THEN  
                IF(lakebasin(lakebasinindex(i))%last)THEN
                  errorw = errorw * lakearea / lakeareain
                ENDIF
              ENDIF
            ENDIF
            arcorr = - errorw
          ELSE  !no observation, use wAR to update wst and q
            arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
          ENDIF
        ELSE  !no observation, use wAR to update wst and q
          arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
        ENDIF
        !>Apply AR correction to waterstage and discharge
        corrWst = corrWst - arcorr
        CALL calculate_flow_from_lake_waterstage(i,2,lakeareain,corrWstold,qoutold,lakestate,branchflow)
        CALL calculate_flow_from_lake_waterstage(i,2,lakeareain,corrWst,qoutnew,lakestate,branchflow)
        corrFlow = 0.5*(qoutold+qoutnew)
        IF(corrFlow<0.) corrFlow=0.
      ENDIF
    ENDIF
  
  END SUBROUTINE apply_warupd
   
  !>Update outlet lake water stage to observed value
  !-----------------------------------------------------------
  SUBROUTINE apply_wendupd(i,isystem,lakeareain,wst,lakestate,slwimax)

    USE MODVAR, ONLY : xobsi,   &
                       xobsindex,  &
                       missing_value,   &
                       wendupdstations, &
                       wobsvar
    USE SURFACEWATER_PROCESSES, ONLY : calculate_olake_waterstage
     
    !Argument declaration
    INTEGER, INTENT(IN)            :: i            !<index of subbasin
    INTEGER, INTENT(IN)            :: isystem      !<local or outlet lake
    REAL, INTENT(IN)               :: lakeareain   !<lake area of subbasin (m2)
    REAL, INTENT(INOUT)            :: wst          !<water in lake (=lwi+slwi) (mm) to be written output
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate !<Lake state
    REAL,OPTIONAL, INTENT(IN)      :: slwimax      !<target value for slowlake (mm)
     
    !Local variables
    REAL fraction
    REAL wstobs   !wstr for Xobs in w-ref system
    REAL wstm     !lake water stage in w-ref system
    REAL w0ref    !waterstage reference level (m)
    REAL deltaw   !change in calculated water stage due to updating (mm)
    REAL lakearea !lake area of outlet lake (m2)
     
    IF(wendupdstations(i)) THEN
      wstobs = xobsi(xobsindex(wobsvar,i))   !m
      IF(wstobs/=missing_value)THEN
        !Calculate change of water stage and apply to output variable
        CALL calculate_olake_waterstage(i,wst,lakeareain,lakearea,wstm,lakestate,w0ref)
        deltaw = (wstobs-w0ref-wstm)*1000.         !mm
        IF(ALLOCATED(lakebasinindex))THEN   !For last lake basin, transform change to local deltaw
          IF(lakebasinindex(i)>0)THEN  
            IF(lakebasin(lakebasinindex(i))%last)THEN
              deltaw = deltaw * lakearea / lakeareain
            ENDIF
          ENDIF
        ENDIF
        wst = wst + deltaw      !wst = wstobs (mm)
        !Apply updating on state variables for lake water
        IF(PRESENT(slwimax))THEN
          fraction = lakestate%water(isystem,i)/(lakestate%water(isystem,i)+lakestate%slowwater(isystem,i))
          lakestate%slowwater(isystem,i)=(1.-fraction)*wst
          lakestate%water(isystem,i)=fraction*wst
          IF(lakestate%slowwater(isystem,i)>slwimax)THEN
            lakestate%water(isystem,i) = lakestate%water(isystem,i) + (lakestate%slowwater(isystem,i)-slwimax)
            lakestate%slowwater(isystem,i)=slwimax
          ENDIF
        ELSEIF(ALLOCATED(lakestate%slowwater))THEN
          lakestate%slowwater(isystem,i) = wst
        ELSE
          lakestate%water(isystem,i) = wst
        ENDIF
      ENDIF
    ENDIF
     
  END SUBROUTINE apply_wendupd
   
  !>Update concentration out of subbasin to fraction of modelled value
  !----------------------------------------------------------------------
  SUBROUTINE apply_nutrientcorr(correction,conc1,conc2)

    !Argument declaration
    REAL, INTENT(IN)    :: correction   !<update correction value
    REAL, INTENT(INOUT) :: conc1        !<simulated concentration (SP,IN) of outflow of subbasin (clakeoutflow)
    REAL, INTENT(INOUT) :: conc2        !<simulated concentration (PP,ON) of outflow of subbasin (clakeoutflow)
     
    REAL  factor    !current correction factor
     
    IF(correction/=0.)THEN
      factor = 1. + correction
      conc1  = conc1*factor
      conc2  = conc2*factor
    ENDIF
     
  END SUBROUTINE apply_nutrientcorr

  !>Reads files with model specific input
  !>For HYPE it is files with different observations (xom0..xom9,xos0..xos9)
  !----------------------------------------------------------------------
  SUBROUTINE load_modeldefined_input(dir,nsmax,ns,indexarray,bdate,edate,status)

    USE HYPE_INDATA, ONLY : load_xoms_files,  &
                            load_data_for_regression_parameter_estimate
    USE LibDate
    
    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(IN)  :: nsmax         !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: ns            !<Number of subbasins, submodel
    INTEGER, INTENT(IN) :: indexarray(ns) !<index for basemodel
    TYPE(DateType), INTENT(IN)  :: bdate  !<Begin simulation date
    TYPE(DateType), INTENT(IN)  :: edate  !<End simulation date
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    status = 0
    IF(conductxoms) CALL load_xoms_files(dir,ns,bdate,edate,status)
    IF(status/=0) RETURN
    IF(conductregest) CALL load_data_for_regression_parameter_estimate(dir,nsmax,ns,indexarray,status)
    IF(status/=0) RETURN
    
  END SUBROUTINE load_modeldefined_input

  !>Reads files with model specific input
  !>For HYPE it is files with different observations (xom0..xom9,xos0..xos9)
  !----------------------------------------------------------------------
  SUBROUTINE reload_modeldefined_observations(dir,status)

    USE HYPE_INDATA, ONLY : close_xoms_files, &
                            reload_xoms_files
    
    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    status = 0
    IF(conductxoms) CALL close_xoms_files()
    IF(conductxoms) CALL reload_xoms_files(dir,status)
    IF(status/=0) RETURN
    
  END SUBROUTINE reload_modeldefined_observations

  !>Opens files for printout
  !--------------------------------------------------------------------
  SUBROUTINE open_modeldefined_outputfiles(dir,n,na,iens,runens)

    USE HYPE_WATERBALANCE, ONLY : prepare_waterbalance_files
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    
    IF(conductwb) CALL prepare_waterbalance_files(dir,n,na,basin%subid)  !These cannot be printed for several emsembles.

  END SUBROUTINE open_modeldefined_outputfiles

  !>Close files for printout
  !--------------------------------------------------------------------
  SUBROUTINE close_modeldefined_outputfiles(n,na,iens)

    USE HYPE_WATERBALANCE, ONLY : close_waterbalance_files
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    
    IF(conductwb) CALL close_waterbalance_files(na)

  END SUBROUTINE close_modeldefined_outputfiles
    
END MODULE

