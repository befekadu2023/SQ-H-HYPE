!> \file general_func.f90
!> Contains module general_functions.

!>General equations that are used by hydrological models
MODULE GENERAL_FUNCTIONS

!Copyright 2012,2014 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tempfactor, &
            halfsatconcfactor, &
            moisturefactor, &
            exponential_moisturefactor, &
            simple_rating_curve

CONTAINS

  !---------------------------------------------------------------
  ! >\brief Calculates a temperature factor 
  !! Based on: Q10=2,reference rate is at 20 degrees and with
  !! thresholds at 0 and 5 degrees (from COUP-model).
  !!
  !! @param[in] temp Current temperature
  ! --------------------------------------------------------------
  REAL FUNCTION tempfactor(temp)  
      
    REAL, INTENT(IN)  :: temp     !current temperature
    !Local variables
    REAL fcn
    
    fcn = 2**((temp - 20.0) / 10.0)
    IF(temp < 5.0) fcn = fcn * (temp / 5.0)
    IF(temp < 0.)  fcn = 0.
    tempfactor = fcn

  END FUNCTION tempfactor
      
  !----------------------------------------------------------------
  !>\brief Calculates a concentration dependence factor 
  !! Based on: half saturation function
  !!
  !! @param[in] conc Current concentration
  !! @param[in] par  Half saturation concentration
  !----------------------------------------------------------------
  REAL FUNCTION halfsatconcfactor(conc,par)  
      
    REAL, INTENT(IN)  :: conc     !current concentration
    REAL, INTENT(IN)  :: par      !half saturation concentration
    !Local variables
    REAL fcn
    
    fcn = conc / (conc + par)
    halfsatconcfactor = fcn
  
  END FUNCTION halfsatconcfactor
      
  !-----------------------------------------------------------------------------
  !> Calculates a soil moisture dependence factor
  !!
  !! @param[in] sm       Soil moisture (mm)
  !! @param[in] wp       Wilting point pore wolume (mm)
  !! @param[in] pw       Total pore wolume (mm)
  !! @param[in] thickm   Thickness of soil layer (m)
  !! @param[in] satsmf   Saturated moisturefactor (satact)
  !! @param[in] exp      Exponent of moisturefactor (thetapow)
  !! @param[in] thetalow Low(?) moisture coefficient (thetalow)
  !! @param[in] thetaupp High(?) moisture coefficient (thetaupp)
  !-----------------------------------------------------------------------------
  REAL FUNCTION moisturefactor(sm,wp,pw,thickm,satsmf,exp,thetalow,thetaupp)
      
    REAL, INTENT(IN)  :: sm     !soil moisture (mm)
    REAL, INTENT(IN)  :: wp     !wilting point pore wolume (mm)
    REAL, INTENT(IN)  :: pw     !total pore wolume (mm)
    REAL, INTENT(IN)  :: thickm !thickness of soil layer (m)
    REAL, INTENT(IN)  :: satsmf !saturated moisturefactor (satact)
    REAL, INTENT(IN)  :: exp    !exponent of moisturefactor (thetapow)
    REAL, INTENT(IN)  :: thetalow !low(?) moisture coefficient (thetalow)
    REAL, INTENT(IN)  :: thetaupp !high(?) moisture coefficient (thetaupp)
    ! Local variables
    REAL thickness  !thickness of soil layer (mm)
    REAL smfcn      !soil moisture dependence factor
    
    !Initiations    
    smfcn = 0.
    thickness = thickm * 1000.      !mm
    
    !Caclulate soilmoisture function value
    IF(thickness>0)THEN
       IF(sm >= pw) THEN   
          smfcn = satsmf
       ELSEIF(sm <= wp) THEN
          smfcn=0.
       ELSE 
          smfcn = min(1., (1-satsmf)*((pw-sm)/((thetaupp/100.)*thickness))**exp + satsmf, ((sm-wp)/((thetalow/100.)*thickness))**exp)    
       ENDIF
    ENDIF
    
    moisturefactor = smfcn

  END FUNCTION moisturefactor
      
  !-----------------------------------------------------------------------------
  !>Calculates an exponential soil moisture dependence factor
  !! 
  !! @param[in] sm     Soil moisture (mm)
  !! @param[in] pw     Total pore wolume (mm)
  !! @param[in] limpar Limitation parameter of moisturefactor (mm) (smfdenitlim)
  !! @param[in] exp    Exponent of moisturefactor (smfdenitpow)
  !-----------------------------------------------------------------------------
  REAL FUNCTION exponential_moisturefactor(sm,pw,limpar,exp)
      
    REAL, INTENT(IN)  :: sm     !soil moisture (mm)
    REAL, INTENT(IN)  :: pw     !total pore wolume (mm)
    REAL, INTENT(IN)  :: limpar !limitation parameter of moisturefactor (mm) (smfdenitlim)
    REAL, INTENT(IN)  :: exp    !exponent of moisturefactor (smfdenitpow)
    ! Local variables
    REAL smfcn      !soil moisture dependence factor
    
    !Initiations    
    smfcn = 0.
    
    !Calculate soilmoisture function value
    IF((sm/pw) > limpar)THEN
       smfcn = (((sm / pw)-limpar)/(1.0-limpar))**exp
    ENDIF
    
    exponential_moisturefactor = smfcn
    
  END FUNCTION exponential_moisturefactor

  !>Calculates momentanous flow from current water level with simple 
  !>rating curve equation: Q = k*(w-w0)**p
  !-----------------------------------------------------------------------------
  REAL FUNCTION simple_rating_curve(wst,k,p,w0)
      
    !Argument declarations
    REAL, INTENT(IN) :: wst    !<current water level (m)
    REAL, INTENT(IN) :: k      !<rating curve coefficient
    REAL, INTENT(IN) :: p      !<rating curve exponent
    REAL, INTENT(IN) :: w0     !<rating curve threshold (m)

    !Local variables
    REAL flow      !discharge (m3/s)
    
    !Calculate discharge
    IF(wst > w0)THEN
      flow = k * (wst-w0)**p
    ELSE
      flow = 0.
    ENDIF
    
    simple_rating_curve = flow
    
  END FUNCTION simple_rating_curve

END MODULE GENERAL_FUNCTIONS
