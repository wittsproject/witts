! This module is for the Coriolis effect.
!
  MODULE coriolis

  USE parameters
  USE field_shared

  CONTAINS

!==========================================================!
!              CALCULATE CORIOLIS FORCES                   !
!==========================================================!
    SUBROUTINE CORIO()
    IMPLICIT NONE
    REAL(KIND=DP):: FRE_COR

    FRE_COR=2.0*OMEGA_EARTH*SIN(LATITUDE)

    IF(IGEOB.EQ.1)THEN    ! ASSUME GEOSTROPHIC BALANCE
      FX=FX+(V-VG)*FRE_COR
      FY=FY-(U-UG)*FRE_COR
    ELSE
      FX=FX+V*FRE_COR
      FY=FY-U*FRE_COR
    END IF

    END SUBROUTINE
  END MODULE

