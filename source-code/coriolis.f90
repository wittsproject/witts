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
    REAL(KIND=DP):: FRE_COR,FCORX,FCORY

    FRE_COR=2.0*OMEGA_EARTH*SIN(LATITUDE*PI/180.0)

    IF(IGEOB.EQ.1)THEN    ! ASSUME GEOSTROPHIC BALANCE
      FX=FX+(V-VG0)*FRE_COR
      FY=FY-(U-UG0)*FRE_COR
    ELSE
      FX=FX+V*FRE_COR
      FY=FY-U*FRE_COR
    END IF

    IF(MYID.EQ.0)THEN
      OPEN(1,FILE='Coriolis.echo')
      WRITE(1,*)'Time =', TIME
      WRITE(1,*)'Latitude =',LATITUDE
      WRITE(1,*)'Frequency =',FRE_COR
      WRITE(1,*)'Rossby number =',1.0/FRE_COR
      WRITE(1,*)'IGEOB =',IGEOB
      IF(IGEOB.EQ.1)THEN
        WRITE(1,*)'UG =',UG0
        WRITE(1,*)'VG =',VG0
      END IF
      CLOSE(1)
    END IF

    END SUBROUTINE
  END MODULE

