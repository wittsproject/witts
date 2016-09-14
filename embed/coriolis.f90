! This module is for the Coriolis effect.
!
  MODULE coriolis

  USE parameters
  USE class_shared

  CONTAINS

!==========================================================!
!              CALCULATE CORIOLIS FORCES                   !
!==========================================================!
    SUBROUTINE CORIO()
    IMPLICIT NONE
    REAL(KIND=DP):: FRE_COR

    FRE_COR=2.0*OMEGA_EARTH*SIN(LATITUDE*PI/180.0)    

    IF(IGEOB.EQ.1)THEN    ! ASSUME GEOSTROPHIC BALANCE
      DO M=1,TOTAL_CELL 
        CELL_FV(M)%CELL_FX=CELL_FV(M)%CELL_FX+(CELL_FV(M)%CELL_VEL(2)-VG0)*FRE_COR
        CELL_FV(M)%CELL_FY=CELL_FV(M)%CELL_FY-(CELL_FV(M)%CELL_VEL(1)-UG0)*FRE_COR
      END DO  
    ELSE
      DO M=1,TOTAL_CELL 
        CELL_FV(M)%CELL_FX=CELL_FV(M)%CELL_FX+CELL_FV(M)%CELL_VEL(2)*FRE_COR
        CELL_FV(M)%CELL_FY=CELL_FV(M)%CELL_FY-CELL_FV(M)%CELL_VEL(1)*FRE_COR
      END DO
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

