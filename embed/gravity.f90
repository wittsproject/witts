! This module is for the gravity (buoyancy) effect
!
  MODULE gravity

  USE parameters
  USE class_shared

  IMPLICIT NONE

  CONTAINS

!===============================================================!
!      GET GRAVITY FORCE USING BOUSSINESQ APPROXIMATION         !
!===============================================================!
    SUBROUTINE BOUSSINESQ()
    IMPLICIT NONE
    
    INTEGER :: M

    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN 
        CELL_FV(M)%CELL_FZ= &
              CELL_FV(M)%CELL_FZ+G*(CELL_FV(M)%CELL_TE-TR)/TR
      END IF   
    END DO
   
    END SUBROUTINE
  
  END MODULE
