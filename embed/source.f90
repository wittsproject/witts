! This module is used for source/sink terms
!
  MODULE source

  USE parameters
  USE class_shared
  USE tools

  CONTAINS

!==========================================================!
!         IMPOSING A CONSTANT PRESSURE GRADIENT            !
!==========================================================! 
    SUBROUTINE SOURCE_MOM()
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE :: UPA,VPA,WPA
    INTEGER :: M,I,J,K
    INTEGER :: ISOURCE,SOURCE_FIX_DIR,SOURCE_SPG_DIR 
    REAL(KIND=DP) :: SOURCE_X,SOURCE_Y,SOURCE_Z,LOC_FIX,LOC_SPG_START,LOC_SPG_END
    REAL(KIND=DP) :: UFIX,VFIX,WFIX,USPG,VSPG,WSPG
    REAL(KIND=DP) :: ALFA,LS,LOC

    OPEN(1,FILE="source.in")
    READ(1,*)
    READ(1,*)ISOURCE
    READ(1,*)
    READ(1,*)SOURCE_X
    READ(1,*)SOURCE_Y
    READ(1,*)SOURCE_Z
    READ(1,*)
    READ(1,*)SOURCE_FIX_DIR
    READ(1,*)LOC_FIX
    READ(1,*)UFIX
    READ(1,*)VFIX
    READ(1,*)WFIX
    READ(1,*)
    READ(1,*)SOURCE_SPG_DIR
    READ(1,*)LOC_SPG_START
    READ(1,*)LOC_SPG_END
    READ(1,*)USPG
    READ(1,*)VSPG
    READ(1,*)WSPG
    CLOSE(1)

    IF(ISOURCE.EQ.1)THEN        ! ISOURCE=1: CONSTANT SOURCE
      DO M=1,TOTAL_CELL
        CELL_FV(M)%CELL_VAR(10)=CELL_FV(M)%CELL_VAR(10)+SOURCE_X
        CELL_FV(M)%CELL_VAR(11)=CELL_FV(M)%CELL_VAR(11)+SOURCE_Y
        CELL_FV(M)%CELL_VAR(12)=CELL_FV(M)%CELL_VAR(12)+SOURCE_Z
      END DO
    ELSE IF(ISOURCE.EQ.3)THEN  ! ISOURCE=3: SET SPONGE LAYER BETWEEN  LOC_SPG_START AND LOC_SPG_END
      LS=LOC_SPG_END-LOC_SPG_START
      DO M=1,TOTAL_CELL
        IF(SOURCE_SPG_DIR.EQ.1)THEN        ! ALONG X DIRECTION
          LOC=CELL_FV(M)%CELL_X
        ELSE IF(SOURCE_SPG_DIR.EQ.2)THEN   ! ALONG Y DIRECTION
          LOC=CELL_FV(M)%CELL_Y
        ELSE                               ! ALONG Z DIRECTION
          LOC=CELL_FV(M)%CELL_Z
        END IF   

        IF(LOC.GE.LOC_SPG_START.AND.LOC.LE.LOC_SPG_END)THEN
          CELL_FV(M)%CELL_VAR(10)=CELL_FV(M)%CELL_VAR(10)-  &
                                  SITA(LOC,LOC_SPG_START,LS,2)*(CELL_FV(M)%CELL_VAR(1)-USPG)
          CELL_FV(M)%CELL_VAR(11)=CELL_FV(M)%CELL_VAR(11)-  &
                                  SITA(LOC,LOC_SPG_START,LS,2)*(CELL_FV(M)%CELL_VAR(2)-VSPG)
          CELL_FV(M)%CELL_VAR(12)=CELL_FV(M)%CELL_VAR(12)-  &
                                  SITA(LOC,LOC_SPG_START,LS,2)*(CELL_FV(M)%CELL_VAR(3)-WSPG)
        END IF
      END DO
    END IF  

    END SUBROUTINE
  
  END MODULE

