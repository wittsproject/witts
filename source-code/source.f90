! This module is used for source/sink terms
!
  MODULE source

  USE parameters
  USE field_shared
  USE tools

  CONTAINS

!==========================================================!
!         IMPOSING A CONSTANT PRESSURE GRADIENT            !
!==========================================================! 
    SUBROUTINE SOURCE_MOM()
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE :: UPA,VPA,WPA
    INTEGER :: I,J,K
    INTEGER :: ISOURCE,SOURCE_FIX_DIR,SOURCE_SPG_DIR 
    REAL(KIND=DP) :: SOURCE_X,SOURCE_Y,SOURCE_Z,LOC_FIX,LOC_SPG_START,LOC_SPG_END
    REAL(KIND=DP) :: UFIX,VFIX,WFIX,USPG,VSPG,WSPG
    REAL(KIND=DP) :: ALFA,LS,XLOC,YLOC,ZLOC

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
      FX=FX+SOURCE_X
      FY=FY+SOURCE_Y
      FZ=FZ+SOURCE_Z
    ELSE IF(ISOURCE.EQ.2)THEN   ! ISOURCE=2: FIX MEAN VELOCITY AT A SECTION (ALONG X,Y,OR Z DIR)
      ALFA=0.5
      IF(SOURCE_FIX_DIR.EQ.1)THEN              ! ALONG X DIRECTION
        ALLOCATE(UPA(NXT),VPA(NXT),WPA(NXT))
        CALL AVE_P_GLOBAL(U,NX1,NY1,NZ1,UPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(V,NX1,NY1,NZ1,VPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(W,NX1,NY1,NZ1,WPA,1,SOURCE_FIX_DIR)
        DO I=1,NXT-1
          IF(XI(I).LE.LOC_FIX.AND.XI(I+1).GT.LOC_FIX)THEN
            SOURCE_X=-ALFA*(UPA(I)-UFIX)*DT
            SOURCE_Y=-ALFA*(VPA(I)-VFIX)*DT
            SOURCE_Z=-ALFA*(WPA(I)-WFIX)*DT
          END IF
        END DO
        DEALLOCATE(UPA,VPA,WPA)
      ELSE IF(SOURCE_FIX_DIR.EQ.2)THEN        ! ALONG Y DIRECTION
        ALLOCATE(UPA(NYT),VPA(NYT),WPA(NYT))
        CALL AVE_P_GLOBAL(U,NX1,NY1,NZ1,UPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(V,NX1,NY1,NZ1,VPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(W,NX1,NY1,NZ1,WPA,1,SOURCE_FIX_DIR)
        DO I=1,NYT-1
          IF(YI(I).LE.LOC_FIX.AND.YI(I+1).GT.LOC_FIX)THEN
            SOURCE_X=-ALFA*(UPA(I)-UFIX)*DT
            SOURCE_Y=-ALFA*(VPA(I)-VFIX)*DT
            SOURCE_Z=-ALFA*(WPA(I)-WFIX)*DT
          END IF
        END DO
        DEALLOCATE(UPA,VPA,WPA)
      ELSE IF(SOURCE_FIX_DIR.EQ.3)THEN       ! ALONG Z DIRECTION
        ALLOCATE(UPA(NZT),VPA(NZT),WPA(NZT))
        CALL AVE_P_GLOBAL(U,NX1,NY1,NZ1,UPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(V,NX1,NY1,NZ1,VPA,1,SOURCE_FIX_DIR)
        CALL AVE_P_GLOBAL(W,NX1,NY1,NZ1,WPA,1,SOURCE_FIX_DIR)
        DO I=1,NZT-1
          IF(ZI(I).LE.LOC_FIX.AND.ZI(I+1).GT.LOC_FIX)THEN
            SOURCE_X=-ALFA*(UPA(I)-UFIX)*DT
            SOURCE_Y=-ALFA*(VPA(I)-VFIX)*DT
            SOURCE_Z=-ALFA*(WPA(I)-WFIX)*DT
          END IF
        END DO
        DEALLOCATE(UPA,VPA,WPA)
      END IF
      FX=FX+SOURCE_X
      FY=FY+SOURCE_Y
      FZ=FZ+SOURCE_Z  
    ELSE IF(ISOURCE.EQ.3)THEN  ! ISOURCE=3: SET SPONGE LAYER BETWEEN  LOC_SPG_START AND LOC_SPG_END
      LS=LOC_SPG_END-LOC_SPG_START
      IF(SOURCE_SPG_DIR.EQ.1)THEN        ! ALONG X DIRECTION
        DO I=1,NX
          IF(ICOLL.EQ.1)THEN
            XLOC=XI(I+MYIDX*NX)
          ELSE
            XLOC=X(I+MYIDX*NX)
          END IF
          IF(XLOC.GE.LOC_SPG_START.AND.XLOC.LE.LOC_SPG_END)THEN
            DO K=1,NZ
              DO J=1,NY
                FX(I,J,K)=FX(I,J,K)-SITA(XLOC,LOC_SPG_START,LS,2)*(U(I,J,K)-USPG)
                FY(I,J,K)=FY(I,J,K)-SITA(XLOC,LOC_SPG_START,LS,2)*(V(I,J,K)-VSPG)
                FZ(I,J,K)=FZ(I,J,K)-SITA(XLOC,LOC_SPG_START,LS,2)*(W(I,J,K)-WSPG)
              END DO
            END DO
          END IF
        END DO
      ELSE IF (SOURCE_SPG_DIR.EQ.2)THEN        ! ALONG Y DIRECTION
        DO J=1,NY
          IF(ICOLL.EQ.1)THEN
            YLOC=YI(J+MYIDY*NY)
          ELSE
            YLOC=Y(J+MYIDY*NY)
          END IF
          IF(YLOC.GE.LOC_SPG_START.AND.YLOC.LE.LOC_SPG_END)THEN
            DO I=1,NX
              DO K=1,NZ
                FX(I,J,K)=FX(I,J,K)-SITA(YLOC,LOC_SPG_START,LS,2)*(U(I,J,K)-USPG)
                FY(I,J,K)=FY(I,J,K)-SITA(YLOC,LOC_SPG_START,LS,2)*(V(I,J,K)-VSPG)
                FZ(I,J,K)=FZ(I,J,K)-SITA(YLOC,LOC_SPG_START,LS,2)*(W(I,J,K)-WSPG)
              END DO
            END DO
          END IF
        END DO  
     ELSE IF (SOURCE_SPG_DIR.EQ.3)THEN        ! ALONG Z DIRECTION
        DO K=1,NZ
          IF(ICOLL.EQ.1)THEN
            ZLOC=ZI(K+MYIDZ*NZ)
          ELSE
            ZLOC=Z(K+MYIDZ*NZ)
          END IF
          IF(ZLOC.GE.LOC_SPG_START.AND.ZLOC.LE.LOC_SPG_END)THEN
            DO I=1,NX
              DO J=1,NY
                FX(I,J,K)=FX(I,J,K)-SITA(ZLOC,LOC_SPG_START,LS,2)*(U(I,J,K)-USPG)                
                FY(I,J,K)=FY(I,J,K)-SITA(ZLOC,LOC_SPG_START,LS,2)*(V(I,J,K)-VSPG)
                FZ(I,J,K)=FZ(I,J,K)-SITA(ZLOC,LOC_SPG_START,LS,2)*(W(I,J,K)-WSPG)
              END DO
            END DO
          END IF
        END DO
      END IF   
    END IF

    END SUBROUTINE
  
  END MODULE

