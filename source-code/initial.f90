! This module is for initialization and restart.
!
  MODULE initial

  USE mpi
  USE parameters
  USE field_shared
  USE boundary
  USE tools
  USE grid_mesh

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!                         MAIN INITIALIZATION                             !
!=========================================================================!
    SUBROUTINE INITIALIZE()

    IMPLICIT NONE

!---GET DOMAIN DECOMPOSITION INFORMATION
    CALL GETMYID()
!---GRID AND MESH
    CALL MESH()

    FX=0.0
    FY=0.0
    FZ=0.0
    MU=MU0
    RHO=RHO0

    IF(IRESTART.EQ.0)THEN         ! START A NEW RUN
      CALL GENINI_BL_S()
    ELSE IF(IRESTART.EQ.1)THEN
      CALL RESTART()
    ELSE IF(IRESTART.EQ.2)THEN    ! START BY MAPING PREVIOUS RESULT
      CALL MAP()
    END IF

    CALL BOUNDARY_VEL(NX,NY,NZ)

    END SUBROUTINE
!=========================================================================!
!           SUBROUTINE OF GENERATE INITIAL FIELD OF BOUNDARY LAYER        !
!=========================================================================!
!   FOR THE STABLE BOUNDARY LAYER
    SUBROUTINE GENINI_BL_S()

    IMPLICIT NONE

    INCLUDE "mpif.h"
    INTEGER :: I,J,K,IPERTURB
    REAL(KIND=DP):: EPSFAC,ZP,NPERIOD

    EPSFAC=0.1
    NPERIOD=4.0
    IPERTURB=2

    TIME=0.0

    U=UG0
    V=VG0
    W=WG0
    PHI=1.0
    PD=0.0

    ZP=0.3
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          IF(IPERTURB.EQ.1)THEN
            IF(ICOLL.EQ.1)THEN
              U(I,J,K)=U(I,J,K)+EPSFAC*LX/2*COS(NPERIOD*2.0*PI*ZI(K+MYIDZ*NZ)/LZ)*  &
                                            COS(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*  &
                                            SIN(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)
              V(I,J,K)=V(I,J,K)+EPSFAC*LY/2*COS(NPERIOD*2.0*PI*ZI(K+MYIDZ*NZ)/LZ)*  &
                                            SIN(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*  &
                                            COS(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)
              W(I,J,K)=W(I,J,K)+EPSFAC*LZ  *SIN(NPERIOD*2.0*PI*ZI(K+MYIDZ*NZ)/LZ)*  &
                                            SIN(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*  &
                                            SIN(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)
            ELSE          
              U(I,J,K)=U(I,J,K)+EPSFAC*LX/2*COS(NPERIOD*2.0*PI*ZI(K+MYIDZ*NZ)/LZ)*  &
                                            COS(NPERIOD*2.0*PI* X(I+MYIDX*NX)/LX)*  &
                                            SIN(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)
              V(I,J,K)=V(I,J,K)+EPSFAC*LY/2*COS(NPERIOD*2.0*PI*ZI(K+MYIDZ*NZ)/LZ)*  &
                                            SIN(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*  &
                                            COS(NPERIOD*2.0*PI* Y(J+MYIDY*NY)/LY)
              W(I,J,K)=W(I,J,K)+EPSFAC*LZ  *SIN(NPERIOD*2.0*PI* Z(K+MYIDZ*NZ)/LZ)*  &
                                            SIN(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*  &
                                            SIN(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)
            END IF
          ELSE IF(IPERTURB.EQ.2)THEN
            U(I,J,K)=U(I,J,K)+EPSFAC*EXP(0.5)*COS(NPERIOD*2.0*PI*YI(J+MYIDY*NY)/LY)*   &
                              ZI(K+MYIDZ*NZ)/(ZP*LZ)*EXP(-0.5*(ZI(K+MYIDZ*NZ)/(ZP*LZ))**2)

            V(I,J,K)=V(I,J,K)+EPSFAC*EXP(0.5)*SIN(NPERIOD*2.0*PI*XI(I+MYIDX*NX)/LX)*   &
                              ZI(K+MYIDZ*NZ)/(ZP*LZ)*EXP(-0.5*(ZI(K+MYIDZ*NZ)/(ZP*LZ))**2)
            W(I,J,K)=0.0
          END IF
        END DO
      END DO
    END DO
    CALL BOUNDARY_VEL(NX,NY,NZ)

!   FOR TEMPERATURE
    DO K=1,NZ
      DO J=1,NY
        DO I=1,NX
          IF(ZI(K+MYIDZ*NZ).LE.100.0)THEN
            TE(I,J,K)=265.0
          ELSE
            TE(I,J,K)=265.0+(ZI(K+MYIDZ*NZ)-100.0)*0.01
          END IF
        END DO
      END DO
    END DO

    CALL GET_BC(NX,NY,NZ,TE,NBX,NBY,NBZ,0,BC(1,4),BV(1,4))

  END SUBROUTINE GENINI_BL_S
!=========================================================================!
!                     MAP INITIAL FIELD FROM PREVIOUS DATA                !
!=========================================================================! 
      SUBROUTINE RESTART()
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER:: I,J,K,DUMI,M
      REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: TRAN
      REAL(KIND=DP):: DUM
      LOGICAL :: FILE_EXIST

      ALLOCATE(TRAN(NXT,NYT,NZT,5))    
      OPEN(UNIT=1,FILE="RESTART_PRIM")
      READ(1,*) NSTART,TIME
      DO K=1,NZT
        DO J=1,NYT
          DO I=1,NXT
            READ(1,*)DUM,DUM,DUM,(TRAN(I,J,K,M),M=1,5)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(1)
   
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            U(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,1)
            V(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,2)
            W(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,3)
            TE(I,J,K)=TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,4)
            RHO(I,J,K)=TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,5)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(TRAN)

      IF(ITYPE.EQ.1.AND.ISGS.EQ.2.OR.ISGS.EQ.3)THEN   ! RESTART FILE IS NEEDED FOR LAGRANGIAN-TYPE SGS MODEL
        INQUIRE(FILE="RESTART_SGS", EXIST=FILE_EXIST)
 
        IF(FILE_EXIST)THEN
          ALLOCATE(TRAN(NXT,NYT,NZT,4))
          OPEN(UNIT=1,FILE="RESTART_SGS")
          READ(1,*)
          DO K=1,NZT
            DO J=1,NYT
              DO I=1,NXT
                READ(1,*)DUM,DUM,DUM,(TRAN(I,J,K,M),M=1,4)
              ENDDO
            ENDDO
          ENDDO
          CLOSE(1)
   
          DO K=1,NZ
            DO J=1,NY
              DO I=1,NX
                PLM(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,1)
                PMM(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,2)
                PQN(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,3)
                PNN(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,4)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(TRAN)

          LAG_START=0
        ELSE           ! If RESTART_SGS DOESN'T EXIST, THEN START THE LAGRANGIAN PROCESS FROM BEGINNING
          LAG_START=1
        END IF
      END IF
      END SUBROUTINE
!=========================================================================!
!             MAP INITIAL FIELD BY EXTENDING RESTART DATA                 !
!=========================================================================! 
      SUBROUTINE MAP_EXTEND()
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER:: I,J,K,DUMI,M,N
      REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: TRAN
      REAL(KIND=DP):: DUM

      ALLOCATE(TRAN(NXT,NYT,NZT,5))
      OPEN(UNIT=1,FILE="RESTART_PRIM")
      READ(1,*) N,TIME
      DO K=1,NZT/IPROLZ
        DO J=1,NYT/IPROLY
          DO I=1,NXT/IPROLX
            READ(1,*)DUM,DUM,DUM,(TRAN(I,J,K,M),M=1,5)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(1)
!-----PROLONGATE THE DOMAIN IN THE X DIRECTION
      DO K=1,NZT/IPROLZ
        DO J=1,NYT/IPROLY
          DO M=1,IPROLX-1
            DO I=NXT/IPROLX*M+1,NXT/IPROLX*(M+1)
              DO N=1,5 
                TRAN(I,J,K,N)=TRAN(I-NXT/IPROLX*M,J,K,N)
              END DO
            ENDDO
          END DO
        ENDDO
      ENDDO     
!-----PROLONGATE THE DOMAIN IN THE Y DIRECTION
      DO I=1,NXT/IPROLX
        DO K=1,NZT/IPROLZ
          DO M=1,IPROLY-1
            DO J=NYT/IPROLY*M+1,NYT/IPROLY*(M+1)
              DO N=1,5
                TRAN(I,J,K,N)=TRAN(I,J-NYT/IPROLY*M,K,N)
              END DO
            ENDDO
          END DO
        ENDDO
      ENDDO
!-----PROLONGATE THE DOMAIN IN THE Z DIRECTION
      DO J=1,NYT
        DO I=1,NXT
          DO M=1,IPROLZ-1
            DO K=NZT/IPROLZ*M+1,NZT/IPROLZ*(M+1)
              DO N=1,5
               TRAN(I,J,K,N)=TRAN(I,J,K-NZT/IPROLZ*M,N)
              END DO
            ENDDO
          END DO
        ENDDO
      ENDDO
   
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            U(I,J,K)  =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,1)
            V(I,J,K)  =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,2)
            W(I,J,K)  =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,3)
            TE(I,J,K) =TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,4)
            RHO(I,J,K)=TRAN(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,5)
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE MAP_EXTEND
!=========================================================================!
!           MAP INITIAL FIELD BY EXTRACTING FROM RESTART DATA             !
!=========================================================================! 
      SUBROUTINE MAP_EXTRACT()
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER:: I,J,K,DUMI,M
      INTEGER::NXT0,NYT0,NZT0
      REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: TRAN
      REAL(KIND=DP):: DUM,X0(NXT0),Y0(YT0),Z(NZT0),XINT,YINT,ZINT
      LOGICAL :: FILE_EXIST

      ALLOCATE(TRAN(NX0,NY0,NZ0,5))    
      OPEN(UNIT=1,FILE="MAP.DAT")
      READ(1,*) NSTART,TIME
      DO K=1,NZT0
        DO J=1,NYT0
          DO I=1,NXT0
            READ(1,*)X0(I),Y0(J),Z0(K),(TRAN(I,J,K,M),M=1,5)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(1)
   
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,1),1,1,1,U(I,J,K))
            CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,2),1,1,1,V(I,J,K))
            CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,3),1,1,1,W(I,J,K))            
            CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,4),1,1,1,TE(I,J,K))
            CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,5),1,1,1,RHO(I,J,K))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(TRAN)

      IF(ITYPE.EQ.1.AND.ISGS.EQ.2.OR.ISGS.EQ.3)THEN   ! RESTART FILE IS NEEDED FOR LAGRANGIAN-TYPE SGS MODEL
        INQUIRE(FILE="RESTART_SGS", EXIST=FILE_EXIST)
 
        IF(FILE_EXIST)THEN
          ALLOCATE(TRAN(NXT0,NYT0,NZT0,4))
          OPEN(UNIT=1,FILE="RESTART_SGS")
          READ(1,*)
          DO K=1,NZT0
            DO J=1,NYT0
              DO I=1,NXT0
                READ(1,*)X0(I),Y0(J),Z0(K),(TRAN(I,J,K,M),M=1,4)
              ENDDO
            ENDDO
          ENDDO
          CLOSE(1)
  
          DO K=1,NZ
            DO J=1,NY
              DO I=1,NX
                CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                                  NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,1),1,1,1,PLM(I,J,K))
                CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                                  NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,2),1,1,1,PMM(I,J,K))
                CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                                  NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,3),1,1,1,PQN(I,J,K))            
                CALL INTER_GLOBAL(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                                  NXT0,NYT0,NZT0,X0,Y0,Z0,1,TRAN(1,1,1,4),1,1,1,PNN(I,J,K))
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(TRAN)

          LAG_START=0
        ELSE           ! If RESTART_SGS DOESN'T EXIST, THEN START THE LAGRANGIAN PROCESS FROM BEGINNING
          LAG_START=1
        END IF
      END IF
      END SUBROUTINE 
  END MODULE
   
