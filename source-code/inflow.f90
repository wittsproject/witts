! This module is for initialization and restart.
!
  MODULE inflow

  USE mpi
  USE parameters
  USE field_shared
  USE boundary
  USE tools
  USE grid_mesh

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!            READ CONTROL PARAMETERS FOR READING INFLOW FILES             !
!=========================================================================!
    SUBROUTINE INFLOW_COE()
    IMPLICIT NONE

    OPEN(1,FILE='inflow.in')
    READ(1,*)INF_INTER
    READ(1,*)INF_NUM_VAR
    READ(1,*)INF_NR_INTERVAL
    READ(1,*)INF_NXT,INF_NYT,INF_NZT
    READ(1,*)INF_DX,INF_DY,INF_DZ
    READ(1,*)INF_TIME_START
    CLOSE(1)

    END SUBROUTINE   
!=========================================================================!
!            SUBROUTINE OF READING INFLOW FILES AT EVERY TIME STEP        !
!=========================================================================!
  SUBROUTINE INFLOW_READ()

  IMPLICIT NONE 
  INCLUDE "mpif.h"
  INTEGER :: N0
  REAL(KIND=DP):: INDATA0(INF_NXT,INF_NYT,INF_NZT,INF_NUM_VAR), &
                  INDATA1(INF_NXT,INF_NYT,INF_NZT,INF_NUM_VAR)
  REAL(KIND=DP):: INF_X(INF_NXT),INF_Y(INF_NYT),INF_Z(INF_NZT)
  REAL(KIND=DP):: VAR0(INF_NXT,INF_NY,INF_NZ,INF_NUM_VAR)
  REAL(KIND=DP):: VAR(INF_NXT,NY,NZ,INF_NUM_VAR)  
  REAL(KIND=DP):: TIME0,TIME1
  REAL(KIND=DP):: DUMX,DUMY,DUMZ
  INTEGER :: I,J,K,M,NT,JUDGE,DUMI
  INTEGER :: MYID,IERR,NZ1,NZ2,NR,NRT
  INTEGER :: NRI,NRTI,REASON
!--ROTATION
  REAL(KIND=DP):: YHUB,UHUB,WHUB,ALFA0,ALFA

  CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYID,IERR)

  OPEN(1,FILE="read_inflow.echo")
  READ(1,*)NR,NRT
  CLOSE(1)

  TIME0=INF_TIME_START
  TIME1=INF_TIME_START

  INDATA0=0.0
  INDATA1=0.0
!---------------------------------------------------------------------
!     IF TIME > TIME1, THEN REPLACE TIME1 BY TIME0, AND READ ONE MORE INFLOW CONDITION AS TIME1
!     MAKE SURE THAT WE ALWAYS HAVE TIME0 < TIME < TIME.
  DO WHILE(TIME.GE.TIME1)
    READ(NR,*,IOSTAT=REASON)DUMI, TIME1 
    DO K=1,INF_NZT
      DO J=1,INF_NYT
        DO I=1,INF_NXT
           READ(NR,*,IOSTAT=REASON)(INDATA1(I,J,K,M),M=1,INF_NUM_VAR)
        ENDDO
        IF(REASON.GT.0)THEN  ! An error is encountered
          IF(MYID.EQ.0)THEN
            PRINT*,'STOP: encounter errors of reading inflow files'
          END IF
          CALL MPI_FINALIZE(IERR)
          STOP
        ELSE IF(REASON.LT.0.AND.NR.LT.NRT)THEN   ! End of file encountered
          NR=NR+INF_NR_INTERVAL
        END IF
      ENDDO
    ENDDO

    IF(TIME.GE.TIME1)THEN
      TIME0=TIME1
      INDATA0=INDATA1
    END IF

    IF(REASON.LT.0.AND.NR.GT.NRT)THEN
      IF(MYID.EQ.0)THEN
        PRINT*,'STOP: run out of inflow files.'
      END IF
      CALL MPI_FINALIZE(IERR)
      STOP
      RETURN
    END IF 
  END DO
! INTERPOLATE IN TIME
  DO K=1,INF_NZT
    DO J=1,INF_NYT
      DO I=1,INF_NXT
        DO M=1,INF_NUM_VAR
          VAR0(I,J,K,M)=(INDATA0(I,J,K,M)*(TIME1-TIME)+ &
                         INDATA1(I,J,K,M)*(TIME-TIME0))/(TIME1-TIME0)
        END DO
      END DO
    ENDDO
  ENDDO
! INTERPOLATE IN SPACE
  IF(INF_INTER.EQ.1)THEN
    DO I=1,INF_NXT
      INF_X(I)=-INF_DX*(INF_NXT-I)
    END DO
    DO J=1,INF_NYT
      INF_Y(J)=INF_DX/2+INF_DX*(J-1)
    END DO
    DO K=1,INF_NZT
      INF_Z(K)=INF_DZ/2+INF_DZ*(K-1)
    END DO   
   
    END DO
    DO K=1,NX
      DO J=1,NY
        DO I=1,INF_NXT
          DO M=1,INF_NUM_VAR
            CALL INTER_GLOBAL(-DX0*(I-1),YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ), &
                              INF_NXT,INF_NYT,INF_NZT,INF_X,INF_Y,INF_Z, &
                              1,VAR0(1,1,1,M),1,1,1,VAR(I,J,K,M))
          END DO
        END DO
      END DO 
    END DO
  ELSE IF(NYT.NE.INF_NYT.OR.NZT.NE.INF_NZT)THEN
    IF(MYID.EQ.0)THEN
      PRINT*,'ERROR: resolution does not match in inflow'
    END IF
    CALL MPI_FINALIZE(IERR)
    STOP
  ELSE
    VAR=VAR0
  END IF    
  
  IF(MYID.EQ.0)THEN
    OPEN(10,FILE="read_inflow_echo")
    WRITE(10,*)NR,NRT
    CLOSE(10)

    PRINT*,TIME0,TIME,TIME1
  END IF
  END SUBROUTINE 
!=========================================================================!
!            SUBROUTINE OF INSTANTANEOUS FLOW FIELD AT ONE TIME STEP      !
!=========================================================================!
      SUBROUTINE GENINFLOW()

      IMPLICIT NONE
      INCLUDE "mpif.h"
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE:: UA,VA,WA,TA
      INTEGER :: I,J,K,MYID,IERR,NEXP

      ALLOCATE(UA(NXT,NYT,NZT),VA(NXT,NYT,NZT), &
               WA(NXT,NYT,NZT),TA(NXT,NYT,NZT))

      CALL ASSEM_ROOT(U,NX1,NY1,NZ1,UA)   
      CALL ASSEM_ROOT(V,NX1,NY1,NZ1,VA)
      CALL ASSEM_ROOT(W,NX1,NY1,NZ1,WA)
      CALL ASSEM_ROOT(TE,NX1,NY1,NZ1,TA)

!      IF(N.EQ.NSTART)THEN
        NEXP=NSTART
!      END IF
!      IF(MOD(N,1000).EQ.0)THEN
!        NEXP=NSTART+1000*INT(N/1000)
!      END IF
!      IF(MYID.EQ.0)THEN
!        PRINT*,NEXP
!      END IF
!     OUTPUT DOMAIN DATA
      IF(MYID.EQ.0)THEN 
        WRITE(NEXP,*)N,TIME
        DO K=1,NZT
          DO J=1,NYT
            DO I=1,3
              WRITE(NEXP,*)UA(I,J,K),VA(I,J,K),WA(I,J,K),TA(I,J,K)
            END DO
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(UA,VA,WA,TA)

      END SUBROUTINE 
  END MODULE
   
