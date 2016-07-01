! This module is used for obtaining time step dynamically
!
  MODULE dynamic_dt
  
  USE mpi
  USE parameters, ONLY: DP,CFL,DT_MAX,NX,NY,NZ,DT,MYID,IERR,SCREEN_LEVEL, &
                        ITURBINE,RADIUS_TURBINE,ROT_SPEED_TURBINE
  USE field_shared, ONLY: U,V,W

  CONTAINS
!****************************************************!
!     GENERATE DT DYNAMICALLY                        !
!****************************************************!
    SUBROUTINE GENDT(DX,DY,DZ)

    IMPLICIT NONE

    REAL(KIND=DP) :: DX,DY,DZ,MIND,MINDR,MAXU,MAXUR
    REAL(KIND=DP) :: DT_WT

    INTEGER :: I,J,K

    MIND=1.0E10
    MIND=MIN(DX,MIND)
    MIND=MIN(DY,MIND)
    MIND=MIN(DZ,MIND)

    CALL MPI_ALLREDUCE(MIND,MINDR,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERR)
    MAXU=0.0
    DO K=1,NZ
      DO J=1,NY
        DO I=1,NX
          MAXU=DMAX1(ABS(U(I,J,K)),MAXU)
          MAXU=DMAX1(ABS(V(I,J,K)),MAXU)
          MAXU=DMAX1(ABS(W(I,J,K)),MAXU)
        ENDDO
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(MAXU,MAXUR,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERR)

!---LIMITATION BY CFL NUMBER
    DT=MINDR/MAXUR*CFL
!---LIMITATION BY PRESCRIBED MAXIMUM VALUE
    DT=DMIN1(DT,DT_MAX)
!---LIMITATION BY THE MOTION OF WIND TURBINE 
    IF(ITURBINE.EQ.1)THEN
      DT_WT=MINDR/(ROT_SPEED_TURBINE*RADIUS_TURBINE)
      DT=DMIN1(DT,DT_WT)
    END IF

    END SUBROUTINE

  END MODULE
