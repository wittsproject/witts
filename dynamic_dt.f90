! This module is used for obtaining time step dynamically
!
  MODULE dynamic_dt
  
  USE mpi
  USE parameters, ONLY: DP,CFL,DT_MAX,NX,NY,NZ,DT,MYID,IERR,SCREEN_LEVEL
  USE field_shared, ONLY: U,V,W

  CONTAINS
!****************************************************!
!     GENERATE DT DYNAMICALLY                        !
!****************************************************!
    SUBROUTINE GENDT(DX,DY,DZ)

    IMPLICIT NONE

    REAL(KIND=DP) :: DX,DY,DZ,MIND,MINDR,MAXU,MAXUR
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

    DT=MINDR/MAXUR*CFL
    IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
      IF(DT.LT.DT_MAX)THEN
        PRINT*,'DT IS LIMITED BY CFL'
      ELSE
        PRINT*,'DT IS LIMITED BY DT_MAX'
      END IF
    END IF
    DT=DMIN1(DT,DT_MAX)

    END SUBROUTINE

  END MODULE
