! This module is for the gravity (buoyancy) effect
!
  MODULE gravity

  USE parameters
  USE field_shared
  USE tools

  IMPLICIT NONE

  CONTAINS

!===============================================================!
!      GET GRAVITY FORCE USING BOUSSINESQ APPROXIMATION         !
!===============================================================!
    SUBROUTINE BOUSSINESQ()
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: TPA
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: TW
    INTEGER :: I,J,K

    ALLOCATE(TPA(NZ))
   
    IF(ICOLL.EQ.0)THEN   ! AT STAGGERED GRID
      ALLOCATE(TW(NX,NY,NZ))
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            TW(I,J,K)=(TE(I,J,K)+TE(I,J,K-1))/2.0
          END DO
        END DO
      END DO
      CALL AVE_P(TW,1,1,1,TPA,1,3)
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            FZ(I,J,K)=FZ(I,J,K)+G*(TW(I,J,K)-TPA(K))/TR
          END DO
        END DO
      END DO
      DEALLOCATE(TW)
    ELSE                ! AT COLLOCATED GRID
      CALL AVE_P(TE,1,1,1,TPA,1,3)
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            FZ(I,J,K)=FZ(I,J,K)+G*(TE(I,J,K)-TPA(K))/TR
          END DO
        END DO
      END DO
    END IF

    DEALLOCATE(TPA)

    END SUBROUTINE
  
  END MODULE
