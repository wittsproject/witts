! This module is used for generating computational grid
!
  MODULE grid_mesh

  USE parameters, ONLY: NX,NY,NZ,NPX,NPY,NPZ,NXT,NYT,NZT, &
                        NX1,NY1,NZ1,NBX,NBY,NBZ,          &
                        MYIDX,MYIDY,MYIDZ,MYID,LX,LY,LZ,  &
                        DX0,DY0,DZ0
  USE field_shared, ONLY: X,Y,Z,XI,YI,ZI

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!                      SUBROUTINE OF GENERATING THE MESH                  !
!=========================================================================!
      SUBROUTINE MESH()   

      INCLUDE "mpif.h"
      INTEGER :: I,J,K,IERR
!-----GET BASE GRID SPACING
      DX0=LX/NXT
      DY0=LY/NYT
      DZ0=LZ/NZT
!-----FOR X
      X(1)=0.0      
      DO I=2,NXT+NBX
        X(I)=X(I-1)+DX0
      ENDDO
      DO I=0,NX1,-1
        X(I)=X(I+1)-DX0
      END DO
!-----FOR Y     
      Y(1)=0.0
      DO J=2,NYT+NBY
        Y(J)=Y(J-1)+DY0
      ENDDO
      DO J=0,NY1,-1
        Y(J)=Y(J+1)-DY0
      END DO
!-----FOR Z
      Z(1)=0.0      
      DO K=2,NZT+NBZ
        Z(K)=Z(K-1)+DZ0
      ENDDO
      DO K=0,NZ1,-1
        Z(K)=Z(K+1)-DZ0
      END DO 
!-----STAGGERED MESH
      DO I=NX1,NXT+NBX-1
        XI(I)=(X(I)+X(I+1))/2.0
      END DO
      XI(NXT+NBX)=XI(NXT+NBX-1)+DX0

      DO J=NY1,NYT+NBY-1
        YI(J)=(Y(J)+Y(J+1))/2.0
      END DO
      YI(NYT+NBY)=YI(NYT+NBY-1)+DY0

      DO K=NZ1,NZT+NBZ-1
        ZI(K)=(Z(K)+Z(K+1))/2.0
      END DO
      ZI(NZT+NBZ)=ZI(NZT+NBZ-1)+DZ0
 
      END SUBROUTINE 

  END MODULE
