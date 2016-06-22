! This module is used for generating computational grid
!
  MODULE grid_mesh

  USE parameters, ONLY: NX,NY,NZ,NPX,NPY,NPZ,NXT,NYT,NZT, &
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
      DO I=2,NXT+1
        X(I)=X(I-1)+DX0
      ENDDO
!-----FOR Y     
      DO J=2,NYT+1
        Y(J)=Y(J-1)+DY0
      ENDDO
!-----FOR Z
      Z(1)=0.0      
      DO K=2,NZT+1
        Z(K)=Z(K-1)+DZ0
      ENDDO 
!-----STAGGERED MESH
      DO I=1,NXT
        XI(I)=(X(I)+X(I+1))/2.0
      END DO
      XI(NXT+1)=2.0*XI(NXT)-XI(NXT-1)
      DO J=1,NYT
        YI(J)=(Y(J)+Y(J+1))/2.0
      END DO
      YI(NYT+1)=2.0*YI(NYT)-YI(NYT-1)
      DO K=1,NZT
        ZI(K)=(Z(K)+Z(K+1))/2.0
      END DO
      ZI(NZT+1)=2.0*ZI(NZT)-ZI(NZT-1) 
      
      END SUBROUTINE 

  END MODULE
