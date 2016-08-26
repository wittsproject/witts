! This module declares the sizes of dynamic global arrays.
! Also, it will deallocate those arrays at the end of simulation.
  MODULE ALLOCATION

  USE parameters
  USE field_shared

  IMPLICIT NONE

  CONTAINS
!-----------------------------------------------------------------!
!                  ALLOCATE GLOBAL ARRAYS                         !
!-----------------------------------------------------------------!
    SUBROUTINE ALLOCATE_GLOBAL()
      
    ALLOCATE(X (NX1:NXT+NBX),Y (NY1:NYT+NBY),Z (NZ1:NZT+NBZ))
    ALLOCATE(XI(NX1:NXT+NBX),YI(NY1:NYT+NBY),ZI(NZ1:NZT+NBZ))
    ALLOCATE(RHO(NX1:NX2,NY1:NY2,NZ1:NZ2),MU(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(U(NX1:NX2,NY1:NY2,NZ1:NZ2),V(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             W(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(PHI(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(TE(NX1:NX2,NY1:NY2,NZ1:NZ2))


    X=0.0
    Y=0.0
    Z=0.0
    XI=0.0
    YI=0.0
    ZI=0.0

    RHO=0.0
    MU=0.0
    U=0.0
    V=0.0
    W=0.0
    PHI=0.0
    TE=0.0

    END SUBROUTINE

!-----------------------------------------------------------------!
!                  DEALLOCATE GLOBAL ARRAYS                       !
!-----------------------------------------------------------------!
    SUBROUTINE DEALLOCATE_GLOBAL()

    DEALLOCATE(X,Y,Z,XI,YI,ZI,RHO,MU,PHI)
    DEALLOCATE(U,V,W,TE)

    END SUBROUTINE

  END MODULE
