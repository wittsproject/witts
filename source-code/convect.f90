! This module is used for calculation of the unlinear convective term
!
  MODULE convect

  USE parameters
  USE field_shared, ONLY: U,V,W,UF,VF,WF,FX,FY,FZ,CHECK,CHECK0
  USE tools
  USE flux

  IMPLICIT NONE
  
  CONTAINS
!=========================================================================!
!                    MAIN SUBROUTINE FOR THE CONVECTIVE TERM              !
!=========================================================================!
    SUBROUTINE CONVECT_WRAP(DX,DY,DZ)
    
    IMPLICIT NONE

    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2
 
    REAL(KIND=DP):: ALFA,DX,DY,DZ  
    INTEGER :: I,J,K

    ALLOCATE(FCX1(NX1:NX2,NY1:NY2,NZ1:NZ2),FCY1(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             FCZ1(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             FCX2(NX1:NX2,NY1:NY2,NZ1:NZ2),FCY2(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             FCZ2(NX1:NX2,NY1:NY2,NZ1:NZ2))

    IF(ISCHEME.EQ.1)THEN     ! FINITE DIFFERENCE SCHEME
      IF(ICOLL.EQ.1)THEN     ! COLLOCATED GRID
        CALL GETCON_DIV_COL(DX,DY,DZ,FCX1,FCY1,FCZ1)
        CALL GETCON_ADV_COL(DX,DY,DZ,FCX2,FCY2,FCZ2)
      ELSE                   ! STAGGERED GRID
        IF(ORDER_CON.EQ.4)THEN  ! 4TH ORDER
          CALL GETCON4_DIV_STA(DX,DY,DZ,FCX1,FCY1,FCZ1)
          CALL GETCON4_ADV_STA(DX,DY,DZ,FCX2,FCY2,FCZ2)
        ELSE IF(ORDER_CON.EQ.2)THEN  ! 2ND ORDER
          CALL GETCON2_DIV_STA(DX,DY,DZ,FCX1,FCY1,FCZ1)
          CALL GETCON2_ADV_STA(DX,DY,DZ,FCX2,FCY2,FCZ2)                
        END IF
      END IF
    ELSE  ! PSEUDO-SPECTRAL SCHEME
!      CALL GETCON_DIV_SPE(DX,DY,DZ,FCX1,FCY1,FCZ1)
!      CALL GETCON_ADV_SPE(DX,DY,DZ,FCX2,FCY2,FCZ2)      
    END IF

    FX=FX-(FCX1*BLEND_CON+FCX2*(1.0-BLEND_CON))
    FY=FY-(FCY1*BLEND_CON+FCY2*(1.0-BLEND_CON))    
    FZ=FZ-(FCZ1*BLEND_CON+FCZ2*(1.0-BLEND_CON))

    DEALLOCATE(FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2)
    END SUBROUTINE    
!=========================================================================!
!    GET CONVECTIVE TERM, 4TH ORDER, DIVERGENCE FORM, COLLOCATED GRID     !
!=========================================================================!
      SUBROUTINE GETCON_DIV_COL(DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE

      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER :: I,J,K

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            FCX(I,J,K)=CON_DIV(U,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
            FCY(I,J,K)=CON_DIV(V,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
            FCZ(I,J,K)=CON_DIV(W,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
          END DO
        END DO
      END DO
      END SUBROUTINE 
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, ADVECTIVE FORM, COLLOCATED GRID     !
!=========================================================================!
!     4TH ORDER ACCURACY, ADVECTIVE FORM (Oleg V. Vasilyev, JCP,  2000)     
      SUBROUTINE GETCON_ADV_COL(DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE

      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER I,J,K

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            FCX(I,J,K)=CON_ADV(U,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
            FCY(I,J,K)=CON_ADV(V,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
            FCZ(I,J,K)=CON_ADV(W,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_CON)
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, ADVECTIVE FORM, STAGGERED GRID      !
!=========================================================================!     
      SUBROUTINE GETCON4_ADV_STA(DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ
      INTEGER:: I,J,K

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------          
          DDX1=(STENCIL_DX(U,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I,  J,K,1,4)+  &
                STENCIL_DX(U,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I+1,J,K,1,4))/2.0
          DDX2=(STENCIL_DX(U,NX1,NY1,NZ1,I-1,J,K,3,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I-1,J,K,1,4)+  &
                STENCIL_DX(U,NX1,NY1,NZ1,I+2,J,K,3,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I+2,J,K,1,4))/2.0
          DDX=(DDX1*9.0/8.0-DDX2/8.0)

          DDY1=(STENCIL_DY(U,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J,  K,1,4)+  &
                STENCIL_DY(U,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J+1,K,1,4))/2.0
          DDY2=(STENCIL_DY(U,NX1,NY1,NZ1,I,J-1,K,3,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J-1,K,1,4)+  &
                STENCIL_DY(U,NX1,NY1,NZ1,I,J+2,K,3,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J+2,K,1,4))/2.0
          DDY=(DDY1*9.0/8.0-DDY2/8.0)

          DDZ1=(STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K,1  ,4)+  &
                STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+1,1,4))/2.0
          DDZ2=(STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K-1,3,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K-1,1,4)+  &
                STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K+2,3,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+2,1,4))/2.0
          DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

          FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------  
          DDX1=(STENCIL_DX(V,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I,  J,K,1,4)+  &
                STENCIL_DX(V,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I+1,J,K,1,4))/2.0
          DDX2=(STENCIL_DX(V,NX1,NY1,NZ1,I-1,J,K,3,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I-1,J,K,1,4)+  &
                STENCIL_DX(V,NX1,NY1,NZ1,I+2,J,K,3,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I+2,J,K,1,4))/2.0
          DDX=(DDX1*9.0/8.0-DDX2/8.0)

          DDY1=(STENCIL_DY(V,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J,  K,1,4)+  &
                STENCIL_DY(V,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+1,K,1,4))/2.0
          DDY2=(STENCIL_DY(V,NX1,NY1,NZ1,I,J-1,K,3,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J-1,K,1,4)+  &
                STENCIL_DY(V,NX1,NY1,NZ1,I,J+2,K,3,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+2,K,1,4))/2.0
          DDY=(DDY1*9.0/8.0-DDY2/8.0)

          DDZ1=(STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K,1  ,4)+  &
                STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+1,1,4))/2.0
          DDZ2=(STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K-1,3,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K-1,1,4)+  &
                STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K+2,3,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+2,1,4))/2.0
          DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

          FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------     
          DDX1=(STENCIL_DX(W,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I,  J,K,1,4)+  &
                STENCIL_DX(W,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I+1,J,K,1,4))/2.0
          DDX2=(STENCIL_DX(W,NX1,NY1,NZ1,I-1,J,K,3,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I-1,J,K,1,4)+  &
                STENCIL_DX(W,NX1,NY1,NZ1,I+2,J,K,3,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I+2,J,K,1,4))/2.0
          DDX=(DDX1*9.0/8.0-DDX2/8.0)

          DDY1=(STENCIL_DY(W,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J,  K,1,4)+  &
                STENCIL_DY(W,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+1,K,1,4))/2.0
          DDY2=(STENCIL_DY(W,NX1,NY1,NZ1,I,J-1,K,3,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J-1,K,1,4)+  &
                STENCIL_DY(W,NX1,NY1,NZ1,I,J+2,K,3,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+2,K,1,4))/2.0
          DDY=(DDY1*9.0/8.0-DDY2/8.0)

          DDZ1=(STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K,1  ,4)+  &
                STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+1,1,4))/2.0
          DDZ2=(STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K-1,3,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K-1,1,4)+  &
                STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K+2,3,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+2,1,4))/2.0
          DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

          FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, DIVERGENCE FORM, STAGGERED GRID     !
!=========================================================================!
      SUBROUTINE GETCON4_DIV_STA(DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: U1A,U1B,U3A,U3B,V1A,V1B,V3A,V3B,W1A,W1B,W3A,W3B
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ
      INTEGER:: I,J,K
            
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------
          DDX1=(VAR_INTER_X(U,NX1,NY1,NZ1,I+1,J,K,1,4)*STENCIL_IX(U,NX1,NY1,NZ1,I+1,J,K,1)-  &
                VAR_INTER_X(U,NX1,NY1,NZ1,I,  J,K,1,4)*STENCIL_IX(U,NX1,NY1,NZ1,I,  J,K,1))/DX
          DDX2=(VAR_INTER_X(U,NX1,NY1,NZ1,I+2,J,K,1,4)*STENCIL_IX(U,NX1,NY1,NZ1,I+2,J,K,3)-  &
                VAR_INTER_X(U,NX1,NY1,NZ1,I-1,J,K,1,4)*STENCIL_IX(U,NX1,NY1,NZ1,I-1,J,K,3))/(DX*3.0)
          DDX=DDX1*9.0/8.0-DDX2/8.0

          DDY1=(VAR_INTER_X(V,NX1,NY1,NZ1,I,J+1,K,1,4)*STENCIL_IY(U,NX1,NY1,NZ1,I,J+1,K,1)-  &
                VAR_INTER_X(V,NX1,NY1,NZ1,I,J,  K,1,4)*STENCIL_IY(U,NX1,NY1,NZ1,I,J,  K,1))/DY
          DDY2=(VAR_INTER_X(V,NX1,NY1,NZ1,I,J+2,K,1,4)*STENCIL_IY(U,NX1,NY1,NZ1,I,J+2,K,3)-  &
                VAR_INTER_X(V,NX1,NY1,NZ1,I,J-1,K,1,4)*STENCIL_IY(U,NX1,NY1,NZ1,I,J-1,K,3))/(DY*3.0)
          DDY=DDY1*9.0/8.0-DDY2/8.0

          DDZ1=(VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+1,1,4)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K+1,1)-  &
                VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K,  1,4)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K,  1))/DZ
          DDZ2=(VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+2,1,4)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K+2,3)-  &
                VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K-1,1,4)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K-1,3))/(DZ*3.0)
          DDZ=DDZ1*9.0/8.0-DDZ2/8.0
            
          FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------
          DDX1=(VAR_INTER_Y(U,NX1,NY1,NZ1,I+1,J,K,1,4)*STENCIL_IX(V,NX1,NY1,NZ1,I+1,J,K,1)-  &
                VAR_INTER_Y(U,NX1,NY1,NZ1,I,  J,K,1,4)*STENCIL_IX(V,NX1,NY1,NZ1,I,  J,K,1))/DX
          DDX2=(VAR_INTER_Y(U,NX1,NY1,NZ1,I+2,J,K,1,4)*STENCIL_IX(V,NX1,NY1,NZ1,I+2,J,K,3)-  &
                VAR_INTER_Y(U,NX1,NY1,NZ1,I-1,J,K,1,4)*STENCIL_IX(V,NX1,NY1,NZ1,I-1,J,K,3))/(DX*3.0)
          DDX=DDX1*9.0/8.0-DDX2/8.0

          DDY1=(VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+1,K,1,4)*STENCIL_IY(V,NX1,NY1,NZ1,I,J+1,K,1)-  &
                VAR_INTER_Y(V,NX1,NY1,NZ1,I,J,  K,1,4)*STENCIL_IY(V,NX1,NY1,NZ1,I,J,  K,1))/DY
          DDY2=(VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+2,K,1,4)*STENCIL_IY(V,NX1,NY1,NZ1,I,J+2,K,3)-  &
                VAR_INTER_Y(V,NX1,NY1,NZ1,I,J-1,K,1,4)*STENCIL_IY(V,NX1,NY1,NZ1,I,J-1,K,3))/(DY*3.0)
          DDY=DDY1*9.0/8.0-DDY2/8.0

          DDZ1=(VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+1,1,4)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K+1,1)-  &
                VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K,  1,4)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K,  1))/DZ
          DDZ2=(VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+2,1,4)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K+2,3)-  &
                VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K-1,1,4)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K-1,3))/(DZ*3.0)
          DDZ=DDZ1*9.0/8.0-DDZ2/8.0    
       
          FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------
          DDX1=(VAR_INTER_Z(U,NX1,NY1,NZ1,I+1,J,K,1,4)*STENCIL_IX(W,NX1,NY1,NZ1,I+1,J,K,1)-  &
                VAR_INTER_Z(U,NX1,NY1,NZ1,I,  J,K,1,4)*STENCIL_IX(W,NX1,NY1,NZ1,I,  J,K,1))/DX
          DDX2=(VAR_INTER_Z(U,NX1,NY1,NZ1,I+2,J,K,1,4)*STENCIL_IX(W,NX1,NY1,NZ1,I+2,J,K,3)-  &
                VAR_INTER_Z(U,NX1,NY1,NZ1,I-1,J,K,1,4)*STENCIL_IX(W,NX1,NY1,NZ1,I-1,J,K,3))/(DX*3.0)
          DDX=DDX1*9.0/8.0-DDX2/8.0

          DDY1=(VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+1,K,1,4)*STENCIL_IY(W,NX1,NY1,NZ1,I,J+1,K,1)-  &
                VAR_INTER_Z(V,NX1,NY1,NZ1,I,J,  K,1,4)*STENCIL_IY(W,NX1,NY1,NZ1,I,J,  K,1))/DY
          DDY2=(VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+2,K,1,4)*STENCIL_IY(W,NX1,NY1,NZ1,I,J+2,K,3)-  &
                VAR_INTER_Z(V,NX1,NY1,NZ1,I,J-1,K,1,4)*STENCIL_IY(W,NX1,NY1,NZ1,I,J-1,K,3))/(DY*3.0)
          DDY=DDY1*9.0/8.0-DDY2/8.0

          DDZ1=(VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+1,1,4)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K+1,1)-  &
                VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K,  1,4)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K,  1))/DZ
          DDZ2=(VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+2,1,4)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K+2,3)-  &
                VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K-1,1,4)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K-1,3))/(DZ*3.0)
          DDZ=DDZ1*9.0/8.0-DDZ2/8.0

          FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 2ND ORDER, ADVECTIVE FORM, STAGGERED GRID      !
!=========================================================================!     
      SUBROUTINE GETCON2_ADV_STA(DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
      INCLUDE "mpif.h"
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DDX,DDY,DDZ
      INTEGER:: I,J,K

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------          
          DDX=(STENCIL_DX(U,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I,  J,K,1,2)+  &
               STENCIL_DX(U,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_X(U,NX1,NY1,NZ1,I+1,J,K,1,2))/2.0

          DDY=(STENCIL_DY(U,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J,  K,1,2)+  &
               STENCIL_DY(U,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_X(V,NX1,NY1,NZ1,I,J+1,K,1,2))/2.0

          DDZ=(STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K,1  ,2)+  &
               STENCIL_DZ(U,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+1,1,2))/2.0

          FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------  
          DDX=(STENCIL_DX(V,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I,  J,K,1,2)+  &
               STENCIL_DX(V,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_Y(U,NX1,NY1,NZ1,I+1,J,K,1,2))/2.0

          DDY=(STENCIL_DY(V,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J,  K,1,2)+  &
               STENCIL_DY(V,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+1,K,1,2))/2.0

          DDZ=(STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K,1  ,2)+  &
               STENCIL_DZ(V,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+1,1,2))/2.0

          FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------     
          DDX=(STENCIL_DX(W,NX1,NY1,NZ1,I,  J,K,1,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I,  J,K,1,2)+  &
               STENCIL_DX(W,NX1,NY1,NZ1,I+1,J,K,1,DX)*VAR_INTER_Z(U,NX1,NY1,NZ1,I+1,J,K,1,2))/2.0

          DDY=(STENCIL_DY(W,NX1,NY1,NZ1,I,J,  K,1,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J,  K,1,2)+  &
               STENCIL_DY(W,NX1,NY1,NZ1,I,J+1,K,1,DY)*VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+1,K,1,2))/2.0

          DDZ=(STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K,  1,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K,1  ,2)+  &
               STENCIL_DZ(W,NX1,NY1,NZ1,I,J,K+1,1,DZ)*VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+1,1,2))/2.0

          FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 2ND ORDER, DIVERGENCE FORM, STAGGERED GRID     !
!=========================================================================!
      SUBROUTINE GETCON2_DIV_STA(DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DDX,DDY,DDZ
      INTEGER:: I,J,K
            
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------
          DDX=(VAR_INTER_X(U,NX1,NY1,NZ1,I+1,J,K,1,2)*STENCIL_IX(U,NX1,NY1,NZ1,I+1,J,K,1)-  &
               VAR_INTER_X(U,NX1,NY1,NZ1,I,  J,K,1,2)*STENCIL_IX(U,NX1,NY1,NZ1,I,  J,K,1))/DX

          DDY=(VAR_INTER_X(V,NX1,NY1,NZ1,I,J+1,K,1,2)*STENCIL_IY(U,NX1,NY1,NZ1,I,J+1,K,1)-  &
               VAR_INTER_X(V,NX1,NY1,NZ1,I,J,  K,1,2)*STENCIL_IY(U,NX1,NY1,NZ1,I,J,  K,1))/DY

          DDZ=(VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K+1,1,2)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K+1,1)-  &
               VAR_INTER_X(W,NX1,NY1,NZ1,I,J,K,  1,2)*STENCIL_IZ(U,NX1,NY1,NZ1,I,J,K,  1))/DZ
            
          FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------
          DDX=(VAR_INTER_Y(U,NX1,NY1,NZ1,I+1,J,K,1,2)*STENCIL_IX(V,NX1,NY1,NZ1,I+1,J,K,1)-  &
               VAR_INTER_Y(U,NX1,NY1,NZ1,I,  J,K,1,2)*STENCIL_IX(V,NX1,NY1,NZ1,I,  J,K,1))/DX

          DDY=(VAR_INTER_Y(V,NX1,NY1,NZ1,I,J+1,K,1,2)*STENCIL_IY(V,NX1,NY1,NZ1,I,J+1,K,1)-  &
               VAR_INTER_Y(V,NX1,NY1,NZ1,I,J,  K,1,2)*STENCIL_IY(V,NX1,NY1,NZ1,I,J,  K,1))/DY

          DDZ=(VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K+1,1,2)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K+1,1)-  &
               VAR_INTER_Y(W,NX1,NY1,NZ1,I,J,K,  1,2)*STENCIL_IZ(V,NX1,NY1,NZ1,I,J,K,  1))/DZ

          FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------
          DDX=(VAR_INTER_Z(U,NX1,NY1,NZ1,I+1,J,K,1,2)*STENCIL_IX(W,NX1,NY1,NZ1,I+1,J,K,1)-  &
               VAR_INTER_Z(U,NX1,NY1,NZ1,I,  J,K,1,2)*STENCIL_IX(W,NX1,NY1,NZ1,I,  J,K,1))/DX

          DDY=(VAR_INTER_Z(V,NX1,NY1,NZ1,I,J+1,K,1,2)*STENCIL_IY(W,NX1,NY1,NZ1,I,J+1,K,1)-  &
               VAR_INTER_Z(V,NX1,NY1,NZ1,I,J,  K,1,2)*STENCIL_IY(W,NX1,NY1,NZ1,I,J,  K,1))/DY

          DDZ=(VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K+1,1,2)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K+1,1)-  &
               VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K,  1,2)*STENCIL_IZ(W,NX1,NY1,NZ1,I,J,K,  1))/DZ

          FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO


      END SUBROUTINE
!=========================================================================!
!              GET CONVECTIVE TERM, PSEUDO-SPEC, ADVECTIVE FORM           !
!=========================================================================!     
      SUBROUTINE GETCON_ADV_SPE(DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE
      INTEGER:: I,J,K
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE::TRAN,UDX,UDZ,VDX,VDZ,WDX,WDZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP):: U1A,U1B,V1A,V1B,W1A,W1B
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ, &
             DDY1A,DDY1B
      REAL(KIND=DP):: FX,FY,FZ,LX,LZ
!-----USE SPECTRAL METHOD TO CALCULATE THE HORIZONTAL GRADIENTS
      ALLOCATE(UDX(NX1:NX2,NY1:NY2,NZ1:NZ2),UDZ(NX1:NX2,NY1:NY2,NZ1:NZ2))
      ALLOCATE(VDX(NX1:NX2,NY1:NY2,NZ1:NZ2),VDZ(NX1:NX2,NY1:NY2,NZ1:NZ2))
      ALLOCATE(WDX(NX1:NX2,NY1:NY2,NZ1:NZ2),WDZ(NX1:NX2,NY1:NY2,NZ1:NZ2))

      CALL GRADIENT_SPEC(U,NX1,NY1,NZ1,UDX,UDZ,1)
      CALL GRADIENT_SPEC(V,NX1,NY1,NZ1,VDX,VDZ,1)
      CALL GRADIENT_SPEC(W,NX1,NY1,NZ1,WDX,WDZ,1)
!-----DEALIASING USING 2/3 RULE
      CALL DEALIAS(UDX,NX1,NY1,NZ1)
      CALL DEALIAS(UDZ,NX1,NY1,NZ1)
      CALL DEALIAS(VDX,NX1,NY1,NZ1)
      CALL DEALIAS(VDZ,NX1,NY1,NZ1)
      CALL DEALIAS(WDX,NX1,NY1,NZ1)
      CALL DEALIAS(WDZ,NX1,NY1,NZ1)

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FCX------------------------------------------------
            FX=U(I,J,K)*UDX(I,J,K)

            V1A=V(I,J,K)
            V1B=V(I,J+1,K)
            DDY1A=(U(I,J,  K)-U(I,J-1,K))/DY
            DDY1B=(U(I,J+1,K)-U(I,J,  K))/DY
            DDY1=(DDY1B*V1B+DDY1A*V1A)/2.0
            FY=DDY1

            FZ=W(I,J,K)*UDZ(I,J,K)

            FCX(I,J,K)=FX+FY+FZ
!-----FCY------------------------------------------------
            FX=(U(I,J,K)+U(I,J-1,K))/2.0*VDX(I,J,K)  

            DDY1A=(V(I,J,  K)-V(I,J-1,K))/DY
            DDY1B=(V(I,J+1,K)-V(I,J,  K))/DY
            DDY1=(DDY1B+DDY1A)/2.0
            FY=V(I,J,K)*DDY1

            FZ=(W(I,J,K)+W(I,J-1,K))/2.0*VDZ(I,J,K) 

            FCY(I,J,K)=FX+FY+FZ
!-----FCZ------------------------------------------------
            FX=U(I,J,K)*WDX(I,J,K)

            V1A=V(I,J,  K)
            V1B=V(I,J+1,K)
            DDY1A=(W(I,J,  K)-W(I,J-1,K))/DY
            DDY1B=(W(I,J+1,K)-W(I,J,  K))/DY
            DDY1=(DDY1B*V1B+DDY1A*V1A)/2.0
            FY=DDY1

            FZ=W(I,J,K)*WDZ(I,J,K)      

            FCZ(I,J,K)=FX+FY+FZ
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(UDX,UDZ,VDX,VDZ,WDX,WDZ)

      END SUBROUTINE
!=========================================================================!
!                            DIVERGENCE FORM                              !
!=========================================================================!
      REAL(KIND=DP) FUNCTION CON_DIV(VAR,SI1,SI2,SI3,DX,DY,DZ,I,J,K,ORDER)

      IMPLICIT NONE

      INTEGER :: SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER :: I,J,K,ORDER

      IF(ORDER.EQ.2)THEN
        CON_DIV=  &   
         (UF(I+1,J,K)*VAR_INTER_X(VAR,SI1,SI2,SI3,I+1,J,K,1,2)-     &
          UF(I,  J,K)*VAR_INTER_X(VAR,SI1,SI2,SI3,I,  J,K,1,2))/DX+ &
         (VF(I,J+1,K)*VAR_INTER_Y(VAR,SI1,SI2,SI3,I,J+1,K,1,2)-     &
          VF(I,J,  K)*VAR_INTER_Y(VAR,SI1,SI2,SI3,I,J,  K,1,2))/DY+ &
         (WF(I,J,K+1)*VAR_INTER_Z(VAR,SI1,SI2,SI3,I,J,K+1,1,2)-     &
          WF(I,J,K  )*VAR_INTER_Z(VAR,SI1,SI2,SI3,I,J,K  ,1,2))/DZ 
      ELSE IF(ORDER.EQ.4)THEN
         CON_DIV=  &
 9.0/8.0*((UF(I+1,J,K)*VAR_INTER_X(VAR,SI1,SI2,SI3,I+1,J,K,1,4)-            &
           UF(I,  J,K)*VAR_INTER_X(VAR,SI1,SI2,SI3,I,  J,K,1,4))/DX+        &
          (VF(I,J+1,K)*VAR_INTER_Y(VAR,SI1,SI2,SI3,I,J+1,K,1,4)-            &
           VF(I,J,  K)*VAR_INTER_Y(VAR,SI1,SI2,SI3,I,J,  K,1,4))/DY+        &
          (WF(I,J,K+1)*VAR_INTER_Z(VAR,SI1,SI2,SI3,I,J,K+1,1,4)-            &
           WF(I,J,K  )*VAR_INTER_Z(VAR,SI1,SI2,SI3,I,J,K,  1,4))/DZ)-       &            
 1.0/8.0*((UF(I+2,J,K)*(STENCIL_IX(VAR,SI1,SI2,SI3,I+2,J,K,1)*9.0/8.0-STENCIL_IX(VAR,SI1,SI2,SI3,I+2,J,K,3)/8.0)-            &
           UF(I-1,J,K)*(STENCIL_IX(VAR,SI1,SI2,SI3,I-1,J,K,1)*9.0/8.0-STENCIL_IX(VAR,SI1,SI2,SI3,I-1,J,K,3)/8.0))/(DX*3.0)+  &
          (VF(I,J+2,K)*(STENCIL_IY(VAR,SI1,SI2,SI3,I,J+2,K,1)*9.0/8.0-STENCIL_IY(VAR,SI1,SI2,SI3,I,J+2,K,3)/8.0)-            &
           VF(I,J-1,K)*(STENCIL_IY(VAR,SI1,SI2,SI3,I,J-1,K,1)*9.0/8.0-STENCIL_IY(VAR,SI1,SI2,SI3,I,J-1,K,3)/8.0))/(DY*3.0)+  &
          (WF(I,J,K+2)*(STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+2,1)*9.0/8.0-STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+2,3)/8.0)-            &
           WF(I,J,K-1)*(STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K-1,1)*9.0/8.0-STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K-1,3)/8.0))/(DZ*3.0))
      END IF

      END FUNCTION
!=========================================================================!
!                     DIVERGENCE FORM USING MUSCL SCHEME                  !
!=========================================================================!
!     SECOND ORDER ONLY
!     ORDER: THE NUMERICAL ORDER ONLY FOR THE CENTRAL SCHEME PART
      REAL(KIND=DP) FUNCTION CON_DIV_MUSCL(VAR,SI1,SI2,SI3,DX,DY,DZ,I,J,K,ORDER)

      IMPLICIT NONE

      INTEGER :: SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP):: VARX1,VARX2,VARY1,VARY2,VARZ1,VARZ2
      INTEGER :: I,J,K,ORDER


      VARX1=FLUX_MUSCL(VAR,SI1,SI2,SI3,UF(I,  J,K),I,  J,K,1,ORDER,DX)
      VARX2=FLUX_MUSCL(VAR,SI1,SI2,SI3,UF(I+1,J,K),I+1,J,K,1,ORDER,DX)

      VARY1=FLUX_MUSCL(VAR,SI1,SI2,SI3,VF(I,J,  K),I,J,  K,2,ORDER,DY)
      VARY2=FLUX_MUSCL(VAR,SI1,SI2,SI3,VF(I,J+1,K),I,J+1,K,2,ORDER,DY)
    
      VARZ1=FLUX_MUSCL(VAR,SI1,SI2,SI3,WF(I,J,K  ),I,J,K,  3,ORDER,DZ)
      VARZ2=FLUX_MUSCL(VAR,SI1,SI2,SI3,WF(I,J,K+1),I,J,K+1,3,ORDER,DZ)

      CON_DIV_MUSCL=(VARX2*UF(I+1,J,K)-VARX1*UF(I,J,K))/DX+ &
                    (VARY2*VF(I,J+1,K)-VARY1*VF(I,J,K))/DY+ &
                    (VARZ2*WF(I,J,K+1)-VARZ1*WF(I,J,K))/DZ

      END FUNCTION
!=========================================================================!
!                          DIVERGENCE FORM                                !
!=========================================================================!
      REAL(KIND=DP) FUNCTION CON_ADV(VAR,SI1,SI2,SI3,DX,DY,DZ,I,J,K,ORDER)

      IMPLICIT NONE

      INTEGER :: SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER :: I,J,K,ORDER

      IF(ORDER.EQ.2)THEN
        CON_ADV=  &
        ((UF(I  ,J,K)*DERIV_X(VAR,SI1,SI2,SI3,I,  J,K,1,2,DX)+      &
          UF(I+1,J,K)*DERIV_X(VAR,SI1,SI2,SI3,I+1,J,K,1,2,DX))/2.0+ &
         (VF(I,J,  K)*DERIV_Y(VAR,SI1,SI2,SI3,I,J,  K,1,2,DY)+      &
          VF(I,J+1,K)*DERIV_Y(VAR,SI1,SI2,SI3,I,J+1,K,1,2,DY))/2.0+ &
         (WF(I,J,K  )*DERIV_Z(VAR,SI1,SI2,SI3,I,J,K,  1,2,DZ)+      &
          WF(I,J,K+1)*DERIV_Z(VAR,SI1,SI2,SI3,I,J,K+1,1,2,DZ))/2.0)
      ELSE IF(ORDER.EQ.4)THEN
        CON_ADV=  &
 9.0/8.0*((UF(I,  J,K)*DERIV_X(VAR,SI1,SI2,SI3,I,  J,K,1,4,DX)+       &
           UF(I+1,J,K)*DERIV_X(VAR,SI1,SI2,SI3,I+1,J,K,1,4,DX))/2.0+  &
          (VF(I,J,  K)*DERIV_Y(VAR,SI1,SI2,SI3,I,J,  K,1,4,DY)+       &
           VF(I,J+1,K)*DERIV_Y(VAR,SI1,SI2,SI3,I,J+1,K,1,4,DY))/2.0+  &
          (WF(I,J,K  )*DERIV_Z(VAR,SI1,SI2,SI3,I,J,K,  1,4,DZ)+       &
           WF(I,J,K+1)*DERIV_Z(VAR,SI1,SI2,SI3,I,J,K+1,1,4,DZ))/2.0)- &
 1.0/8.0*(((STENCIL_DX(VAR,SI1,SI2,SI3,I+2,J,K,1,DX)*9.0/8.0-STENCIL_DX(VAR,SI1,SI2,SI3,I+2,J,K,3,DX)/8.0)*UF(I+2,J,K)+ &
           (STENCIL_DX(VAR,SI1,SI2,SI3,I-1,J,K,1,DX)*9.0/8.0-STENCIL_DX(VAR,SI1,SI2,SI3,I-1,J,K,3,DX)/8.0)*UF(I-1,J,K))/2.0+ &
          ((STENCIL_DY(VAR,SI1,SI2,SI3,I,J+2,K,1,DY)*9.0/8.0-STENCIL_DY(VAR,SI1,SI2,SI3,I,J+2,K,3,DY)/8.0)*VF(I,J+2,K)+ &
           (STENCIL_DY(VAR,SI1,SI2,SI3,I,J-1,K,1,DY)*9.0/8.0-STENCIL_DY(VAR,SI1,SI2,SI3,I,J-1,K,3,DY)/8.0)*VF(I,J-1,K))/2.0+ &
          ((STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+2,1,DZ)*9.0/8.0-STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+2,3,DZ)/8.0)*WF(I,J,K+2)+ &
           (STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K-1,1,DZ)*9.0/8.0-STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K-1,3,DZ)/8.0)*WF(I,J,K-1))/2.0)
      END IF

      END FUNCTION

!-------------------------------------------------------------------
      SUBROUTINE GETCONVEC4(DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
 
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: U1A,U1B,U3A,U3B,V1A,V1B,V3A,V3B,W1A,W1B,W3A,W3B
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ
      INTEGER:: I,J,K
            
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------
            U1A=(U(I,J,K)+U(I-1,J,K))/2.0
            U1B=(U(I,J,K)+U(I+1,J,K))/2.0
            U3A=(U(I+1,J,K)+U(I-2,J,K))/2.0
            U3B=(U(I-1,J,K)+U(I+2,J,K))/2.0
            DDX1=((U1B*9.0/8.0-U3B/8.0)*U1B-(U1A*9.0/8.0-U3A/8.0)*U1A)/DX
            U1A=(U(I-1,J,K)+U(I-2,J,K))/2.0
            U1B=(U(I+1,J,K)+U(I+2,J,K))/2.0
            U3A=(U(I,J,K)+U(I-3,J,K))/2.0
            U3B=(U(I,J,K)+U(I+3,J,K))/2.0
            DDX2=((U1B*9.0/8.0-U3B/8.0)*U3B-(U1A*9.0/8.0-U3A/8.0)*U3A)/(DX*3.0)
            DDX=DDX1*9.0/8.0-DDX2/8.0

            V1A=(V(I,J  ,K)+V(I-1,J  ,K))/2.0
            V1B=(V(I,J+1,K)+V(I-1,J+1,K))/2.0
            V3A=(V(I+1,J  ,K)+V(I-2,J  ,K))/2.0
            V3B=(V(I+1,J+1,K)+V(I-2,J+1,K))/2.0
            U1A=(U(I,J,K)+U(I,J-1,K))/2.0
            U1B=(U(I,J,K)+U(I,J+1,K))/2.0            
            DDY1=((V1B*9.0/8.0-V3B/8.0)*U1B-(V1A*9.0/8.0-V3A/8.0)*U1A)/DY
            V1A=(V(I,J-1,K)+V(I-1,J-1,K))/2.0
            V1B=(V(I,J+2,K)+V(I-1,J+2,K))/2.0
            V3A=(V(I+1,J-1,K)+V(I-2,J-1,K))/2.0
            V3B=(V(I+1,J+2,K)+V(I-2,J+2,K))/2.0
            U3A=(U(I,J,K)+U(I,J-3,K))/2.0
            U3B=(U(I,J,K)+U(I,J+3,K))/2.0            
            DDY2=((V1B*9.0/8.0-V3B/8.0)*U3B-(V1A*9.0/8.0-V3A/8.0)*U3A)/(DY*3.0)
            DDY=DDY1*9.0/8.0-DDY2/8.0

            W1A=(W(I,J,K  )+W(I-1,J,K  ))/2.0
            W1B=(W(I,J,K+1)+W(I-1,J,K+1))/2.0
            W3A=(W(I+1,J,K  )+W(I-2,J,K  ))/2.0
            W3B=(W(I+1,J,K+1)+W(I-2,J,K+1))/2.0
            U1A=(U(I,J,K)+U(I,J,K-1))/2.0
            U1B=(U(I,J,K)+U(I,J,K+1))/2.0            
            DDZ1=((W1B*9.0/8.0-W3B/8.0)*U1B-(W1A*9.0/8.0-W3A/8.0)*U1A)/DZ
            W1A=(W(I,J,K-1)+W(I-1,J,K-1))/2.0
            W1B=(W(I,J,K+2)+W(I-1,J,K+2))/2.0
            W3A=(W(I+1,J,K-1)+W(I-2,J,K-1))/2.0
            W3B=(W(I+1,J,K+2)+W(I-2,J,K+2))/2.0
            U3A=(U(I,J,K)+U(I,J,K-3))/2.0
            U3B=(U(I,J,K)+U(I,J,K+3))/2.0            
            DDZ2=((W1B*9.0/8.0-W3B/8.0)*U3B-(W1A*9.0/8.0-W3A/8.0)*U3A)/(DZ*3.0)
            DDZ=DDZ1*9.0/8.0-DDZ2/8.0
            
            FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------
            U1A=(U(I  ,J,K)+U(I  ,J-1,K))/2.0
            U1B=(U(I+1,J,K)+U(I+1,J-1,K))/2.0
            U3A=(U(I  ,J+1,K)+U(I  ,J-2,K))/2.0
            U3B=(U(I+1,J+1,K)+U(I+1,J-2,K))/2.0
            V1A=(V(I,J,K)+V(I-1,J,K))/2.0
            V1B=(V(I,J,K)+V(I+1,J,K))/2.0            
            DDX1=((U1B*9.0/8.0-U3B/8.0)*V1B-(U1A*9.0/8.0-U3A/8.0)*V1A)/DX
            U1A=(U(I-1,J,K)+U(I-1,J-1,K))/2.0
            U1B=(U(I+2,J,K)+U(I+2,J-1,K))/2.0
            U3A=(U(I-1,J+1,K)+U(I-1,J-2,K))/2.0
            U3B=(U(I+2,J+1,K)+U(I+2,J-2,K))/2.0
            V3A=(V(I,J,K)+V(I-3,J,K))/2.0
            V3B=(V(I,J,K)+V(I+3,J,K))/2.0            
            DDX2=((U1B*9.0/8.0-U3B/8.0)*V3B-(U1A*9.0/8.0-U3A/8.0)*V3A)/(DX*3.0)
            DDX=DDX1*9.0/8.0-DDX2/8.0

            V1A=(V(I,J,K)+V(I,J-1,K))/2.0
            V1B=(V(I,J,K)+V(I,J+1,K))/2.0
            V3A=(V(I,J+1,K)+V(I,J-2,K))/2.0
            V3B=(V(I,J-1,K)+V(I,J+2,K))/2.0
            DDY1=((V1B*9.0/8.0-V3B/8.0)*V1B-(V1A*9.0/8.0-V3A/8.0)*V1A)/DY
            V1A=(V(I,J-1,K)+V(I,J-2,K))/2.0
            V1B=(V(I,J+1,K)+V(I,J+2,K))/2.0
            V3A=(V(I,J,K)+V(I,J-3,K))/2.0
            V3B=(V(I,J,K)+V(I,J+3,K))/2.0
            DDY2=((V1B*9.0/8.0-V3B/8.0)*V3B-(V1A*9.0/8.0-V3A/8.0)*V3A)/(DY*3.0)
            DDY=DDY1*9.0/8.0-DDY2/8.0

            W1A=(W(I,J,K  )+W(I,J-1,K  ))/2.0
            W1B=(W(I,J,K+1)+W(I,J-1,K+1))/2.0
            W3A=(W(I,J+1,K  )+W(I,J-2,K  ))/2.0
            W3B=(W(I,J+1,K+1)+W(I,J-2,K+1))/2.0
            V1A=(V(I,J,K)+V(I,J,K-1))/2.0
            V1B=(V(I,J,K)+V(I,J,K+1))/2.0            
            DDZ1=((W1B*9.0/8.0-W3B/8.0)*V1B-(W1A*9.0/8.0-W3A/8.0)*V1A)/DZ
            W1A=(W(I,J,K-1)+W(I,J-1,K-1))/2.0
            W1B=(W(I,J,K+2)+W(I,J-1,K+2))/2.0
            W3A=(W(I,J+1,K-1)+W(I,J-2,K-1))/2.0
            W3B=(W(I,J+1,K+2)+W(I,J-2,K+2))/2.0
            V3A=(V(I,J,K)+V(I,J,K-3))/2.0
            V3B=(V(I,J,K)+V(I,J,K+3))/2.0            
            DDZ2=((W1B*9.0/8.0-W3B/8.0)*V3B-(W1A*9.0/8.0-W3A/8.0)*V3A)/(DZ*3.0)
            DDZ=DDZ1*9.0/8.0-DDZ2/8.0
            
            FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------
            U1A=(U(I  ,J,K)+U(I  ,J,K-1))/2.0
            U1B=(U(I+1,J,K)+U(I+1,J,K-1))/2.0
            U3A=(U(I  ,J,K+1)+U(I  ,J,K-2))/2.0
            U3B=(U(I+1,J,K+1)+U(I+1,J,K-2))/2.0
            W1A=(W(I,J,K)+W(I-1,J,K))/2.0
            W1B=(W(I,J,K)+W(I+1,J,K))/2.0            
            DDX1=((U1B*9.0/8.0-U3B/8.0)*W1B-(U1A*9.0/8.0-U3A/8.0)*W1A)/DX
            U1A=(U(I-1,J,K)+U(I-1,J,K-1))/2.0
            U1B=(U(I+2,J,K)+U(I+2,J,K-1))/2.0
            U3A=(U(I-1,J,K+1)+U(I-1,J,K-2))/2.0
            U3B=(U(I+2,J,K+1)+U(I+2,J,K-2))/2.0
            W3A=(W(I,J,K)+W(I-3,J,K))/2.0
            W3B=(W(I,J,K)+W(I+3,J,K))/2.0            
            DDX2=((U1B*9.0/8.0-U3B/8.0)*W3B-(U1A*9.0/8.0-U3A/8.0)*W3A)/(DX*3.0)
            DDX=DDX1*9.0/8.0-DDX2/8.0

            V1A=(V(I,J  ,K)+V(I,J  ,K-1))/2.0
            V1B=(V(I,J+1,K)+V(I,J+1,K-1))/2.0
            V3A=(V(I,J  ,K+1)+V(I,J  ,K-2))/2.0
            V3B=(V(I,J+1,K+1)+V(I,J+1,K-2))/2.0
            W1A=(W(I,J,K)+W(I,J-1,K))/2.0
            W1B=(W(I,J,K)+W(I,J+1,K))/2.0            
            DDY1=((V1B*9.0/8.0-V3B/8.0)*W1B-(V1A*9.0/8.0-V3A/8.0)*W1A)/DY
            V1A=(V(I,J-1,K)+V(I,J-1,K-1))/2.0
            V1B=(V(I,J+2,K)+V(I,J+2,K-1))/2.0
            V3A=(V(I,J-1,K+1)+V(I,J-1,K-2))/2.0
            V3B=(V(I,J+2,K+1)+V(I,J+2,K-2))/2.0
            W3A=(W(I,J,K)+W(I,J-3,K))/2.0
            W3B=(W(I,J,K)+W(I,J+3,K))/2.0            
            DDY2=((V1B*9.0/8.0-V3B/8.0)*W3B-(V1A*9.0/8.0-V3A/8.0)*W3A)/(DY*3.0)
            DDY=DDY1*9.0/8.0-DDY2/8.0

            W1A=(W(I,J,K)+W(I,J,K-1))/2.0
            W1B=(W(I,J,K)+W(I,J,K+1))/2.0
            W3A=(W(I,J,K+1)+W(I,J,K-2))/2.0
            W3B=(W(I,J,K-1)+W(I,J,K+2))/2.0
            DDZ1=((W1B*9.0/8.0-W3B/8.0)*W1B-(W1A*9.0/8.0-W3A/8.0)*W1A)/DZ
            W1A=(W(I,J,K-1)+W(I,J,K-2))/2.0
            W1B=(W(I,J,K+1)+W(I,J,K+2))/2.0
            W3A=(W(I,J,K)+W(I,J,K-3))/2.0
            W3B=(W(I,J,K)+W(I,J,K+3))/2.0
            DDZ2=((W1B*9.0/8.0-W3B/8.0)*W3B-(W1A*9.0/8.0-W3A/8.0)*W3A)/(DZ*3.0)
            DDZ=DDZ1*9.0/8.0-DDZ2/8.0

            FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!*************************************************************************! 
!            SUBROUTINE OF CALCULATING THE CONVECTIVE TERMS               ! 
!*************************************************************************! 
!     4TH ORDER ACCURACY, ADVECTIVE FORM (Oleg V. Vasilyev, JCP,  2000)     
      SUBROUTINE GETCON4_AD(DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
 
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: FCX,FCY,FCZ
      REAL(KIND=DP):: U1A,U1B,U3A,U3B,V1A,V1B,V3A,V3B,W1A,W1B,W3A,W3B
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ
      REAL(KIND=DP):: DDX1A,DDX1B,DDX3A,DDX3B
      REAL(KIND=DP):: DDY1A,DDY1B,DDY3A,DDY3B
      REAL(KIND=DP):: DDZ1A,DDZ1B,DDZ3A,DDZ3B
      INTEGER:: I,J,K

!-----GET GRID SPACING IN UNIFORM COMPUTATIONAL DOMAIN                                                                         
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----FDX---------------------------------------------          
            DDX1A=(U(I,J,K)-U(I-1,J,K))/DX*          &
                 ((U(I,J,K)+U(I-1,J,K))/2.0*9.0/8.0-  &
                  (U(I+1,J,K)+U(I-2,J,K))/2.0/8.0)
            DDX1B=(U(I+1,J,K)-U(I,J,K))/DX*          &
                 ((U(I+1,J,K)+U(I,J,K))/2.0*9.0/8.0-  &
                  (U(I+2,J,K)+U(I-1,J,K))/2.0/8.0)
            DDX1=(DDX1B+DDX1A)/2.0
            DDX3A=(U(I,J,K)-U(I-3,J,K))/(DX*3.0)*    &
                 ((U(I-1,J,K)+U(I-2,J,K))/2.0*9.0/8.0-&
                  (U(I,J,K)+U(I-3,J,K))/2.0/8.0)
            DDX3B=(U(I+3,J,K)-U(I,J,K))/(DX*3.0)*    &
                 ((U(I+2,J,K)+U(I+1,J,K))/2.0*9.0/8.0-&
                  (U(I+3,J,K)+U(I,J,K))/2.0/8.0)
            DDX2=(DDX3B+DDX3A)/2.0
            DDX=(DDX1*9.0/8.0-DDX2/8.0)

            V1A=(V(I,J,K)+V(I-1,J,K))/2.0
            V3A=(V(I+1,J,K)+V(I-2,J,K))/2.0
            V1B=(V(I,J+1,K)+V(I-1,J+1,K))/2.0
            V3B=(V(I+1,J+1,K)+V(I-2,J+1,K))/2.0
            DDY1A=(U(I,J,K)-U(I,J-1,K))/DY
            DDY1B=(U(I,J+1,K)-U(I,J,K))/DY
            DDY1=(DDY1B*(V1B*9.0/8.0-V3B/8.0)+ &
                  DDY1A*(V1A*9.0/8.0-V3A/8.0))/2.0
            V1A=(V(I,J-1,K)+V(I-1,J-1,K))/2.0
            V3A=(V(I+1,J-1,K)+V(I-2,J-1,K))/2.0
            V1B=(V(I,J+2,K)+V(I-1,J+2,K))/2.0
            V3B=(V(I+1,J+2,K)+V(I-2,J+2,K))/2.0
            DDY3A=(U(I,J,K)-U(I,J-3,K))/(DY*3.0)
            DDY3B=(U(I,J+3,K)-U(I,J,K))/(DY*3.0)
            DDY2=(DDY3B*(V1B*9.0/8.0-V3B/8.0)+ &
                  DDY3A*(V1A*9.0/8.0-V3A/8.0))/2.0
            DDY=(DDY1*9.0/8.0-DDY2/8.0)

            W1A=(W(I,J,K)+W(I-1,J,K))/2.0
            W3A=(W(I+1,J,K)+W(I-2,J,K))/2.0
            W1B=(W(I,J,K+1)+W(I-1,J,K+1))/2.0
            W3B=(W(I+1,J,K+1)+W(I-2,J,K+1))/2.0
            DDZ1A=(U(I,J,K)-U(I,J,K-1))/DZ
            DDZ1B=(U(I,J,K+1)-U(I,J,K))/DZ
            DDZ1=(DDZ1B*(W1B*9.0/8.0-W3B/8.0)+ &
                  DDZ1A*(W1A*9.0/8.0-W3A/8.0))/2.0
            W1A=(W(I,J,K-1)+W(I-1,J,K-1))/2.0
            W3A=(W(I+1,J,K-1)+W(I-2,J,K-1))/2.0
            W1B=(W(I,J,K+2)+W(I-1,J,K+2))/2.0
            W3B=(W(I+1,J,K+2)+W(I-2,J,K+2))/2.0
            DDZ3A=(U(I,J,K)-U(I,J,K-3))/(DZ*3.0)
            DDZ3B=(U(I,J,K+3)-U(I,J,K))/(DZ*3.0)
            DDZ2=(DDZ3B*(W1B*9.0/8.0-W3B/8.0)+ &
                  DDZ3A*(W1A*9.0/8.0-W3A/8.0))/2.0
            DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

            FCX(I,J,K)=DDX+DDY+DDZ
!-----FDY---------------------------------------------  
            U1A=(U(I,J,K)+U(I,J-1,K))/2.0
            U3A=(U(I,J+1,K)+U(I,J-2,K))/2.0
            U1B=(U(I+1,J,K)+U(I+1,J-1,K))/2.0
            U3B=(U(I+1,J+1,K)+U(I+1,J-2,K))/2.0
            DDX1A=(V(I,J,K)-V(I-1,J,K))/DX
            DDX1B=(V(I+1,J,K)-V(I,J,K))/DX
            DDX1=(DDX1B*(U1B*9.0/8.0-U3B/8.0)+ &
                  DDX1A*(U1A*9.0/8.0-U3A/8.0))/2.0
            U1A=(U(I-1,J,K)+U(I-1,J-1,K))/2.0
            U3A=(U(I-1,J+1,K)+U(I-1,J-2,K))/2.0
            U1B=(U(I+2,J,K)+U(I+2,J-1,K))/2.0
            U3B=(U(I+2,J+1,K)+U(I+2,J-2,K))/2.0
            DDX3A=(V(I,J,K)-V(I-3,J,K))/(DX*3.0)
            DDX3B=(V(I+3,J,K)-V(I,J,K))/(DX*3.0)
            DDX2=(DDX3B*(U1B*9.0/8.0-U3B/8.0)+ &
                  DDX3A*(U1A*9.0/8.0-U3A/8.0))/2.0
            DDX=(DDX1*9.0/8.0-DDX2/8.0)


            DDY1A=(V(I,J,K)-V(I,J-1,K))/DY*           &
                 ((V(I,J,K)+V(I,J-1,K))/2.0*9.0/8.0-   &
                  (V(I,J+1,K)+V(I,J-2,K))/2.0/8.0)
            DDY1B=(V(I,J+1,K)-V(I,J,K))/DY*           &
                 ((V(I,J+1,K)+V(I,J,K))/2.0*9.0/8.0-   &
                  (V(I,J+2,K)+V(I,J-1,K))/2.0/8.0)
            DDY1=(DDY1B+DDY1A)/2.0
            DDY3A=(V(I,J,K)-V(I,J-3,K))/(DY*3.0)*     &
                 ((V(I,J-1,K)+V(I,J-2,K))/2.0*9.0/8.0- &
                  (V(I,J,K)+V(I,J-3,K))/2.0/8.0)
            DDY3B=(V(I,J+3,K)-V(I,J,K))/(DY*3.0)*     &
                 ((V(I,J+2,K)+V(I,J+1,K))/2.0*9.0/8.0- &
                  (V(I,J+3,K)+V(I,J,K))/2.0/8.0)
            DDY2=(DDY3B+DDY3A)/2.0
            DDY=(DDY1*9.0/8.0-DDY2/8.0)


            W1A=(W(I,J,K)+W(I,J-1,K))/2.0
            W3A=(W(I,J+1,K)+W(I,J-2,K))/2.0
            W1B=(W(I,J,K+1)+W(I,J-1,K+1))/2.0
            W3B=(W(I,J+1,K+1)+W(I,J-2,K+1))/2.0
            DDZ1A=(V(I,J,K)-V(I,J,K-1))/DZ
            DDZ1B=(V(I,J,K+1)-V(I,J,K))/DZ
            DDZ1=(DDZ1B*(W1B*9.0/8.0-W3B/8.0)+ &
                  DDZ1A*(W1A*9.0/8.0-W3A/8.0))/2.0
            W1A=(W(I,J,K-1)+W(I,J-1,K-1))/2.0
            W3A=(W(I,J+1,K-1)+W(I,J-2,K-1))/2.0
            W1B=(W(I,J,K+2)+W(I,J-1,K+2))/2.0
            W3B=(W(I,J+1,K+2)+W(I,J-2,K+2))/2.0
            DDZ3A=(V(I,J,K)-V(I,J,K-3))/(DZ*3.0)
            DDZ3B=(V(I,J,K+3)-V(I,J,K))/(DZ*3.0)
            DDZ2=(DDZ3B*(W1B*9.0/8.0-W3B/8.0)+ &
                  DDZ3A*(W1A*9.0/8.0-W3A/8.0))/2.0
            DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

            FCY(I,J,K)=DDX+DDY+DDZ
!-----FDZ---------------------------------------------     
            U1A=(U(I,J,K)+U(I,J,K-1))/2.0
            U3A=(U(I,J,K+1)+U(I,J,K-2))/2.0
            U1B=(U(I+1,J,K)+U(I+1,J,K-1))/2.0
            U3B=(U(I+1,J,K+1)+U(I+1,J,K-2))/2.0
            DDX1A=(W(I,J,K)-W(I-1,J,K))/DX
            DDX1B=(W(I+1,J,K)-W(I,J,K))/DX
            DDX1=(DDX1B*(U1B*9.0/8.0-U3B/8.0)+ &
                  DDX1A*(U1A*9.0/8.0-U3A/8.0))/2.0
            U1A=(U(I-1,J,K)+U(I-1,J,K-1))/2.0
            U3A=(U(I-1,J,K+1)+U(I-1,J,K-2))/2.0
            U1B=(U(I+2,J,K)+U(I+2,J,K-1))/2.0
            U3B=(U(I+2,J,K+1)+U(I+2,J,K-2))/2.0
            DDX3A=(W(I,J,K)-W(I-3,J,K))/(DX*3.0)
            DDX3B=(W(I+3,J,K)-W(I,J,K))/(DX*3.0)
            DDX2=(DDX3B*(U1B*9.0/8.0-U3B/8.0)+ &
                  DDX3A*(U1A*9.0/8.0-U3A/8.0))/2.0
            DDX=(DDX1*9.0/8.0-DDX2/8.0)

            V1A=(V(I,J,K)+V(I,J,K-1))/2.0
            V3A=(V(I,J,K+1)+V(I,J,K-2))/2.0
            V1B=(V(I,J+1,K)+V(I,J+1,K-1))/2.0
            V3B=(V(I,J+1,K+1)+V(I,J+1,K-2))/2.0
            DDY1A=(W(I,J,K)-W(I,J-1,K))/DY
            DDY1B=(W(I,J+1,K)-W(I,J,K))/DY
            DDY1=(DDY1B*(V1B*9.0/8.0-V3B/8.0)+ &
                  DDY1A*(V1A*9.0/8.0-V3A/8.0))/2.0
            V1A=(V(I,J-1,K)+V(I,J-1,K-1))/2.0
            V3A=(V(I,J-1,K+1)+V(I,J-1,K-2))/2.0
            V1B=(V(I,J+2,K)+V(I,J+2,K-1))/2.0
            V3B=(V(I,J+2,K+1)+V(I,J+2,K-2))/2.0
            DDY3A=(W(I,J,K)-W(I,J-3,K))/(DY*3.0)
            DDY3B=(W(I,J+3,K)-W(I,J,K))/(DY*3.0)
            DDY2=(DDY3B*(V1B*9.0/8.0-V3B/8.0)+         &
                  DDY3A*(V1A*9.0/8.0-V3A/8.0))/2.0
            DDY=(DDY1*9.0/8.0-DDY2/8.0)


            DDZ1A=(W(I,J,K)-W(I,J,K-1))/DZ*           &
                 ((W(I,J,K)+W(I,J,K-1))/2.0*9.0/8.0-   &
                  (W(I,J,K+1)+W(I,J,K-2))/2.0/8.0)
            DDZ1B=(W(I,J,K+1)-W(I,J,K))/DZ*           &
                 ((W(I,J,K+1)+W(I,J,K))/2.0*9.0/8.0-   &
                  (W(I,J,K+2)+W(I,J,K-1))/2.0/8.0)
            DDZ1=(DDZ1B+DDZ1A)/2.0
            DDZ3A=(W(I,J,K)-W(I,J,K-3))/(DZ*3.0)*     &
                 ((W(I,J,K-1)+W(I,J,K-2))/2.0*9.0/8.0- &
                  (W(I,J,K)+W(I,J,K-3))/2.0/8.0)
            DDZ3B=(W(I,J,K+3)-W(I,J,K))/(DZ*3.0)*     &
                 ((W(I,J,K+2)+W(I,J,K+1))/2.0*9.0/8.0- &
                  (W(I,J,K+3)+W(I,J,K))/2.0/8.0)
            DDZ2=(DDZ3B+DDZ3A)/2.0
            DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

            FCZ(I,J,K)=DDX+DDY+DDZ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

  END MODULE
