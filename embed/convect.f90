! This module is used for calculation of the unlinear convective term
!
  MODULE convect

  USE parameters
  USE field_shared, ONLY: U,V,W,UF,VF,WF,FX,FY,FZ,CHECK,CHECK0
  USE tools
  USE class_shared
  USE class_cell

  IMPLICIT NONE
  
  CONTAINS
!=========================================================================!
!                    MAIN SUBROUTINE FOR THE CONVECTIVE TERM              !
!=========================================================================!
    SUBROUTINE CONVECT_WRAP(DX,DY,DZ)
    
    IMPLICIT NONE

    REAL(KIND=DP):: FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2
 
    REAL(KIND=DP):: DX,DY,DZ  
    INTEGER :: I,J,K


    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
!---------COLLOCATED GRID---------------------------------------            
          IF(ICOLL.EQ.1)THEN     
             CALL GETCON_DIV_COL(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                 DX,DY,DZ,FCX1,FCY1,FCZ1)
             CALL GETCON_ADV_COL(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                 DX,DY,DZ,FCX2,FCY2,FCZ2)
!---------STAGGERED GRID----------------------------------------             
          ELSE                   
            IF(ORDER_CON.EQ.4)THEN  ! 4TH ORDER
              CALL GETCON4_DIV_STA(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                   DX,DY,DZ,FCX1,FCY1,FCZ1)
              CALL GETCON4_ADV_STA(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                   DX,DY,DZ,FCX2,FCY2,FCZ2)
            ELSE IF(ORDER_CON.EQ.2)THEN  ! 2ND ORDER
              CALL GETCON2_DIV_STA(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                   DX,DY,DZ,FCX1,FCY1,FCZ1)
              CALL GETCON2_ADV_STA(U,V,W,NX1,NY1,NZ1,I,J,K, &
                                   DX,DY,DZ,FCX2,FCY2,FCZ2)                
            END IF
          END IF
 
          FX(I,J,K)=FX(I,J,K)-(FCX1*BLEND_CON+FCX2*(1.0-BLEND_CON))
          FY(I,J,K)=FY(I,J,K)-(FCY1*BLEND_CON+FCY2*(1.0-BLEND_CON))    
          FZ(I,J,K)=FZ(I,J,K)-(FCZ1*BLEND_CON+FCZ2*(1.0-BLEND_CON))
        END DO
      END DO 
    END DO
   
    END SUBROUTINE CONVECT_WRAP
!=========================================================================!
!    GET CONVECTIVE TERM, 4TH ORDER, DIVERGENCE FORM, COLLOCATED GRID     !
!=========================================================================!
      SUBROUTINE GETCON_DIV_COL(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE

      INTEGER I,J,K,NB1,NB2,NB3
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: DX,DY,DZ,FCX,FCY,FCZ      

      FCX=CON_DIV(VAR1,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)
      FCY=CON_DIV(VAR2,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)
      FCZ=CON_DIV(VAR3,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)

      END SUBROUTINE 
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, ADVECTIVE FORM, COLLOCATED GRID     !
!=========================================================================!
!     4TH ORDER ACCURACY, ADVECTIVE FORM (Oleg V. Vasilyev, JCP,  2000)     
      SUBROUTINE GETCON_ADV_COL(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE
        
      INTEGER I,J,K,NB1,NB2,NB3
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: DX,DY,DZ,FCX,FCY,FCZ


      FCX=CON_ADV(VAR1,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)
      FCY=CON_ADV(VAR2,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)
      FCZ=CON_ADV(VAR3,NB1,NB2,NB3,DX,DY,DZ,I,J,K,ORDER_CON)

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, ADVECTIVE FORM, STAGGERED GRID      !
!=========================================================================!     
      SUBROUTINE GETCON4_ADV_STA(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                 DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
      INTEGER:: NB1,NB2,NB3,I,J,K        
      REAL(KIND=DP):: DX,DY,DZ,FCX,FCY,FCZ
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ


!-----FDX---------------------------------------------          
      DDX1=(STENCIL_DX(VAR1,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I,  J,K,1,4)+  &
            STENCIL_DX(VAR1,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I+1,J,K,1,4))/2.0
      DDX2=(STENCIL_DX(VAR1,NB1,NB2,NB3,I-1,J,K,3,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)+  &
            STENCIL_DX(VAR1,NB1,NB2,NB3,I+2,J,K,3,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I+2,J,K,1,4))/2.0
      DDX=(DDX1*9.0/8.0-DDX2/8.0)

      DDY1=(STENCIL_DY(VAR1,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J,  K,1,4)+  &
            STENCIL_DY(VAR1,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+1,K,1,4))/2.0
      DDY2=(STENCIL_DY(VAR1,NB1,NB2,NB3,I,J-1,K,3,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)+  &
            STENCIL_DY(VAR1,NB1,NB2,NB3,I,J+2,K,3,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+2,K,1,4))/2.0
      DDY=(DDY1*9.0/8.0-DDY2/8.0)

      DDZ1=(STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K,  1,4)+  &
            STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+1,1,4))/2.0
      DDZ2=(STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K-1,3,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)+  &
            STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K+2,3,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+2,1,4))/2.0
      DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

      FCX=DDX+DDY+DDZ
!-----FDY---------------------------------------------  
      DDX1=(STENCIL_DX(VAR2,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I,  J,K,1,4)+  &
            STENCIL_DX(VAR2,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+1,J,K,1,4))/2.0
      DDX2=(STENCIL_DX(VAR2,NB1,NB2,NB3,I-1,J,K,3,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)+  &
            STENCIL_DX(VAR2,NB1,NB2,NB3,I+2,J,K,3,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+2,J,K,1,4))/2.0
      DDX=(DDX1*9.0/8.0-DDX2/8.0)

      DDY1=(STENCIL_DY(VAR2,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J,  K,1,4)+  &
            STENCIL_DY(VAR2,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+1,K,1,4))/2.0
      DDY2=(STENCIL_DY(VAR2,NB1,NB2,NB3,I,J-1,K,3,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)+  &
            STENCIL_DY(VAR2,NB1,NB2,NB3,I,J+2,K,3,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+2,K,1,4))/2.0
      DDY=(DDY1*9.0/8.0-DDY2/8.0)

      DDZ1=(STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K,  1,4)+  &
            STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+1,1,4))/2.0
      DDZ2=(STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K-1,3,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)+  &
            STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K+2,3,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+2,1,4))/2.0
      DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

      FCY=DDX+DDY+DDZ
!-----FDZ---------------------------------------------     
      DDX1=(STENCIL_DX(VAR3,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I,  J,K,1,4)+  &
            STENCIL_DX(VAR3,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+1,J,K,1,4))/2.0
      DDX2=(STENCIL_DX(VAR3,NB1,NB2,NB3,I-1,J,K,3,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)+  &
            STENCIL_DX(VAR3,NB1,NB2,NB3,I+2,J,K,3,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+2,J,K,1,4))/2.0
      DDX=(DDX1*9.0/8.0-DDX2/8.0)

      DDY1=(STENCIL_DY(VAR3,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J,  K,1,4)+  &
            STENCIL_DY(VAR3,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+1,K,1,4))/2.0
      DDY2=(STENCIL_DY(VAR3,NB1,NB2,NB3,I,J-1,K,3,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)+  &
            STENCIL_DY(VAR3,NB1,NB2,NB3,I,J+2,K,3,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+2,K,1,4))/2.0
      DDY=(DDY1*9.0/8.0-DDY2/8.0)

      DDZ1=(STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K,  1,4)+  &
            STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+1,1,4))/2.0
      DDZ2=(STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K-1,3,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)+  &
            STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K+2,3,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+2,1,4))/2.0
      DDZ=(DDZ1*9.0/8.0-DDZ2/8.0)

      FCZ=DDX+DDY+DDZ

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, DIVERGENCE FORM, STAGGERED GRID     !
!=========================================================================!
      SUBROUTINE GETCON4_DIV_STA(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                 DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE
        
      INTEGER:: NB1,NB2,NB3,I,J,K
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: U1A,U1B,U3A,U3B,V1A,V1B,V3A,V3B,W1A,W1B,W3A,W3B
      REAL(KIND=DP):: DDX1,DDX2,DDX,DDY1,DDY2,DDY,DDZ1,DDZ2,DDZ
      REAL(KIND=DP):: FCX,FCY,FCZ

!-----FDX---------------------------------------------
      DDX1=(VAR_INTER_X(VAR1,NB1,NB2,NB3,I+1,J,K,1,4)*STENCIL_IX(VAR1,NB1,NB2,NB3,I+1,J,K,1)-  &
            VAR_INTER_X(VAR1,NB1,NB2,NB3,I,  J,K,1,4)*STENCIL_IX(VAR1,NB1,NB2,NB3,I,  J,K,1))/DX
      DDX2=(VAR_INTER_X(VAR1,NB1,NB2,NB3,I+2,J,K,1,4)*STENCIL_IX(VAR1,NB1,NB2,NB3,I+2,J,K,3)-  &
            VAR_INTER_X(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)*STENCIL_IX(VAR1,NB1,NB2,NB3,I-1,J,K,3))/(DX*3.0)
      DDX=DDX1*9.0/8.0-DDX2/8.0

      DDY1=(VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+1,K,1,4)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J+1,K,1)-  &
            VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J,  K,1,4)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J,  K,1))/DY
      DDY2=(VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+2,K,1,4)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J+2,K,3)-  &
            VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J-1,K,3))/(DY*3.0)
      DDY=DDY1*9.0/8.0-DDY2/8.0

      DDZ1=(VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+1,1,4)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K+1,1)-  &
            VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K,  1,4)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K,  1))/DZ
      DDZ2=(VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+2,1,4)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K+2,3)-  &
            VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K-1,3))/(DZ*3.0)
      DDZ=DDZ1*9.0/8.0-DDZ2/8.0
            
      FCX=DDX+DDY+DDZ
!-----FDY---------------------------------------------
      DDX1=(VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+1,J,K,1,4)*STENCIL_IX(VAR2,NB1,NB2,NB3,I+1,J,K,1)-  &
            VAR_INTER_Y(VAR1,NB1,NB2,NB3,I,  J,K,1,4)*STENCIL_IX(VAR2,NB1,NB2,NB3,I,  J,K,1))/DX
      DDX2=(VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+2,J,K,1,4)*STENCIL_IX(VAR2,NB1,NB2,NB3,I+2,J,K,3)-  &
            VAR_INTER_Y(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)*STENCIL_IX(VAR2,NB1,NB2,NB3,I-1,J,K,3))/(DX*3.0)
      DDX=DDX1*9.0/8.0-DDX2/8.0

      DDY1=(VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+1,K,1,4)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J+1,K,1)-  &
            VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J,  K,1,4)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J,  K,1))/DY
      DDY2=(VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+2,K,1,4)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J+2,K,3)-  &
            VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J-1,K,3))/(DY*3.0)
      DDY=DDY1*9.0/8.0-DDY2/8.0

      DDZ1=(VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+1,1,4)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K+1,1)-  &
            VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K,  1,4)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K,  1))/DZ
      DDZ2=(VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+2,1,4)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K+2,3)-  &
            VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K-1,3))/(DZ*3.0)
      DDZ=DDZ1*9.0/8.0-DDZ2/8.0    
       
      FCY=DDX+DDY+DDZ
!-----FDZ---------------------------------------------
      DDX1=(VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+1,J,K,1,4)*STENCIL_IX(VAR3,NB1,NB2,NB3,I+1,J,K,1)-  &
            VAR_INTER_Z(VAR1,NB1,NB2,NB3,I,  J,K,1,4)*STENCIL_IX(VAR3,NB1,NB2,NB3,I,  J,K,1))/DX
      DDX2=(VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+2,J,K,1,4)*STENCIL_IX(VAR3,NB1,NB2,NB3,I+2,J,K,3)-  &
            VAR_INTER_Z(VAR1,NB1,NB2,NB3,I-1,J,K,1,4)*STENCIL_IX(VAR3,NB1,NB2,NB3,I-1,J,K,3))/(DX*3.0)
      DDX=DDX1*9.0/8.0-DDX2/8.0

      DDY1=(VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+1,K,1,4)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J+1,K,1)-  &
            VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J,  K,1,4)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J,  K,1))/DY
      DDY2=(VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+2,K,1,4)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J+2,K,3)-  &
            VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J-1,K,1,4)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J-1,K,3))/(DY*3.0)
      DDY=DDY1*9.0/8.0-DDY2/8.0

      DDZ1=(VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+1,1,4)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K+1,1)-  &
            VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K,  1,4)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K,  1))/DZ
      DDZ2=(VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+2,1,4)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K+2,3)-  &
            VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K-1,1,4)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K-1,3))/(DZ*3.0)
      DDZ=DDZ1*9.0/8.0-DDZ2/8.0

      FCZ=DDX+DDY+DDZ

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 2ND ORDER, ADVECTIVE FORM, STAGGERED GRID      !
!=========================================================================!     
      SUBROUTINE GETCON2_ADV_STA(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                 DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
        
      INTEGER:: I,J,K,NB1,NB2,NB3
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: DDX,DDY,DDZ,FCX,FCY,FCZ

!-----FDX---------------------------------------------          
      DDX=(STENCIL_DX(VAR1,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I,  J,K,1,2)+  &
           STENCIL_DX(VAR1,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_X(VAR1,NB1,NB2,NB3,I+1,J,K,1,2))/2.0

      DDY=(STENCIL_DY(VAR1,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J,  K,1,2)+  &
           STENCIL_DY(VAR1,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+1,K,1,2))/2.0

      DDZ=(STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K,1  ,2)+  &
           STENCIL_DZ(VAR1,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+1,1,2))/2.0

      FCX=DDX+DDY+DDZ
!-----FDY---------------------------------------------  
      DDX=(STENCIL_DX(VAR2,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I,  J,K,1,2)+  &
           STENCIL_DX(VAR2,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+1,J,K,1,2))/2.0

      DDY=(STENCIL_DY(VAR2,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J,  K,1,2)+  &
           STENCIL_DY(VAR2,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+1,K,1,2))/2.0

      DDZ=(STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K,1  ,2)+  &
           STENCIL_DZ(VAR2,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+1,1,2))/2.0

      FCY=DDX+DDY+DDZ
!-----FDZ---------------------------------------------     
      DDX=(STENCIL_DX(VAR3,NB1,NB2,NB3,I,  J,K,1,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I,  J,K,1,2)+  &
           STENCIL_DX(VAR3,NB1,NB2,NB3,I+1,J,K,1,DX)*VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+1,J,K,1,2))/2.0

      DDY=(STENCIL_DY(VAR3,NB1,NB2,NB3,I,J,  K,1,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J,  K,1,2)+  &
           STENCIL_DY(VAR3,NB1,NB2,NB3,I,J+1,K,1,DY)*VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+1,K,1,2))/2.0

      DDZ=(STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K,  1,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K,1  ,2)+  &
           STENCIL_DZ(VAR3,NB1,NB2,NB3,I,J,K+1,1,DZ)*VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+1,1,2))/2.0

      FCZ=DDX+DDY+DDZ

      END SUBROUTINE
!=========================================================================!
!     GET CONVECTIVE TERM, 2ND ORDER, DIVERGENCE FORM, STAGGERED GRID     !
!=========================================================================!
      SUBROUTINE GETCON2_DIV_STA(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                 DX,DY,DZ,FCX,FCY,FCZ)

      IMPLICIT NONE
        
      INTEGER:: I,J,K,NB1,NB2,NB3
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(NB1:,NB2:,NB3:):: VAR1,VAR2,VAR3
      REAL(KIND=DP):: DDX,DDY,DDZ,FCX,FCY,FCZ

!-----FDX---------------------------------------------
      DDX=(VAR_INTER_X(VAR1,NB1,NB2,NB3,I+1,J,K,1,2)*STENCIL_IX(VAR1,NB1,NB2,NB3,I+1,J,K,1)-  &
           VAR_INTER_X(VAR1,NB1,NB2,NB3,I,  J,K,1,2)*STENCIL_IX(VAR1,NB1,NB2,NB3,I,  J,K,1))/DX

      DDY=(VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J+1,K,1,2)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J+1,K,1)-  &
           VAR_INTER_X(VAR2,NB1,NB2,NB3,I,J,  K,1,2)*STENCIL_IY(VAR1,NB1,NB2,NB3,I,J,  K,1))/DY

      DDZ=(VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K+1,1,2)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K+1,1)-  &
           VAR_INTER_X(VAR3,NB1,NB2,NB3,I,J,K,  1,2)*STENCIL_IZ(VAR1,NB1,NB2,NB3,I,J,K,  1))/DZ
            
      FCX=DDX+DDY+DDZ
!-----FDY---------------------------------------------
      DDX=(VAR_INTER_Y(VAR1,NB1,NB2,NB3,I+1,J,K,1,2)*STENCIL_IX(VAR2,NB1,NB2,NB3,I+1,J,K,1)-  &
           VAR_INTER_Y(VAR1,NB1,NB2,NB3,I,  J,K,1,2)*STENCIL_IX(VAR2,NB1,NB2,NB3,I,  J,K,1))/DX
      
      DDY=(VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J+1,K,1,2)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J+1,K,1)-  &
           VAR_INTER_Y(VAR2,NB1,NB2,NB3,I,J,  K,1,2)*STENCIL_IY(VAR2,NB1,NB2,NB3,I,J,  K,1))/DY

      DDZ=(VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K+1,1,2)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K+1,1)-  &
           VAR_INTER_Y(VAR3,NB1,NB2,NB3,I,J,K,  1,2)*STENCIL_IZ(VAR2,NB1,NB2,NB3,I,J,K,  1))/DZ

      FCY=DDX+DDY+DDZ
!-----FDZ---------------------------------------------
      DDX=(VAR_INTER_Z(VAR1,NB1,NB2,NB3,I+1,J,K,1,2)*STENCIL_IX(VAR3,NB1,NB2,NB3,I+1,J,K,1)-  &
           VAR_INTER_Z(VAR1,NB1,NB2,NB3,I,  J,K,1,2)*STENCIL_IX(VAR3,NB1,NB2,NB3,I,  J,K,1))/DX

      DDY=(VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J+1,K,1,2)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J+1,K,1)-  &
           VAR_INTER_Z(VAR2,NB1,NB2,NB3,I,J,  K,1,2)*STENCIL_IY(VAR3,NB1,NB2,NB3,I,J,  K,1))/DY

      DDZ=(VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K+1,1,2)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K+1,1)-  &
           VAR_INTER_Z(VAR3,NB1,NB2,NB3,I,J,K,  1,2)*STENCIL_IZ(VAR3,NB1,NB2,NB3,I,J,K,  1))/DZ

      FCZ=DDX+DDY+DDZ

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
         (UF(I+1,J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I+1,J,K,1)-     &
          UF(I,  J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I,  J,K,1))/DX+ &
         (VF(I,J+1,K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+1,K,1)-     &
          VF(I,J,  K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J,  K,1))/DY+ &
         (WF(I,J,K+1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+1,1)-     &
          WF(I,J,K  )*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K  ,1))/DZ 
      ELSE IF(ORDER.EQ.4)THEN
        CON_DIV=  &
    9.0/8.0*((UF(I+1,J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I+1,J,K,1)-            &
              UF(I,  J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I,  J,K,1))/DX+        &
             (VF(I,J+1,K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+1,K,1)-            &
              VF(I,  J,K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J,  K,1))/DY+        &
             (WF(I,J,K+1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+1,1)-            &
              WF(I,  J,K)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,  1))/DZ)-       &
    1.0/8.0*((UF(I+2,J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I+2,J,K,3)-            &
              UF(I-1,J,K)*STENCIL_IX(VAR,SI1,SI2,SI3,I-1,J,K,3))/(DX*3.0)+  &
             (VF(I,J+2,K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+2,K,3)-            &
              VF(I,J-1,K)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J-1,K,3))/(DY*3.0)+  &
             (WF(I,J,K+2)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+2,3)-            &
              WF(I,J,K-1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K-1,3))/(DZ*3.0))
      END IF

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
        ((UF(I  ,J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I,  J,K,1,DX)+      &
          UF(I+1,J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I+1,J,K,1,DX))/2.0+ &
         (VF(I,J,  K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J,  K,1,DY)+      &
          VF(I,J+1,K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+1,K,1,DY))/2.0+ &
         (WF(I,J,K  )*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,  1,DZ)+      &
          WF(I,J,K+1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+1,1,DZ))/2.0)
      ELSE IF(ORDER.EQ.4)THEN
        CON_ADV=  &
    9.0/8.0*((UF(I,  J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I,  J,K,1,DX)+       &
              UF(I+1,J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I+1,J,K,1,DX))/2.0+  &
             (VF(I,J,  K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J,  K,1,DY)+       &
              VF(I,J+1,K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+1,K,1,DY))/2.0+  &
             (WF(I,J,K  )*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,  1,DZ)+       &
              WF(I,J,K+1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+1,1,DZ))/2.0)- &
    1.0/8.0*((UF(I-1,J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I-1,J,K,3,DX)+       &
              UF(I+2,J,K)*STENCIL_DX(VAR,SI1,SI2,SI3,I+2,J,K,3,DX))/2.0+  &
             (VF(I,J-1,K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J-1,K,3,DY)+       &
              VF(I,J+2,K)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+2,K,3,DY))/2.0+  &
             (WF(I,J,K-1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K-1,3,DZ)+       &
              WF(I,J,K+2)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+2,3,DZ))/2.0)
      END IF

      END FUNCTION
!******************************************************************************!
!               THIS LINE BELOW IS FOR THE FV CELL APPROACH                    !
!******************************************************************************!
      
!=========================================================================!
!            WRAPPING FOR THE CONVECTIVE TERM BASED ON FV CELL            !
!=========================================================================!
    SUBROUTINE CONVECT_CELL_WRAP()
    IMPLICIT NONE
      
    REAL(KIND=DP),DIMENSION(:,:,:) ALLOCATABLE::VAR1,VAR2,VAR3
    REAL(KIND=DP):: FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2
    REAL(KIND=DP):: DXC,DYC,DZC 
    INTEGER:: M,NB

    NB=ORDER_CON-1
    
    ALLOCATE(VAR1(-NB:NB,-NB:NB,-NB:NB),VAR2(-NB:NB,-NB:NB,-NB:NB), &
             VAR3(-NB:NB,-NB:NB,-NB:NB))
    
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND. &
         CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
         
        CALL CELL_TO_STRUCT(M,NB,1,VAR1)
        CALL CELL_TO_STRUCT(M,NB,2,VAR2)
        CALL CELL_TO_STRUCT(M,NB,3,VAR3)

        DXC=CELL_FV(M)%CELL_DX
        DYC=CELL_FV(M)%CELL_DY
        DZC=CELL_FV(M)%CELL_DZ        
!---------COLLOCATED GRID---------------------------------------
        IF(ICOLL.EQ.1)THEN     ! COLLOCATED GRID
          CALL GETCON_DIV_COL(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                              DXC,DYC,DZC,FCX1,FCY1,FCZ1)
          CALL GETCON_ADV_COL(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                              DXC,DYC,DZC,FCX2,FCY2,FCZ2)
!---------STAGGERED GRID----------------------------------------         
        ELSE                   ! STAGGERED GRID
          IF(ORDER_CON.EQ.4)THEN  ! 4TH ORDER    
            CALL GETCON4_DIV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                                 DXC,DYC,DZC,FCX1,FCY1,FCZ1)
            CALL GETCON4_ADV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                                 DXC,DYC,DZC,FCX2,FCY2,FCZ2)
          ELSE IF(ORDER_CON.EQ.2)THEN  ! 2ND ORDER
            CALL GETCON2_DIV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                                 DXC,DYC,DZC,FCX1,FCY1,FCZ1)
            CALL GETCON2_ADV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,I,J,K, &
                                 DXC,DYC,DZC,FCX2,FCY2,FCZ2)               
          END IF
        END IF

        CELL_FV(INDEX)%CELL_VAR(10)=CELL_FV(INDEX)%CELL_VAR(10)- &
                                    (FCX1*BLEND_CON+FCX2*(1.0-BLEND_CON))
        CELL_FV(INDEX)%CELL_VAR(11)=CELL_FV(INDEX)%CELL_VAR(11)- &
                                    (FCY1*BLEND_CON+FCY2*(1.0-BLEND_CON))    
        CELL_FV(INDEX)%CELL_VAR(12)=CELL_FV(INDEX)%CELL_VAR(12)- &
                                    (FCZ1*BLEND_CON+FCZ2*(1.0-BLEND_CON))
      END IF
    END DO
   
    DEALLOCATE(VAR1,VAR2,VAR3)
    
    END SUBROUTINE CONVECT_CELL

  END MODULE
