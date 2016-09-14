! This module is used for calculation of the unlinear convective term
!
  MODULE convect

  USE parameters
  USE tools
  USE class_shared
  USE class_cell

  IMPLICIT NONE
  
  CONTAINS
!=========================================================================!
!            WRAPPING FOR THE CONVECTIVE TERM BASED ON FV CELL            !
!=========================================================================!
    SUBROUTINE CONVECT_CELL_WRAP()
    IMPLICIT NONE
      
    REAL(KIND=DP),DIMENSION(:,:,:), ALLOCATABLE::VAR1,VAR2,VAR3
    REAL(KIND=DP):: FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2
    REAL(KIND=DP):: DXC,DYC,DZC 
    INTEGER:: M,NB

    NB=ORDER_CON-1
    
    ALLOCATE(VAR1(-NB:NB,-NB:NB,-NB:NB),VAR2(-NB:NB,-NB:NB,-NB:NB), &
             VAR3(-NB:NB,-NB:NB,-NB:NB))
    
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
!---------COLLOCATED GRID--------------------------------------
        IF(ICOLL.EQ.1)THEN
           CELL_FV(M)%CELL_FX=CELL_FV(M)%CELL_FX+ &
                              CONVECT(CELL_FV(:)%CELL_VEL(1),M,BLEND_CON,ORDER_CON)
           CELL_FV(M)%CELL_FY=CELL_FV(M)%CELL_FY+ &
                              CONVECT(CELL_FV(:)%CELL_VEL(2),M,BLEND_CON,ORDER_CON)
           CELL_FV(M)%CELL_FZ=CELL_FV(M)%CELL_FZ+ &
                              CONVECT(CELL_FV(:)%CELL_VEL(3),M,BLEND_CON,ORDER_CON)       
!---------STAGGERED GRID----------------------------------------
        ELSE           
          CALL CELL_TO_STRUCT(CELL_FV(:)%CELL_VEL(1),M,NB,VAR1)
          CALL CELL_TO_STRUCT(CELL_FV(:)%CELL_VEL(2),M,NB,VAR2)
          CALL CELL_TO_STRUCT(CELL_FV(:)%CELL_VEL(3),M,NB,VAR3)

          DXC=CELL_FV(M)%CELL_DX
          DYC=CELL_FV(M)%CELL_DY
          DZC=CELL_FV(M)%CELL_DZ        
        
          IF(ORDER_CON.EQ.4)THEN  ! 4TH ORDER    
            CALL GETCON4_DIV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,0,0,0, &
                                 DXC,DYC,DZC,FCX1,FCY1,FCZ1)
            CALL GETCON4_ADV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,0,0,0, &
                                 DXC,DYC,DZC,FCX2,FCY2,FCZ2)
          ELSE IF(ORDER_CON.EQ.2)THEN  ! 2ND ORDER
            CALL GETCON2_DIV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,0,0,0, &
                                 DXC,DYC,DZC,FCX1,FCY1,FCZ1)
            CALL GETCON2_ADV_STA(VAR1,VAR2,VAR3,-NB,-NB,-NB,0,0,0, &
                                 DXC,DYC,DZC,FCX2,FCY2,FCZ2)               
          END IF
        END IF
        
        CELL_FV(M)%CELL_FX=CELL_FV(M)%CELL_FX- &
                                (FCX1*BLEND_CON+FCX2*(1.0-BLEND_CON))
        CELL_FV(M)%CELL_FY=CELL_FV(M)%CELL_FY- &
                                (FCY1*BLEND_CON+FCY2*(1.0-BLEND_CON))    
        CELL_FV(M)%CELL_FZ=CELL_FV(M)%CELL_FZ- &
                                (FCZ1*BLEND_CON+FCZ2*(1.0-BLEND_CON))
      END IF
    END DO
   
    DEALLOCATE(VAR1,VAR2,VAR3)
    
    END SUBROUTINE CONVECT_CELL_WRAP

  
!=========================================================================!
!             CALCULATE THE CONVECTION TERM AT THE CELL CENTER            !
!=========================================================================!
!   adv: u*dQ/dx+v*dQ/dy+w*dQ/dz
!   div: d(u*Q)/dx+d(v*Q)/dy+d(w*Q)/dz
!   return: -[div*BLEND+adv*(1-BLEND)]
!   Q is VAR_CELL here;
!   u,v,w are velocity fluxes on cell faces    
    REAL(KIND=DP) FUNCTION CONVECT(VAR_CELL,INDEX,BLEND,ORDER)
    IMPLICIT NONE

    INTEGER:: INDEX,M,ORDER,NB
    REAL(KIND=DP):: BLEND
    REAL(KIND=DP),DIMENSION(:):: VAR_CELL
    REAL(KIND=DP),DIMENSION(:,:,:), ALLOCATABLE::VAR  
    REAL(KIND=DP):: FCX1,FCY1,FCZ1,FCX2,FCY2,FCZ2
    
    NB=ORDER-1
    
    ALLOCATE(VAR(-NB:NB,-NB:NB,-NB:NB))
         
    CALL CELL_TO_STRUCT(VAR_CELL,INDEX,NB,VAR)

    DX=CELL_FV(INDEX)%CELL_DX
    DY=CELL_FV(INDEX)%CELL_DY
    DZ=CELL_FV(INDEX)%CELL_DZ         

    DIV=CON_DIV(INDEX,VAR,-NB,-NB,-NB,DX,DY,DZ,0,0,0,ORDER)
    ADV=CON_ADV(INDEX,VAR,-NB,-NB,-NB,DX,DY,DZ,0,0,0,ORDER)

    CONVECT=-(DIV*BLEND+ADV*(1.0-BLEND))
     
    END FUNCTION CONVECT  
!=========================================================================!
!     GET CONVECTIVE TERM, 4TH ORDER, ADVECTIVE FORM, STAGGERED GRID      !
!=========================================================================!     
      SUBROUTINE GETCON4_ADV_STA(VAR1,VAR2,VAR3,NB1,NB2,NB3,I,J,K, &
                                 DX,DY,DZ,FCX,FCY,FCZ)
      IMPLICIT NONE
      INTEGER:: NB1,NB2,NB3,I,J,K,INDEX        
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

      END SUBROUTINE GETCON4_ADV_STA    
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

      END SUBROUTINE GETCON4_DIV_STA    
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

      END SUBROUTINE GETCON2_ADV_STA    
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

      END SUBROUTINE GETCON2_DIV_STA    
!=========================================================================!
!                            DIVERGENCE FORM                              !
!=========================================================================!
      REAL(KIND=DP) FUNCTION CON_DIV(INDEX,VAR,SI1,SI2,SI3,DX,DY,DZ,I,J,K,ORDER)

      IMPLICIT NONE

      INTEGER :: INDEX,SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER :: I,J,K,ORDER
      INTEGER :: IN_X_P1,IN_X_M1,IN_Y_P1,IN_Y_M1,IN_Z_P1,IN_Z_M1

      IF(ORDER.EQ.2)THEN
        CON_DIV=  &
         (CELL_FV(INDEX)%CELL_UF(1,2)*STENCIL_IX(VAR,SI1,SI2,SI3,I+1,J,K,1)-     &
          CELL_FV(INDEX)%CELL_UF(1,1)*STENCIL_IX(VAR,SI1,SI2,SI3,I,  J,K,1))/DX+ &
         (CELL_FV(INDEX)%CELL_UF(2,2)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+1,K,1)-     &
          CELL_FV(INDEX)%CELL_UF(2,1)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J,  K,1))/DY+ &
         (CELL_FV(INDEX)%CELL_UF(3,2)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+1,1)-     &
          CELL_FV(INDEX)%CELL_UF(3,1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K  ,1))/DZ 
      ELSE IF(ORDER.EQ.4)THEN
        IN_X_P1=LOOKUP_NEI(INDEX,1, 1)  
        IN_X_M1=LOOKUP_NEI(INDEX,1,-1)
        IN_Y_P1=LOOKUP_NEI(INDEX,2, 1)  
        IN_Y_M1=LOOKUP_NEI(INDEX,2,-1)
        IN_Z_P1=LOOKUP_NEI(INDEX,3, 1)  
        IN_Z_M1=LOOKUP_NEI(INDEX,3,-1)
         
        CON_DIV=  &
    9.0/8.0*((CELL_FV(INDEX)%CELL_UF(1,2)*STENCIL_IX(VAR,SI1,SI2,SI3,I+1,J,K,1)-            &
              CELL_FV(INDEX)%CELL_UF(1,1)*STENCIL_IX(VAR,SI1,SI2,SI3,I,  J,K,1))/DX+        &
             (CELL_FV(INDEX)%CELL_UF(2,2)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+1,K,1)-            &
              CELL_FV(INDEX)%CELL_UF(2,1)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J,  K,1))/DY+        &
             (CELL_FV(INDEX)%CELL_UF(3,2)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+1,1)-            &
              CELL_FV(INDEX)%CELL_UF(3,1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,  1))/DZ)-       &
    1.0/8.0*((CELL_FV(IN_X_P1)%CELL_UF(1,2)*STENCIL_IX(VAR,SI1,SI2,SI3,I+2,J,K,3)-            &
              CELL_FV(IN_X_M1)%CELL_UF(1,1)*STENCIL_IX(VAR,SI1,SI2,SI3,I-1,J,K,3))/(DX*3.0)+  &
             (CELL_FV(IN_Y_P1)%CELL_UF(2,2)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J+2,K,3)-            &
              CELL_FV(IN_Y_M1)%CELL_UF(2,1)*STENCIL_IY(VAR,SI1,SI2,SI3,I,J-1,K,3))/(DY*3.0)+  &
             (CELL_FV(IN_Z_P1)%CELL_UF(3,2)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K+2,3)-            &
              CELL_FV(IN_Z_M1)%CELL_UF(3,1)*STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K-1,3))/(DZ*3.0))
      END IF

      END FUNCTION CON_DIV    
!=========================================================================!
!                          DIVERGENCE FORM                                !
!=========================================================================!
      REAL(KIND=DP) FUNCTION CON_ADV(INDEX,VAR,SI1,SI2,SI3,DX,DY,DZ,I,J,K,ORDER)

      IMPLICIT NONE

      INTEGER :: INDEX,SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER :: I,J,K,ORDER
      INTEGER :: IN_X_P1,IN_X_M1,IN_Y_P1,IN_Y_M1,IN_Z_P1,IN_Z_M1
      
      IF(ORDER.EQ.2)THEN
        CON_ADV=  &
        ((CELL_FV(INDEX)%CELL_UF(1,1)*STENCIL_DX(VAR,SI1,SI2,SI3,I,  J,K,1,DX)+      &
          CELL_FV(INDEX)%CELL_UF(1,2)*STENCIL_DX(VAR,SI1,SI2,SI3,I+1,J,K,1,DX))/2.0+ &
         (CELL_FV(INDEX)%CELL_UF(2,1)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J,  K,1,DY)+      &
          CELL_FV(INDEX)%CELL_UF(2,2)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+1,K,1,DY))/2.0+ &
         (CELL_FV(INDEX)%CELL_UF(3,1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,  1,DZ)+      &
          CELL_FV(INDEX)%CELL_UF(3,2)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+1,1,DZ))/2.0)
      ELSE IF(ORDER.EQ.4)THEN
        IN_X_P1=LOOKUP_NEI(INDEX,1, 1)  
        IN_X_M1=LOOKUP_NEI(INDEX,1,-1)
        IN_Y_P1=LOOKUP_NEI(INDEX,2, 1)  
        IN_Y_M1=LOOKUP_NEI(INDEX,2,-1)
        IN_Z_P1=LOOKUP_NEI(INDEX,3, 1)  
        IN_Z_M1=LOOKUP_NEI(INDEX,3,-1)
         
        CON_ADV=  &
    9.0/8.0*((CELL_FV(INDEX)%CELL_UF(1,1)*STENCIL_DX(VAR,SI1,SI2,SI3,I,  J,K,1,DX)+       &
              CELL_FV(INDEX)%CELL_UF(1,2)*STENCIL_DX(VAR,SI1,SI2,SI3,I+1,J,K,1,DX))/2.0+  &
             (CELL_FV(INDEX)%CELL_UF(2,1)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J,  K,1,DY)+       &
              CELL_FV(INDEX)%CELL_UF(2,2)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+1,K,1,DY))/2.0+  &
             (CELL_FV(INDEX)%CELL_UF(3,1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,  1,DZ)+       &
              CELL_FV(INDEX)%CELL_UF(3,2)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+1,1,DZ))/2.0)- &
    1.0/8.0*((CELL_FV(IN_X_M1)%CELL_UF(1,1)*STENCIL_DX(VAR,SI1,SI2,SI3,I-1,J,K,3,DX)+       &
              CELL_FV(IN_X_P1)%CELL_UF(1,2)*STENCIL_DX(VAR,SI1,SI2,SI3,I+2,J,K,3,DX))/2.0+  &
             (CELL_FV(IN_X_M1)%CELL_UF(2,1)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J-1,K,3,DY)+       &
              CELL_FV(IN_X_P1)%CELL_UF(2,2)*STENCIL_DY(VAR,SI1,SI2,SI3,I,J+2,K,3,DY))/2.0+  &
             (CELL_FV(IN_X_M1)%CELL_UF(3,1)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K-1,3,DZ)+       &
              CELL_FV(IN_X_P1)%CELL_UF(3,2)*STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K+2,3,DZ))/2.0)
      END IF

      END FUNCTION CON_ADV
    
  END MODULE
