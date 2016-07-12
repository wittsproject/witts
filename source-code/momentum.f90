! This module is for solving the momentum N-S equations
!
  MODULE momentum

  USE parameters
  USE field_shared
  USE tools
  USE boundary
  USE flux
  USE convect
  USE viscous
  USE wall_model
  USE source
  USE turbine
  USE coriolis
  USE gravity
  USE poisson

  CONTAINS
!=============================================================================!
!	         UPDATE MOMENTUM EQUATION BY SEMI A-B SCHEME                  !
!=============================================================================!
!       THE SEMI A-B SCHEME: THE PROJECTION STEP IS OUTSIDE THE LOOP
        SUBROUTINE MOMENTUM_AB(DX,DY,DZ)

        IMPLICIT NONE 
        INTEGER :: M,MM,I,J,K
        REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE:: A
        REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: FXI,FYI,FZI
        REAL(KIND=DP) :: DX,DY,DZ
 
        ALLOCATE(A(ORDER_TIM,ORDER_TIM))
!-------PARAMETERS FOR THE A-B SCHEME
        IF(ORDER_TIM.EQ.1)THEN
          A(1,1)=1.0
        ELSE IF(ORDER_TIM.EQ.2)THEN
          A(1,1)=1.0
          A(1,2)=-0.5
          A(2,2)=1.5
        ELSE IF(ORDER_TIM.EQ.3)THEN
          A(1,1)=1.0
          A(1,2)=-0.5
          A(2,2)=1.5          
          A(1,3)=5.0/12.0
          A(2,3)=-4.0/3.0
          A(3,3)=23.0/12.0
        ELSE IF(ORDER_TIM.EQ.4)THEN
          A(1,1)=1.0
          A(1,2)=-0.5
          A(2,2)=1.5          
          A(1,3)=5.0/12.0
          A(2,3)=-4.0/3.0
          A(3,3)=23.0/12.0
          A(1,4)=-3.0/8.0
          A(2,4)=37.0/24.0
          A(3,4)=-59.0/24.0
          A(4,4)=55.0/24.0
        END IF

        ALLOCATE(FXI(NX1:NX2,NX1:NY2,NX1:NZ2,ORDER_TIM),FYI(NX1:NX2,NX1:NY2,NX1:NZ2,ORDER_TIM), &
                 FZI(NX1:NX2,NX1:NY2,NX1:NZ2,ORDER_TIM))
        DO M=1,ORDER_TIM 
          FX=0.0
          FY=0.0
          FZ=0.0
!---------GET WIND TURBINE BODY FORCE TERM
          IF(ITURBINE.EQ.1)THEN
            CALL TURBINE_WRAP(DX,DY,DZ,DT/ORDER_TIM)
          END IF
!---------GET WALL SHEAR STRESS AND HEAT FLUX BY USING WALL MODEL
          CALL WALL_MODEL_WRAP(DX,DY,DZ)
!---------GET VISCOUS FORCES  
          CALL VISCOUS_WRAP(DX,DY,DZ)
!---------GET CONVECTIVE FORCES
          CALL CONVECT_WRAP(DX,DY,DZ)
!---------GET GRAVITY FORCE USING BOUSSINESQ APPROXIMATION
          IF(IBOUS.EQ.1)THEN
            CALL BOUSSINESQ()
          END IF 
!---------GET CORIOLIS FORCES
          IF(ICORI.EQ.1)THEN
            CALL CORIO()
          END IF
!---------SET SOURCE TERM
          CALL SOURCE_MOM()
!---------SAVE RHS
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                FXI(I,J,K,M)=FX(I,J,K)
                FYI(I,J,K,M)=FY(I,J,K)
                FZI(I,J,K,M)=FZ(I,J,K)
!---------CORRECT THE VELOCITY (AND THE FLUXES)
                DO MM=1,M
                  U(I,J,K)=U(I,J,K)+DT/ORDER_TIM*A(MM,M)*FXI(I,J,K,MM)
                  V(I,J,K)=V(I,J,K)+DT/ORDER_TIM*A(MM,M)*FYI(I,J,K,MM)
                  W(I,J,K)=W(I,J,K)+DT/ORDER_TIM*A(MM,M)*FZI(I,J,K,MM)
                END DO
              END DO
            END DO
          END DO
!---------UPDATE VELOCITY BOUNDARY CONDITIONS
          CALL BOUNDARY_VEL(NX,NY,NZ)
        END DO
!      UPDATE INTERMEDIATE VELOCITY FLUX
        CALL FLUX_VELOCITY(DX,DY,DZ,ORDER_CON)
!
!-------PROJECTION STEP----------------------------------------------
!       SOLVE THE POISSON EQUATION TO GET THE PRESSURE GRADIENT TERM  
!       NOTE, THE VELOCITY FLUX IS USED TO CALCULATE THE DIVERGENCE  
        CALL PRESSURE_WRAP(DX,DY,DZ,DT) 
!       CORRECT THE VELOCITY FLUXES AT STAGGERED LOCATIONS. NOTE, U=UF,V=VF,W=WF WHEN ICOLL=0
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ   
              UF(I,J,K)=UF(I,J,K)-DERIV_X(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DX)*DT
              VF(I,J,K)=VF(I,J,K)-DERIV_Y(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DY)*DT
              WF(I,J,K)=WF(I,J,K)-DERIV_Z(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DZ)*DT   
!              UF(I,J,K)=U(I,J,K)-((PD(I,  J,K)-PD(I-1,J,K))/DX*9.0/8.0- &
!                                  (PD(I+1,J,K)-PD(I-2,J,K))/(DX*3.0)/8.0)*DT
!              VF(I,J,K)=V(I,J,K)-((PD(I,J,  K)-PD(I,J-1,K))/DY*9.0/8.0- &
!                                  (PD(I,J+1,K)-PD(I,J-2,K))/(DY*3.0)/8.0)*DT
!              WF(I,J,K)=W(I,J,K)-((PD(I,J,K)  -PD(I,J,K-1))/DZ*9.0/8.0- &
!                                  (PD(I,J,K+1)-PD(I,J,K-2))/(DZ*3.0)/8.0)*DT      
            END DO
          END DO
        END DO

        IF(ICOLL.EQ.0)THEN
          U=UF
          V=VF
          W=WF
        ELSE   !  IF ICOLL=1, THEN USE INTERPOLATION TO GET VELOCITY AT THE CELL CENTER
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                U(I,J,K)=VAR_INTER_X(UF,NX1,NY1,NZ1,I+1,J,K,1,ORDER_CON)
                V(I,J,K)=VAR_INTER_Y(VF,NX1,NY1,NZ1,I,J+1,K,1,ORDER_CON)
                W(I,J,K)=VAR_INTER_Z(WF,NX1,NY1,NZ1,I,J,K+1,1,ORDER_CON)
              END DO
            END DO
          END DO
        END IF

!       UPDATE VELOCITY BOUNDARY CONDITIONS
        CALL BOUNDARY_VEL(NX,NY,NZ)

        DEALLOCATE(A,FXI,FYI,FZI)

        END SUBROUTINE
!=============================================================================!
!	            UPDATE MOMENTUM EQUATION BY A-B SCHEME                    !
!=============================================================================!
!       UPDATE MOMENTUM BY USING PROJECTION METHOD:
!       A 4th-ORDER Runge-Kutta METHOD IS USED
        SUBROUTINE MOMENTUM_RK(DX,DY,DZ)

        IMPLICIT NONE 
        INTEGER :: M,MM,I,J,K
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: C,B
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: FX0,FY0,FZ0,UN,VN,WN,UI,VI,WI
        REAL(KIND=DP) :: DX,DY,DZ
    
        ALLOCATE(C(ORDER_TIM),B(ORDER_TIM))
!-------PARAMETERS FOR THE R-K SCHEME
        IF(ORDER_TIM.EQ.1)THEN
          B(1)=1.0
          C(1)=0.0
        ELSE IF(ORDER_TIM.EQ.2)THEN
          B(1)=0.0
          B(2)=1.0
          C(1)=0.0
          C(2)=0.5
        ELSE IF(ORDER_TIM.EQ.4)THEN
          B(1)=1.0/6.0
          B(2)=1.0/3.0
          B(3)=1.0/3.0
          B(4)=1.0/6.0
          C(1)=0.0
          C(2)=0.5
          C(3)=0.5
          C(4)=1.0
        END IF

        ALLOCATE(UN(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 VN(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 WN(NX1:NX2,NY1:NY2,NZ1:NZ2))
        ALLOCATE(UI(NX1:NX2,NX1:NY2,NX1:NZ2),VI(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 WI(NX1:NX2,NX1:NY2,NX1:NZ2)) 
        FX=0.0
        FY=0.0
        FZ=0.0       
        UN=U
        VN=V
        WN=W
        UI=U
        VI=V
        WI=W

        DO MM=1,ORDER_TIM  
          FX=0.0
          FY=0.0
          FZ=0.0
!-------GET WIND TURBINE BODY FORCE TERM    
          IF(ITURBINE.EQ.1)THEN
            CALL TURBINE_WRAP(DX,DY,DZ,DT*B(MM))
          END IF
!---------GET WALL SHEAR STRESS AND HEAT FLUX BY USING WALL MODEL
          CALL WALL_MODEL_WRAP(DX,DY,DZ)
!---------GET VISCOUS FORCES  
          CALL VISCOUS_WRAP(DX,DY,DZ)
!---------GET CONVECTIVE FORCES
          CALL CONVECT_WRAP(DX,DY,DZ)
!---------GET GRAVITY FORCE USING BOUSSINESQ APPROXIMATION
          IF(IBOUS.EQ.1)THEN
            CALL BOUSSINESQ()
          END IF 
!---------GET CORIOLIS FORCES
          IF(ICORI.EQ.1)THEN
            CALL CORIO()
          END IF
!---------SET SOURCE TERM
          CALL SOURCE_MOM()
!---------CALCULATE THE INTERMEDIATE VELOCITY 
          U=UN+DT*FX
          V=VN+DT*FY
          W=WN+DT*FZ
!---------EXPAND VELOCITY FIELD
          CALL BOUNDARY_VEL(NX,NY,NZ)
          CALL FLUX_VELOCITY(DX,DY,DZ,ORDER_CON)
!---------SOLVE THE POISSON EQUATION TO GET THE PRESSURE GRADIENT TERM    
          CALL PRESSURE_WRAP(DX,DY,DZ,DT) 
!---------CORRECT THE VELOCITY (AND THE FLUXES)
          IF(MM.LT.ORDER_TIM)THEN
            U=UN+C(MM+1)*DT*FX
            V=VN+C(MM+1)*DT*FY
            W=WN+C(MM+1)*DT*FZ
!---------UPDATE VELOCITY BOUNDARY CONDITIONS
            CALL BOUNDARY_VEL(NX,NY,NZ)
!---------FOR COLLOCATED GRID, USE RHIE-CHOW SCHEME TO UPDATE VELOCITY FLUXES
            IF(ICOLL.EQ.1)THEN  
              DO I=NX1+1,NX2-1
                DO J=NY1+1,NY2-1
                  DO K=NZ1+1,NZ2-1   
                    UF(I,J,K)=UF(I,J,K)-DERIV_X(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DX)*DT*C(MM+1)
                    VF(I,J,K)=VF(I,J,K)-DERIV_Y(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DY)*DT*C(MM+1)                                  
                    WF(I,J,K)=WF(I,J,K)-DERIV_Z(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DZ)*DT*C(MM+1)                               
                  END DO
                END DO
              END DO  
            ELSE
              UF=U
              VF=V
              WF=W
            END IF
          END IF

          UI=UI+B(MM)*DT*FX
          VI=VI+B(MM)*DT*FY
          WI=WI+B(MM)*DT*FZ
        END DO

        U=UI
        V=VI
        W=WI
!-------EXPAND VELOCITY FIELD
        CALL BOUNDARY_VEL(NX,NY,NZ)
        CALL FLUX_VELOCITY(DX,DY,DZ,ORDER_CON)
       
        DEALLOCATE(B,C,UN,VN,WN,UI,VI,WI)

        END SUBROUTINE  
!=============================================================================!
!	         UPDATE MOMENTUM EQUATION BY SEMI R-K SCHEME                  !
!=============================================================================!
!       THE SEMI R-K SCHEME: THE PROJECTION STEP IS OUTSIDE THE LOOP
        SUBROUTINE MOMENTUM_SEMI_RK(DX,DY,DZ)

        IMPLICIT NONE 
        INTEGER :: M,MM,I,J,K
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: C,B
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: FX0,FY0,FZ0,UN,VN,WN,UI,VI,WI
        REAL(KIND=DP) :: DX,DY,DZ
    
        ALLOCATE(C(ORDER_TIM),B(ORDER_TIM))
!-------PARAMETERS FOR THE R-K SCHEME
        IF(ORDER_TIM.EQ.1)THEN
          B(1)=1.0
          C(1)=0.0
        ELSE IF(ORDER_TIM.EQ.2)THEN
          B(1)=0.0
          B(2)=1.0
          C(1)=0.0
          C(2)=0.5
        ELSE IF(ORDER_TIM.EQ.4)THEN
          B(1)=1.0/6.0
          B(2)=1.0/3.0
          B(3)=1.0/3.0
          B(4)=1.0/6.0
          C(1)=0.0
          C(2)=0.5
          C(3)=0.5
          C(4)=1.0
        END IF

        ALLOCATE(UN(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 VN(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 WN(NX1:NX2,NY1:NY2,NZ1:NZ2))
        ALLOCATE(UI(NX1:NX2,NX1:NY2,NX1:NZ2),VI(NX1:NX2,NX1:NY2,NX1:NZ2), &
                 WI(NX1:NX2,NX1:NY2,NX1:NZ2)) 
           
        UN=U
        VN=V
        WN=W
        UI=U
        VI=V
        WI=W

        DO MM=1,ORDER_TIM
          FX=0.0
          FY=0.0
          FZ=0.0  
!-------GET WIND TURBINE BODY FORCE TERM                                                                                                                                                         
          IF(ITURBINE.EQ.1)THEN
            CALL TURBINE_WRAP(DX,DY,DZ,DT*B(MM))
          END IF
!---------GET WALL SHEAR STRESS AND HEAT FLUX BY USING WALL MODEL
          CALL WALL_MODEL_WRAP(DX,DY,DZ)
!---------GET VISCOUS FORCES  
          CALL VISCOUS_WRAP(DX,DY,DZ)
!---------GET CONVECTIVE FORCES
          CALL CONVECT_WRAP(DX,DY,DZ)
!---------GET GRAVITY FORCE USING BOUSSINESQ APPROXIMATION
          IF(IBOUS.EQ.1)THEN
            CALL BOUSSINESQ()
          END IF 
!---------GET CORIOLIS FORCES
          IF(ICORI.EQ.1)THEN
            CALL CORIO()
          END IF
!---------SET SOURCE TERM
          CALL SOURCE_MOM()
!---------CALCULATE THE INTERMEDIATE VELOCITY 
          IF(MM.LT.ORDER_TIM)THEN
            U=UN+C(MM+1)*DT*FX
            V=VN+C(MM+1)*DT*FY
            W=WN+C(MM+1)*DT*FZ
!---------EXPAND VELOCITY FIELD
            CALL BOUNDARY_VEL(NX,NY,NZ)
            CALL FLUX_VELOCITY(DX,DY,DZ,ORDER_CON)
          END IF

          UI=UI+B(MM)*DT*FX
          VI=VI+B(MM)*DT*FY
          WI=WI+B(MM)*DT*FZ
        END DO

        U=UI
        V=VI
        W=WI
!-------EXPAND VELOCITY FIELD
        CALL BOUNDARY_VEL(NX,NY,NZ)
        CALL FLUX_VELOCITY(DX,DY,DZ,ORDER_CON)

!-------PROJECTION STEP----------------------------------------------
!       SOLVE THE POISSON EQUATION TO GET THE PRESSURE GRADIENT TERM  
!       NOTE, THE VELOCITY FLUX IS USED TO CALCULATE THE DIVERGENCE  
        CALL PRESSURE_WRAP(DX,DY,DZ,DT) 
!       CORRECT THE VELOCITY FLUXES AT STAGGERED LOCATIONS. NOTE, U=UF,V=VF,W=WF WHEN ICOLL=0
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ   
              UF(I,J,K)=UF(I,J,K)-DERIV_X(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DX)*DT
              VF(I,J,K)=VF(I,J,K)-DERIV_Y(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DY)*DT                                  
              WF(I,J,K)=WF(I,J,K)-DERIV_Z(PD,NX1,NY1,NZ1,I,J,K,1,ORDER_POI,DZ)*DT                               
            END DO
          END DO
        END DO
        IF(ICOLL.EQ.0)THEN
          U=UF
          V=VF
          W=VF
        ELSE   !  IF ICOLL=1, THEN USE INTERPOLATION TO GET VELOCITY AT THE CELL CENTER
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                U(I,J,K)=VAR_INTER_X(UF,NX1,NY1,NZ1,I+1,J,K,1,ORDER_CON)
                V(I,J,K)=VAR_INTER_Y(VF,NX1,NY1,NZ1,I,J+1,K,1,ORDER_CON)
                W(I,J,K)=VAR_INTER_Z(WF,NX1,NY1,NZ1,I,J,K+1,1,ORDER_CON)
              END DO
            END DO
          END DO
        END IF
!       UPDATE VELOCITY BOUNDARY CONDITIONS
        CALL BOUNDARY_VEL(NX,NY,NZ)
        
        DEALLOCATE(B,C,FX0,FY0,FZ0,UN,VN,WN,UI,VI,WI)

        END SUBROUTINE
  END MODULE
