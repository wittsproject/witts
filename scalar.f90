! This module is for solving the scalar transport equation.
!
  MODULE scalar 

  USE parameters
  USE field_shared
  USE convect

  CONTAINS

!=============================================================================!
!	             SUBROUTINE OF UPDATING TEMPERATURE 		      !
!=============================================================================!
!       UPDATE TEMPERATURE BY SOLVING A TRANSPORT EQUATION:
!       D (T0)/ DT = - DIV (Q), WHERE D()/ DT IS THE TOTAL DERIVATIVE, AND Q IS THE 
!       HEAT FLUX  (Lu and Fernando Porte-Agel, 2011)
!       A 4th-ORDER Runge-Kutta METHOD IS USED
        SUBROUTINE SCALAR_MAIN(VAR,SI1,SI2,SI3,DX,DY,DZ)

        IMPLICIT NONE 
        INTEGER :: N,I,J,K,N_RK_SCA,SI1,SI2,SI3,NB
        REAL(KIND=DP):: DX,DY,DZ
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: C,B
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: VAR0,VARI,DIV
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: F
        REAL(KIND=DP):: CONTE,SUM,SUMT,TB,SH,BLEND_SCA
        REAL(KIND=DP):: SCHMIDT
        INTEGER :: M,ORDER_SCA

        OPEN(1,FILE="scalar.in")
        READ(1,*)N_RK_SCA
        READ(1,*)BLEND_SCA
        READ(1,*)ORDER_SCA
        READ(1,*)SCHMIDT
        CLOSE(1)

        ALLOCATE(C(N_RK_SCA),B(N_RK_SCA))

        IF(N_RK_SCA.EQ.1)THEN
          B(1)=0.0
          C(1)=0.0
        ELSE IF(N_RK_SCA.EQ.2)THEN
          B(2)=1.0
          B(1)=1.0
          C(1)=0.0
          C(2)=0.5
        ELSE IF(N_RK_SCA.EQ.4)THEN
          B(1)=1.0/6.0
          B(2)=1.0/3.0
          B(3)=1.0/3.0
          B(4)=1.0/6.0
          C(1)=0.0
          C(2)=0.5
          C(3)=0.5
          C(4)=1.0
        END IF

        IF(ORDER_SCA.EQ.2)THEN
          NB=1
        ELSE IF(ORDER_SCA.EQ.4)THEN
          NB=3
        END IF

        ALLOCATE(VAR0(1-NB:NX+NB,1-NB:NY+NB,1-NB:NZ+NB),VARI(1-NB:NX+NB,1-NB:NY+NB,1-NB:NZ+NB), &
                 F(NX,NY,NZ),DIV(NX,NY,NZ))

        DO I=1-NB,NX+NB
          DO J=1-NB,NY+NB
            DO K=1-NB,NZ+NB
              VAR0(I,J,K)=VAR(I,J,K)
              VARI(I,J,K)=VAR(I,J,K)
            END DO
          END DO
        END DO

        DO N=1,N_RK_SCA
          CALL SCALAR_VIS(VARI,1-NB,1-NB,1-NB,SCHMIDT,DIV)
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                F(I,J,K)=-(CON_DIV(VARI,1-NB,1-NB,1-NB,DX,DY,DZ,I,J,K,ORDER_SCA)*BLEND_SCA+         &
                           CON_ADV(VARI,1-NB,1-NB,1-NB,DX,DY,DZ,I,J,K,ORDER_SCA)*(1.0-BLEND_SCA))+  &
                          DIV(I,J,K)
              END DO
            END DO
          END DO

          VARI=VAR0+C(N)*DT*F
          VAR =VAR +B(N)*DT*F
        END DO     
        
        DEALLOCATE(B,C,VAR0,VARI,F,DIV)

        END SUBROUTINE
!=========================================================================!
!                          GET THE SCALAR DIFFUSION                       !
!=========================================================================!
      SUBROUTINE SCALAR_VIS(VAR,SI1,SI2,SI3,SCHMIDT,DIV)
      IMPLICIT NONE
      INTEGER :: I,J,K,SI1,SI2,SI3
      REAL(KIND=DP) :: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: VAR
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: DIV
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: Q1,Q2,Q3
      REAL(KIND=DP) :: SCHMIDT

      ALLOCATE(Q1(0:NX+1,0:NY+1,0:NY+1),Q2(0:NX+1,0:NY+1,0:NY+1), &
               Q3(0:NX+1,0:NY+1,0:NY+1))

      DO I=0,NX+1
        DO J=0,NY+1
          DO K=0,NZ+1
            Q1(I,J,K)=NU(I,J,K)/SCHMIDT*DERIV_X(VAR,SI1,SI2,SI3,I,J,K,1,2,DX)
            Q2(I,J,K)=NU(I,J,K)/SCHMIDT*DERIV_Y(VAR,SI1,SI2,SI3,I,J,K,1,2,DY)
            Q3(I,J,K)=NU(I,J,K)/SCHMIDT*DERIV_Z(VAR,SI1,SI2,SI3,I,J,K,1,2,DZ)
          END DO
        END DO
      END DO
!-----WALL MODEL AT THE FIRST POINT OFF THE WALL
      IF(IWALL(1).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            Q1(0,J,K)=-Q1W_1(J,K)*2.0-Q1(1,J,K)
          END DO
        END DO   
      END IF     

      IF(IWALL(2).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            Q1(NX+1,J,K)=-Q1W_2(J,K)*2.0-Q1(NX,J,K)
          END DO
        END DO   
      END IF    

      IF(IWALL(3).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            Q2(I,0,K)=-Q2W_3(I,K)*2.0-Q2(I,1,K)
          END DO
        END DO  
      END IF     

      IF(IWALL(4).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            Q2(I,NY+1,K)=-Q2W_4(I,K)*2.0-Q2(I,NY,K)
          END DO
        END DO  
      END IF  

      IF(IWALL(5).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            Q3(I,J,0)=-Q3W_5(I,J)*2.0-Q3(I,J,1)
          END DO
        END DO 
      END IF     

      IF(IWALL(6).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            Q3(I,J,NZ+1)=-Q3W_6(I,J)*2.0-Q3(I,J,NZ)
          END DO
        END DO 
      END IF      
!------------------------------------------------------------
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            DIV(I,J,K)=DERIV_X(Q1,0,0,0,I+1,J,K,1,2,DX)+ &
                       DERIV_Y(Q2,0,0,0,I,J+1,K,1,2,DY)+ &
                       DERIV_Z(Q3,0,0,0,I,J,K+1,1,2,DZ)
          END DO
        END DO
      END DO

      DEALLOCATE(Q1,Q2,Q3)
      
      END SUBROUTINE
    
  END MODULE  

      
