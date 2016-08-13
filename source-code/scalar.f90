! This module is for solving the scalar transport equation.
!
  MODULE scalar 

  USE mpi
  USE parameters
  USE field_shared
  USE boundary
  USE convect
  USE tools

  CONTAINS

!=============================================================================!
!	             SUBROUTINE OF UPDATING TEMPERATURE 		      !
!=============================================================================!
!       UPDATE TEMPERATURE BY SOLVING A TRANSPORT EQUATION:
!       D (T0)/ DT = - DIV (Q), WHERE D()/ DT IS THE TOTAL DERIVATIVE, AND Q IS THE 
!       HEAT FLUX  (Lu and Fernando Porte-Agel, 2011)
!       A 4th-ORDER Runge-Kutta METHOD IS USED
        SUBROUTINE SCALAR_WRAP(VAR,SI1,SI2,SI3,DX,DY,DZ)

        IMPLICIT NONE 
        INTEGER :: N,I,J,K,N_RK_SCA,MUSCL_SCA,SI1,SI2,SI3,NB
        REAL(KIND=DP):: DX,DY,DZ
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: C,B
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: VAR0,VARI,DIV
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: F
        REAL(KIND=DP):: CONTE,SUM,SUMT,TB,SH,BLEND_SCA
        REAL(KIND=DP):: SCHMIDT,F_CON_DIV,F_CON_ADV
        INTEGER :: M,ORDER_SCA

        OPEN(1,FILE="scalar.in")
        READ(1,*)N_RK_SCA
        READ(1,*)MUSCL_SCA
        READ(1,*)BLEND_SCA
        READ(1,*)ORDER_SCA
        READ(1,*)SCHMIDT
        CLOSE(1)
             
        ALLOCATE(C(N_RK_SCA),B(N_RK_SCA))

        IF(N_RK_SCA.EQ.1)THEN
          B(1)=1.0
          C(1)=0.0
        ELSE IF(N_RK_SCA.EQ.2)THEN
          B(1)=0.0
          B(2)=1.0
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

        ALLOCATE(VAR0(NX1:NX2,NY1:NY2,NZ1:NZ2),VARI(NX1:NX2,NY1:NY2,NZ1:NZ2), &
                 F(NX1:NX2,NY1:NY2,NZ1:NZ2),DIV(NX,NY,NZ))

        VAR0=VAR
        VARI=VAR

        DO N=1,N_RK_SCA
          CALL SCALAR_VIS(VARI,NX1,NY1,NZ1,DX,DY,DZ,SCHMIDT,DIV)
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                IF(MUSCL_SCA.EQ.1)THEN
                  F_CON_DIV=CON_DIV_MUSCL(VARI,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_SCA)
                ELSE
                  F_CON_DIV=CON_DIV(VARI,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_SCA)   
                END IF

                F_CON_ADV=CON_ADV(VARI,NX1,NY1,NZ1,DX,DY,DZ,I,J,K,ORDER_SCA)     
          
                F(I,J,K)=-(F_CON_DIV*BLEND_SCA+F_CON_ADV*(1.0-BLEND_SCA))+  &
                           DIV(I,J,K)
              END DO
            END DO
          END DO

          IF(N.LT.N_RK_SCA)THEN
             VARI=VAR0+C(N+1)*DT*F
            CALL GET_BC(NX,NY,NZ,VARI,NB,NB,NB,0,  &
                        BC(1,4),BV(1,4)) 
          END IF

          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ                  
                VAR(I,J,K)=VAR(I,J,K)+B(N)*DT*F(I,J,K)
              END DO
            END DO
          END DO             
       END DO     

       CALL GET_BC(NX,NY,NZ,VAR,NB,NB,NB,0,  &
                   BC(1,4),BV(1,4))       
        
       DEALLOCATE(B,C,VAR0,VARI,F,DIV)

       END SUBROUTINE
!=========================================================================!
!                          GET THE SCALAR DIFFUSION                       !
!=========================================================================!
      SUBROUTINE SCALAR_VIS(VAR,SI1,SI2,SI3,DX,DY,DZ,SCHMIDT,DIV)
      IMPLICIT NONE
      INTEGER :: I,J,K,SI1,SI2,SI3
      REAL(KIND=DP) :: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: VAR
      REAL(KIND=DP),DIMENSION(:,:,:) :: DIV
      REAL(KIND=DP) :: SCHMIDT

      DO I=0,NX+1
        DO J=0,NY+1
           DO K=0,NZ+1
            Q1(I,J,K)=-NU(I,J,K)/SCHMIDT*DERIV_X(VAR,NX1,NY1,NZ1,I,J,K,1,2,DX)
            Q2(I,J,K)=-NU(I,J,K)/SCHMIDT*DERIV_Y(VAR,NX1,NY1,NZ1,I,J,K,1,2,DY)
            Q3(I,J,K)=-NU(I,J,K)/SCHMIDT*DERIV_Z(VAR,NX1,NY1,NZ1,I,J,K,1,2,DZ)
          END DO
        END DO
     END DO
!-----WALL MODEL AT THE FIRST POINT OFF THE WALL
     IF(IWALL(1).NE.0.AND.MYIDX.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            Q1(0,J,K)=Q1W_1(J,K)*2.0-Q1(1,J,K)
          END DO
        END DO   
      END IF     

      IF(IWALL(2).NE.0.AND.MYIDX.EQ.NPX-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            Q1(NX+1,J,K)=Q1W_2(J,K)*2.0-Q1(NX,J,K)
          END DO
        END DO   
      END IF    

      IF(IWALL(3).NE.0.AND.MYIDY.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            Q2(I,0,K)=Q2W_3(I,K)*2.0-Q2(I,1,K)
          END DO
        END DO  
      END IF     

      IF(IWALL(4).NE.0.AND.MYIDY.EQ.NPY-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            Q2(I,NY+1,K)=Q2W_4(I,K)*2.0-Q2(I,NY,K)
          END DO
        END DO  
      END IF  

      IF(IWALL(5).NE.0.AND.MYIDZ.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            Q3(I,J,0)=Q3W_5(I,J)*2.0-Q3(I,J,1)
          END DO
        END DO 
      END IF     

      IF(IWALL(6).NE.0.AND.MYIDZ.EQ.NPZ-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            Q3(I,J,NZ+1)=Q3W_6(I,J)*2.0-Q3(I,J,NZ)
          END DO
        END DO 
      END IF      
 !------------------------------------------------------------
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            DIV(I,J,K)=-(DERIV_X(Q1,NX1,NY1,NZ1,I+1,J,K,1,2,DX)+ &
                         DERIV_Y(Q2,NX1,NY1,NZ1,I,J+1,K,1,2,DY)+ &
                         DERIV_Z(Q3,NX1,NY1,NZ1,I,J,K+1,1,2,DZ))
          END DO
        END DO
      END DO

      END SUBROUTINE
    
  END MODULE  

      
