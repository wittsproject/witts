! This module is used for calculating viscous forces.
!
  MODULE viscous

  USE parameters
  USE field_shared   
  USE tools
  USE sgs

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!                 CALCULATE THE VISCOUS FORCE (MOLECULAR+SGS)             !
!=========================================================================!
     SUBROUTINE VISCOUS_WRAP(DX,DY,DZ)
     IMPLICIT NONE
     INTEGER :: I,J,K
     REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: SFX,SFY,SFZ
     REAL(KIND=DP):: DX,DY,DZ

     ALLOCATE(S12(NX1:NX2,NY1:NY2,NZ1:NZ2),S13(NX1:NX2,NY1:NY2,NZ1:NZ2), &
              S23(NX1:NX2,NY1:NY2,NZ1:NZ2))
     ALLOCATE(S11(NX1:NX2,NY1:NY2,NZ1:NZ2),S22(NX1:NX2,NY1:NY2,NZ1:NZ2), &
              S33(NX1:NX2,NY1:NY2,NZ1:NZ2))
     ALLOCATE(S12C(NX1:NX2,NY1:NY2,NZ1:NZ2),S13C(NX1:NX2,NY1:NY2,NZ1:NZ2), &
              S23C(NX1:NX2,NY1:NY2,NZ1:NZ2),S(NX1:NX2,NY1:NY2,NZ1:NZ2))  
 
     NU=0.0
!----CALCUALTE THE EDDY VISCOSITY
     IF(ITYPE.EQ.1)THEN        ! IF LES IS USED
       IF(ISCHEME.EQ.1)THEN     ! FINITE-DIFFERENCE METHOD
         IF(ICOLL.EQ.1)THEN       ! COLLOCATED GRID
           CALL STRAIN_COL(DX,DY,DZ)
         ELSE                     ! STAGGERED GRID
           CALL STRAIN_STA(DX,DY,DZ)
         END IF
       ELSE                     ! PSEUDO-SPECTRAL METHOD
         CALL STRAIN_SPEC(DY)
       END IF

       IF(MOD(N,N_SGS_SKIP).EQ.0.OR.N.EQ.NSTART)THEN
         IF(N.EQ.0)THEN
           LAG_START=1
         END IF
         IF(ISGS.EQ.3)THEN        ! LASD MODEL 
           CALL SGS_LASD(DX,DY,DZ)
         ELSE IF(ISGS.EQ.2)THEN   ! LASI MODEL
           CALL SGS_LASI(DX,DY,DZ)
         ELSE IF(ISGS.EQ.1)THEN   ! CONSTANT SMAGORINSKY MODEL
          CALL SGS_C(DX,DY,DZ)
         END IF
       END IF
     END IF
!----GET TOTAL VISCOSITY (SGS+MOLECULAR)
     NU=NU+MU/RHO
!----CALCUATE THE STRESS FORCING TERM
     ALLOCATE(SFX(NX1:NX2,NY1:NY2,NZ1:NZ2), &
              SFY(NX1:NX2,NY1:NY2,NZ1:NZ2), &
              SFZ(NX1:NX2,NY1:NY2,NZ1:NZ2))
     SFX=0.0
     SFY=0.0
     SFZ=0.0
     IF(ICOLL.EQ.1)THEN    ! COLLOCATED GRID
       CALL GETSF_COL(DX,DY,DZ,SFX,SFY,SFZ)
     ELSE                  ! STAGGERED GRID
       CALL GETSF_STA(DX,DY,DZ,SFX,SFY,SFZ)
     END IF

     FX=FX+SFX
     FY=FY+SFY
     FZ=FZ+SFZ

     DEALLOCATE (S11,S22,S33,S12,S13,S23,S12C,S13C,S23C,S, &
                 SFX,SFY,SFZ)

     END SUBROUTINE
!=========================================================================!
!             CALCULATING SHEAR FORCE AT COLLOCATED GRID                  !
!=========================================================================!
      SUBROUTINE GETSF_COL(DX,DY,DZ,SFX,SFY,SFZ)

      IMPLICIT NONE
  
      REAL(KIND=DP) :: DX,DY,DZ,TRACE     
      INTEGER :: I,J,K
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: SFX,SFY,SFZ
      REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE:: TAU12T,TAU23T

      DO I=0,NX+1
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU11(I,J,K)=NU(I,J,K)*S11(I,J,K)*2.0
            TAU22(I,J,K)=NU(I,J,K)*S22(I,J,K)*2.0
            TAU33(I,J,K)=NU(I,J,K)*S33(I,J,K)*2.0
            TAU12(I,J,K)=NU(I,J,K)*S12(I,J,K)*2.0
            TAU13(I,J,K)=NU(I,J,K)*S13(I,J,K)*2.0
            TAU23(I,J,K)=NU(I,J,K)*S23(I,J,K)*2.0

            TRACE=(TAU11(I,J,K)+TAU22(I,J,K)+TAU33(I,J,K))/3.0       
            TAU11(I,J,K)=TAU11(I,J,K)-TRACE
            TAU22(I,J,K)=TAU22(I,J,K)-TRACE
            TAU33(I,J,K)=TAU33(I,J,K)-TRACE
          ENDDO
        ENDDO
      ENDDO
!-----WALL MODEL AT THE FIRST POINT OFF THE WALL
      IF(IWALL(1).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU12(0,J,K)=-TAU12W_1(J,K)*2.0-TAU12(1,J,K)
            TAU13(0,J,K)=-TAU13W_1(J,K)*2.0-TAU13(1,J,K)
          END DO
        END DO   
      END IF     

      IF(IWALL(2).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU12(NX+1,J,K)=-TAU12W_2(J,K)*2.0-TAU12(NX,J,K)
            TAU13(NX+1,J,K)=-TAU13W_2(J,K)*2.0-TAU13(NX,J,K)
          END DO
        END DO   
      END IF    

      IF(IWALL(3).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            TAU12(I,0,K)=-TAU12W_3(I,K)*2.0-TAU12(I,1,K)
            TAU23(I,0,K)=-TAU23W_3(I,K)*2.0-TAU23(I,1,K)
          END DO
        END DO  
      END IF     

      IF(IWALL(4).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            TAU12(I,NY+1,K)=-TAU12W_4(I,K)*2.0-TAU12(I,NY,K)
            TAU23(I,NY+1,K)=-TAU23W_4(I,K)*2.0-TAU23(I,NY,K)
          END DO
        END DO  
      END IF  

      IF(IWALL(5).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            TAU13(I,J,0)=-TAU13W_5(I,J)*2.0-TAU13(I,J,1)
            TAU23(I,J,0)=-TAU23W_5(I,J)*2.0-TAU23(I,J,1)
          END DO
        END DO 
      END IF     

      IF(IWALL(6).NE.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            TAU13(I,J,NZ+1)=-TAU13W_6(I,J)*2.0-TAU13(I,J,NZ)
            TAU23(I,J,NZ+1)=-TAU23W_6(I,J)*2.0-TAU23(I,J,NZ)
          END DO
        END DO 
      END IF
!-----CALCULATE THE STRESS FORCES
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            SFX(I,J,K)=  &
            (TAU11(I+1,J,K)-TAU11(I-1,J,K))/(DX*2.0)+ &
            (TAU12(I,J+1,K)-TAU12(I,J-1,K))/(DY*2.0)+ &
            (TAU13(I,J,K+1)-TAU13(I,J,K-1))/(DZ*2.0)

            SFZ(I,J,K)=  &
            (TAU12(I+1,J,K)-TAU12(I-1,J,K))/(DX*2.0)+ &
            (TAU22(I,J+1,K)-TAU22(I,J-1,K))/(DY*2.0)+ &
            (TAU23(I,J,K+1)-TAU23(I,J,K-1))/(DZ*2.0)
 
            SFZ(I,J,K)=  &
            (TAU13(I+1,J,K)-TAU13(I-1,J,K))/(DX*2.0)+ &
            (TAU23(I,J+1,K)-TAU23(I,J-1,K))/(DY*2.0)+ &
            (TAU33(I,J,K+1)-TAU33(I,J,K-1))/(DZ*2.0)
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!             CALCULATING SHEAR FORCE AT STAGGERED GRID                   !
!=========================================================================!
      SUBROUTINE GETSF_STA(DX,DY,DZ,SFX,SFY,SFZ)
      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ,TRACE
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: SFX,SFY,SFZ
      INTEGER:: I,J,K

      DO I=0,NX+1
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU11(I,J,K)=2*NU(I,J,K)*S11(I,J,K)*2.0
            TAU22(I,J,K)=2*NU(I,J,K)*S22(I,J,K)*2.0
            TAU33(I,J,K)=2*NU(I,J,K)*S33(I,J,K)*2.0


            TAU12(I,J,K)=(NU(I,J-1,K)+NU(I-1,J-1,K)+NU(I,J,K)+NU(I-1,J,K))/4.0* &
                         S12(I,J,K)*2.0
      
            TAU13(I,J,K)=(NU(I-1,J,K-1)+NU(I-1,J,K)+NU(I,J,K-1)+NU(I,J,K))/4.0* &
                          S13(I,J,K)*2.0

            TAU23(I,J,K)=(NU(I,J-1,K-1)+NU(I,J-1,K)+NU(I,J,K-1)+NU(I,J,K))/4.0* &
                         S23(I,J,K)*2.0

            TRACE=(TAU11(I,J,K)+TAU22(I,J,K)+TAU33(I,J,K))/3.0       
            TAU11(I,J,K)=TAU11(I,J,K)-TRACE
            TAU22(I,J,K)=TAU22(I,J,K)-TRACE
            TAU33(I,J,K)=TAU33(I,J,K)-TRACE
          ENDDO
        ENDDO
      ENDDO
!-----WALL MODEL AT THE FIRST POINT OFF THE WALL 
! Note: if the shear stress makes the flow parcel rotate counter-clockwisely (from the along-axis point of view, then it is "+", otherwise, it is "-".
      IF(IWALL(1).NE.0.AND.MYIDX.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU12(1,J,K)=-TAU12W_1(J,K)
            TAU13(1,J,K)=-TAU13W_1(J,K)
          END DO
        END DO
      END IF   

      IF(IWALL(2).NE.0.AND.MYIDX.EQ.NPX-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO J=0,NY+1
          DO K=0,NZ+1
            TAU12(NX+1,J,K)= TAU12W_2(J,K)
            TAU13(NX+1,J,K)= TAU13W_2(J,K)
          END DO
        END DO
      END IF    

      IF(IWALL(3).NE.0.AND.MYIDY.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            TAU12(I,1,K)=-TAU12W_3(I,K)
            TAU23(I,1,K)=-TAU23W_3(I,K)
          END DO
        END DO  
      END IF     

      IF(IWALL(4).NE.0.AND.MYIDY.EQ.NPY-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO K=0,NZ+1
            TAU12(I,NY+1,K)= TAU12W_4(I,K)
            TAU23(I,NY+1,K)= TAU23W_4(I,K)
          END DO
        END DO  
      END IF  

      IF(IWALL(5).NE.0.AND.MYIDZ.EQ.0)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            TAU13(I,J,1)=-TAU13W_5(I,J)
            TAU23(I,J,1)=-TAU23W_5(I,J)
          END DO
        END DO   
      END IF     

      IF(IWALL(6).NE.0.AND.MYIDZ.EQ.NPZ-1)THEN           ! USE WALL MODEL SHEAR STRESS
        DO I=0,NX+1
          DO J=0,NY+1
            TAU13(I,J,NZ+1)= TAU13W_6(I,J)
            TAU23(I,J,NZ+1)= TAU23W_6(I,J)
          END DO
        END DO  
      END IF
!-----CALCULATE THE STRESS FORCES
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            SFX(I,J,K)=(TAU11(I,J,  K  )-TAU11(I-1,J,K))/DX+ &
                       (TAU12(I,J+1,K  )-TAU12(I,  J,K))/DY+ &
                       (TAU13(I,J,  K+1)-TAU13(I,  J,K))/DZ

            SFY(I,J,K)=(TAU12(I+1,J,K  )-TAU12(I,J,  K))/DX+ &
                       (TAU22(I,  J,K  )-TAU22(I,J-1,K))/DY+ &
                       (TAU23(I,  J,K+1)-TAU23(I,J,  K))/DZ
          

            SFZ(I,J,K)=(TAU13(I+1,J,  K)-TAU13(I,J,K  ))/DX+ &
                       (TAU23(I,  J+1,K)-TAU23(I,J,K  ))/DY+ &
                       (TAU33(I,  J,  K)-TAU33(I,J,K-1))/DZ  
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
!=========================================================================!
!                 STRAIN RATE TENSOR AT COLLOCATED GRID                   !
!=========================================================================!
      SUBROUTINE STRAIN_COL(DX,DY,DZ)

      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ   
      INTEGER :: I,J,K

      DO I=0,NX+1
        DO J=0,NY+1
          DO K=0,NZ+1
!-----------NORMAL STRAIN AT THE CENTER
            S11(I,J,K)=(U(I+1,J,K)-U(I-1,J,K))/(DX*2.0)         
            S22(I,J,K)=(V(I,J+1,K)-V(I,J-1,K))/(DY*2.0)
            S33(I,J,K)=(W(I,J,K+1)-W(I,J,K-1))/(DZ*2.0)
!-----------SHEAR STRAIN AT THE FACES 
            S12(I,J,K)=((U(I,J+1,K)-U(I,J-1,K))/(DY*2.0)+  &
                        (V(I+1,J,K)-V(I-1,J,K))/(DX*2.0))/2.0

            S13(I,J,K)=((U(I,J,K+1)-U(I,J,K-1))/(DZ*2.0)+  &
                        (W(I+1,J,K)-W(I-1,J,K))/(DX*2.0))/2.0

            S23(I,J,K)=((V(I,J,K+1)-V(I,J,K-1))/(DZ*2.0)+  &
                        (W(I,J+1,K)-W(I,J-1,K))/(DY*2.0))/2.0 
          END DO
        END DO
      END DO         

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
          S(I,J,K)=SQRT((S11(I,J,K)**2+S12(I,J,K)**2+S13(I,J,K)**2+  &
                         S12(I,J,K)**2+S22(I,J,K)**2+S23(I,J,K)**2+  &
                         S13(I,J,K)**2+S23(I,J,K)**2+S33(I,J,K)**2)*2.0)
          ENDDO
        ENDDO
      ENDDO
      END SUBROUTINE
!=========================================================================!
!                 STRAIN RATE TENSOR AT STAGGERED GRID                    !
!=========================================================================!
      SUBROUTINE STRAIN_STA(DX,DY,DZ)
      IMPLICIT NONE
      REAL(KIND=DP):: DX,DY,DZ
      INTEGER I,J,K
      REAL(KIND=DP):: U0,U1,U2,V0,V1,V2,W0,W1,W2,UDY
      REAL(KIND=DP):: X1,X2,Y1,Y2,Z1,Z2

      DO I=0,NX+1
        DO J=0,NY+1
          DO K=0,NZ+1
!-----------NORMAL STRAIN AT THE CENTER
            S11(I,J,K)=(U(I+1,J,K)-U(I,J,K))/DX         
            S22(I,J,K)=(V(I,J+1,K)-V(I,J,K))/DY
            S33(I,J,K)=(W(I,J,K+1)-W(I,J,K))/DZ
!-----------SHEAR STRAIN AT THE FACES 
            S12(I,J,K)=((U(I,J,K)-U(I,J-1,K))/DY+ &
                        (V(I,J,K)-V(I-1,J,K))/DX)/2.0

            S13(I,J,K)=((U(I,J,K)-U(I,J,K-1))/DZ+ &
                        (W(I,J,K)-W(I-1,J,K))/DX)/2.0

            S23(I,J,K)=((V(I,J,K)-V(I,J,K-1))/DZ+ &
                        (W(I,J,K)-W(I,J-1,K))/DY)/2.0
          END DO
        END DO
      END DO         

      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
!-----------SHEAR STRAIN AT THE CENTER
            S12C(I,J,K)=(S12(I,J,K)+S12(I+1,J,K)+  &
                         S12(I,J+1,K)+S12(I+1,J+1,K))/4.0
            S13C(I,J,K)=(S13(I,J,K)+S13(I+1,J,K)+  &
                         S13(I,J,K+1)+S13(I+1,J,K+1))/4.0
            S23C(I,J,K)=(S23(I,J,K)+S23(I,J+1,K)+  &
                         S23(I,J,K+1)+S23(I,J+1,K+1))/4.0
!-----------MODULUS OF THE STRAIN RATE TENSOR
          S(I,J,K)=SQRT((S11 (I,J,K)**2+S12C(I,J,K)**2+S13C(I,J,K)**2+ &
                         S12C(I,J,K)**2+S22 (I,J,K)**2+S23C(I,J,K)**2+ &
                         S13C(I,J,K)**2+S23C(I,J,K)**2+S33 (I,J,K)**2)*2.0)
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE 
!=========================================================================!
!            STRAIN RATE TENSOR USING PSEUDO-SPECTRAL METHOD              !
!=========================================================================!
      SUBROUTINE STRAIN_SPEC(DY)

      IMPLICIT NONE
      REAL(KIND=DP):: DY     
      INTEGER :: I,J,K
      REAL(KIND=DP),DIMENSION(:,:,:), ALLOCATABLE:: UDX,UDZ,VDX,VDZ,WDX,WDZ
      REAL(KIND=DP):: UDY,WDY,DDY1,DDY2
      REAL(KIND=DP):: X1,X2,Y1,Y2,Z1,Z2
!-----CALCULATE VELOCITY GRADIENTS IN THE HORIZONTAL PLANE USING SPECTRAL METHOD
      ALLOCATE(UDX(NX,NY,NZ),UDZ(NX,NY,NZ), &
               VDX(NX,NY,NZ),VDZ(NX,NY,NZ), &
               WDX(NX,NY,NZ),WDZ(NX,NY,NZ)) 

      CALL GRADIENT_SPEC(U, NX1,NY1,NZ1,UDX,UDZ,1)
      CALL GRADIENT_SPEC(VC,NX1,NY1,NZ1,VDX,VDZ,1)
      CALL GRADIENT_SPEC(W, NX1,NY1,NZ1,WDX,WDZ,1)
!-----CALCULATE THE STRAIN RATE TENSOR
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ  
            S11(I,J,K)=UDX(I,J,K)
            S22(I,J,K)=(V(I,J+1,K)-V(I,J,K))/DY
            S33(I,J,K)=WDZ(I,J,K)
 
            UDY=(U(I,J+1,K)-U(I,J-1,K))/(DY*2.0)   
            S12(I,J,K)=(UDY+VDX(I,J,K))/2.0

            S13(I,J,K)=(UDZ(I,J,K)+WDX(I,J,K))/2.0

            WDY=(W(I,J+1,K)-W(I,J-1,K))/(DY*2.0)            
            S23(I,J,K)=(VDZ(I,J,K)+WDY)/2.0
!-----------MODULUS OF THE STRAIN RATE TENSOR
          S(I,J,K)=SQRT((S11(I,J,K)**2+S12(I,J,K)**2+S13(I,J,K)**2+  &
                         S12(I,J,K)**2+S22(I,J,K)**2+S23(I,J,K)**2+  &
                         S13(I,J,K)**2+S23(I,J,K)**2+S33(I,J,K)**2)*2.0)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(UDX,UDZ,VDX,VDZ,WDX,WDZ)
      END SUBROUTINE

  END MODULE
 
