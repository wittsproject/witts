! This module is used for calculating viscous forces.
!
  MODULE viscous

  USE parameters
  USE field_shared
  USE class_shared
  USE class_cell
  USE tools
  USE sgs

  IMPLICIT NONE

  CONTAINS
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
!     RANGE: 20-25      
      SUBROUTINE STRESS_CELL(DX,DY,DZ,SFX,SFY,SFZ)
      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ,TRACE
      REAL(KIND=DP),DIMENSION(NX1:,NY1:,NZ1:):: SFX,SFY,SFZ
      INTEGER:: I,J,K

      DO I=0,2
        CELL_FV(INDEX)%CELL_VAR(20+I)=2*NU_TOTAL*CELL_FV(INDEX)%CELL_VAR(13+I)*2.0
      END DO

      IF(ICOLL.EQ.1)THEN
        DO I=3,5
          CELL_FV(INDEX)%CELL_VAR(20+I)=2*NU_TOTAL*CELL_FV(INDEX)%CELL_VAR(13+I)*2.0
        END DO
      ELSE 
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
!                 STRAIN RATE TENSOR AT STAGGERED GRID                    !
!=========================================================================!
      SUBROUTINE STRAIN(DX,DY,DZ,I,J,K, &
                        VEL1,VEL2,VEL3,SI1,SI2,SI3, &
                        SR,S)
      IMPLICIT NONE

      INTEGER:: SI1,SI2,SI3,I,J,K,II,JJ,KK
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VEL1,VEL2,VEL3      
      REAL(KIND=DP):: DX,DY,DZ,SR(6),S      
      REAL(KIND=DP):: S12T(0:1,0:1),S13T(0:1,0:1),S23T(0:1,0:1)

!-----------NORMAL STRAIN AT THE CENTER
      SR(1)=(VEL1(I+1,J,K)-VEL1(I,J,K))/DX         
      SR(2)=(VEL2(I,J+1,K)-VEL2(I,J,K))/DY
      SR(3)=(VEL3(I,J,K+1)-VEL3(I,J,K))/DZ
!-----------COLLOCATED GRID      
      IF(ICOLL.EQ.1)THEN
        SR(4)=((VEL1(I,J+1,K)-VEL1(I,J-1,K))/(DY*2.0)+  &
               (VEL2(I+1,J,K)-VEL2(I-1,J,K))/(DX*2.0))/2.0

        SR(5)=((VEL1(I,J,K+1)-VEL1(I,J,K-1))/(DZ*2.0)+  &
               (VEL3(I+1,J,K)-VEL3(I-1,J,K))/(DX*2.0))/2.0
      
        SR(6)=((VEL2(I,J,K+1)-VEL2(I,J,K-1))/(DZ*2.0)+  &
             (VEL3(I,J+1,K)-VEL3(I,J-1,K))/(DY*2.0))/2.0
!-----------STAGGERED GRID        
      ELSE  
        DO II=0,1
          DO JJ=0,1 
            S12T(II,JJ)=((VEL1(I+II,J+JJ,K)-VEL1(I+II,J-1+JJ,K))/DY+ &
                         (VEL2(I+II,J+JJ,K)-VEL2(I-1+II,J+JJ,K))/DX)/2.0
          END DO
        END DO
     
        DO II=0,1
          DO KK=0,1
            S13T(II,KK)=((VEL1(I+II,J,K+KK)-VEL1(I+II,J,K-1+KK))/DZ+ &
                         (VEL3(I+II,J,K+KK)-VEL3(I-1+II,J,K+KK))/DX)/2.0
          END DO
        END DO 

        DO JJ=0,1
          DO KK=0,1 
            S23T(JJ,KK)=((VEL2(I,J+JJ,K+KK)-VEL2(I,J+JJ,K-1+KK))/DZ+ &
                         (VEL3(I,J+JJ,K+KK)-VEL3(I,J-1+JJ,K+KK))/DY)/2.0
          END DO
        END DO  
!       SHEAR STRAIN AT THE CENTER
        SR(4)=0.0      
        DO II=0,1
          DO JJ=0,1
            SR(4)=SR(4)+S12T(II,JJ)
          END DO
        END DO
     
        SR(5)=0.0      
        DO II=0,1
          DO KK=0,1
            SR(5)=SR(5)+S13T(II,KK)
          END DO
        END DO

        SR(6)=0.0      
        DO KK=0,1
          DO JJ=0,1
            SR(6)=SR(6)+S23T(JJ,KK)
          END DO
        END DO
      END IF 
!-----------MODULUS OF THE STRAIN RATE TENSOR
      S=SQRT((SR(1)**2+SR(4)**2+SR(5)**2+  &
              SR(4)**2+SR(2)**2+SR(6)**2+  &
              SR(5)**2+SR(6)**2+SR(3)**2)*2.0)      
      END SUBROUTINE 
!=========================================================================!
!                 CALCULATE THE VISCOUS FORCE (MOLECULAR+SGS)             !
!=========================================================================!
     SUBROUTINE VISCOUS_WRAP_CELL()
     IMPLICIT NONE
     INTEGER :: M,NB,I,J,K
     REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: VEL1,VEL2,VEL3
     REAL(KIND=DP):: DX,DY,DZ,SFX,SFY,SFZ,DXC,DYC,DZC
     REAL(KIND=DP):: SR(6),S,NU_TOTAL

     NB=2

     ALLOCATE(VEL1(-NB:NB),VEL2(-NB:NB),VEL3(-NB:NB))

     DO M=1,TOTAL_CELL
       IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
         CALL CELL_TO_STRUCT(M,2,1,VEL1)
         CALL CELL_TO_STRUCT(M,2,2,VEL2)
         CALL CELL_TO_STRUCT(M,2,3,VEL3)

         DXC=CELL_FV(M)%CELL_DX
         DYC=CELL_FV(M)%CELL_DY
         DZC=CELL_FV(M)%CELL_DZ     
!----CALCUALTE THE EDDY VISCOSITY
         IF(ITYPE.EQ.1)THEN        ! IF LES IS USED
           CALL STRAIN(DXC,DYC,DZC,0,0,0, &  ! CALCULATE STRAIN RATE TENSOR
                       VEL1,VEL2,VEL2,-NB,-NB,-NB, &
                       SR,S)
           DO I=1,6
             CELL_FV(M)%CELL_VAR(12+I)=SR(I)  ! RANGE: 13-18
           END DO
           CELL_FV(M)%CELL_VAR(19)=S      

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
     END DO 
!----GET TOTAL VISCOSITY (SGS+MOLECULAR)
     NU_TOTAL=CELL_FV(INDEX)%CELL_VAR(6)+ &
              CELL_FV(INDEX)%CELL_VAR(7)/CELL_FV(INDEX)%CELL_VAR(8)
!----CALCUATE THE STRESS FORCING TERM

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

     DEALLOCATE (SFX,SFY,SFZ)

     END SUBROUTINE
  END MODULE
 
