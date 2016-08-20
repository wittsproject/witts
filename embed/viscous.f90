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
!                 CALCULATE THE VISCOUS FORCE (MOLECULAR+SGS)             !
!=========================================================================!
     SUBROUTINE VISCOUS_CELL_WRAP()
     IMPLICIT NONE
     INTEGER :: M,NB,I,J,K
     REAL(KIND=DP):: FV(3)

     DO M=1,TOTAL_CELL
       IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
         CALL STRAIN_CELL(M) 
       END IF
    END DO

    DO I=13,19
      CALL GHOST_BOUNDARY(I)
    END DO   
!---FOR SGS MODELING---------------------------------
    IF(MOD(N,N_SGS_SKIP).EQ.0.OR.N.EQ.NSTART)THEN
      IF(N.EQ.0)THEN
        LAG_START=1
      END IF
      IF(ISGS.EQ.3)THEN        ! LASD MODEL 
        CALL SGS_LASD()
      ELSE IF(ISGS.EQ.2)THEN   ! LASI MODEL
!        CALL SGS_LASI()
      ELSE IF(ISGS.EQ.1)THEN   ! CONSTANT SMAGORINSKY MODEL
        CALL SGS_C()
      END IF
    END IF
!---CALCUATE THE STRESS TENSOR
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
        CALL STRESS_CELL(M)
      END IF
    END DO
!---WALL SHEAR STRESS
    CALL WALL_MODEL_WRAP()
!---GET THE STRESS TENSOR ON GHOST CELLS
    DO I=20,25        
      CALL GHOST_BOUNDARY(I)  
    END DO 
!---CALCUATE THE VISCOUS FORCES
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
        CALL VISCOUS_FORCE_CELL(M,VF)
        DO I=0,2
          CELL_FV(M)%CELL_VAR(10+I)=CELL_FV(M)%CELL_VAR(10+I)+VF(I+1)
        END DO
      END IF
    END DO
          
    END SUBROUTINE VISCOUS_CELL_WRAP  
!=========================================================================!
!                 CALCULATE STRESS TENSOR ON A FV CELL                    !
!=========================================================================!
!     RANGE: 20-25      
      SUBROUTINE STRESS_CELL(INDEX)
      IMPLICIT NONE

      REAL(KIND=DP):: DX,DY,DZ,TRACE
      REAL(KIND=DP):: TAU(6)
      REAL(KIND=DP),DIMENSION(-1:1,-1:1,-1:1):: NU,VEL1,VEL2,VEL3, &
                                                NU_STR,MU_STR,RHO_STR
      INTEGER:: I,J,K,INDEX,M

!-----For collocated grid
      IF(ICOLL.EQ.1)THEN
        NU(0,0,0)=CELL_FV(INDEX)%CELL_VAR(6)+    &    
                  CELL_FV(INDEX)%CELL_VAR(7)/CELL_FV(INDEX)%CELL_VAR(8)
        DO I=0,6 
          CELL_FV(INDEX)%CELL_VAR(20+I)=NU(0,0,0)*CELL_FV(INDEX)%CELL_VAR(13+I)*2.0
        END DO  
!-----For staggered grid
      ELSE
        CALL CELL_TO_STRUCT(INDEX,1,1,VEL1)
        CALL CELL_TO_STRUCT(INDEX,1,2,VEL2)
        CALL CELL_TO_STRUCT(INDEX,1,3,VEL3)
        CALL CELL_TO_STRUCT(INDEX,1,6,NU_STR)
        CALL CELL_TO_STRUCT(INDEX,1,7,MU_STR)
        CALL CELL_TO_STRUCT(INDEX,1,8,RHO_STR)

        DO I=-1,1
          DO J=-1,1
            DO K=-1,1
              NU(I,J,K)=NU_STR(I,J,K)+MU_STR(I,J,K)/RHO_STR(I,J,K)
            END DO
          END DO
        END DO

        DX=CELL_FV(%)CELL_DX
        DY=CELL_FV(%)CELL_DY
        DZ=CELL_FV(%)CELL_DZ

        TAU(1)=NU(0,0,0)*(VEL1(1,0,0)-VEL1(0,0,0))/DX 
        TAU(2)=NU(0,0,0)*(VEL2(0,1,0)-VEL2(0,0,0))/DY
        TAU(3)=NU(0,0,0)*(VEL3(0,0,1)-VEL3(0,0,0))/DZ 
        TAU(4)=(NU(0,-1,0)+NU(-1,-1,0)+NU(0,0,0)+NU(-1,0,0))/4.0* &
               ((VEL1(0,0,0)-VEL1(0,-1,0))/DY+(VEL2(0,0,0)-VEL2(-1,0,0))/DX)
      
        TAU(5)=(NU(-1,0,0-1)+NU(0-1,0,0)+NU(0,0,-1)+NU(0,0,0))/4.0* &
               ((VEL1(0,0,0)-VEL1(0,0,0-1))/DZ+(VEL3(0,0,0)-VEL3(-1,0,0))/DX)       

        TAU(6)=(NU(0,-1,-1)+NU(0,0-1,0)+NU(0,0,-1)+NU(0,0,0))/4.0* &
               ((VEL2(0,0,0)-VEL2(0,0,-1))/DZ+(VEL3(0,0,0)-VEL3(0,-1,0))/DY)

        TRACE=(TAU(1)+TAU(2)+TAU(3))/3.0 
        DO I=1,3      
          TAU(I)=TAU11(I)-TRACE
        END DO
!-------TRANSFER THE STRUCTURED ARRAY BACK TO THE FV-CELL ARRAY
        DO M=0,5
          CELL_FV(INDEX)%CELL_VAR(20+M)=TAU(M+1)
        END DO
      END IF

      END SUBROUTINE STRESS_STR
!=========================================================================!
!                CALCULATE VISCOUS FORCES ON A FV CELL                    !
!=========================================================================!
!     RANGE: 20-25      
      SUBROUTINE VISCOUS_FORCE_CELL(INDEX,VF)
      IMPLICIT NONE

      INTEGER:: INDEX
      REAL(KIND=DP),DIMENSION(-1:1,-1:1,-1:1):: TAU11,TAU22,TAU33,TAU12,TAU13,TAU23
      REAL(KIND=DP):: DX,DY,DZ,VF(3)


      CALL CELL_TO_STRUCT(INDEX,1,20,TAU11)
      CALL CELL_TO_STRUCT(INDEX,1,21,TAU22)
      CALL CELL_TO_STRUCT(INDEX,1,22,TAU33)
      CALL CELL_TO_STRUCT(INDEX,1,23,TAU12)
      CALL CELL_TO_STRUCT(INDEX,1,24,TAU13)
      CALL CELL_TO_STRUCT(INDEX,1,25,TAU23)

      DX=CELL_FV(%)CELL_DX
      DY=CELL_FV(%)CELL_DY
      DZ=CELL_FV(%)CELL_DZ
!-----CALCULATE THE STRESS FORCES
      IF(ICOLL.EQ.1)THEN
        VF(1)=(TAU11(1,0,0)-TAU11(-1,0,0))/(DX*2.0)+ &
              (TAU12(0,1,0)-TAU12(0,-1,0))/(DY*2.0)+ &
              (TAU13(0,0,1)-TAU13(0,0,-1))/(DZ*2.0)

        VF(2)=(TAU12(1,0,0)-TAU12(-1,0,0))/(DX*2.0)+ &
              (TAU22(0,1,0)-TAU22(0,-1,0))/(DY*2.0)+ &
              (TAU23(0,0,1)-TAU23(0,0,-1))/(DZ*2.0)
 
        VF(3)=(TAU13(1,0,0)-TAU13(-1,0,0))/(DX*2.0)+ &
              (TAU23(0,1,0)-TAU23(0,-1,0))/(DY*2.0)+ &
              (TAU33(0,0,1)-TAU33(0,0,-1))/(DZ*2.0)
      ELSE
        VF(1)=(TAU11(0,0,0)-TAU11(-1,0,0))/DX+ &
              (TAU12(0,1,0)-TAU12( 0,0,0))/DY+ &
              (TAU13(0,0,1)-TAU13( 0,0,0))/DZ

        VF(2)=(TAU12(1,0,0)-TAU12(0, 0,0))/DX+ &
              (TAU22(0,0,0)-TAU22(0,-1,0))/DY+ &
              (TAU23(0,0,1)-TAU23(0, 0,0))/DZ
          
        VF(3)=(TAU13(1,0,0)-TAU13(0,0, 0))/DX+ &
              (TAU23(0,1,0)-TAU23(0,0, 0))/DY+ &
              (TAU33(0,0,0)-TAU33(0,0,-1))/DZ  
      END IF

      END SUBROUTINE
!=========================================================================!
!                 STRAIN RATE TENSOR AT STAGGERED GRID                    !
!=========================================================================!
      SUBROUTINE STRAIN_CELL(INDEX)
      IMPLICIT NONE

      INTEGER:: INDEX,I,J,K,II,JJ,KK
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: VEL1,VEL2,VEL3   
      REAL(KIND=DP):: DX,DY,DZ,SR(6),S      
      REAL(KIND=DP):: S12T(0:1,0:1),S13T(0:1,0:1),S23T(0:1,0:1)

      NB=1

      ALLOCATE(VEL1(-NB:NB,-NB:NB,-NB:NB), &
               VEL2(-NB:NB,-NB:NB,-NB:NB), &
               VEL3(-NB:NB,-NB:NB,-NB:NB))

      CALL CELL_TO_STRUCT(INDEX,NB,1,VEL1)
      CALL CELL_TO_STRUCT(INDEX,NB,2,VEL2)
      CALL CELL_TO_STRUCT(INDEX,NB,3,VEL3)

      DX=CELL_FV(M)%CELL_DX
      DY=CELL_FV(M)%CELL_DY
      DZ=CELL_FV(M)%CELL_DZ
!-----NORMAL STRAIN AT THE CENTER
      SR(1)=(VEL1(1,0,0)-VEL1(0,0,0))/DX         
      SR(2)=(VEL2(0,1,0)-VEL2(0,0,0))/DY
      SR(3)=(VEL3(0,0,1)-VEL3(0,0,0))/DZ
!-----COLLOCATED GRID      
      IF(ICOLL.EQ.1)THEN
        SR(4)=((VEL1(0,1,0)-VEL1(0,-1,0))/(DY*2.0)+  &
               (VEL2(1,0,0)-VEL2(-1,0,0))/(DX*2.0))/2.0

        SR(5)=((VEL1(0,0,1)-VEL1(0,0,-1))/(DZ*2.0)+  &
               (VEL3(1,0,0)-VEL3(-1,0,0))/(DX*2.0))/2.0
      
        SR(6)=((VEL2(0,0,1)-VEL2(0,0,-1))/(DZ*2.0)+  &
               (VEL3(0,1,0)-VEL3(0,-1,))/(DY*2.0))/2.0
!-----STAGGERED GRID        
      ELSE  
        DO II=0,1
          DO JJ=0,1 
            S12T(II,JJ)=((VEL1(II,JJ,0)-VEL1(II,-1+JJ,0))/DY+ &
                         (VEL2(II,JJ,0)-VEL2(-1+II,JJ,0))/DX)/2.0
          END DO
        END DO
     
        DO II=0,1
          DO KK=0,1
            S13T(II,KK)=((VEL1(II,0,KK)-VEL1(II,0,-1+KK))/DZ+ &
                         (VEL3(II,0,KK)-VEL3(-1+II,0,KK))/DX)/2.0
          END DO
        END DO 

        DO JJ=0,1
          DO KK=0,1 
            S23T(JJ,KK)=((VEL2(0,JJ,KK)-VEL2(0,JJ,-1+KK))/DZ+ &
                         (VEL3(0,JJ,KK)-VEL3(0,-1+JJ,KK))/DY)/2.0
          END DO
        END DO  
!       SHEAR STRAIN AT THE CENTER
        SR(4)=0.0      
        DO II=0,1
          DO JJ=0,1
            SR(4)=SR(4)+S12T(II,JJ)
          END DO
        END DO
        SR(4)=SR(4)/4.0
     
        SR(5)=0.0      
        DO II=0,1
          DO KK=0,1
            SR(5)=SR(5)+S13T(II,KK)
          END DO
        END DO
        SR(5)=SR(5)/4.0

        SR(6)=0.0      
        DO KK=0,1
          DO JJ=0,1
            SR(6)=SR(6)+S23T(JJ,KK)
          END DO
        END DO
      END IF
      SR(6)=SR(6)/4.0
!-----MODULUS OF THE STRAIN RATE TENSOR
      S=SQRT((SR(1)**2+SR(4)**2+SR(5)**2+  &
              SR(4)**2+SR(2)**2+SR(6)**2+  &
              SR(5)**2+SR(6)**2+SR(3)**2)*2.0)      
!-----TRANSFER THE STRUCTURED ARRAY BACK TO THE FV-CELL ARRAY
      DO I=1,6
        CELL_FV(M)%CELL_VAR(12+I)=SR(I)  ! RANGE: 13-18
      END DO
      CELL_FV(M)%CELL_VAR(19)=S 

      DEALLOCATE(VEL1,VEL2,VEL3)

      END SUBROUTINE 

  END MODULE
 
