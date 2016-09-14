! This module is for solving the scalar transport equation.
!
  MODULE scalar 

  USE mpi
  USE parameters
  USE class_shared
  USE call_cell
  USE convect
  USE tools

  CONTAINS

!=============================================================================!
!	             SUBROUTINE OF SOLVING SCALAR TRANSPORT EQ 		      !
!=============================================================================!
!       Update scalar field by solving:
!       D (s)/ Dt = - div (q), where D()/ Dt is the total derivative, and q is the 
!       scalar flux
!       NUM=4: s is temperature
        SUBROUTINE SCALAR_WRAP(NUM)

        IMPLICIT NONE 
        INTEGER :: N,NUM,N_RK_SCA,SI1,SI2,SI3,NB,SCALAR_FLAG
        REAL(KIND=DP):: DX,DY,DZ
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: C,B,Q1,Q2,Q3
        REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: VAR0,VARI,VIS,F
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
       
        ALLOCATE(Q1(TOTAL_CELL),Q2(TOTAL_CELL),Q3(TOTAL_CELL))
        ALLOCATE(VAR0(TOTAL_CELL),VARI(TOTAL_CELL),F(TOTAL_CELL))

        DO I=1,TOTAL_CELL
          VAR0(M)=CELL_FV(M)%CELL_VAR(NUM)
          VARI(M)=CELL_FV(M)%CELL_VAR(NUM)
        END DO
!-------GET SCALAR FLUX ON CELL FACES
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.CELL_FV(M)%CELL_WALL.EQ.0)THEN
            Q1(M)=SCALAR_FLUX(VAR0,M,SCHMIDT,1)
            Q2(M)=SCALAR_FLUX(VAR0,M,SCHMIDT,2)
            Q3(M)=SCALAR_FLUX(VAR0,M,SCHMIDT,3)
          END IF
        END IF
        CALL GHOST_BOUNDARY(Q1)
        CALL GHOST_BOUNDARY(Q2)
        CALL GHOST_BOUNDARY(Q3)

        IF(NUM.EQ.4)THEN
          DO M=1,TOTAL_CELL
            CELL_FV(M)%CELL_VAR(26)=Q1(M)
            CELL_FV(M)%CELL_VAR(27)=Q2(M)
            CELL_FV(M)%CELL_VAR(28)=Q3(M)
          END DO
        END IF 
!-------R-K METHOD TO SOLVE THE SCALAR TRANSPORT EQUATION        
        DO N=1,N_RK_SCA
          DO M=1,TOTAL_CELL
            IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
              VIS=DERIV_CELL(Q1,M,1,1,10,2)+ &
                  DERIV_CELL(Q2,M,2,1,10,2)+ &
                  DERIV_CELL(Q3,M,3,1,10,2)  
              F(M)=VIS+CONVECT(VARI,M,BLEND_SCA,ORDER_SCA)
            END IF
          END DO

          IF(N.LT.N_RK_SCA)THEN
            DO M=1,TOTAL_CELL
              IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN 
                VARI(M)=VAR0(M)+C(N+1)*DT*F(M)
              END IF
            END DO              
          END IF

          DO M=1,TOTAL_CELL
            IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN                  
              CELL_FV(M)%CELL_VAR(NUM)=CELL_FV(M)%CELL_VAR(NUM)+ &
                                       B(N)*DT*F(M)
            END IF
          END DO             
       END DO     

       CALL GHOST_BOUNDARY(CELL_FV(:)%CELL_VAR(NUM))      
        
       DEALLOCATE(B,C,VAR0,VARI,F,Q1,Q2,Q3)

       END SUBROUTINE
!=========================================================================!
!                            GET SCALAR FLUX                              !
!=========================================================================!
!   This subroutine is used to claclulate scalar flux at cell faces:
!   ID=1: q1=-c*dQ/dx at lower x face,
!   ID=2: q2=-c*dQ/dy at lower y face,
!   ID=3: q3=-c*dQ/dz at lower z face.
!   Here, c is scalar diffusivity (or conductivity), Q is the scalar (VAR is used)
!   c is modeled as: c=nu/Sc, Sc is the Schmidt number (Prandtl numer for heat flux)       
!   The results are stored in CELL_VAR(0) temporarily.       
    REAL(KIND=DP) FUNCTION SCALAR_FLUX(VAR,INDEX,SCHMIDT,ID)
    IMPLICIT NONE
    INTEGER :: INDEX,ID,SCHEME,ORDER,ISTAG
    REAL(KIND=DP),DIMEMSION(:):: VAR
    REAL(KIND=DP) :: C,SCHMIDT

    ID=1

    SCHEME=1
    ISTAG=1
    ORDER=2

    C=CELL_FV(INDEX)%CELL_VAR(6)/SCHMIDT
       
    SCALAR_FLUX=-C*DERIV_CELL(VAR,INDEX,ID,SCHEME,ISTAG,ORDER)

    END FUNCTION SCALAR_FLUX  
      
  END MODULE  

      
