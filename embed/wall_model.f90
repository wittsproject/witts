! This module is used for wall models.
!
  MODULE wall_model

  USE mpi  
  USE parameters
  USE class_shared

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!                              WALL MODEL                                 !
!=========================================================================!
    SUBROUTINE WALL_MODEL_WRAP()

    IMPLICIT NONE
    INTEGER :: ID,M
    REAL(KIND=DP):: TAUW1,TAUW2,QWALL

    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_WALL.NE.0)THEN
        ID=CELL_FV(M)%CELL_WALL
!-------AT BOUNDARY #1 & #2----------------------------------------------------------    
        IF(IWALL(1).EQ.1.AND.MYIDX.EQ.0.OR.IWALL(2).EQ.1.AND.MYIDX.EQ.NPX-1)THEN  
          CALL WALLMODEL(CELL_FV(M)%CELL_DX/2.0,CELL_FV(M)%CELL_VEL(2),CELL_FV(M)%CELL_VEL(3), &
                         CELL_FV(M)%CELL_TE,BV(ID,4),Z0(ID),BC(ID,4),QS(ID),TAUW1,TAUW2,QWALL)
          CELL_FV(M)%CELL_TAU(4)=TAUW1    ! TAU12
          CELL_FV(M)%CELL_TAU(5)=TAUW2    ! TAU13
          CELL_FV(M)%CELL_HF(1)=QWALL    ! Q1
        END IF   
!-------AT BOUNDARY #3 & #4----------------------------------------------------------    
        IF(IWALL(3).EQ.1.AND.MYIDY.EQ.0.OR.IWALL(4).EQ.1.AND.MYIDY.EQ.NPY-1)THEN  
          CALL WALLMODEL(CELL_FV(M)%CELL_DY/2.0,CELL_FV(M)%CELL_VEL(1),CELL_FV(M)%CELL_VEL(3), &
                         CELL_FV(M)%CELL_TE,BV(ID,4),Z0(ID),BC(ID,4),QS(ID),TAUW1,TAUW2,QWALL)
          CELL_FV(M)%CELL_TAU(4)=TAUW1    ! TAU12
          CELL_FV(M)%CELL_TAU(6)=TAUW2    ! TAU23
          CELL_FV(M)%CELL_HF(2)=QWALL    ! Q2
        END IF   
!-------AT BOUNDARY #5 & #6----------------------------------------------------------    
        IF(IWALL(5).EQ.1.AND.MYIDZ.EQ.0.OR.IWALL(6).EQ.1.AND.MYIDZ.EQ.NPZ-1)THEN  
          CALL WALLMODEL(CELL_FV(M)%CELL_DZ/2.0,CELL_FV(M)%CELL_VEL(1),CELL_FV(M)%CELL_VEL(2), &
                         CELL_FV(M)%CELL_TE,BV(ID,4),Z0(ID),BC(ID,4),QS(ID),TAUW1,TAUW2,QWALL)
          CELL_FV(M)%CELL_TAU(5)=TAUW1    ! TAU13
          CELL_FV(M)%CELL_TAU(6)=TAUW2    ! TAU23
          CELL_FV(M)%CELL_HF(3)=QWALL    ! Q3
        END IF   
      END IF
    END DO
   
    END SUBROUTINE WALL_MODEL_WRAP
!========================================================================!
!                 MODELLING WALL SHEAR STRESS AND HEAT FLUX              !
!========================================================================!
!  INPUT: D,U1_1,U1_1,T1,TB,Z0,BC_T,QS 
!  OUTPUT: USTAR,TAUW1,TAUW2,QWALL
      SUBROUTINE WALLMODEL(D,U1_1,U1_2,T1,TB,Z0I,BC_T,QS0, &
                           TAUW1,TAUW2,QWALL)   

      IMPLICIT NONE
      REAL(KIND=DP) :: D,U1_1,U1_2,T1,TB,Z0I,QS0
      REAL(KIND=DP) :: USTAR,TAUW1,TAUW2,QWALL
      REAL(KIND=DP) :: ZERO,GAMMAM,GAMMAH,TOL,D1,DELTAT,TAUW,UH
      REAL(KIND=DP) :: USTAR_OLD,QWALL_OLD,L,PSIM,PSIH,X
      INTEGER :: IT,ITMAX,BC_T

      ZERO=1.E-12
      GAMMAM = 4.8
      GAMMAH = 7.8
      TOL=1.0E-6
      ITMAX=30

      D1=D/2.0

      UH=SQRT(U1_1**2+U1_2**2)

      DELTAT=T1-TB

      ! Set initial guesse of u* and qw (they are also used for the neutral condition)
      USTAR = (KAPPA*UH)/LOG(D1/Z0I)
      USTAR = DMAX1(USTAR,0.0D0)     ! Limit u* to always be zero or positive
      QWALL = (-DELTAT*USTAR*KAPPA) /(LOG(D1/Z0I))

      ! stable-----------------------------------------------------------------------
      ! see book "Modelling of Atmospheric Flow Field", editors D. Lalas and C. Ratto
      ! chapter 2--Modelling the Vertical ABL Structure, D. Etling, pp. 56--57
      IF(DELTAT.GT.ZERO)THEN
        IT=1
                
        ! limit qw0 to always be negative and non-zero  
 1      QWALL = DMIN1(QWALL,-1.0D-10)
        USTAR_OLD = USTAR
        QWALL_OLD = QWALL
        ! set initial guesses at L
        L = -(USTAR**3.0)*T1/(KAPPA*G*QWALL)
        ! limit L to always be positive and finite
        L = DMAX1(L,1.0D-10)
        ! update u*
        PSIM = -GAMMAM*D1/L
        USTAR = (KAPPA*UH)/(LOG(D1/Z0I)-PSIM)
        ! update qw
        IF(BC_T.EQ.2)THEN
          QWALL=QS0
        ELSE
          PSIH = -GAMMAH*D1/L
          QWALL =-(KAPPA*USTAR*DELTAT)/(LOG(D1/Z0I)-PSIH)
        END IF

        IF((ABS(USTAR-USTAR_OLD).GT.TOL.OR.ABS(QWALL-QWALL_OLD).GT.TOL).AND.IT.LT.ITMAX)THEN
          IT=IT+1
        GOTO 1
        ELSE IF(IT.EQ.ITMAX)THEN
          IF(MYID.EQ.0)THEN
            PRINT*,'Max qw, uStar iterations reached!!!'
          END IF
        END IF
        ! unstable-----------------------------------------------------------------------
        ! see book "Modelling of Atmospheric Flow Field", editors D. Lalas and C. Ratto
        ! chapter 2--Modelling the Vertical ABL Structure, D. Etling, pp. 56--57 and
        ! article "The Mathematical Representation of Wind Speed and Temperature Profiles
        ! in the Unstable Atmospheric Surface Layer", C. Paulson, Journal of Applied Meteorology, Vol 9, 1970, pp. 857--861.
      ELSE IF(DELTAT.LT.-ZERO)THEN
        IT=1

 2      USTAR_OLD=USTAR
        QWALL_OLD = QWALL

        L = -USTAR**3*T1/(KAPPA*G*QWALL)     ! LOCAL OBUKHOV LENGTH

        X = (1.0-15.0*D1/L)**0.25
        ! STABILITY CORRECTION (Stull (1988) and Arya (2001)):
        PSIM = 2.0*LOG((1.0+X)/2.0)+LOG((1.0+X**2)/2.0)-2.0*ATAN(X)+PI/2.0               
        USTAR = (KAPPA*UH)/(LOG(D1/Z0I)-PSIM)

         IF(BC_T.EQ.2)THEN
          QWALL=QS0
        ELSE         
          PSIH = 2.0*LOG((1.0+X**2)/2.0)
          QWALL = -(KAPPA*USTAR*DELTAT)/(LOG(D1/Z0I)-PSIH)
        END IF
   
        IF((ABS(USTAR-USTAR_OLD).GT.TOL.OR.ABS(QWALL-QWALL_OLD).GT.TOL).AND.IT.LT.ITMAX)THEN
          IT=IT+1
          GOTO 2
        ELSE IF(IT.EQ.ITMAX)THEN
           PRINT*,'Max qw, uStar iterations reached!!!'
        END IF
      END IF

      TAUW=-USTAR**2
      TAUW1=TAUW*(U1_1/UH)
      TAUW2=TAUW*(U1_2/UH)

      END SUBROUTINE
 
  END MODULE
