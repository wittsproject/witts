! This module is used for wall models.
!
  MODULE wall_model

  USE parameters
  USE field_shared, ONLY: U,V,W,TE,TAU12W_1,TAU13W_1,TAU12W_2,TAU13W_2, &
                          TAU12W_3,TAU23W_3,TAU12W_4,TAU23W_4,TAU13W_5, &
                          TAU23W_5,TAU13W_6,TAU23W_6,Q1W_1,Q1W_2,Q2W_3, &
                          Q2W_4,Q3W_5,Q3W_6
  USE tools, ONLY: AVE_P

  IMPLICIT NONE

  CONTAINS
!=========================================================================!
!                              WALL MODEL                                 !
!=========================================================================!
    SUBROUTINE WALL_MODEL_WRAP(DX,DY,DZ)

    IMPLICIT NONE
    REAL(KIND=DP) :: DX,DY,DZ
    INTEGER :: I,J,K

    IF(IWALL(1).EQ.1)THEN
      IF(MYIDX.EQ.0)THEN
        DO J=0,NY+1
           DO K=0,NZ+1
            CALL WALLMODEL(DX/2.0,V(1,J,K),W(1,J,K),TE(1,J,K),BV(1,4),Z0(1),BC(1,4),QS(1),    &
                            TAU12W_1(J,K),TAU13W_1(J,K),Q1W_1(J,K))
          END DO
        END DO
      END IF
    END IF
       
    IF(IWALL(2).EQ.1)THEN
      IF(MYIDX.EQ.NPX-1)THEN
        DO J=0,NY+1
          DO K=0,NZ+1
            CALL WALLMODEL(DX/2.0,V(NX,J,K),W(NX,J,K),TE(NX,J,K),BV(2,4),Z0(2),BC(2,4),QS(2), &
                            TAU12W_2(J,K),TAU13W_2(J,K),Q1W_2(J,K))
          END DO
        END DO
      END IF
    END IF

    IF(IWALL(3).EQ.1)THEN
      IF(MYIDY.EQ.0)THEN
        DO I=0,NX+1
          DO K=0,NZ+1
            CALL WALLMODEL(DY/2.0,U(I,1,K),W(I,1,K),TE(I,1,K),BV(3,4),Z0(3),BC(3,4),QS(3),    &
                            TAU12W_3(I,K),TAU23W_3(I,K),Q2W_3(I,K))
          END DO
        END DO
      END IF
    END IF

    IF(IWALL(4).EQ.1)THEN
      IF(MYIDY.EQ.NPY-1)THEN
        DO I=0,NX+1
          DO K=0,NZ+1
            CALL WALLMODEL(DY/2.0,U(I,NY,K),W(I,NY,K),TE(I,NY,K),BV(4,4),Z0(4),BC(4,4),QS(4), &
                            TAU12W_4(I,K),TAU23W_4(I,K),Q2W_4(I,K))
          END DO
        END DO
      END IF
    END IF

    IF(IWALL(5).EQ.1)THEN
      IF(MYIDZ.EQ.0)THEN
        DO I=0,NX+1
          DO J=0,NY+1
            CALL WALLMODEL(DZ/2.0,U(I,J,1),W(I,J,1),TE(I,J,1),BV(5,4),Z0(5),BC(5,4),QS(5),    &
                           TAU13W_5(I,J),TAU23W_5(I,J),Q3W_5(I,J))
          END DO
        END DO
      END IF
    END IF
 
    IF(IWALL(6).EQ.1)THEN
      IF(MYIDZ.EQ.NPZ-1)THEN
        DO I=0,NX+1
          DO J=0,NY+1
            CALL WALLMODEL(DZ/2.0,U(I,J,NZ),W(I,J,NZ),TE(I,J,NZ),BV(6,4),Z0(6),BC(6,4),QS(6), &
                            TAU13W_6(I,J),TAU23W_6(I,J),Q3W_6(I,J))
          END DO
        END DO
      END IF
    END IF
    END SUBROUTINE 
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
