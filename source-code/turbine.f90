! This module is use for modeling of wind turbines.
!
  MODULE turbine

  USE mpi
  USE parameters, ONLY: DP,N,NSTART,NX,NY,NZ,NXT,NYT,NZT,NX1,NY1,NZ1, &
                        MYIDX,MYIDY,MYIDZ,IERR,MYID, &
                        PI,TIME,I_END,ICOLL,ROT_SPEED_TURBINE,RADIUS_TURBINE
  USE field_shared, ONLY: X,Y,Z,XI,YI,ZI,U,V,W,FX,FY,FZ
  USE tools

  IMPLICIT NONE

  INTEGER:: NUM,NUMB,I_WT_SCH,IBEM,NS,NS_TOW,I_WT_NAC,I_WT_TOW,I_WT_TSR, &
            I_WT_OUT,I_WT_OUT_FORCE,N_OUT_FORCE
  REAL(KIND=DP):: UP_FACTOR
  REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: WT_ANGLE
 
  CONTAINS
!=============================================================================!
!                            WIND TURBINE MODELLING                           !
!=============================================================================!
    SUBROUTINE TURBINE_WRAP(DX,DY,DZ,DT)

    IMPLICIT NONE
    INTEGER :: I,J,K,DUMI,STAT
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: U_G,V_G,W_G
    REAL(KIND=DP) :: TIME_WT_START,DX,DY,DZ,DT,SUMM
    LOGICAL :: FILE_EXIST
   
    OPEN(1,FILE="turbine.in")
    READ(1,*)TIME_WT_START
    READ(1,*)NUM
    READ(1,*)NUMB
    READ(1,*)I_WT_SCH
    READ(1,*)IBEM
    READ(1,*)NS
    READ(1,*)I_WT_NAC
    READ(1,*)I_WT_TOW
    READ(1,*)NS_TOW
    READ(1,*)I_WT_TSR
    READ(1,*)UP_FACTOR
    READ(1,*)I_WT_OUT
    READ(1,*)I_WT_OUT_FORCE
    READ(1,*)N_OUT_FORCE
    CLOSE(1)

!---DETERMINE IF TIME > TIME_WT_START. IF NO, THEN DELETE .echo FILE AND RETURN
    IF(TIME.LE.TIME_WT_START)THEN
      INQUIRE(FILE="turbine_angle.echo", EXIST=FILE_EXIST)
      IF(FILE_EXIST)THEN      
        OPEN(1, IOSTAT=STAT, FILE="turbine_angle.echo", STATUS='old')    ! Delete existing *.echo file
        IF (stat == 0) CLOSE(1, STATUS='delete')
      END IF
      RETURN
    END IF
!---READ WIND TURBINE ANGLES
    ALLOCATE(WT_ANGLE(NUM))

    INQUIRE(FILE="turbine_angle.echo", EXIST=FILE_EXIST) ! Determine if the file exists

    IF(FILE_EXIST)THEN
      OPEN(1,FILE='turbine_angle.echo')
      DO I=1,NUM
        READ(1,*)DUMI,WT_ANGLE(I)
      END DO
      CLOSE(1)
    ELSE
      DO I=1,NUM
        WT_ANGLE(I)=0.0
      END DO
    END IF

    ALLOCATE(U_G(NXT,NYT,NZT),V_G(NXT,NYT,NZT),W_G(NXT,NYT,NZT))

    CALL ASSEM_ALL(U,NX1,NY1,NZ1,U_G)
    CALL ASSEM_ALL(V,NX1,NY1,NZ1,V_G)
    CALL ASSEM_ALL(W,NX1,NY1,NZ1,W_G)
    IF(I_WT_SCH.EQ.1)THEN   
      CALL HAWT_ALM(U_G,V_G,W_G,1,1,1,DX,DY,DZ,DT)
    END IF
!---SAVE WIND TURBINE ANGLES
    OPEN(1,FILE='turbine_angle.echo')
    DO I=1,NUM
      WRITE(1,*)I,WT_ANGLE(I)
    END DO
    CLOSE(1)

    DEALLOCATE(WT_ANGLE)
    DEALLOCATE(U_G,V_G,W_G)

    END SUBROUTINE
!=============================================================================!
!            WIND TURBINE BLADES SIMULATED BY ACTUATOR LINE MODEL             !
!=============================================================================!
!       INPUT:
!          NX,NY,NZ: NUMBER OF POINTS
!          RA:       RADIUS OF THE TURBINE DISK
!          X,Y,Z:    COORDINATES
!          U:        STREAMWISE VELOCITY
!          DT:       TIME STEP
!       OUTPUT:
!          FX,FY,FZ: BODY FORCE COMPONENT
!          RADIUS_TURBINE:      MAXIMUM RADIUS OF THE ROTOR
!          ROT_SPEED_TURBINE:   MAXIMUM ANGULAR VELOCITY
!          WT_ANGLE: ROTATIONAL ANGLE OF ONE OF THE BLADES FOR EACH TURBINE
        SUBROUTINE HAWT_ALM(U_G,V_G,W_G,SI1,SI2,SI3,DX,DY,DZ,DT)

 	IMPLICIT NONE  
  
        INTEGER:: M,I,J,K,SI1,SI2,SI3
        REAL(KIND=DP):: DX,DY,DZ,DT
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: U_G,V_G,W_G
        REAL(KIND=DP):: RA(NUM),D(NUM),X0(NUM),Y0(NUM),Z0(NUM),SR(NUM)
        REAL(KIND=DP):: RN(NUM),NR(NUM),TOR(NUM)
        REAL(KIND=DP):: U_RO,UA_RO(NUM),U_UP,UA_UP(NUM),CDN(NUM)
        REAL(KIND=DP):: OMEGA(NUM),SFX(NUM),POWER(NUM)
        REAL(KIND=DP):: CT(NUM),CP(NUM)
	REAL(KIND=DP):: ALFA,ANGLE(NUMB),POWERT
	REAL(KIND=DP):: TW,DV,DUM,DIS
        REAL(KIND=DP):: AF,T
        CHARACTER:: DUMC
        INTEGER:: ID(NUM),DUMI
!-------TOWER
        REAL(KIND=DP):: CDT(NUM),TS,FT(NS_TOW),RT(NUM)
        REAL(KIND=DP):: TX,TY(NS_TOW),TZ(NS_TOW),UB,EPS,DY0,DZ0
        LOGICAL :: DIR_E

        RADIUS_TURBINE=0.0
        ROT_SPEED_TURBINE=0.0
!-------READ THE TURBINE INFORMATION
        OPEN(1,FILE='turbine.dat')
        READ(1,*)DUMC,DUMC,DUMC,DUMC,DUMC,DUMC,DUMC,DUMC,DUMC,DUMC,DUMC
        DO I=1,NUM
          READ(1,*)DUMI,RA(I),X0(I),Y0(I),Z0(I),SR(I),RN(I),RT(I),CDN(I),CDT(I)
          RADIUS_TURBINE=MAX(RADIUS_TURBINE,RA(I))  !  GET MAXIMUM ROTOR RADIUS
        END DO
        CLOSE(1)

        D=RA*2  
!-------NACELLE
        IF(I_WT_NAC.EQ.1)THEN
          CALL NACELLE(X0,Y0,Z0,RN,DX,CDN)
        END IF
!-------ROTOR 
!       GET ROTATIONAL SPEED OF THE TURBINES (RAD/S)
        CALL ROTATION_SPEED(X0,Y0,Z0,SR,RA,OMEGA)

        DO M=1,NUM
          CALL ROTOR_ALM(U_G,V_G,W_G,SI1,SI2,SI3,DX,DY,DZ,OMEGA(M),RA(M),&
                         X0(M),Y0(M),Z0(M),WT_ANGLE(M),RN(M),SFX(M),TOR(M))
        END DO
!-------TOWER
        IF(I_WT_TOW.EQ.1)THEN
          CALL TOWER(U_G,SI1,SI2,SI3,DX,DY,DZ,NS_TOW,D,CDT,RT,X0,Y0,Z0)
        END IF
!-------UPDATE BLADE POSITIONS           
        DO I=1,NUM               
          WT_ANGLE(I)=WT_ANGLE(I)+DT*OMEGA(I)
        END DO
!-------OUTPUT OF SOME PERFORMANCE INFORMATION
!       CALCULATE EXTRACTED POWER AND THRUST AND POWER COEFFICIENTS
        IF(I_WT_OUT.EQ.1)THEN
          DO M=1,NUM
            UA_UP(M)=AVE_WIND(X0(M)-UP_FACTOR*D(M),Y0(M),Z0(M),RA(M))
            UA_RO(M)=AVE_WIND(X0(M),Y0(M),Z0(M),RA(M))
          END DO

          DO M=1,NUM
            POWER(M)=UA_RO(M)*SFX(M)
            CT(M)=SFX(M)/(0.5*PI*RA(M)**2*UA_UP(M)**2)
            CP(M)=POWER(M)/(0.5*PI*RA(M)**2*UA_UP(M)**3)
          END DO

          INQUIRE(FILE='OUTPUT_WT', EXIST=DIR_E)
          IF(MYID.EQ.0)THEN
            IF(DIR_E)THEN
            ELSE
              CALL system('mkdir OUTPUT_WT')
            END IF  

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_THRUST.OUT')
            ELSE  
              OPEN(1,FILE='OUTPUT_WT/WT_THRUST.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(SFX(M),M=1,NUM)
            CLOSE(1) 
  
            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_TORQUE.OUT')
            ELSE    
              OPEN(1,FILE='OUTPUT_WT/WT_TORQUE.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(TOR(M),M=1,NUM)
            CLOSE(1)

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_WIND_SPEED.OUT')
            ELSE
              OPEN(1,FILE='OUTPUT_WT/WT_WIND_SPEED.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(UA_UP(M),M=1,NUM),(UA_RO(M),M=1,NUM)
            CLOSE(1)

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_POWER.OUT')
            ELSE            
              OPEN(1,FILE='OUTPUT_WT/WT_POWER.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(POWER(M)/1000.0D0,M=1,NUM)
            CLOSE(1)

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_ROT_SPEED.OUT')
            ELSE
              OPEN(1,FILE='OUTPUT_WT/WT_ROT_SPEED.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(OMEGA(M),M=1,NUM)
            CLOSE(1)

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_THRUST_COE.OUT')
            ELSE
              OPEN(1,FILE='OUTPUT_WT/WT_THRUST_COE.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(CT(M),M=1,NUM)
            CLOSE(1)

            IF(N.EQ.1)THEN
              OPEN(1,FILE='OUTPUT_WT/WT_POWER_COE.OUT')
            ELSE
              OPEN(1,FILE='OUTPUT_WT/WT_POWER_COE.OUT',ACCESS='APPEND')
            END IF
            WRITE(1,*)TIME,(CP(M),M=1,NUM)
            CLOSE(1)
          END IF
        END IF

	END SUBROUTINE
!--------------------------------------------------------------!
!               DRAG DISK MODEL FOR THE NACELLE PART           !
!--------------------------------------------------------------!
!  INPUT:
!       X0,Y0,Z0: LOCATIONS OF THE TURBINE 
!       RN:       RADIUS OF THE HUB
!       DX:       GRID SPACING IN THE X DIRECTION
!       CDN:      DRAG COEFFICIENT OF THE NACELLE
!  OUTPUT:
!       FX:       IT IS A GLOBAL ARRAY
        SUBROUTINE NACELLE(X0,Y0,Z0,RN,DX,CDN)
        IMPLICIT NONE
        INTEGER :: J,K,M,ID(NUM)
	REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: GAMMA
        REAL(KIND=DP) :: DX,X0(NUM),Y0(NUM),Z0(NUM),RN(NUM),CDN(NUM) 
        REAL(KIND=DP) :: UA_RO
             
        ALLOCATE(GAMMA(NX,NY,NZ))
        GAMMA=0.0
!-------DETERMINE THE X INDEX OF THE ROTOR: ID
        DO M=1,NUM
          IF(ICOLL.EQ.0)THEN
            ID(NUM)=INDEX_X(X0(M),X)
          ELSE
            ID(NUM)=INDEX_X(X0(M),XI)
          END IF
        END DO   
!-------DRAG DISK MODEL     
        DO M=1,NUM
          UA_RO=AVE_WIND(X0(M),Y0(M),Z0(M),RN(M))         
          IF(ID(NUM).NE.0)THEN
            CALL CAL_GAMMA(ID(NUM),Y0(M),Z0(M),RN(M),GAMMA)
          END IF
  	  DO J=1,NY
	    DO K=1,NZ
              IF(ID(NUM).NE.0)THEN
                FX(ID(NUM),J,K)=FX(ID(NUM),J,K)+ &
               UA_RO**2/2.0*CDN(M)/DX*GAMMA(ID(NUM),J,K)
              END IF
	    END DO
	  END DO
        END DO 
        DEALLOCATE(GAMMA)

        END SUBROUTINE        
!--------------------------------------------------------------!
!                 ROTATIONAL SPEED OF THE ROTOR                !
!--------------------------------------------------------------!
!  INPUT:
!       X0,Y0,Z0:  LOCATIONS OF THE TURBINE 
!       SR:        TIP SPEED RATIO, IF I_WT_TSR=1; RMS, IF I_WT_TSR=0
!       RA:        RADIUS OF THE ROTOR
!  OUTPUT:
!       OMEGA:     ROTATIONAL SPEED OF THE ROTOR (RAD)
!       ROT_SPEED_TURBINE:    MAXIMUM ROTATIONAL SPEED AMONG ALL TURBINES
        SUBROUTINE ROTATION_SPEED(X0,Y0,Z0,SR,RA,OMEGA)
        IMPLICIT NONE
        INTEGER :: M  
        REAL(KIND=DP) :: X0(NUM),Y0(NUM),Z0(NUM)   
        REAL(KIND=DP) :: SR(NUM),RA(NUM),OMEGA(NUM)
        REAL(KIND=DP) :: UA_UP

        ROT_SPEED_TURBINE=0.0
        DO M=1,NUM
          IF(I_WT_TSR.EQ.1)THEN          
            UA_UP=AVE_WIND(X0(M)-UP_FACTOR*RA(M)*2.0,Y0(M),Z0(M),RA(M))
            OMEGA(M)=SR(M)*UA_UP/RA(M) ! FIX THE TIP SPEED RATIO
          ELSE
            OMEGA(M)=SR(M)*2.0*PI/60.0    !FIX THE ANGULAR VELOCITY
          END IF
          ROT_SPEED_TURBINE=MAX(ROT_SPEED_TURBINE,ABS(OMEGA(M)))
        END DO
 
        END SUBROUTINE
!--------------------------------------------------------------!
!                 ALM + BEM MODEL FOR THE ROTOR                !
!--------------------------------------------------------------!
!  INPUT:
!       DX,DY,DZ:  GRIP SPACING 
!       OMEGA:     ROTATIONAL SPEED OF THE ROTOR (RAD)
!       RA:        RADIUS OF THE ROTOR
!       X0,Y0,Z0:  LOCATIONS OF THE TURBINE 
!       ANGLES0:   ROTATIONAL ANGLE OF ONE OF THE BLADES
!       RN:        RADIUS OF THE HUB
!  OUTPUT:
!       FX:        IT IS A GLOBAL ARRAY
!       SFX:       TOTAL THRUST FORCE
!       TOR:       TOTAL TORQUE
!  NOTE, THE ROTOR ROTATES WITH X AXIS IN THE COUNTER CLOCKWISE DIRECTION
!  THE PITCH ANGLE IS THE ANGLE BETWEEN THE BLADE SECTION (AIRFOIL) AND THE ROTOR PLANE
!  THE ANGLE_WIND IS THE ANGLE BETWEEN THE INCIDENT WIND AND THE ROTOR PLANE
!  THE ANGLE OF ATTACK IS ANGLE_WIND MINUS PITCH ANGLE
        SUBROUTINE ROTOR_ALM(U_G,V_G,W_G,SI1,SI2,SI3,DX,DY,DZ,OMEGA,RA, &
                             X0,Y0,Z0,ANGLE0,RN,SFX,TOR)

        IMPLICIT NONE
        INTEGER :: I,J,K,M
        INTEGER :: ND,IO,NSEC,SI1,SI2,SI3
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: U_G,V_G,W_G
        REAL(KIND=DP),DIMENSION(:), ALLOCATABLE :: RADIUS,PITCH,CHORD
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: BFX,BFY,BFZ
        REAL(KIND=DP):: DX,DY,DZ,RA,X0,Y0,Z0,SR,ANGLE0,RN
        REAL(KIND=DP):: OMEGA,RHOA,R,DIS,ALFA
        REAL(KIND=DP):: XREF,XS(NS,NUMB),YS(NS,NUMB),ZS(NS,NUMB)
        REAL(KIND=DP):: UB,VB,WB,XLOC,YLOC,ZLOC
        REAL(KIND=DP):: C(NS),PA(NS)
        REAL(KIND=DP):: VREL(2)
        REAL(KIND=DP):: DS,DSS,N1,N2,N3
        REAL(KIND=DP):: CL,CD,L(NS,NUMB),D(NS,NUMB),F2D(2)
        REAL(KIND=DP):: F2DX(NS,NUMB),F2DY(NS,NUMB),F2DZ(NS,NUMB)
        REAL(KIND=DP):: ANGLE(NUMB),ANGLE_WIND,ATTACK
        REAL(KIND=DP):: EPS,IANGLE,T,TOR
        REAL(KIND=DP):: SFX,SFY,SFZ,SUM,SUMT
        CHARACTER*20 :: DUMC
        CHARACTER*10,DIMENSION(:),ALLOCATABLE:: FOILTYPE
        CHARACTER*10,DIMENSION(NS):: FOIL

!-------READ TURBINE BLADE INFORMATION------------------
        OPEN(1,FILE="turbine_blade.dat")
        READ(1,*)
        NSEC = 0
        DO 
          READ(1,*,IOSTAT=IO) 
          IF (IO /= 0) EXIT
          NSEC = NSEC + 1
        END DO

        REWIND(1)
    
        ALLOCATE(RADIUS(NSEC),PITCH(NSEC),CHORD(NSEC),FOILTYPE(NSEC))
        READ(1,*)
        DO I=1,NSEC
          READ(1,*)RADIUS(I),PITCH(I),CHORD(I),FOILTYPE(I)
        END DO
        CLOSE(1)

        DS=(RA-RN)/(NS-1)

        DO I=1,NS
          R=DS*(I-1)+RN
          IF(R.LE.RADIUS(1))THEN
            PA(I)=PITCH(1)
            C(I)=CHORD(1)
            FOIL(I)=FOILTYPE(I)
          ELSE IF(R.GE.RADIUS(NSEC))THEN  
            PA(I)=PITCH(NSEC)
            C(I)=CHORD(NSEC)
            FOIL(I)=FOILTYPE(NSEC)
          ELSE   
            J=1      
            DO 
              IF(R.GT.RADIUS(J).AND.R.LE.RADIUS(J+1))THEN
                PA(I)=PITCH(J+1)
                C(I)=CHORD(J+1)
                FOIL(I)=FOILTYPE(J+1)
                EXIT
              END IF
              J=J+1
            END DO
          END IF   
        END DO
        DEALLOCATE(RADIUS,PITCH,CHORD,FOILTYPE)

        DO I=1,NS
          PA(I)=PA(I)*PI/180.0
        END DO
!-------ROTATIONAL ANGLES OF THE BLADES
        DO I=1,NUMB
          ANGLE(I)=ANGLE0+(I-1)*360.0/NUMB*PI/180.0
        END DO
!-------POSITIONS AT BLADE SECTIONS
        DO I=1,NUMB
          XS(1,I)=X0
          YS(1,I)=Y0+RN*SIN(ANGLE(I))
          ZS(1,I)=Z0+RN*COS(ANGLE(I))
        END DO 

        DO J=1,NUMB
          DO I=2,NS
            XS(I,J)=X0
            YS(I,J)=YS(I-1,J)+DS*SIN(ANGLE(J))
            ZS(I,J)=ZS(I-1,J)+DS*COS(ANGLE(J))
          END DO
        END DO
!-------GET FORCING TERMS AT EACH ELEMENT--------------------------
        DO J=1,NUMB
          DO I=1,NS
            F2DX(I,J)=0.0
            F2DY(I,J)=0.0
            F2DZ(I,J)=0.0
          END DO
        END DO

        TOR=0.0          
  
        DO J=1,NUMB
          DO I=1,NS
            R=DMIN1(RA,RN+DS*(I-1))
            UB=0.0
            VB=0.0
            WB=0.0
!-----------DIRECTLY CALCULATE AXIAL AND TANGENTIAL VELOCITIES-------------------
            IF(IBEM.EQ.0)THEN
              IF(ICOLL.EQ.0)THEN
                XREF=DMAX1(X0-RA,X(1))
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,X,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
                XREF=DMAX1(X0-RA,XI(1))
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,XI,Y,ZI,1,V_G,SI1,SI2,SI3,VB)
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,XI,Y,ZI,1,W_G,SI1,SI2,SI3,WB)
              ELSE
                XREF=DMAX1(X0-RA,XI(1))
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,XI,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,XI,YI,ZI,1,V_G,SI1,SI2,SI3,VB)
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J),NXT,NYT,NZT,XI,YI,ZI,1,W_G,SI1,SI2,SI3,WB)
              END IF   

              N1=0.0                       ! UNIT TANGENTIAL VECTOR (N1,N2,N3)
              N2=-DSIN(ANGLE(J))
              N3= DCOS(ANGLE(J))
              VREL(1)=UB                   ! VELOCITY IN THE AXIAL DIRECTION
              VREL(2)=VB*N2+WB*N3-OMEGA*R  ! VELOCITY IN THE TANGENTIAL DIRECTION   
!-----------USE BEM METHOD TO CALCULATE AXIAL AND TANGENTIAL VELOCITIES------------
            ELSE
              IF(ICOLL.EQ.0)THEN
                XREF=DMAX1(X0-RA*2.0*UP_FACTOR,X(1))
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J), &
                                  NXT,NYT,NZT,X,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
              ELSE
                XREF=DMAX1(X0-RA*2.0*UP_FACTOR,XI(1))
                CALL INTER_GLOBAL(XREF,YS(I,J),ZS(I,J), &
                                  NXT,NYT,NZT,XI,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
              END IF  
             
              CALL BEM(UB,R,RA,C(I),PA(I),OMEGA,VREL(1),VREL(2),FOIL(I))     

            END IF 
!-----------CALCULATE FORCES-------------------------------------------------------
!       ANGLE BETWEEN APPARENT WIND AND ROTOR PLANE   
            IF(VREL(2).GT.0.0)THEN
              ANGLE_WIND=PI+ATAN(VREL(1)/(-VREL(2)))
            ELSE
              ANGLE_WIND=ATAN(VREL(1)/(-VREL(2)))
            END IF 

            ATTACK=ANGLE_WIND-PA(I)               ! ANGLE OF ATTACK

            CALL LIFT_DRAG(ATTACK,FOIL(I),CL,CD)  !  LIFT AND DRAG COEFFS FROM TABULATED DATA     
!       CALCULATE THE LIFT AND DRAG
            L(I,J)=(VREL(1)**2+VREL(2)**2)*C(I)*CL/2.0
            D(I,J)=(VREL(1)**2+VREL(2)**2)*C(I)*CD/2.0
!       TRANSFER THE FORCE TO X-Y-Z COORDINATE
            F2D(1)=L(I,J)*DCOS(ANGLE_WIND)+D(I,J)*DSIN(ANGLE_WIND)
            F2D(2)=L(I,J)*DSIN(ANGLE_WIND)-D(I,J)*DCOS(ANGLE_WIND)
            F2DX(I,J)=F2D(1)
            F2DY(I,J)=-F2D(2)*DCOS(ANGLE(J))
            F2DZ(I,J)= F2D(2)*DSIN(ANGLE(J))
!       CALCUALTE TORQUE
            TOR=TOR+F2D(2)*DS*R
          END DO
        END DO  
!-------DISTRIBUTE THE FORCE TO GRID POINTS
!       PARAMETER OF CONCENTRATION: EPS
        EPS=(DX*DZ*DY)**(1.0/3.0)

        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              IF(ICOLL.EQ.0)THEN
                XLOC=X (I+MYIDX*NX)
                YLOC=Y (J+MYIDY*NY)
                ZLOC=Z (K+MYIDZ*NZ)
              ELSE
                XLOC=XI(I+MYIDX*NX)
                YLOC=YI(J+MYIDY*NY)
                ZLOC=ZI(K+MYIDZ*NZ)
              END IF

              DIS=SQRT((XLOC-X0)**2+ &
                       (YI(J+MYIDY*NY)-Y0)**2+ &
                       (ZI(K+MYIDZ*NZ)-Z0)**2)
              IF(DIS.LT.RA*1.2)THEN
                DO M=1,NUMB
                  FX(I,J,K)=FX(I,J,K)-  &
                            GAUSSD(XLOC,YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ),DS, &
                                   XS(1:,M),YS(1:,M),ZS(1:,M),F2DX(1:,M),NS,EPS)
                END DO
              END IF
             
              DIS=SQRT((XI(I+MYIDX*NX)-X0)**2+ &
                       (YLOC-Y0)**2+ &
                       (ZI(K+MYIDZ*NZ)-Z0)**2)
              IF(DIS.LT.RA*1.2)THEN
                DO M=1,NUMB
                  FY(I,J,K)=FY(I,J,K)-  &
                            GAUSSD(XI(I+MYIDX*NX),YLOC,ZI(K+MYIDZ*NZ),DS, &
                                   XS(1:,M),YS(1:,M),ZS(1:,M),F2DY(1:,M),NS,EPS)
                END DO
              END IF

              DIS=SQRT((XI(I+MYIDX*NX)-X0)**2+ &
                       (YI(J+MYIDY*NY)-Y0)**2+ &
                       (ZLOC-Z0)**2)
              IF(DIS.LT.RA*1.2)THEN
                DO M=1,NUMB
                  FZ(I,J,K)=FZ(I,J,K)-  &
                            GAUSSD(XI(I+MYIDX*NX),YI(J+MYIDY*NY),ZLOC,DS, &
                                   XS(1:,M),YS(1:,M),ZS(1:,M),F2DZ(1:,M),NS,EPS)
                END DO
              END IF
            END DO
          END DO
        END DO
!-------GET TOTAL THRUST FORCE
        SFX=0.0
        SFY=0.0
        SFZ=0.0
        DO J=1,NUMB
          DO I=1,NS
            IF(I.EQ.1.OR.I.EQ.NS)THEN
              DSS=DS/2.0
            ELSE
              DSS=DS
            END IF
            SFX=SFX+F2DX(I,J)*DSS
            SFY=SFY+F2DY(I,J)*DSS
            SFZ=SFZ+F2DZ(I,J)*DSS
          END DO
        END DO

!       SUM=0.0
!       DO I=1,NX
!         DO J=1,NY
!           DO K=1,NZ
!             SUM=SUM+FX(I,J,K)*DX*DY*DZ
!           END DO
!         END DO
!       END DO
!       CALL MPI_ALLREDUCE(SUM,SUMT,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                          MPI_COMM_WORLD,IERR)
!       IF(MYID.EQ.0)THEN
!         PRINT*,SFX,SUMT
!       END IF
!-------EXPORT LIFT AND DRAG FORCES ON THE BLADES
        IF(I_WT_OUT_FORCE.EQ.1)THEN
          IF(MOD(N,N_OUT_FORCE).EQ.0)THEN
            CALL OUTPUT_WT_FORCE(L,D,RN,RA)
          END IF
        END IF

      END SUBROUTINE
!--------------------------------------------------------------!
!                   LINE MODEL FOR THE TOWER                   !
!--------------------------------------------------------------!
!  INPUT:
!       DX,DY,DZ:  GRIP SPACING 
!       NS_TOW:    NUMBER OF SECTIONS IN EACH TOWER
!       D:         DIAMETER OF THE ROTOR
!       CDT:       DRAG COEFFICIENT OF THE TOWER
!       RT:        RADIUS OF THE TOWER
!       X0,Y0,Z0:  LOCATIONS OF THE TURBINE 
!  OUTPUT:
!       FX:        IT IS A GLOBAL ARRAY
      SUBROUTINE TOWER(U_G,SI1,SI2,SI3,DX,DY,DZ,NS_TOW,D,CDT,RT,X0,Y0,Z0)

      IMPLICIT NONE

      INTEGER :: I,J,K,M,N,NS_TOW,SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: U_G
      REAL(KIND=DP) :: X0(NUM),Y0(NUM),Z0(NUM)      
      REAL(KIND=DP) :: D(NUM),CDT(NUM),RT(NUM)
      REAL(KIND=DP) :: TX(NS_TOW),TY(NS_TOW),TZ(NS_TOW),FT(NS_TOW)
      REAL(KIND=DP) :: EPS,DX,DY,DZ,TS,UB,XLOC

      EPS=SQRT(DY**2+DZ**2)
      DO M=1,NUM
!-----CALCULATE LINE FORCES ATTACHED ON THE TOWER SECTIONS
        TS=Z0(M)/NS_TOW

        TX=X0(M)+5.0
        TY=Y0(M)
        TZ(1)=TS/2
        DO N=2,NS_TOW
          TZ(N)=TZ(N-1)+TS
        END DO


        DO N=1,NS_TOW
          IF(ICOLL.EQ.0)THEN
            CALL INTER_GLOBAL(TX(N)-UP_FACTOR*D(M),TY(N),TZ(N), &
                              NXT,NYT,NZT,X,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
          ELSE
            CALL INTER_GLOBAL(TX(N)-UP_FACTOR*D(M),TY(N),TZ(N), &
                              NXT,NYT,NZT,XI,YI,ZI,1,U_G,SI1,SI2,SI3,UB)
          END IF
          FT(N)=-UB**2*CDT(M)*RT(M)    ! DRAG FORCE PER UNIT LENGTH
        END DO
 !----DISTRIBUTE THE LINE FORCES TO NEIGHBORING GRID POINTS USING GUSSIAN FUNCTION
        DO I=1,NX
	  DO J=1,NY
	    DO K=1,NZ
              IF(ICOLL.EQ.0)THEN
                XLOC=X(I+MYIDX*NX)
              ELSE
                XLOC=XI(I+MYIDX*NX)
              END IF

	      FX(I,J,K)=FX(I,J,K)-  &
                        GAUSSD(XLOC,YI(J+MYIDY*NY),ZI(K+MYIDZ*NZ),TS, &
                               TX,TY,TZ,FT,NS_TOW,EPS)
	    END DO
	  END DO
  	END DO
      END DO

      END SUBROUTINE
!--------------------------------------------------------------!
!                    BLADE ELEMENT MOMENTUM METHOD             !
!--------------------------------------------------------------!
!     Grant Ingram (2011)
      SUBROUTINE BEM(UB,R,RA,C,PA0,OMEGA,VX,VSITA,FOIL)

      IMPLICIT NONE
      INTEGER :: COUNT,NT,I
      REAL(KIND=DP):: UB,R,RA,OMEGA,PA0,VX,VSITA,PA
      REAL(KIND=DP):: SR,BETA,BETA0,A,AP,B,C,SIGMA,Q
      REAL(KIND=DP):: CL,CD,IANGLE
      REAL(KIND=DP):: RHS,ERROR,ZERO
      CHARACTER*10 :: FOIL

      DATA B /3.0/
      DATA ERROR,NT/1.E-4,200/
      DATA ZERO /1.E-8/

!-----GET AEROFOIL INLET ANGLE: PA
      PA=PI/2.0-PA0
!-----GET THE LOCAL SPEED RATIO: SR
      SR=OMEGA*R/UB
!-----LOCAL SOLIDITY: SIGMA
      SIGMA=B*C/(R*2.0*PI)
!-----INITIAL GUESS OF RELATIVE FLOW ANGLE: BETA       
      BETA=PI/2.0-ATAN(1.0/SR)*2.0/3.0
      COUNT=0
!-----INTERPOLATE LIFT AND DRAG COEFFICIENTS: CL,CD
2     IANGLE=PA-BETA   !THIS IS THE INCIDENT ANGLE
      CALL LIFT_DRAG(IANGLE,FOIL,CL,CD)                   
!-----GET A AND AP         
      IF(COUNT.EQ.0)THEN   
!       INITIAL VALUE OF AXIAL INDUCTION FACTOR: A
!       INITIAL VALUE OF ANGULAR INDUCTION FACTOR: AP    
        A=(1.0+4.0*DCOS(BETA)**2/(SIGMA*CL*DSIN(BETA)))**(-1)       
        AP=(1.0-3.0*A)/(4.0*A-1.0)
      ELSE
! THE LOSS CORRECTION FACTOR DUE TO TIP VORTEX: Q 
! Q varies from 0 to 1 where the magnitude determines the reduction in forces along the blade
        Q=2.0/PI*ACOS(EXP(-(B/2.0*(1.0-R/RA)/DMAX1(ZERO,(R/RA*DCOS(BETA))))))+ZERO
!       UPDAT A & AP
        RHS=SIGMA*(CL*DSIN(BETA)+CD*DCOS(BETA))/(4.0*Q*DCOS(BETA)**2)
        A=RHS/(RHS+1.0)
        RHS=SIGMA*(CL*DCOS(BETA)-CD*DSIN(BETA))/(4.0*Q*SR*DCOS(BETA)**2)
        AP=RHS*(1.0-A)
      END IF
!-----UPDATE THE RELATIVE FLOW ANGLE: BETA
      BETA0=BETA
      BETA=ATAN(SR*(1+AP)/(1.0-A))
      BETA=(BETA*3.0+BETA0)/4.0
!-----CHECK CONVERGENCY
      IF(ABS(BETA-BETA0).GT.ERROR.AND.COUNT.LT.NT)THEN
        COUNT=COUNT+1
        GOTO 2
      ELSE IF(ABS(BETA-BETA0).GT.ERROR.AND.COUNT.EQ.NT)THEN
        IF(MYID.EQ.0)THEN
          PRINT*,'BEM DOES NOT CONVERGE.'
        END IF
      END IF
!-----GET AXIAL AND TANGENTIAL VELOCITIES: VX & VSITA
 10   CONTINUE

      VX=UB*(1.0-A)
      VSITA=OMEGA*R*(1.0+AP)

      END SUBROUTINE  
!*********************************************************************C
!             GAUSSIAN DISTRIBUTION OF A LINE FUNCTION FB 	      C
!*********************************************************************C
!     INPUT: 
!     X,Y,Z: POSITION OF THE POINT
!     XS,YS,ZS: COORDINATES OF ELEMENTS
!     FB: FUNCTION AT THE ELEMENTS
!     NUMBER: NUMBER OF THE ELEMENTS
!     EPS: PARAMETER OF CONCENTRATION
!     OUTPUT:
!     F: DISTRIBUTED FUNCTION VALUE AT THE POINT  
      REAL(KIND=DP) FUNCTION GAUSSD(X0,Y0,Z0,DS,XS,YS,ZS,FB,NUMBER,EPS)

      IMPLICIT NONE
      INTEGER :: M,N,NUMBER
      REAL(KIND=DP):: X0,Y0,Z0
      REAL(KIND=DP),DIMENSION(:):: XS,YS,ZS,FB
      REAL(KIND=DP):: DIS,ETA,EPS,DS

      GAUSSD=0.0

      DO M=1,NUMBER
        DIS=SQRT((X0-XS(M))**2+(Y0-YS(M))**2+(Z0-ZS(M))**2)
        ETA=EXP(-(DIS/EPS)**2)/(EPS**3*PI**1.5)
        GAUSSD=GAUSSD+FB(M)*ETA*DS
      END DO

      END FUNCTION
!***********************************************************************!
!               CALCULATE THE OVERLAPED AREA FACTOR: GAMMA              !
!***********************************************************************!
        SUBROUTINE CAL_GAMMA(ID,Y0,Z0,R,GAMMA)

	IMPLICIT NONE

	INTEGER :: J,K,ID
	REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE:: GAMMA
	REAL(KIND=DP):: X0,Y0,Z0,R
	REAL(KIND=DP):: D1,D2,D3,D4
	REAL(KIND=DP):: DZ0,DY0,RY,RZ,RX1,RX2,RY1,RY2,RZ1,RZ2,A

	DO J=1,NY
	  DO K=1,NZ
            A=0.0
	    DZ0=Z(K+MYIDZ*NZ+1)-Z(K+MYIDZ*NZ)
	    DY0=Y(J+MYIDY*NY+1)-Y(J+MYIDY*NY)
	    D1=SQRT((Y(J+MYIDY*NY  )-Y0)**2+(Z(K+MYIDZ*NZ  )-Z0)**2)-R
	    D2=SQRT((Y(J+MYIDY*NY  )-Y0)**2+(Z(K+MYIDZ*NZ+1)-Z0)**2)-R
	    D3=SQRT((Y(J+MYIDY*NY+1)-Y0)**2+(Z(K+MYIDZ*NZ  )-Z0)**2)-R
	    D4=SQRT((Y(J+MYIDY*NY+1)-Y0)**2+(Z(K+MYIDZ*NZ+1)-Z0)**2)-R   
!    FOR CELL INSIDE THE DISK
 	    IF(D1.LE.0.0.AND.D2.LE.0.0.AND.D3.LE.0.0.AND.D4.LE.0.0)THEN
              A=DZ0*DY0
!    FOR CELL WITH THREE POINTS INSIDE THE DISK
            ELSE IF(D1.GE.0.0.AND.D2.LE.0.0.AND.D3.LE.0.0.AND.D4.LE.0.0)THEN
	      RY=D1*DY0/(ABS(D3)+D1)
	      RZ=D1*DZ0/(ABS(D2)+D1)
	      A=DZ0*DY0-RZ*RY/2.0
	    ELSE IF(D1.LE.0.0.AND.D2.GE.0.0.AND.D3.LE.0.0.AND.D4.LE.0.0)THEN
	      RY=D2*DY0/(ABS(D4)+D2)
	      RZ=D2*DZ0/(ABS(D1)+D2)
	      A=DZ0*DY0-RZ*RY/2.0
	    ELSE IF(D1.LE.0.0.AND.D2.LE.0.0.AND.D3.GE.0.0.AND.D4.LE.0.0)THEN
	      RY=D3*DY0/(ABS(D1)+D3)
	      RZ=D3*DZ0/(ABS(D4)+D3)
	      A=DZ0*DY0-RZ*RY/2.0
	    ELSE IF(D1.LE.0.0.AND.D2.LE.0.0.AND.D3.LE.0.0.AND.D4.GE.0.0)THEN
	      RY=D4*DY0/(ABS(D2)+D4)
	      RZ=D4*DZ0/(ABS(D3)+D4)
	      A=DZ0*DY0-RZ*RY/2.0
!    FOR CELL WITH TWO POINTS INSIDE THE DISK
	    ELSE IF(D1.GE.0.0.AND.D2.GE.0.0.AND.D3.LE.0.0.AND.D4.LE.0.0)THEN
	      RY1=D1*DY0/(ABS(D3)+D1)
	      RY2=D2*DY0/(ABS(D4)+D2)
	      A=DZ0*DY0-(RY1+RY2)*DZ0/2.0
	    ELSE IF(D1.GE.0.0.AND.D2.LE.0.0.AND.D3.GE.0.0.AND.D4.LE.0.0)THEN
	      RZ1=D1*DZ0/(ABS(D2)+D1)
	      RZ2=D3*DZ0/(ABS(D4)+D3)
	      A=DZ0*DY0-(RZ1+RZ2)*DY0/2.0
	    ELSE IF(D1.LE.0.0.AND.D2.GE.0.0.AND.D3.LE.0.0.AND.D4.GE.0.0)THEN
	      RZ1=D2*DZ0/(ABS(D1)+D2)
	      RZ2=D4*DZ0/(ABS(D3)+D4)
	      A=DZ0*DY0-(RZ1+RZ2)*DY0/2.0
	    ELSE IF(D1.LE.0.0.AND.D2.LE.0.0.AND.D3.GE.0.0.AND.D4.GE.0.0)THEN
	      RY1=D3*DY0/(ABS(D1)+D3)
	      RY2=D4*DY0/(ABS(D2)+D4)
	      A=DZ0*DY0-(RY1+RY2)*DZ0/2.0
!    FOR CELL WITH ONE POINTS INSIDE THE DISK
	    ELSE IF(D1.LE.0.0.AND.D2.GE.0.0.AND.D3.GE.0.0.AND.D4.GE.0.0)THEN
	      RY=ABS(D1)*DY0/(D3+ABS(D1))
	      RZ=ABS(D1)*DZ0/(D2+ABS(D1))
	      A=RZ*RY/2.0
	    ELSE IF(D1.GE.0.0.AND.D2.LE.0.0.AND.D3.GE.0.0.AND.D4.GE.0.0)THEN
	      RY=ABS(D2)*DY0/(D4+ABS(D2))
	      RZ=ABS(D2)*DZ0/(D1+ABS(D2))
	      A=RZ*RY/2.0
	    ELSE IF(D1.GE.0.0.AND.D2.GE.0.0.AND.D3.LE.0.0.AND.D4.GE.0.0)THEN
	      RY=ABS(D3)*DY0/(D1+ABS(D3))
	      RZ=ABS(D3)*DZ0/(D4+ABS(D3))
	      A=RZ*RY/2.0
	    ELSE IF(D1.GE.0.0.AND.D2.GE.0.0.AND.D3.GE.0.0.AND.D4.LE.0.0)THEN
	      RY=ABS(D4)*DY0/(D2+ABS(D4))
	      RZ=ABS(D4)*DZ0/(D3+ABS(D4))
	      A=RZ*RY/2.0
	    END IF
	    GAMMA(ID+1,J,K)=A/(DZ0*DY0)
	  END DO
	END DO
      END SUBROUTINE 
!--------------------------------------------------------------!
!               AVERAGED WIND SPEED OVER ROTOR DISK            !
!--------------------------------------------------------------!
      REAL(KIND=DP) FUNCTION AVE_WIND(XC,YC,ZC,RAD)

      IMPLICIT NONE
      REAL(KIND=DP):: XC,YC,ZC,RAD
      REAL(KIND=DP):: U_AVE,DIS
      INTEGER:: COUNTI,COUNT,I,J,K,ID
      
      ID=0    

      IF(ICOLL.EQ.0)THEN
        XC=DMIN1(DMAX1(XC,X(1)),X(NXT+1))
        ID=INDEX_X(XC,X)
      ELSE
        XC=DMIN1(DMAX1(XC,XI(1)),XI(NXT+1))
        ID=INDEX_X(XC,XI)
      END IF 

      U_AVE=0.0
      COUNTI=0
      IF(ID.NE.0)THEN
        DO J=1,NY
          DO K=1,NZ
            DIS=SQRT((YI(J+MYIDY*NY)-YC)**2+(ZI(K+MYIDZ*NZ)-ZC)**2) 
            IF(DIS.LT.RAD)THEN
              U_AVE=U_AVE+U(ID,J,K)
              COUNTI=COUNTI+1
            END IF     
          END DO        
        END DO
      END IF 

      CALL MPI_ALLREDUCE(U_AVE,AVE_WIND,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(COUNTI,COUNT,1,MPI_INTEGER,MPI_SUM, &
                         MPI_COMM_WORLD,IERR)

      AVE_WIND=AVE_WIND/(COUNT+1.0E-12)
 
      END FUNCTION
!--------------------------------------------------------------!
!         GET LIFT AND DRAG FROM TABULATED AIRFOIL DATA        !
!--------------------------------------------------------------!
!  INPUT: 
!    ATTACK: ATTACK ANGLE
!    FOIL:   AIRFOIL TYPE
!  OUTPUT:
!    CL, CD: LIFT AND DRAG COEFFICIENT
      SUBROUTINE LIFT_DRAG(ATTACK,FOIL,CL,CD)

      IMPLICIT NONE

      INTEGER :: ND,I,N,IO
      REAL(KIND=DP) :: ATTACK,CL,CD
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE:: ALFA_T,CL_T,CD_T
      CHARACTER*10 :: FOIL,STR1
      CHARACTER*20 :: FILENAME

!-----READ LIFT AND DRAG COEFFICIENTS FROM TABULATED FILE
      STR1='.dat'
      FILENAME=TRIM(FOIL) // TRIM(STR1)

      OPEN(1,FILE=FILENAME)
      READ(1,*)
      ND = 0
      DO 
        READ(1,*,IOSTAT=IO) 
        IF (IO /= 0) EXIT
        ND = ND + 1
      END DO

      REWIND(1)
    
      ALLOCATE(ALFA_T(ND),CL_T(ND),CD_T(ND))
      READ(1,*)
      DO I=1,ND
        READ(1,*)ALFA_T(I),CL_T(I),CD_T(I)
      END DO
      CLOSE(1)
      DO I=1,ND
        ALFA_T(I)=ALFA_T(I)*PI/180.0
      END DO
!-----INTERPOLATE LIFT AND DRAG COEFFICIENTS: CL,CD
      DO I=1,ND-1
        IF(ATTACK.GE.ALFA_T(I).AND.ATTACK.LT.ALFA_T(I+1))THEN
          CL=(CL_T(I+1)-CL_T(I))/(ALFA_T(I+1)-ALFA_T(I))*  &
             (ATTACK-ALFA_T(I))+CL_T(I)
          CD=(CD_T(I+1)-CD_T(I))/(ALFA_T(I+1)-ALFA_T(I))*  &
             (ATTACK-ALFA_T(I))+CD_T(I)
          EXIT
        END IF
      END DO
      IF(ATTACK.LT.ALFA_T(1))THEN
        CL=0.45
        CD=1.3
      ELSE IF(ATTACK.GE.ALFA_T(ND))THEN
        CL=0.0
        CD=1.3
      END IF
   
      DEALLOCATE(ALFA_T,CL_T,CD_T)

      END SUBROUTINE
!=========================================================================!
!           EXPORT LIFT AND DRAG FORCES ON THE WIND TURBINE BLADES        !
!=========================================================================!
      SUBROUTINE OUTPUT_WT_FORCE(LIFT,DRAG,RN,RA)
      IMPLICIT NONE
      INTEGER :: I,J
      REAL(KIND=DP) :: UB(NS,3),LIFT(NS,3),DRAG(NS,3),R(NS)
      REAL(KIND=DP) :: RN,RA,DS
      CHARACTER*10 :: DIR,SR1,SR2,SR3,SR4
      CHARACTER*30 :: FILENAME
      LOGICAL :: DIR_E
      
      DS=(RA-RN)/(NS-1)
      DO I=1,NS
        R(I)=MIN(RA,RN+DS*(I-1))/RA
      END DO

      DIR="OUTPUT_WT/"
      SR1="WT_"
      WRITE(SR2,'(I6.6)')N
      SR3="_NUM"
      WRITE(SR4,'(I1)')NUM
      FILENAME=TRIM(DIR)//TRIM(SR1)//TRIM(SR2)//TRIM(SR3)//TRIM(SR4)
      
      INQUIRE(FILE='OUTPUT_WT', EXIST=DIR_E)

      IF(MYID.EQ.0)THEN
        IF(DIR_E)THEN
        ELSE
          CALL system('mkdir OUTPUT_WT')
        END IF

        OPEN(UNIT=1,FILE=FILENAME)
        WRITE(1,*)'VARIABLES="R/R0","LIFT","DRAG"'
        DO J=1,NUMB
          WRITE(1,*)'ZONE T="',J,'"'
          DO I=1,NS
            WRITE(1,*)R(I),LIFT(I,J),DRAG(I,J)
          END DO
        END DO
        CLOSE(1)
      END IF

      END SUBROUTINE    
  END MODULE turbine
