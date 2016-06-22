! This module contains the Poisson solvers
!
  MODULE poisson

  USE mpi
  USE parameters, ONLY: N,DP,NX,NY,NZ,NXT,NYT,NZT,NX1,NY1,NZ1,NPX,NPY,NPZ,MYID,MYIDX,MYIDY,MYIDZ,  &
                        ICOLL,ISCHEME,ISCHE_POI,ORDER_POI,NMUL_POI,NITE_POI,TOLE_POI,BC,SCREEN_LEVEL
  USE field_shared, ONLY: XI,YI,ZI,UF,VF,WF,PD,CHECK,CHECK0
  USE boundary
  USE tools

  IMPLICIT NONE
! INCLUDE "fftw3.f"
  
  CONTAINS
!--------------------------------------------------------------!
!            Define and solve the poisson equation             !
!--------------------------------------------------------------!
!   ALFA: THE SCALING FACTOR OF THE TIME STEP
    SUBROUTINE PRESSURE_WRAP(DX,DY,DZ,DT)

    IMPLICIT NONE
    REAL(KIND=DP) :: DX,DY,DZ,DT
    REAL(KIND=DP), DIMENSION (:,:,:), ALLOCATABLE:: SIG
    INTEGER :: I,J,K,ID_COL,SI

    ALLOCATE(SIG(NX,NY,NZ))

    IF(ISCHEME.EQ.1)THEN
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            SIG(I,J,K)=(DERIV_X(UF,NX1,NY1,NZ1,I+1,J,K,1,ORDER_POI,DX)+ &
                        DERIV_Y(VF,NX1,NY1,NZ1,I,J+1,K,1,ORDER_POI,DY)+ &                          
                        DERIV_Z(WF,NX1,NY1,NZ1,I,J,K+1,1,ORDER_POI,DZ))/DT
!              SIG(I,J,K)=((UF(I+1,J,K)-UF(I,J,K))/DX*9.0/8.0-      &
!                          (UF(I+2,J,K)-UF(I-1,J,K))/(DX*3.0)/8.0+  &
!                          (VF(I,J+1,K)-VF(I,J,K))/DY*9.0/8.0-      &
!                          (VF(I,J+2,K)-VF(I,J-1,K))/(DY*3.0)/8.0+  &
!                          (WF(I,J,K+1)-WF(I,J,K))/DZ*9.0/8.0-      &
!                          (WF(I,J,K+2)-WF(I,J,K-1))/(DZ*3.0)/8.0)/DT
          END DO
        END DO
      END DO

      IF(ISCHE_POI.EQ.1)THEN   ! MULTIGRID METHOD WITH DYNAMIC SOR SMOOTHER
        CALL MULTIGRID(NX,NY,NZ,NPX,NPY,NPZ, &
                       DX,DY,DZ,NMUL_POI,NITE_POI,TOLE_POI,SIG,PD,NX1,NY1,NZ1, &
                       BC(:,5)) 
      ELSE IF(ISCHE_POI.EQ.2)THEN ! DYNAMIC SOR METHOD
        CALL DSOR(NX,NY,NZ,DX,DY,DZ,BC(:,5),NITE_POI,SIG,PD,NX1,NY1,NZ1,ORDER_POI)
      END IF       
    END IF

    DEALLOCATE(SIG)
!---UPDATE VELOCITY FIELD BY PRESSURE CORRECTION
    ID_COL=ABS(ICOLL-1)  ! COLLOCATED: ICOLL=1, THEN ID_COL=0; STAGGERED: ICOLL=0, THEN ID_COL=1
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          FX(I,J,K)=FX(I,J,K)-DERIV_X(PD,NX1,NY1,NZ1,I,J,K,ID_COL,ORDER_POI,DX)
          FY(I,J,K)=FY(I,J,K)-DERIV_Y(PD,NX1,NY1,NZ1,I,J,K,ID_COL,ORDER_POI,DY)                                  
          FZ(I,J,K)=FZ(I,J,K)-DERIV_Z(PD,NX1,NY1,NZ1,I,J,K,ID_COL,ORDER_POI,DZ)
        END DO
      END DO
    END DO
 
    END SUBROUTINE
!--------------------------------------------------------------!
!          MULTIGRID METHOD FOR SOLVING THE POISSON EQ         !
!--------------------------------------------------------------!
    SUBROUTINE MULTIGRID(IM,JM,KM,NPX,NPY,NPZ,            &
                         DX,DY,DZ,NHT0,NTOTAL,TOLE,F0,P0,SI1,SI2,SI3,  &
                         BC0)    ! OPTIONAL

    IMPLICIT NONE 
    INTEGER :: IM,JM,KM,NPX,NPY,NPZ
    INTEGER :: NH,NHT0,NHT,NB,SI1,SI2,SI3
    INTEGER :: IMC,JMC,KMC,IMF,JMF,KMF,N_SM
    INTEGER :: I,J,K
    INTEGER, OPTIONAL :: BC0(6)
    REAL(KIND=DP):: DX,DY,DZ
    REAL(KIND=DP):: DXH(NHT0),DYH(NHT0),DZH(NHT0)
    REAL(KIND=DP),DIMENSION(:,:,:)::F0
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: P0
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE::R,PI
    REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE::PH,FH         
    REAL(KIND=DP):: SUMR,SUMF,L2,ERROR,ZERO,TOLE
    REAL(KIND=DP):: DX2,DX1,DY2,DY1,DZ2,DZ1,SUMI,SUM,SUM2,SUM3
    REAL(KIND=DP):: VREF,VREFT
    INTEGER :: NI,NTOTAL,IT,I_BC(6),CONVERGE_FLAG

    DATA ZERO /1.0E-8/

    IF(ORDER_POI.EQ.2)THEN
      NB=1
    ELSE IF(ORDER_POI.EQ.4)THEN
      NB=3
    END IF

!---CHECK THE MAXIMUM LEVEL: IF THE CELL NUMBER AT THE HIGHEST LEVEL IS LESS THAN NB, THEN REDUCE THE NUMBER OF LEVELS
    NHT=NHT0
    DO NH=1,NHT0
      IF(MIN(NZ,MIN(NX,NY))/2**(NHT-1).LT.NB)THEN
        NHT=NHT-1
      ELSE
        EXIT
      END IF
    END DO
    NHT=MAX(NHT,1)    
 
    P0=0.0

    CONVERGE_FLAG=0
!---SET BOUNDARY CONDITION
    IF(PRESENT(BC0))THEN
      DO I=1,6
        I_BC(I)=BC0(I)
      END DO
    ELSE
      I_BC=3    ! DEFAULT: NEUMANN BC
    END IF
!------------------------------
    L2=RMS(F0,1,1,1,IM,JM,KM)

    IF(MYID.EQ.0)THEN
      PRINT*,'CONTINUITY ERROR:',L2
    END IF
    IF(L2.LT.TOLE)THEN
      RETURN
    END IF
!---HIERARCHY OF GRID
!   BASE GRID
    DXH(1)=DX
    DYH(1)=DY
    DZH(1)=DZ
!   COARSER GRIDS
    DO NH=2,NHT
      DXH(NH)=DXH(NH-1)*2.0
      DYH(NH)=DYH(NH-1)*2.0
      DZH(NH)=DZH(NH-1)*2.0
    END DO

    ALLOCATE(PH(SI1:IM+NB,SI2:JM+NB,SI3:KM+NB,NHT))
    PH=0.0
    ALLOCATE(FH(IM,JM,KM,NHT))
    FH=0.0
    DO I=1,IM
      DO J=1,JM
        DO K=1,KM
          FH(I,J,K,1)=F0(I,J,K)
        END DO
      END DO
    END DO

    DO NI=1,NTOTAL
!-----FROM FINE TO COARSE LEVELS---------------------------------------
      DO NH=1,NHT-1
        IMF=IM/(2**(NH-1))
        JMF=JM/(2**(NH-1))
        KMF=KM/(2**(NH-1))
        IMC=IM/(2**NH)
        JMC=JM/(2**NH)
        KMC=KM/(2**NH)

        CALL DSOR(IMF,JMF,KMF,DXH(NH),DYH(NH),DZH(NH),I_BC,20, &
                  FH(:,:,:,NH),PH(:,:,:,NH),SI1,SI2,SI3,ORDER_POI)

        ALLOCATE(R(0:IMF+1,0:JMF+1,0:KMF+1))
        DO I=1,IMF
          DO J=1,JMF
            DO K=1,KMF
              R(I,J,K)=FH(I,J,K,NH)-                                 &
                       LAPLACE(DXH(NH),DYH(NH),DZH(NH),IMF,JMF,KMF,  &
                               PH(:,:,:,NH),SI1,SI2,SI3,I,J,K,I_BC,ORDER_POI)
            END DO
          END DO
        END DO
        CALL GET_BC(IMF,JMF,KMF,R,1,1,1,0,I_BC)    

        DO I=1,IMC
          DO J=1,JMC
            DO K=1,KMC
              FH(I,J,K,NH+1)=REST(R,0,0,0,I,J,K) !  GET FH AT COARSER LEVEL
            END DO
           END DO
        END DO

        DEALLOCATE(R)

        DO I=-NB+1,IMC+NB
          DO J=-NB+1,JMC+NB
            DO K=-NB+1,KMC+NB
              PH(I,J,K,NH+1)=0.0
            END DO
          END DO
        END DO
      END DO   
!-----AT THE COARSEST GRID--------------------------------------------
      IMC=IM/(2**(NHT-1))
      JMC=JM/(2**(NHT-1))
      KMC=KM/(2**(NHT-1))

      CALL DSOR(IMC,JMC,KMC,DXH(NHT),DYH(NHT),DZH(NHT),I_BC,100, &
                FH(:,:,:,NHT),PH(:,:,:,NHT),SI1,SI2,SI3,ORDER_POI)
!-----FROM COARSE TO FINE LEVELS--------------------------------------- 
      DO NH=NHT-1,1,-1
        IMC=IM/(2**NH)
        JMC=JM/(2**NH)
        KMC=KM/(2**NH)           
        IMF=IM/(2**(NH-1))
        JMF=JM/(2**(NH-1))
        KMF=KM/(2**(NH-1))         
 
        DO I=1,IMF
          DO J=1,JMF
            DO K=1,KMF
              PH(I,J,K,NH)=PH(I,J,K,NH)+PROL(PH(:,:,:,NH+1),SI1,SI2,SI3,IMF,JMF,KMF,IMC,JMC,KMC,I,J,K)
            END DO
          END DO
        END DO 
        CALL GET_BC(IMF,JMF,KMF,PH(:,:,:,NH),NB,NB,NB,0,I_BC)

        CALL DSOR(IMF,JMF,KMF,DXH(NH),DYH(NH),DZH(NH),I_BC,20,     &
                  FH(:,:,:,NH),PH(:,:,:,NH),SI1,SI2,SI3,ORDER_POI)
      END DO 
!-----CHECK CONVERGENCE  
      ALLOCATE(R(IM,JM,KM))
      DO I=-NB+1,IM+NB
        DO J=-NB+1,JM+NB
          DO K=-NB+1,KM+NB
            P0(I,J,K)=PH(I,J,K,1)
          END DO
        END DO
      END DO
      DO I=1,IM
        DO J=1,JM
          DO K=1,KM
            R(I,J,K)=F0(I,J,K)-LAPLACE(DX,DY,DZ,IM,JM,KM,P0,SI1,SI2,SI3,I,J,K, &
                                       I_BC,ORDER_POI)
          END DO
        END DO
      END DO
      ERROR=RMS(R,1,1,1,IM,JM,KM)

      IF(MYID.EQ.0.AND.SCREEN_LEVEL.EQ.2)THEN
        PRINT*,'MULTIGRID RESIDUAL:',ERROR,' AT STEP ', NI
      END IF

      DEALLOCATE(R)
      IF(ERROR.LT.TOLE)THEN
        CONVERGE_FLAG=1
        EXIT
      END IF
    END DO

    IF(CONVERGE_FLAG.EQ.0)THEN
      IF(MYID.EQ.0)THEN
        PRINT*,'WARNING: mulgrid solver does not converge at error = ',ERROR
      END IF

    END IF
          
    DEALLOCATE(PH,FH)

    END SUBROUTINE 
!*****************************************************************************!
!	     	          DYNAMIC SOR POISSON SOLVER                          !
!*****************************************************************************!   
!   INPUT: INITIAL VALUES OF PI;
!          COEFFICIENTS: C1,C2,C3,C4,C5,C6,C0
!          SOURCE TERM: F
!          BOUNADRY VALUES(DERIVATIVES) BC: BV
!          I_BC=1: DIRICHLET BC; 
!          I_BC=2: NEUMANN BC; 
!          I_BC=3: PERIODIC BC
!   OUTPUT: COMPUTED VALUES OF P
        SUBROUTINE DSOR(IM,JM,KM,DX,DY,DZ,I_BC,NT,F,P,SI1,SI2,SI3,ORDER_LAPLACE)

        IMPLICIT NONE 
  
        INTEGER :: IM,JM,KM,I,J,K,NI,NT,ORDER_LAPLACE,SI1,SI2,SI3
        REAL(KIND=DP),DIMENSION(:,:,:), ALLOCATABLE:: R,CHECK,CHECK0
        REAL(KIND=DP),DIMENSION(:,:,:):: F
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: P
        REAL(KIND=DP):: DX,DY,DZ,RP,W
        REAL(KIND=DP):: VREF,VREFT
	REAL(KIND=DP):: ERROR,TOLE,ERRORP,ERRORPP
        INTEGER :: I_BC(6)

        DATA TOLE /1.D-5/

        ALLOCATE(R(IM,JM,KM))

        ERROR=1.0E3
        ERRORP=1.0E3
        W=1.0

        DO NI=0,NT
          CALL SOR_JACOBI(IM,JM,KM,DX,DY,DZ,F,P,SI1,SI2,SI3, &
                          I_BC,ORDER_LAPLACE,W)
!-------FIX THE ABSOLUTE VALUE OF P AT A POINT
          IF(I_BC(1).NE.1.AND.I_BC(2).NE.1.AND.I_BC(3).NE.1.AND. &
             I_BC(4).NE.1.AND.I_BC(5).NE.1.AND.I_BC(6).NE.1)THEN
            IF(MYID.EQ.0)THEN
              VREF=P(1,1,1)
            ELSE
              VREF=0.0
            END IF
            CALL MPI_ALLREDUCE(VREF,VREFT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR)     

            P=P-VREFT
          END IF             
!-------CHECK THE CONVERGENCE
          ERRORPP=ERRORP
          ERRORP=ERROR
          DO I=1,IM
            DO J=1,JM
              DO K=1,KM
                R(I,J,K)=ABS(F(I,J,K)-LAPLACE(DX,DY,DZ,IM,JM,KM,P,SI1,SI2,SI3,I,J,K, &
                                              I_BC,ORDER_LAPLACE))
              END DO
            END DO
          END DO
          ERROR=RMS(R,1,1,1,IM,JM,KM)          
!--------ADJUST THE SOR COEFFICIENT
          IF(ERROR.GT.ERRORP)THEN
            W=W-0.1
          ELSE 
            IF(ERROR+ERRORPP-ERRORP*2.0.GT.0.0)THEN
              W=W+0.05
            END IF
          END IF     
          W=DMIN1(2.0D0,DMAX1(W,0.2D0))

         IF(MYID.EQ.0.AND.ISCHE_POI.EQ.2.AND.SCREEN_LEVEL.EQ.2)THEN
           PRINT*,'DSOR RESIDUAL=',ERROR,', ITERATION=', NI,', RF=',W
          END IF

          IF(ABS(ERROR).LT.TOLE)THEN
            EXIT
          END IF     
         
        END DO

        DEALLOCATE(R)

        END SUBROUTINE 
!*****************************************************************************!
!		      SOUBROUTINE OF SOR JACOBI ITERATION             	      !
!*****************************************************************************!
!   FOR A POISSON EQUATION AS BELOW
!   CX(-3)*P(I-3,J,K)+CX(-2)*P(I-2,J,K)+CX(-1)*P(I-1,J,K)+
!   CX( 1)*P(I+1,J,K)+CX( 2)*P(I+2,J,K)+CX( 3)*P(I+3,J,K)+
!   CY(-3)*P(I,J-3,K)+CY(-2)*P(I,J-2,K)+CY(-1)*P(I,J-1,K)+
!   CY( 1)*P(I,J+1,K)+CY( 2)*P(I,J+2,K)+CY( 3)*P(I,J+3,K)+
!   CX(-3)*P(I,J,K-3)+CX(-2)*P(I,J,K-2)+CX(-1)*P(I,J,K-1)+
!   CX( 1)*P(I,J,K+1)+CX( 2)*P(I,J,K+2)+CX( 3)*P(I,J,K+3)+
!   C0*P(I,J,K)=F(I,J,K)
!
!   SUCCESSIVE OVER RELAXATION (SOR):
!   PR IS THE INITIAL GUESS;
!   P(I,J,K)=PR(I,J,K)*(1-W)-
!            W/C0*(CX(-3)*PR(I-3,J,K)+CX(-2)*PR(I-2,J,K)+CX(-1)*PR(I-1,J,K)+
!                  CX( 1)*PR(I+1,J,K)+CX( 2)*PR(I+2,J,K)+CX( 3)*PR(I+3,J,K)+
!                  CY(-3)*PR(I,J-3,K)+CY(-2)*PR(I,J-2,K)+CY(-1)*PR(I,J-1,K)+
!                  CY( 1)*PR(I,J+1,K)+CY( 2)*PR(I,J+2,K)+CY( 3)*PR(I,J+3,K)+
!                  CZ(-3)*PR(I,J,K-3)+CZ(-2)*PR(I,J,K-2)+CZ(-1)*PR(I,J,K-1)+
!                  CZ( 1)*PR(I,J,K+1)+CZ( 2)*PR(I,J,K+2)+CZ( 3)*PR(I,J,K+3)-F(I,J,K))
!     W: COE. OF SOR (0 <= W <= 1)
      SUBROUTINE SOR_JACOBI(IM,JM,KM,DX,DY,DZ,F,P,SI1,SI2,SI3, &
                            I_BC,ORDER_LAPLACE,W)
      IMPLICIT NONE
      INTEGER:: IM,JM,KM,I,J,K,M,I_BC(6),NB,ID,ORDER_LAPLACE,SI1,SI2,SI3
      REAL(KIND=DP):: DX,DY,DZ,W
      REAL(KIND=DP):: C0
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE:: CX,CY,CZ
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE:: PR,CHECK,CHECK0
      REAL(KIND=DP), DIMENSION(:,:,:):: F
      REAL(KIND=DP), DIMENSION(SI1:,SI2:,SI3:):: P

      IF(ORDER_LAPLACE.EQ.2)THEN
        NB=1
      ELSE IF(ORDER_LAPLACE.EQ.4)THEN
        NB=3
      END IF

      ALLOCATE(PR(1-NB:IM+NB,1-NB:JM+NB,1-NB:KM+NB))
      ALLOCATE(CX(-NB:NB),CY(-NB:NB),CZ(-NB:NB))

      DO I=-NB+1,IM+NB
        DO J=-NB+1,JM+NB
          DO K=-NB+1,KM+NB
            PR(I,J,K)=P(I,J,K)
          END DO
         END DO
      END DO

      DO J=1,JM
        DO I=1,IM
          DO K=1,KM
            CALL PARA_LAPLACE(CX,CY,CZ,NB,DX,DY,DZ,IM,JM,KM,I,J,K,I_BC)  

            C0=CX(0)+CY(0)+CZ(0)
 
            P(I,J,K)=PR(I,J,K)*(1.0-W)
            DO M=-NB,NB
              IF(M.NE.0)THEN
                P(I,J,K)=P(I,J,K)-W/C0*(CX(M)*PR(I+M,J,K)+ &
                                        CY(M)*PR(I,J+M,K)+ &
                                        CZ(M)*PR(I,J,K+M))
              END IF
            END DO
            P(I,J,K)=P(I,J,K)+F(I,J,K)*W/C0
         END DO
        END DO
      END DO  

      CALL GET_BC(IM,JM,KM,P,NB,NB,NB,0,I_BC)
 
      DEALLOCATE(PR,CX,CY,CZ)
 
      END SUBROUTINE 
!*****************************************************************************!
!		      FUNCTION OF LAPLACE OPERATOR                	      !
!*****************************************************************************!
!  INPUTS:
!     DX,DY,DZ: GRID SPACING
!     I_BC: 
      REAL(KIND=DP) FUNCTION LAPLACE(DX,DY,DZ,IM,JM,KM,P,SI1,SI2,SI3,I,J,K, &
                                     BC_LAPLACE,ORDER_LAPLACE)  ! OPTIONAL
      IMPLICIT NONE
      INTEGER :: NB,IM,JM,KM,I,J,K,ORDER,I_BC(6),M,SI1,SI2,SI3
      INTEGER, OPTIONAL:: ORDER_LAPLACE,BC_LAPLACE(6)
      REAL(KIND=DP) :: DX,DY,DZ
      REAL(KIND=DP), DIMENSION(SI1:,SI2:,SI3:):: P
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE::CX,CY,CZ

      IF(PRESENT(ORDER_LAPLACE))THEN
        ORDER=ORDER_LAPLACE
      ELSE
        ORDER=4     ! SET DEFAULT ORDER
      END IF

      IF(PRESENT(BC_LAPLACE))THEN
        I_BC=BC_LAPLACE
      ELSE
        I_BC=2   ! DEFAULT BC IS NEUMANN
      END IF   

      IF(ORDER.EQ.2)THEN
        NB=1
      ELSE IF(ORDER.EQ.4)THEN
        NB=3
      END IF

      ALLOCATE(CX(-NB:NB),CY(-NB:NB),CZ(-NB:NB))

      CALL PARA_LAPLACE(CX,CY,CZ,NB,DX,DY,DZ,IM,JM,KM,I,J,K,I_BC)  

      LAPLACE=0.0
      DO M=-NB,NB
          LAPLACE=LAPLACE+CX(M)*P(I+M,J,K)+  &
                          CY(M)*P(I,J+M,K)+  &
                          CZ(M)*P(I,J,K+M)
      END DO  

      DEALLOCATE(CX,CY,CZ)

      END FUNCTION
!*****************************************************************************!
!	                Parameters of the Laplace operator                    !
!*****************************************************************************!
        SUBROUTINE PARA_LAPLACE(CX,CY,CZ,NB,DX,DY,DZ,IM,JM,KM,I,J,K, &
                                BC_LAPLACE)  ! OPTIONAL

        IMPLICIT NONE
        INTEGER :: IM,JM,KM,I,J,K,NB,M
        INTEGER :: I_BC(6)
        INTEGER, OPTIONAL:: BC_LAPLACE(6)
        REAL(KIND=DP), DIMENSION(:), ALLOCATABLE:: CX,CY,CZ
        REAL(KIND=DP) :: DX,DY,DZ
!-------INNER POINT-------------------------------------------------
        IF(NB.EQ.3)THEN           ! 4TH-ORDER CENTRAL SCHEME
          CX(-3) = 1.0/(9.0*64.0)/DX**2         !  i - 3
          CX(-2) = -6.0/(64.0)/DX**2            !  i - 2
          CX(-1) = 87.0/(64.0)/DX**2            !  i - 1
          CX(1)  = 87.0/(64.0)/DX**2            !  i + 1
          CX(2)  = -6.0/(64.0)/DX**2            !  i + 2
          CX(3)  = 1.0/(9.0*64.0)/DX**2         !  i + 3
 
          CY(-3) = 1.0/(9.0*64.0)/DY**2         !  j - 3
          CY(-2) = -6.0/(64.0)/DY**2            !  j - 2
          CY(-1) = 87.0/(64.0)/DY**2            !  j - 1
          CY(1)  = 87.0/(64.0)/DY**2            !  j + 1
          CY(2)  = -6.0/(64.0)/DY**2            !  j + 2
          CY(3)  = 1.0/(9.0*64.0)/DY**2         !  j + 3

          CZ(-3) = 1.0/(9.0*64.0)/DZ**2         !  k - 3
          CZ(-2) = -6.0/(64.0)/DZ**2            !  k - 2
          CZ(-1) = 87.0/(64.0)/DZ**2            !  k - 1
          CZ(1)  = 87.0/(64.0)/DZ**2            !  k + 1
          CZ(2)  = -6.0/(64.0)/DZ**2            !  k + 2
          CZ(3)  = 1.0/(9.0*64.0)/DZ**2         !  k + 3
        ELSE IF(NB.EQ.1)THEN      ! 2ND-ORDER CENTRAL SCHEME
          CX(-1) = 1.0/DX**2            !  i - 1
          CX(1)  = 1.0/DX**2            !  i + 1

          CY(-1) = 1.0/DY**2            !  j - 1
          CY(1)  = 1.0/DY**2            !  j + 1

          CZ(-1) = 1.0/DZ**2            !  k - 1
          CZ(1)  = 1.0/DZ**2            !  k + 1 
        END IF
!-------CORRECT BOUNDARY CONDITIONS----------------------------------
        IF(PRESENT(BC_LAPLACE))THEN
          I_BC=BC_LAPLACE
        ELSE
          I_BC=2   ! DEFAULT BC IS NEUMANN
        END IF

        IF(NB.EQ.3)THEN  ! 4TH-ORDER SCHEME 
!---------INFLOW BC
          IF(I_BC(1).EQ.2)THEN
            IF(MYIDX.EQ.0.AND.I.EQ.1)THEN
              CX(-3) = 0.0
              CX(-2) = 0.0
              CX(-1) = 0.0
              CX(1)  = 81.0/64.0/DX**2
              CX(2)  = -53.0/(64.0*9.0)/DX**2
            END IF

            IF(MYIDX.EQ.0.AND.I.EQ.2)THEN
              CX(-3) = 0.0
              CX(-2) = 0.0
              CX(-1) = 81.0/64.0/DX**2
            END IF

            IF(MYIDX.EQ.0.AND.I.EQ.3)THEN
              CX(-3) = 0.0
              CX(-2) = -53.0/(64.0*9.0)/DX**2
            END IF
          END IF
!---------OUTFLOW BC
          IF(I_BC(2).EQ.2)THEN
            IF(MYIDX.EQ.NPX-1.AND.I.EQ.IM)THEN
              CX(-2) = -53.0/(64.0*9.0)/DX**2
              CX(-1) = 81.0/64.0/DX**2
              CX(1)  = 0.0
              CX(2)  = 0.0
              CX(3)  = 0.0
            END IF

            IF(MYIDX.EQ.NPX-1.AND.I.EQ.IM-1)THEN
              CX(1)  = 81.0/64.0/DX**2
              CX(2)  = 0.0
              CX(3)  = 0.0
            END IF

            IF(MYIDX.EQ.NPX-1.AND.I.EQ.IM-2)THEN
              CX(2)  = -53.0/(64.0*9.0)/DX**2
              CX(3)  = 0.0
            END IF
          END IF
!---------LOWER BC
          IF(I_BC(3).EQ.2)THEN
            IF(MYIDY.EQ.0.AND.J.EQ.1)THEN
              CY(-3) = 0.0
              CY(-2) = 0.0
              CY(-1) = 0.0
              CY(1)  = 81.0/64.0/DY**2
              CY(2)  = -53.0/(64.0*9.0)/DY**2
            END IF

            IF(MYIDY.EQ.0.AND.J.EQ.2)THEN
              CY(-3) = 0.0
              CY(-2) = 0.0
              CY(-1) = 81.0/64.0/DY**2
            END IF

            IF(MYIDY.EQ.0.AND.J.EQ.3)THEN
              CY(-3) = 0.0
              CY(-2) = -53.0/(64.0*9.0)/DY**2
            END IF
          END IF
!---------UPPER BC
          IF(I_BC(4).EQ.2)THEN
            IF(MYIDY.EQ.NPY-1.AND.J.EQ.JM)THEN
              CY(-2) = -53.0/(64.0*9.0)/DY**2
              CY(-1) = 81.0/64.0/DY**2
              CY(1)  = 0.0
              CY(2)  = 0.0
              CY(3)  = 0.0
            END IF

            IF(MYIDY.EQ.NPY-1.AND.J.EQ.JM-1)THEN
              CY(1)  = 81.0/64.0/DY**2
              CY(2)  = 0.0
              CY(3)  = 0.0
            END IF

            IF(MYIDY.EQ.NPY-1.AND.J.EQ.JM-2)THEN
              CY(2)  = -53.0/(64.0*9.0)/DY**2
              CY(3)  = 0.0
            END IF
          END IF
!---------LEFT BC
          IF(I_BC(5).EQ.2)THEN
            IF(MYIDZ.EQ.0.AND.K.EQ.1)THEN
              CZ(-3) = 0.0
              CZ(-2) = 0.0
              CZ(-1) = 0.0
              CZ(1)  = 81.0/64.0/DZ**2
              CZ(2)  = -53.0/(64.0*9.0)/DZ**2
            END IF

            IF(MYIDZ.EQ.0.AND.K.EQ.2)THEN
              CZ(-3) = 0.0
              CZ(-2) = 0.0
              CZ(-1) = 81.0/64.0/DZ**2
            END IF

            IF(MYIDZ.EQ.0.AND.K.EQ.3)THEN
              CZ(-3) = 0.0
              CZ(-2) = -53.0/(64.0*9.0)/DZ**2
            END IF
          END IF
!---------RIGHT BC
          IF(I_BC(6).EQ.2)THEN
            IF(MYIDZ.EQ.NPZ-1.AND.K.EQ.KM)THEN
              CZ(-2) = -53.0/(64.0*9.0)/DZ**2
              CZ(-1) = 81.0/64.0/DZ**2
              CZ(1)  = 0.0
              CZ(2)  = 0.0
              CZ(3)  = 0.0
            END IF

            IF(MYIDZ.EQ.NPZ-1.AND.K.EQ.KM-1)THEN
              CZ(1)  = 81.0/64.0/DZ**2
              CZ(2)  = 0.0
              CZ(3)  = 0.0
            END IF

            IF(MYIDZ.EQ.NPZ-1.AND.K.EQ.KM-2)THEN
              CZ(2)  = -53.0/(64.0*9.0)/DZ**2
              CZ(3)  = 0.0
            END IF
          END IF
        ELSE IF(NB.EQ.1)THEN   ! 2ND-ORDER SCHEME   
!---------INFLOW BC
          IF(I_BC(1).EQ.2)THEN
            IF(MYIDX.EQ.0.AND.I.EQ.1)THEN
              CX(-1) = 0.0
            END IF
          END IF
!---------OUTFLOW BC
          IF(I_BC(2).EQ.2)THEN
            IF(MYIDX.EQ.NPX-1.AND.I.EQ.IM)THEN
              CX(1)  = 0.0
            END IF
          END IF
!---------LOWER BC
          IF(I_BC(3).EQ.2)THEN
            IF(MYIDY.EQ.0.AND.J.EQ.1)THEN
              CY(-1) = 0.0
            END IF
          END IF
!---------UPPER BC
          IF(I_BC(4).EQ.2)THEN
            IF(MYIDY.EQ.NPY-1.AND.J.EQ.JM)THEN
              CY(1)  = 0.0
            END IF
          END IF
!---------LEFT BC
          IF(I_BC(5).EQ.2)THEN
            IF(MYIDZ.EQ.0.AND.K.EQ.1)THEN
              CZ(-1) = 0.0
            END IF
          END IF
!---------RIGHT BC
          IF(I_BC(6).EQ.2)THEN
            IF(MYIDZ.EQ.NPZ-1.AND.K.EQ.KM)THEN
              CZ(1)  = 0.0
            END IF
          END IF
        END IF

        CX(0)=0.0
        CY(0)=0.0
        CZ(0)=0.0
        DO M=-NB,NB
          IF(M.NE.0)THEN
            CX(0)=CX(0)-CX(M)
            CY(0)=CY(0)-CY(M)
            CZ(0)=CZ(0)-CZ(M)
          END IF
        END DO

        END SUBROUTINE
!*****************************************************************************!
!		PROLONGATION FROM COASE LEVEL TO FINE LEVEL     	      !
!*****************************************************************************!
!       COARSE LEVEL: PC
!       FINE LEVEL: PF
!       THE ASPECT RATIO IS 2
        REAL(KIND=DP) FUNCTION PROL(PC,SI1,SI2,SI3,NXF,NYF,NZF,NXC,NYC,NZC,I,J,K)

        IMPLICIT NONE
        INTEGER :: I,J,K,M,N,L,NXF,NYF,NZF,NXC,NYC,NZC
        INTEGER :: M1,M2,N1,N2,L1,L2,SI1,SI2,SI3
	REAL(KIND=DP):: LX(2),LY(2),LZ(2)
        REAL(KIND=DP), DIMENSION(SI1:,SI2:,SI3:) :: PC
        REAL(KIND=DP) :: DF,DC
        REAL(KIND=DP) :: XF(NXF),YF(NYF),ZF(NZF),XC(0:NXC+1),YC(0:NYC+1),ZC(0:NZC+1)

        DF=1.0
        DC=2.0

        XF(1)=DF/2.0
        DO M=2,NXF
          XF(M)=XF(M-1)+DF
        END DO
        YF(1)=DF/2.0
        DO M=2,NYF
          YF(M)=YF(M-1)+DF
        END DO
        ZF(1)=DF/2.0
        DO M=2,NZF
          ZF(M)=ZF(M-1)+DF
        END DO

        XC(0)=-DC/2.0
        DO M=1,NXC+1
          XC(M)=XC(M-1)+DC
        END DO

        YC(0)=-DC/2.0
        DO M=1,NYC+1
          YC(M)=YC(M-1)+DC
        END DO

        ZC(0)=-DC/2.0
        DO M=1,NZC+1
          ZC(M)=ZC(M-1)+DC
        END DO           

        DO M=1,2
          LX(M)=1.
	  DO N=1,2
	    IF(M.NE.N)THEN
	      LX(M)=LX(M)*(XF(I)-XC(I/2+N-1))/(XC(I/2+M-1)-XC(I/2+N-1))
	    END IF
	  END DO
        END DO

        DO M=1,2
          LY(M)=1.
          DO N=1,2
            IF(M.NE.N)THEN
              LY(M)=LY(M)*(YF(J)-YC(J/2+N-1))/(YC(J/2+M-1)-YC(J/2+N-1))
            END IF
          END DO
        END DO

        DO M=1,2
          LZ(M)=1.
          DO N=1,2
            IF(M.NE.N)THEN
              LZ(M)=LZ(M)*(ZF(K)-ZC(K/2+N-1))/(ZC(K/2+M-1)-ZC(K/2+N-1))
            END IF
          END DO
        END DO

        PROL=0.0
        DO M=1,2
          DO N=1,2
            DO L=1,2
              PROL=PROL+LX(M)*LY(N)*LZ(L)*PC(I/2+M-1,J/2+N-1,K/2+L-1)
            END DO
          END DO
        END DO

        END FUNCTION
!*****************************************************************************!
!		RESTRICTION FROM FINE LEVEL TO COARSE LEVEL     	      !
!*****************************************************************************!
!     SI: START INDEX
!     FINE LEVEL: PF
!     COARSE LEVEL: PC
!     THE ASPECT RATIO IS 2
      REAL(KIND=DP) FUNCTION REST(PF,SI1,SI2,SI3,I,J,K)
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,N,SI1,SI2,SI3
      REAL(KIND=DP), DIMENSION(SI1:,SI2:,SI3:):: PF
      REAL(KIND=DP):: CX,CY,CZ   

      REST=0.0
      DO L=-2,1
        DO M=-2,1
          DO N=-2,1
            IF(L.EQ.-2.OR.L.EQ.1)THEN
              CX=1.0/8.0
            ELSE
              CX=3.0/8.0
            END IF
            IF(M.EQ.-2.OR.M.EQ.1)THEN
              CY=1.0/8.0
            ELSE
              CY=3.0/8.0
            END IF
            IF(N.EQ.-2.OR.N.EQ.1)THEN
              CZ=1.0/8.0
            ELSE
              CZ=3.0/8.0
            END IF
            REST=REST+CX*CY*CZ*PF(I*2+L,J*2+M,K*2+N)
          END DO
        END DO
      END DO      

      END FUNCTION

  END MODULE
