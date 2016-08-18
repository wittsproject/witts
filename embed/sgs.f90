! This module contains the subgrid-scale (SGS) models
!
  MODULE SGS

  USE mpi
  USE parameters
  USE class_shared   
  USE tools
  USE boundary

  CONTAINS

!=========================================================================!
!                    LAGRANGIAN SCALE DEPENDENT SGS MODEL                 !
!=========================================================================!
!   SCALE-DEPENDENT VERSION (Elie Bou-Zeid et al., 2005) 
!   INPUT: DX,DY,DZ
!   OUTPUT: NU,CS2,DISSIP (DEFINED IN THE MODULE FIELD_SHARED)
    SUBROUTINE SGS_LASD(DX,DY,DZ)
    IMPLICIT NONE
!    INCLUDE "mpif.h"
    INTEGER :: I,J,K,N1,N2,IFILTER
    INTEGER :: I_OUT_SGS,I_AVE_SGS,I_START_AVE,N_TAVE_SKIP,N_FLD_OUT,N_PRF_OUT
    INTEGER :: I_BC(6),NB
    REAL(KIND=DP),DIMENSION(-NB:NB,-NB:NB,-NB:NB):: U_BAR,V_BAR,W_BAR,U_HAT,V_HAT,W_HAT,      &
                                L11,L22,L33,L12,L13,L23,Q11,Q22,Q33,Q12,Q13,Q23, &
                                S11_BAR,S22_BAR,S33_BAR,S12_BAR,S13_BAR,S23_BAR, &
                                S11_HAT,S22_HAT,S33_HAT,S12_HAT,S13_HAT,S23_HAT, &
                                S_BAR,S_HAT,M11,M22,M33,M12,M13,M23,N11,N22,N33,N12,N13,N23, &
                                MM,LM,NN,QN,VART1,VART2
    REAL(KIND=DP)::DX,DY,DZ
    REAL(KIND=DP)::CS20,CS2_2,CS2_4,U0,U1,U2,V0,V1,V2,W0,W1,W2,X1,X2,Y1,Y2,Z1,Z2,DEL,CONST
    REAL(KIND=DP)::UC,VC,WC,XP,YP,ZP,PLMP,PMMP,PQNP,PNNP,EPSI,DUM,FDX,FDY,FDZ,T,TRACE
    REAL(KIND=DP)::ZERO,TF1,TF2,TF1_2,TF2_2,POWCOEFF,LAGRAN_DT
    REAL(KIND=DP)::NU_SA,CS2_SA,BETA_SA,DISSIP_SA
    CHARACTER:: DUMC


    DATA ZERO /1.E-16/

    NB=2
!---SET THE BC FOR THE EDDY VISOCOSITY. THE ASSUMED NEUMANN BC IS FOR TEMPORAL USE
    DO I=1,6
      I_BC(I)=2
    END DO
  
    CS20=CS0**2

    TF1=2.0   ! test filter 1 / grid size
    TF2=4.0	  ! test filter 2 / grid size 
    TF1_2=TF1**2
    TF2_2=TF2**2

    POWCOEFF = -1.0e0/8.0e0   ! cut off coefficient for Beta=cs(2Delta)/Cs(Delta)

    LAGRAN_DT=DT*N_SGS_SKIP      ! cs update every fifth time step,i.e. cs_count=5
!---CALCULATE THE RATE OF STRAIN TENSOR (Sij)
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
 
        DX=CELL_FV(M)%CELL_DX
        DY=CELL_FV(M)%CELL_DY
        DZ=CELL_FV(M)%CELL_DZ

        CALL CELL_TO_STRUCT(M,NB,1,U_BAR)
        CALL CELL_TO_STRUCT(M,NB,2,V_BAR)
        CALL CELL_TO_STRUCT(M,NB,3,W_BAR)

        U_HAT=U_BAR
        V_HAT=V_BAR
        W_HAT=W_BAR

        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB      
              L11(I,J,K)=U_BAR(I,J,K)*U_BAR(I,J,K)  ! temporary variable
              L12(I,J,K)=U_BAR(I,J,K)*V_BAR(I,J,K)  ! temporary variable
              L13(I,J,K)=U_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
              L23(I,J,K)=V_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
              L22(I,J,K)=V_BAR(I,J,K)*V_BAR(I,J,K)	! temporary variable
              L33(I,J,K)=W_BAR(I,J,K)*W_BAR(I,J,K)	! temporary variable

              Q11(I,J,K)=U_BAR(I,J,K)*U_BAR(I,J,K)  ! temporary variable
              Q12(I,J,K)=U_BAR(I,J,K)*V_BAR(I,J,K)  ! temporary variable
              Q13(I,J,K)=U_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
              Q23(I,J,K)=V_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
              Q22(I,J,K)=V_BAR(I,J,K)*V_BAR(I,J,K)	! temporary variable
              Q33(I,J,K)=W_BAR(I,J,K)*W_BAR(I,J,K)	! temporary variable
            END DO
          END DO
        END DO
!  Here we filter the variable (e.g. u_bar at a scale =tf1*gridsize)
        CALL FILTER(DX,DY,DZ,U_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,V_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,W_BAR,-NB,-NB,-NB,TF1)
! computing Lij as defined in Bou-Zeid et al 2005            
        CALL FILTER(DX,DY,DZ,L11,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,L12,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,L13,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,L22,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,L23,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,L33,-NB,-NB,-NB,TF1)


        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB
              L11(I,J,K)=L11(I,J,K)-U_BAR(I,J,K)*U_BAR(I,J,K)
              L12(I,J,K)=L12(I,J,K)-U_BAR(I,J,K)*V_BAR(I,J,K)
              L13(I,J,K)=L13(I,J,K)-U_BAR(I,J,K)*W_BAR(I,J,K)
              L22(I,J,K)=L22(I,J,K)-V_BAR(I,J,K)*V_BAR(I,J,K)
              L23(I,J,K)=L23(I,J,K)-V_BAR(I,J,K)*W_BAR(I,J,K)
              L33(I,J,K)=L33(I,J,K)-W_BAR(I,J,K)*W_BAR(I,J,K)
!       Deviatoric part of the stress
              TRACE=(L11(I,J,K)+L22(I,J,K)+L33(I,J,K))/3.0
              L11(I,J,K)=L11(I,J,K)-TRACE
              L22(I,J,K)=L22(I,J,K)-TRACE
              L33(I,J,K)=L33(I,J,K)-TRACE
            END DO
          END DO
        END DO
!  Here we test filter the variable (e.g. u_bar at a scale =tf2*gridsize)
        CALL FILTER(DX,DY,DZ,U_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,V_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,W_HAT,-NB,-NB,-NB,TF2)
! computing Qij as defined in Bou-Zeid et al 2005
        CALL FILTER(DX,DY,DZ,Q11,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,Q12,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,Q13,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,Q22,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,Q23,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,Q33,-NB,-NB,-NB,TF2)

        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB
              Q11(I,J,K)=Q11(I,J,K)-U_HAT(I,J,K)*U_HAT(I,J,K)
              Q12(I,J,K)=Q12(I,J,K)-U_HAT(I,J,K)*V_HAT(I,J,K)
              Q13(I,J,K)=Q13(I,J,K)-U_HAT(I,J,K)*W_HAT(I,J,K)
              Q22(I,J,K)=Q22(I,J,K)-V_HAT(I,J,K)*V_HAT(I,J,K)
              Q23(I,J,K)=Q23(I,J,K)-V_HAT(I,J,K)*W_HAT(I,J,K)
              Q33(I,J,K)=Q33(I,J,K)-W_HAT(I,J,K)*W_HAT(I,J,K)

              TRACE=(Q11(I,J,K)+Q22(I,J,K)+Q33(I,J,K))/3.0
              Q11(I,J,K)=Q11(I,J,K)-TRACE
              Q22(I,J,K)=Q22(I,J,K)-TRACE
              Q33(I,J,K)=Q33(I,J,K)-TRACE
            END DO
          END DO
        END DO

        CALL CELL_TO_STRUCT(M,NB,20,S11_BAR)
        CALL CELL_TO_STRUCT(M,NB,21,S22_BAR)
        CALL CELL_TO_STRUCT(M,NB,22,S33_BAR)
        CALL CELL_TO_STRUCT(M,NB,20,S12_BAR)
        CALL CELL_TO_STRUCT(M,NB,21,S13_BAR)
        CALL CELL_TO_STRUCT(M,NB,22,S23_BAR)

        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB
              S_BAR(I,J,K)=SQRT(2.*(S11_BAR(I,J,K)**2+S22_BAR(I,J,K)**2+ &
                                    S33_BAR(I,J,K)**2+                   &
                                2.*(S12_BAR(I,J,K)**2+S13_BAR(I,J,K)**2+ &
                                    S23_BAR(I,J,K)**2)))
            END DO
          END DO
        END DO

        S11_HAT=S11_BAR
        S12_HAT=S12_BAR
        S13_HAT=S13_BAR
        S22_HAT=S22_BAR
        S23_HAT=S23_BAR
        S33_HAT=S33_BAR
        S_HAT=S_BAR

        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB
              M11(I,J,K)=S_BAR(I,J,K)*S11_BAR(I,J,K)
              M12(I,J,K)=S_BAR(I,J,K)*S12_BAR(I,J,K)
              M13(I,J,K)=S_BAR(I,J,K)*S13_BAR(I,J,K)
              M22(I,J,K)=S_BAR(I,J,K)*S22_BAR(I,J,K)
              M23(I,J,K)=S_BAR(I,J,K)*S23_BAR(I,J,K)
              M33(I,J,K)=S_BAR(I,J,K)*S33_BAR(I,J,K)

              N11(I,J,K)=M11(I,J,K)
              N12(I,J,K)=M12(I,J,K)
              N13(I,J,K)=M13(I,J,K)
              N22(I,J,K)=M22(I,J,K)
              N23(I,J,K)=M23(I,J,K)
              N33(I,J,K)=M33(I,J,K)
            END DO
          END DO
        END DO
! filtered strain rate tensors 
! _bar refers to the first test filter = 2*gridsize
! _hat refers to the second test filter = 4*gridsize
        CALL FILTER(DX,DY,DZ,S11_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S12_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S13_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S22_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S23_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S33_BAR,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,S11_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S12_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S13_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S22_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S23_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S33_HAT,-NB,-NB,-NB,TF2)

        DO I=-NB,NB
          DO J=-NB,NB
            DO K=-NB,NB
              S_BAR(I,J,K)=SQRT(2.*(S11_BAR(I,J,K)**2+S22_BAR(I,J,K)**2+ &
                                    S33_BAR(I,J,K)**2+                   &
                                2.*(S12_BAR(I,J,K)**2+S13_BAR(I,J,K)**2+ &
                                    S23_BAR(I,J,K)**2)))
              S_HAT(I,J,K)=SQRT(2.*(S11_HAT(I,J,K)**2+S22_HAT(I,J,K)**2+ &
                                    S33_HAT(I,J,K)**2+                   &
                                2.*(S12_HAT(I,J,K)**2+S13_HAT(I,J,K)**2+ &
                                   S23_HAT(I,J,K)**2)))
            END DO
          END DO
        END DO

        CALL FILTER(DX,DY,DZ,M11,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,M12,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,M13,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,M22,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,M23,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,M33,-NB,-NB,-NB,TF1)
        CALL FILTER(DX,DY,DZ,N11,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,N12,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,N13,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,N22,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,N23,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,N33,-NB,-NB,-NB,TF2)
! computing Mij and Nij as defined in Bou-Zeid et al 2005
        DEL=SQRT(DX*DY*DZ)
        CONST=2.0*(DEL**2) 
   	
        M11(0,0,0)=CONST*(M11(0,0,0)-TF1_2*S_BAR(0,0,0)*S11_BAR(0,0,0))
        M12(0,0,0)=CONST*(M12(0,0,0)-TF1_2*S_BAR(0,0,0)*S12_BAR(0,0,0))
        M13(0,0,0)=CONST*(M13(0,0,0)-TF1_2*S_BAR(0,0,0)*S13_BAR(0,0,0))
        M22(0,0,0)=CONST*(M22(0,0,0)-TF1_2*S_BAR(0,0,0)*S22_BAR(0,0,0))
        M23(0,0,0)=CONST*(M23(0,0,0)-TF1_2*S_BAR(0,0,0)*S23_BAR(0,0,0))
        M33(0,0,0)=CONST*(M33(0,0,0)-TF1_2*S_BAR(0,0,0)*S33_BAR(0,0,0))
!       Deviatoric part of the stress
        TRACE=(M11(0,0,0)+M22(0,0,0)+M33(0,0,0))/3.0
        M11(0,0,0)=M11(0,0,0)-TRACE
        M22(0,0,0)=M22(0,0,0)-TRACE
        M33(0,0,0)=M33(0,0,0)-TRACE
	 
        N11(0,0,0)=CONST*(N11(0,0,0)-TF2_2*S_HAT(0,0,0)*S11_HAT(0,0,0))
        N12(0,0,0)=CONST*(N12(0,0,0)-TF2_2*S_HAT(0,0,0)*S12_HAT(0,0,0))
        N13(0,0,0)=CONST*(N13(0,0,0)-TF2_2*S_HAT(0,0,0)*S13_HAT(0,0,0))
        N22(0,0,0)=CONST*(N22(0,0,0)-TF2_2*S_HAT(0,0,0)*S22_HAT(0,0,0))
        N23(0,0,0)=CONST*(N23(0,0,0)-TF2_2*S_HAT(0,0,0)*S23_HAT(0,0,0))
        N33(0,0,0)=CONST*(N33(0,0,0)-TF2_2*S_HAT(0,0,0)*S33_HAT(0,0,0))
!       Deviatoric part of the stress
        TRACE=(N11(0,0,0)+N22(0,0,0)+N33(0,0,0))/3.0
        N11(0,0,0)=N11(0,0,0)-TRACE
        N22(0,0,0)=N22(0,0,0)-TRACE
        N33(0,0,0)=N33(0,0,0)-TRACE
! Contracting the tensors
        LM(0,0,0)=L11(0,0,0)*M11(0,0,0)+L22(0,0,0)*M22(0,0,0)+L33(0,0,0)*M33(0,0,0)+   &
              2.*(L12(0,0,0)*M12(0,0,0)+L13(0,0,0)*M13(0,0,0)+L23(0,0,0)*M23(0,0,0)) 
        MM(0,0,0)=M11(0,0,0)**2+M22(0,0,0)**2+M33(0,0,0)**2+   &
              2.*(M12(0,0,0)**2+M13(0,0,0)**2+M23(0,0,0)**2)   

        QN(0,0,0)=Q11(0,0,0)*N11(0,0,0)+Q22(0,0,0)*N22(0,0,0)+Q33(0,0,0)*N33(0,0,0)+   &
              2.*(Q12(0,0,0)*N12(0,0,0)+Q13(0,0,0)*N13(0,0,0)+Q23(0,0,0)*N23(0,0,0))
        NN(0,0,0)=N11(0,0,0)**2+N22(0,0,0)**2+N33(0,0,0)**2+   &
              2.*(N12(0,0,0)**2+N13(0,0,0)**2+N23(0,0,0)**2)
! for first test filter
! the lines below check if this is the first time step at which we compute P (P is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ratio = cs^2 =CS2
        IF(LAG_START.EQ.1)THEN
          PMM=DMAX1(ZERO,MM(0,0,0))
          PLM=PMM*CS20
        ELSE   ! CALCULATE PLM & PMM USING LAGRANGIAN AVERAGING 
          DEL=(DX*DY*DZ)**(1.0/3.0)
          T=1.5*DEL*ABS(PLM*PMM)**(-1.0/8.0)
          T= DMAX1(1D-24,T)	    ! clip to avoid numerical problem is tend to zero
          EPSI=LAGRAN_DT/T/(LAGRAN_DT/T+1.0)
          XP=CELL_FV(M)%CELL_X-CELL_FV(M)%CELL_VAR(1)*LAGRAN_DT
          YP=CELL_FV(M)%CELL_Y-CELL_FV(M)%CELL_VAR(2)*LAGRAN_DT
          ZP=CELL_FV(M)%CELL_Z-CELL_FV(M)%CELL_VAR(3)*LAGRAN_DT
          XP=DMIN1(DMAX1(XP,XI(1)),XI(NXT))
          YP=DMIN1(DMAX1(YP,YI(1)),YI(NYT))
          ZP=DMIN1(DMAX1(ZP,ZI(1)),ZI(NZT))
          CALL INTER_CELL(XP,YP,ZP,29,PLMP)
          CALL INTER_CELL(XP,YP,ZP,30,PMMP)
          PLM=DMAX1(EPSI*LM(0,0,0)+(1.0D0-EPSI)*PLMP,ZERO)
          PMM=EPSI*MM(0,0,0)+(1.0e0-EPSI)*PMMP  
        END IF
! for second test filter
! the lines below check if this is the first time step at which we compute F (F is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ration = cs^2 =CS20
        IF(LAG_START.EQ.1)THEN
          PNN=DMAX1(ZERO,NN(0,0,0))
          PQN=PNN*CS20               
        ELSE  ! CALCULATE PQN & PNN USING LAGRANGIAN AVERAGING    
          DEL=(DX*DY*DZ)**(1.0/3.0)
          T=1.5d0*DEL*ABS(PQN*PNN)**(-1.0e0/8.0e0)
          T=DMAX1(1D-24,T)	    ! clip to avoid numerical problem is tend to zero
          EPSI=LAGRAN_DT/T/(LAGRAN_DT/T+1.0)
          XP=CELL_FV(M)%CELL_X-CELL_FV(M)%CELL_VAR(1)*LAGRAN_DT
          YP=CELL_FV(M)%CELL_Y-CELL_FV(M)%CELL_VAR(2)*LAGRAN_DT
          ZP=CELL_FV(M)%CELL_Z-CELL_FV(M)%CELL_VAR(3)*LAGRAN_DT
          XP=DMIN1(DMAX1(XP,XI(1)),XI(NXT))
          YP=DMIN1(DMAX1(YP,YI(1)),YI(NYT))
          ZP=DMIN1(DMAX1(ZP,ZI(1)),ZI(NZT))
          CALL INTER_CELL(XP,YP,ZP,31,PQNP)
          CALL INTER_CELL(XP,YP,ZP,32,PNNP)
          PQN=DMAX1(EPSI*QN(0,0,0)+(1.0D0-EPSI)*PQNP,ZERO)
          PNN=EPSI*NN(0,0,0)+(1.0e0-EPSI)*PNNP
        END IF
!---CALCULATE SMAGORINSKY CONSTANT AND EDDY VISCOSITY------------------ 
        DEL=(DX*DY*DZ)**(1.0/3.0)
        CS2_2=PLM/PMM
        CS2_2=DMAX1(ZERO,CS2_2)    ! clip to avoid numerical problem is tend to zero
        CS2_4=PQN/PNN
        CS2_4=DMAX1(ZERO,CS2_4)
! Line with logs is not needed if test filter =2delta and second test filter =4delta
        BETA=(CS2_4/CS2_2)**(log(tf1)/(log(tf2)-log(tf1)))
        BETA=DMAX1(BETA,1.0D0/(tf1*tf2))    ! clipping Beta at 1/8
        CS2=CS2_2/BETA
     
        CELL_FV(%)CELL_VAR(6)=CS2*(DEL**2)*CELL_FV(%)CELL_VAR(19) ! SGS eddy viscosity
        DISSIP=-CS2*(DEL**2)*(CELL_FV(%)CELL_VAR(19)**3)        ! SGS dissipation

        LAG_START=0
      END IF
    END DO
!---OUTPUT SPATIAL AVERAGED NU,CS AND DISSIP
    IF(SCREEN_LEVEL.EQ.2)THEN
      CALL AVE_S_GLOBAL(NU,-NB,-NB,-NB,NU_SA,1)
      CALL AVE_S_GLOBAL(CS2,-NB,-NB,-NB,CS2_SA,1)
      CALL AVE_S_GLOBAL(BETA,-NB,-NB,-NB,BETA_SA,2)
      CALL AVE_S_GLOBAL(DISSIP,-NB,-NB,-NB,DISSIP_SA,1)
      IF(MYID.EQ.0)THEN
        PRINT*,'SPATIAL AVERAGED NU=',NU_SA,' EPS=',DISSIP_SA
        PRINT*,'SPATIAL AVERAGED CS2=',CS2_SA,' BETA=',BETA_SA
      END IF
    END IF   

    END SUBROUTINE 
!=========================================================================!
!                     LAGRANGIAN SCALE INVARIANT SGS MODEL                !
!=========================================================================!
!   SCALE-DEPENDENT VERSION (Elie Bou-Zeid et al., 2005) 
    SUBROUTINE SGS_LASI(DX,DY,DZ)
    IMPLICIT NONE
    INTEGER :: I,J,K,N1,N2,IFILTER
    INTEGER:: I_OUT_SGS,I_AVE_SGS,I_START_AVE,N_TAVE_SKIP,N_FLD_OUT,N_PRF_OUT
    INTEGER :: I_BC(6)
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE::U_BAR,V_BAR,W_BAR,                        &
                                L11,L22,L33,L12,L13,L23,                         &
                                S11_BAR,S22_BAR,S33_BAR,S12_BAR,S13_BAR,S23_BAR, &
                                S_BAR,M11,M22,M33,M12,M13,M23,                   &
                                MM,LM,VART1,VART2
    REAL(KIND=DP)::DX,DY,DZ
    REAL(KIND=DP)::CS20,CS2_2,CS2_4,U0,U1,U2,V0,V1,V2,W0,W1,W2,X1,X2,Y1,Y2,Z1,Z2,DEL,CONST
    REAL(KIND=DP)::UC,VC,WC,XP,YP,ZP,PLMP,PMMP,PQNP,PNNP,EPSI,DUM,DT,FDX,FDY,FDZ,T,TRACE
    REAL(KIND=DP)::ZERO,TF1,TF2,TF1_2,TF2_2,POWCOEFF,LAGRAN_DT
    REAL(KIND=DP)::NU_SA,CS2_SA,DISSIP_SA
    CHARACTER:: DUMC


    DATA ZERO /1.E-16/
!---SET THE BC FOR THE EDDY VISOCOSITY. THE ASSUMED NEUMANN BC IS FOR TEMPORAL USE
    DO I=1,6
      I_BC(I)=2
    END DO

    CS20=CS0**2
!   IFILTER=1: GAUSSIAN FUNCTION
!   IFILTER=2: SHARP SPECTRAL FUNCTION
!   IFILTER=3: BOX FUNCTION
    IF(ISCHEME.EQ.1)THEN
      IFILTER=3
    ELSE
      IFILTER=2
    END IF

    TF1=2.0   ! test filter 1 / grid size
    TF1_2=TF1**2

    POWCOEFF = -1.0e0/8.0e0   ! cut off coefficient for Beta=cs(2Delta)/Cs(Delta)

    LAGRAN_DT=DT*N_SGS_SKIP      ! cs update every fifth time step,i.e. cs_count=5
!---CALCULATE THE RATE OF STRAIN TENSOR (Sij)
    IF(IFILTER.EQ.2)THEN
      CALL DEALIAS(U,-NB,-NB,-NB) 
      CALL DEALIAS(V,-NB,-NB,-NB)
      CALL DEALIAS(W,-NB,-NB,-NB)
    END IF

    ALLOCATE(U_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),V_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),W_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(L11(NX1:NX2,NY1:NY2,NZ1:NZ2),L22(NX1:NX2,NY1:NY2,NZ1:NZ2),L33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             L12(NX1:NX2,NY1:NY2,NZ1:NZ2),L13(NX1:NX2,NY1:NY2,NZ1:NZ2),L23(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(S11_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S22_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S33_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S12_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S13_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S23_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(M11(NX1:NX2,NY1:NY2,NZ1:NZ2),M22(NX1:NX2,NY1:NY2,NZ1:NZ2),M33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             M12(NX1:NX2,NY1:NY2,NZ1:NZ2),M13(NX1:NX2,NY1:NY2,NZ1:NZ2),M23(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(MM(NX1:NX2,NY1:NY2,NZ1:NZ2),LM(NX1:NX2,NY1:NY2,NZ1:NZ2))

    IF(ICOLL.EQ.1)THEN
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            U_BAR(I,J,K)=U(I,J,K)
            V_BAR(I,J,K)=V(I,J,K)
            W_BAR(I,J,K)=W(I,J,K)
          END DO
        END DO
      END DO
    ELSE
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            U_BAR(I,J,K)=(U(I,J,K)+U(I+1,J,K))/2.0
            V_BAR(I,J,K)=(V(I,J,K)+V(I,J+1,K))/2.0
            W_BAR(I,J,K)=(W(I,J,K)+W(I,J,K+1))/2.0
          END DO
        END DO
      END DO
    END IF
 
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ      
          L11(I,J,K)=U_BAR(I,J,K)*U_BAR(I,J,K)  ! temporary variable
          L12(I,J,K)=U_BAR(I,J,K)*V_BAR(I,J,K)  ! temporary variable
          L13(I,J,K)=U_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
          L23(I,J,K)=V_BAR(I,J,K)*W_BAR(I,J,K) 	! temporary variable
          L22(I,J,K)=V_BAR(I,J,K)*V_BAR(I,J,K)	! temporary variable
          L33(I,J,K)=W_BAR(I,J,K)*W_BAR(I,J,K)	! temporary variable
        END DO
      END DO
    END DO
!  Here we filter the variable (e.g. u_bar at a scale =tf1*gridsize)
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,U_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,V_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,W_BAR,-NB,-NB,-NB,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,U_BAR,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,V_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,W_BAR,-NB,-NB,-NB,TF1)
    END IF
! computing Lij as defined in Bou-Zeid et al 2005 
    IF(IFILTER.NE.2)THEN            
      CALL FILTER(DX,DY,DZ,L11,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,L12,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,L13,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,L22,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,L23,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,L33,-NB,-NB,-NB,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,L11,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,L12,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,L13,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,L22,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,L23,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,L33,-NB,-NB,-NB,TF1)
    END IF

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          L11(I,J,K)=L11(I,J,K)-U_BAR(I,J,K)*U_BAR(I,J,K)
          L12(I,J,K)=L12(I,J,K)-U_BAR(I,J,K)*V_BAR(I,J,K)
          L13(I,J,K)=L13(I,J,K)-U_BAR(I,J,K)*W_BAR(I,J,K)
          L22(I,J,K)=L22(I,J,K)-V_BAR(I,J,K)*V_BAR(I,J,K)
          L23(I,J,K)=L23(I,J,K)-V_BAR(I,J,K)*W_BAR(I,J,K)
          L33(I,J,K)=L33(I,J,K)-W_BAR(I,J,K)*W_BAR(I,J,K)
!       Deviatoric part of the stress
          TRACE=(L11(I,J,K)+L22(I,J,K)+L33(I,J,K))/3.0
          L11(I,J,K)=L11(I,J,K)-TRACE
          L22(I,J,K)=L22(I,J,K)-TRACE
          L33(I,J,K)=L33(I,J,K)-TRACE
        END DO
      END DO
    END DO
 
    IF(ICOLL.EQ.1)THEN
      S12_BAR=S12  
      S13_BAR=S13  
      S23_BAR=S23  
    ELSE
      S12_BAR=S12C  
      S13_BAR=S13C  
      S23_BAR=S23C
    END IF
    S11_BAR=S11  
    S22_BAR=S22  
    S33_BAR=S33
! filtered strain rate tensors 
! _bar refers to the first test filter = 2*gridsize
! _hat refers to the second test filter = 4*gridsize
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,S11_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,S12_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,S13_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,S22_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,S23_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,S33_BAR,-NB,-NB,-NB,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,S11_BAR,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,S12_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,S13_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,S22_BAR,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,S23_BAR,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,S33_BAR,-NB,-NB,-NB,TF1)
    END IF

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          S_BAR(I,J,K)=SQRT(2.*(S11_BAR(I,J,K)**2+S22_BAR(I,J,K)**2+ &
                                S33_BAR(I,J,K)**2+                   &
                            2.*(S12_BAR(I,J,K)**2+S13_BAR(I,J,K)**2+ &
                                S23_BAR(I,J,K)**2)))
        END DO
      END DO
    END DO

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          M11(I,J,K)=S(I,J,K)*S11(I,J,K)
          M12(I,J,K)=S(I,J,K)*S12(I,J,K)
          M13(I,J,K)=S(I,J,K)*S13(I,J,K)
          M22(I,J,K)=S(I,J,K)*S22(I,J,K)
          M23(I,J,K)=S(I,J,K)*S23(I,J,K)
          M33(I,J,K)=S(I,J,K)*S33(I,J,K)
        END DO
      END DO
    END DO

    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,M11,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,M12,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,M13,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,M22,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,M23,-NB,-NB,-NB,TF1)
      CALL FILTER(DX,DY,DZ,M33,-NB,-NB,-NB,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,M11,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,M12,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,M13,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,M22,-NB,-NB,-NB,TF1) 
      CALL FILTER_SPEC(DX,DY,M23,-NB,-NB,-NB,TF1)
      CALL FILTER_SPEC(DX,DY,M33,-NB,-NB,-NB,TF1)
    END IF
! computing Mij and Nij as defined in Bou-Zeid et al 2005
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ     
          DEL=SQRT(DX*DY)
          CONST=2.0*(DEL**2) 
   	
          M11(I,J,K)=CONST*(M11(I,J,K) - TF1_2*S_BAR(I,J,K)*S11_BAR(I,J,K))
          M12(I,J,K)=CONST*(M12(I,J,K) - TF1_2*S_BAR(I,J,K)*S12_BAR(I,J,K))
          M13(I,J,K)=CONST*(M13(I,J,K) - TF1_2*S_BAR(I,J,K)*S13_BAR(I,J,K))
          M22(I,J,K)=CONST*(M22(I,J,K) - TF1_2*S_BAR(I,J,K)*S22_BAR(I,J,K))
          M23(I,J,K)=CONST*(M23(I,J,K) - TF1_2*S_BAR(I,J,K)*S23_BAR(I,J,K))
          M33(I,J,K)=CONST*(M33(I,J,K) - TF1_2*S_BAR(I,J,K)*S33_BAR(I,J,K))
!       Deviatoric part of the stress
          TRACE=(M11(I,J,K)+M22(I,J,K)+M33(I,J,K))/3.0
          M11(I,J,K)=M11(I,J,K)-TRACE
          M22(I,J,K)=M22(I,J,K)-TRACE
          M33(I,J,K)=M33(I,J,K)-TRACE
! Contracting the tensors
          LM(I,J,K)=L11(I,J,K)*M11(I,J,K)+L22(I,J,K)*M22(I,J,K)+ &
                    L33(I,J,K)*M33(I,J,K)+                       &
                2.*(L12(I,J,K)*M12(I,J,K)+L13(I,J,K)*M13(I,J,K)+ &
                    L23(I,J,K)*M23(I,J,K)) 
          MM(I,J,K)=M11(I,J,K)**2+M22(I,J,K)**2+M33(I,J,K)**2+   &
                2.*(M12(I,J,K)**2+M13(I,J,K)**2+M23(I,J,K)**2)   
        END DO
      END DO
    END DO

    DEALLOCATE(U_BAR,V_BAR,W_BAR)
    DEALLOCATE(L11,L22,L33,L12,L13,L23)
    DEALLOCATE(S11_BAR,S22_BAR,S33_BAR,S12_BAR,S13_BAR,S23_BAR,S_BAR)
    DEALLOCATE(M11,M22,M33,M12,M13,M23)
! for first test filter
! the lines below check if this is the first time step at which we compute P (P is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ratio = cs^2 =CS20
    IF(LAG_START.EQ.1)THEN
      DO J=1,NY
        DO K=1,NZ
          DO I=1,NX
            PMM(I,J,K)=DMAX1(ZERO,MM(I,J,K))
            PLM(I,J,K)=PMM(I,J,K)*CS20
          END DO
        END DO
      END DO
    ELSE   ! CALCULATE PLM & PMM USING LAGRANGIAN AVERAGING 
      ALLOCATE(VART1(NXT,NYT,NZT),VART2(NXT,NYT,NZT))
      CALL ASSEM_ALL(PLM,-NB,-NB,-NB,VART1)
      CALL ASSEM_ALL(PMM,-NB,-NB,-NB,VART2)
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            DEL=SQRT(DX*DY)
            T=1.5*DEL*ABS(PLM(I,J,K)*PMM(I,J,K))**(-1.0/8.0)
            T= DMAX1(1D-24,T)	    ! clip to avoid numerical problem is tend to zero
            EPSI=LAGRAN_DT/T/(LAGRAN_DT/T+1.0)
            XP=XI(I+MYIDX*NX)-U(I,J,K)*LAGRAN_DT
            YP=YI(J+MYIDY*NY)-V(I,J,K)*LAGRAN_DT
            ZP=ZI(K+MYIDZ*NZ)-W(I,J,K)*LAGRAN_DT
            IF(XP.GT.XI(NXT))THEN
              XP=XI(1)+ABS(XP-XI(NXT))
            ELSE IF(XP.LT.XI(1))THEN
              XP=XI(NXT)-ABS(XI(1)-XP)
            END IF
            IF(ZP.GT.ZI(NZT))THEN
              ZP=ZI(1)+ABS(ZP-ZI(NZT))
            ELSE IF(ZP.LT.ZI(1))THEN
              ZP=ZI(NZT)-ABS(ZI(1)-ZP)
            END IF
            YP=DMIN1(DMAX1(YP,YI(1)),YI(NYT))
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART1,1,1,1,PLMP)
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART2,1,1,1,PMMP)
            PLM(I,J,K)=DMAX1(EPSI*LM(I,J,K)+(1.0D0-EPSI)*PLMP,ZERO)
            PMM(I,J,K)=EPSI*MM(I,J,K)+(1.0e0-EPSI)*PMMP  
          END DO
        END DO
      END DO
      DEALLOCATE(VART1,VART2)
    END IF
    DEALLOCATE(LM,MM)
!---CALCULATE SMAGORINSKY CONSTANT AND EDDY VISCOSITY------------------ 
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          DEL=(DX*DY)**(1.0/2.0)
          CS2_2=PLM(I,J,K)/PMM(I,J,K)
          CS2(I,J,K)=DMAX1(ZERO,CS2_2)    ! clip to avoid numerical problem is tend to zero
     
          NU(I,J,K)=CS2(I,J,K)*(DEL**2)*S(I,J,K)                ! SGS eddy viscosity
          DISSIP(I,J,K)=-CS2(I,J,K)*(DEL**2)*(S(I,J,K)**3)        ! SGS dissipation
        END DO
      END DO
    END DO
    CALL GET_BC(NX,NY,NZ,NU,NBX,NBY,NBZ,0,I_BC)
   
    LAG_START=0
!---OUTPUT SPATIAL AVERAGED NU,CS AND DISSIP
    IF(SCREEN_LEVEL.EQ.2)THEN
      CALL AVE_S_GLOBAL(NU,-NB,-NB,-NB,NU_SA,1)
      CALL AVE_S_GLOBAL(CS2,-NB,-NB,-NB,CS2_SA,1)
      CALL AVE_S_GLOBAL(DISSIP,-NB,-NB,-NB,DISSIP_SA,1)
      IF(MYID.EQ.0)THEN
        PRINT*,'SPATIAL AVERAGED NU=',NU_SA,' EPS=',DISSIP_SA
        PRINT*,'SPATIAL AVERAGED CS2=',CS2_SA
      END IF
    END IF

    END SUBROUTINE
!*************************************************************************!
!              SUBROUTINE OF CONSTANT SMAGORINSKY SGS MODEL               !
!*************************************************************************!
!   SCALE-DEPENDENT VERSION (Elie Bou-Zeid et al., 2005)
    SUBROUTINE SGS_C(DX,DY,DZ) 

    IMPLICIT NONE
    INTEGER:: I,J,K
    INTEGER :: I_BC(6)
    REAL(KIND=DP):: DX,DY,DZ,DEL,CS20,NU_SA,DISSIP_SA
!---SET THE BC FOR THE EDDY VISOCOSITY. THE ASSUMED NEUMANN BC IS FOR TEMPORAL USE
    DO I=1,6
      I_BC(I)=2
    END DO
!---CALCULATE SMAGORINSKY CONSTANT AND EDDY VISCOSITY------------------     
    CS20=CS0**2

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          DEL=(DX*DY*DZ)**(1.0/3.0)
          NU(I,J,K)=CS20*(DEL**2)*S(I,J,K)                ! SGS eddy viscosity
          DISSIP(I,J,K)=-CS20*(DEL**2)*(S(I,J,K)**3)        ! SGS dissipation
        END DO
      END DO
    END DO
    CALL GET_BC(NX,NY,NZ,NU,NBX,NBY,NBZ,0,I_BC)
!---OUTPUT SPATIAL AVERAGED NU,CS AND DISSIP
    IF(SCREEN_LEVEL.EQ.2)THEN
      CALL AVE_S_GLOBAL(NU,-NB,-NB,-NB,NU_SA,1)
      CALL AVE_S_GLOBAL(DISSIP,-NB,-NB,-NB,DISSIP_SA,1)
      IF(MYID.EQ.0)THEN
        PRINT*,'SPATIAL AVERAGED NU=',NU_SA,' EPS=',DISSIP_SA
      END IF
    END IF

    END SUBROUTINE
  
  END MODULE
