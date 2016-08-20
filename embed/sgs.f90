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
    SUBROUTINE SGS_LASD()
    IMPLICIT NONE
!    INCLUDE "mpif.h"
    INTEGER :: I,J,K,M,N1,N2,IFILTER
    INTEGER :: I_OUT_SGS,I_AVE_SGS,I_START_AVE,N_TAVE_SKIP,N_FLD_OUT,N_PRF_OUT
    INTEGER :: NB
    PARAMETER(NB=2)
    REAL(KIND=DP),DIMENSION(-NB:NB,-NB:NB,-NB:NB):: &
                  U_BAR,V_BAR,W_BAR,U_HAT,V_HAT,W_HAT,      &
                  L11,L22,L33,L12,L13,L23,Q11,Q22,Q33,Q12,Q13,Q23, &
                  S11_BAR,S22_BAR,S33_BAR,S12_BAR,S13_BAR,S23_BAR, &
                  S11_HAT,S22_HAT,S33_HAT,S12_HAT,S13_HAT,S23_HAT, &
                  S_BAR,S_HAT,M11,M22,M33,M12,M13,M23,N11,N22,N33,N12,N13,N23, &
                                
    REAL(KIND=DP)::DX,DY,DZ
    REAL(KIND=DP)::CS20,CS2_2,CS2_4,U0,U1,U2,V0,V1,V2,W0,W1,W2,X1,X2,Y1,Y2,Z1,Z2,DEL,CONST
    REAL(KIND=DP)::XP,YP,ZP,MM,LM,NN,QN,PLMP,PMMP,PQNP,PNNP,EPSI,DUM,T,TRACE
    REAL(KIND=DP)::ZERO,TF1,TF2,TF1_2,TF2_2,POWCOEFF,LAGRAN_DT
    REAL(KIND=DP)::NU_SA,CS2_SA,BETA_SA,DISSIP_SA
    CHARACTER:: DUMC


    DATA ZERO /1.E-16/

!---SET THE BC FOR THE EDDY VISOCOSITY. THE ASSUMED NEUMANN BC IS FOR TEMPORAL USE 
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
              L13(I,J,K)=U_BAR(I,J,K)*W_BAR(I,J,K)  ! temporary variable
              L23(I,J,K)=V_BAR(I,J,K)*W_BAR(I,J,K)  ! temporary variable
              L22(I,J,K)=V_BAR(I,J,K)*V_BAR(I,J,K)  ! temporary variable
              L33(I,J,K)=W_BAR(I,J,K)*W_BAR(I,J,K)  ! temporary variable

              Q11(I,J,K)=L11(I,J,K)  ! temporary variable
              Q12(I,J,K)=L12(I,J,K)  ! temporary variable
              Q13(I,J,K)=L13(I,J,K)  ! temporary variable
              Q23(I,J,K)=L23(I,J,K)  ! temporary variable
              Q22(I,J,K)=L22(I,J,K)  ! temporary variable
              Q33(I,J,K)=L33(I,J,K)  ! temporary variable
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

        L11(0,0,0)=L11(0,0,0)-U_BAR(0,0,0)*U_BAR(0,0,0)
        L12(0,0,0)=L12(0,0,0)-U_BAR(0,0,0)*V_BAR(0,0,0)
        L13(0,0,0)=L13(0,0,0)-U_BAR(0,0,0)*W_BAR(0,0,0)
        L22(0,0,0)=L22(0,0,0)-V_BAR(0,0,0)*V_BAR(0,0,0)
        L23(0,0,0)=L23(0,0,0)-V_BAR(0,0,0)*W_BAR(0,0,0)
        L33(0,0,0)=L33(0,0,0)-W_BAR(0,0,0)*W_BAR(0,0,0)
!       Deviatoric part of the stress
        TRACE=(L11(0,0,0)+L22(0,0,0)+L33(0,0,0))/3.0
        L11(0,0,0)=L11(0,0,0)-TRACE
        L22(0,0,0)=L22(0,0,0)-TRACE
        L33(0,0,0)=L33(0,0,0)-TRACE
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

        Q11(0,0,0)=Q11(0,0,0)-U_HAT(0,0,0)*U_HAT(0,0,0)
        Q12(0,0,0)=Q12(0,0,0)-U_HAT(0,0,0)*V_HAT(0,0,0)
        Q13(0,0,0)=Q13(0,0,0)-U_HAT(0,0,0)*W_HAT(0,0,0)
        Q22(0,0,0)=Q22(0,0,0)-V_HAT(0,0,0)*V_HAT(0,0,0)
        Q23(0,0,0)=Q23(0,0,0)-V_HAT(0,0,0)*W_HAT(0,0,0)
        Q33(0,0,0)=Q33(0,0,0)-W_HAT(0,0,0)*W_HAT(0,0,0)

        TRACE=(Q11(0,0,0)+Q22(0,0,0)+Q33(0,0,0))/3.0
        Q11(0,0,0)=Q11(0,0,0)-TRACE
        Q22(0,0,0)=Q22(0,0,0)-TRACE
        Q33(0,0,0)=Q33(0,0,0)-TRACE

        CALL CELL_TO_STRUCT(M,NB,13,S11_BAR)
        CALL CELL_TO_STRUCT(M,NB,14,S22_BAR)
        CALL CELL_TO_STRUCT(M,NB,15,S33_BAR)
        CALL CELL_TO_STRUCT(M,NB,16,S12_BAR)
        CALL CELL_TO_STRUCT(M,NB,17,S13_BAR)
        CALL CELL_TO_STRUCT(M,NB,18,S23_BAR)
        CALL CELL_TO_STRUCT(M,NB,19,S_BAR)        

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
        CALL FILTER(DX,DY,DZ,S_BAR,-NB,-NB,-NB,TF1)
        
        CALL FILTER(DX,DY,DZ,S11_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S12_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S13_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S22_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S23_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S33_HAT,-NB,-NB,-NB,TF2)
        CALL FILTER(DX,DY,DZ,S_HAT,-NB,-NB,-NB,TF2)

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
        DEL=(DX*DY*DZ)**(1.0/3.0)
        CONST=2.0*(DEL**3) 
   	
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
        LM=L11(0,0,0)*M11(0,0,0)+L22(0,0,0)*M22(0,0,0)+L33(0,0,0)*M33(0,0,0)+   &
           2.*(L12(0,0,0)*M12(0,0,0)+L13(0,0,0)*M13(0,0,0)+L23(0,0,0)*M23(0,0,0)) 
        MM=M11(0,0,0)**2+M22(0,0,0)**2+M33(0,0,0)**2+   &
           2.*(M12(0,0,0)**2+M13(0,0,0)**2+M23(0,0,0)**2)   

        QN=Q11(0,0,0)*N11(0,0,0)+Q22(0,0,0)*N22(0,0,0)+Q33(0,0,0)*N33(0,0,0)+   &
           2.*(Q12(0,0,0)*N12(0,0,0)+Q13(0,0,0)*N13(0,0,0)+Q23(0,0,0)*N23(0,0,0))
        NN=N11(0,0,0)**2+N22(0,0,0)**2+N33(0,0,0)**2+   &
           2.*(N12(0,0,0)**2+N13(0,0,0)**2+N23(0,0,0)**2)
! for first test filter
! the lines below check if this is the first time step at which we compute P (P is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ratio = cs^2 =CS2
        IF(LAG_START.EQ.1)THEN
          PMM=DMAX1(ZERO,MM)
          PLM=PMM*CS20
        ELSE   ! CALCULATE PLM & PMM USING LAGRANGIAN AVERAGING
          PLM=CELL_FV(M)%CELL_VAR(29)
          PMM=CELL_FV(M)%CELL_VAR(30)
          
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
          
          PLM=DMAX1(EPSI*LM+(1.0D0-EPSI)*PLMP,ZERO)
          PMM=EPSI*MM+(1.0e0-EPSI)*PMMP  
        END IF
! for second test filter
! the lines below check if this is the first time step at which we compute F (F is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ration = cs^2 =CS20
        IF(LAG_START.EQ.1)THEN
          PNN=DMAX1(ZERO,NN)
          PQN=PNN*CS20               
        ELSE  ! CALCULATE PQN & PNN USING LAGRANGIAN AVERAGING
          PQN=CELL_FV(M)%CELL_VAR(31)
          PNN=CELL_FV(M)%CELL_VAR(32)
           
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
          
          PQN=DMAX1(EPSI*QN+(1.0D0-EPSI)*PQNP,ZERO)
          PNN=EPSI*NN+(1.0e0-EPSI)*PNNP
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
     
        CELL_FV(M)%CELL_VAR(6)=CS2*(DEL**2)*CELL_FV(M)%CELL_VAR(19) ! SGS eddy viscosity
        DISSIP=-CS2*(DEL**2)*(CELL_FV(%)CELL_VAR(19)**3)        ! SGS dissipation

        CELL_FV(M)%CELL_VAR(29)=PLM
        CELL_FV(M)%CELL_VAR(30)=PMM
        CELL_FV(M)%CELL_VAR(31)=PQN
        CELL_FV(M)%CELL_VAR(32)=PNN        

        LAG_START=0        
      END IF    
    END DO

    CALL GHOST_BOUNDARY(19)
    
    END SUBROUTINE 
!*************************************************************************!
!              SUBROUTINE OF CONSTANT SMAGORINSKY SGS MODEL               !
!*************************************************************************!
!   SCALE-DEPENDENT VERSION (Elie Bou-Zeid et al., 2005)
    SUBROUTINE SGS_C() 

    IMPLICIT NONE
    INTEGER:: M
    REAL(KIND=DP):: DX,DY,DZ,DEL,CS20

!---CALCULATE SMAGORINSKY CONSTANT AND EDDY VISCOSITY------------------     
    CS20=CS0**2

    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
        DX=CELL_FV(M)%CELL_DX
        DY=CELL_FV(M)%CELL_DY
        DZ=CELL_FV(M)%CELL_DZ
        
        DEL=(DX*DY*DZ)**(1.0/3.0)
        CELL_FV(M)%CELL_VAR(6)=CS20*(DEL**2)*CELL_FV(M)%CELL_VAR(19)                ! SGS eddy viscosity
        DISSIP=-CS2*(DEL**2)*(CELL_FV(M)%CELL_VAR(19)**3)        ! SGS dissipation        
      END IF
    END DO
    CALL GET_BC(NX,NY,NZ,NU,NBX,NBY,NBZ,0,I_BC)

    END SUBROUTINE
  
  END MODULE
