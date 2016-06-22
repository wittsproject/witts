! This module contains the subgrid-scale (SGS) models
!
  MODULE SGS

  USE mpi
  USE parameters
  USE field_shared   
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
    INTEGER :: I_BC(6)
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE::U_BAR,V_BAR,W_BAR,U_HAT,V_HAT,W_HAT,      &
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
    ALLOCATE(U_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),V_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),W_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(U_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2),V_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2),W_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(L11(NX1:NX2,NY1:NY2,NZ1:NZ2),L22(NX1:NX2,NY1:NY2,NZ1:NZ2),L33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             L12(NX1:NX2,NY1:NY2,NZ1:NZ2),L13(NX1:NX2,NY1:NY2,NZ1:NZ2),L23(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             Q11(NX1:NX2,NY1:NY2,NZ1:NZ2),Q22(NX1:NX2,NY1:NY2,NZ1:NZ2),Q33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             Q12(NX1:NX2,NY1:NY2,NZ1:NZ2),Q13(NX1:NX2,NY1:NY2,NZ1:NZ2),Q23(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(S11_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S22_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S33_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S12_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S13_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S23_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S11_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2),S22_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S33_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2),S12_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S13_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2),S23_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             S_BAR(NX1:NX2,NY1:NY2,NZ1:NZ2),S_HAT(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(M11(NX1:NX2,NY1:NY2,NZ1:NZ2),M22(NX1:NX2,NY1:NY2,NZ1:NZ2),M33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             M12(NX1:NX2,NY1:NY2,NZ1:NZ2),M13(NX1:NX2,NY1:NY2,NZ1:NZ2),M23(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             N11(NX1:NX2,NY1:NY2,NZ1:NZ2),N22(NX1:NX2,NY1:NY2,NZ1:NZ2),N33(NX1:NX2,NY1:NY2,NZ1:NZ2), &
             N12(NX1:NX2,NY1:NY2,NZ1:NZ2),N13(NX1:NX2,NY1:NY2,NZ1:NZ2),N23(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(MM(NX1:NX2,NY1:NY2,NZ1:NZ2),LM(NX1:NX2,NY1:NY2,NZ1:NZ2),NN(NX1:NX2,NY1:NY2,NZ1:NZ2),QN(NX1:NX2,NY1:NY2,NZ1:NZ2))

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

    U_HAT=U_BAR
    V_HAT=V_BAR
    W_HAT=W_BAR

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ      
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
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,U_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,V_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,W_BAR,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,U_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,V_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,W_BAR,NX1,NY1,NZ1,TF1)
    END IF
! computing Lij as defined in Bou-Zeid et al 2005 
    IF(IFILTER.NE.2)THEN            
      CALL FILTER(DX,DY,DZ,L11,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L12,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L13,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L22,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L23,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L33,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,L11,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,L12,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L13,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L22,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,L23,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L33,NX1,NY1,NZ1,TF1)
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
!  Here we test filter the variable (e.g. u_bar at a scale =tf2*gridsize)
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,U_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,V_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,W_HAT,NX1,NY1,NZ1,TF2)
    ELSE
      CALL FILTER_SPEC(DX,DY,U_HAT,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,V_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,W_HAT,NX1,NY1,NZ1,TF2)	
    END IF
! computing Qij as defined in Bou-Zeid et al 2005
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,Q11,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,Q12,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,Q13,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,Q22,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,Q23,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,Q33,NX1,NY1,NZ1,TF2)
    ELSE
      CALL FILTER_SPEC(DX,DY,Q11,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,Q12,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,Q13,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,Q22,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,Q23,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,Q33,NX1,NY1,NZ1,TF2)
    END IF

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          Q11(I,J,K)=Q11(I,J,K)-U_HAT(I,J,K)*U_HAT(I,J,K)
          Q12(I,J,K)=Q12(I,J,K)-U_HAT(I,J,K)*V_HAT(I,J,K)
          Q13(I,J,K)=Q13(I,J,K)-U_HAT(I,J,K)*W_HAT(I,J,K)
          Q22(I,J,K)=Q22(I,J,K)-V_HAT(I,J,K)*V_HAT(I,J,K)
          Q23(I,J,K)=Q23(I,J,K)-V_HAT(I,J,K)*W_HAT(I,J,K)
          Q33(I,J,K)=Q33(I,J,K)-W_HAT(I,J,K)*W_HAT(I,J,K)
!       Deviatoric part of the stress
          TRACE=(Q11(I,J,K)+Q22(I,J,K)+Q33(I,J,K))/3.0
          Q11(I,J,K)=Q11(I,J,K)-TRACE
          Q22(I,J,K)=Q22(I,J,K)-TRACE
          Q33(I,J,K)=Q33(I,J,K)-TRACE
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

    S11_HAT=S11_BAR
    S12_HAT=S12_BAR
    S13_HAT=S13_BAR
    S22_HAT=S22_BAR
    S23_HAT=S23_BAR
    S33_HAT=S33_BAR
! filtered strain rate tensors 
! _bar refers to the first test filter = 2*gridsize
! _hat refers to the second test filter = 4*gridsize
    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,S11_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S12_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S13_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S22_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S23_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S33_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S11_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,S12_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,S13_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,S22_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,S23_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,S33_HAT,NX1,NY1,NZ1,TF2)
    ELSE
      CALL FILTER_SPEC(DX,DY,S11_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,S12_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S13_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S22_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,S23_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S33_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S11_HAT,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,S12_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,S13_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,S22_HAT,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,S23_HAT,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,S33_HAT,NX1,NY1,NZ1,TF2)
    END IF

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
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

    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          M11(I,J,K)=S(I,J,K)*S11(I,J,K)
          M12(I,J,K)=S(I,J,K)*S12(I,J,K)
          M13(I,J,K)=S(I,J,K)*S13(I,J,K)
          M22(I,J,K)=S(I,J,K)*S22(I,J,K)
          M23(I,J,K)=S(I,J,K)*S23(I,J,K)
          M33(I,J,K)=S(I,J,K)*S33(I,J,K)

          N11(I,J,K)=M11(I,J,K)
          N12(I,J,K)=M12(I,J,K)
          N13(I,J,K)=M13(I,J,K)
          N22(I,J,K)=M22(I,J,K)
          N23(I,J,K)=M23(I,J,K)
          N33(I,J,K)=M33(I,J,K)
        END DO
      END DO
    END DO

    IF(IFILTER.NE.2)THEN
      CALL FILTER(DX,DY,DZ,M11,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M12,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M13,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M22,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M23,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M33,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,N11,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,N12,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,N13,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,N22,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,N23,NX1,NY1,NZ1,TF2)
      CALL FILTER(DX,DY,DZ,N33,NX1,NY1,NZ1,TF2)
    ELSE
      CALL FILTER_SPEC(DX,DY,M11,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,M12,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M13,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M22,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,M23,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M33,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,N11,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,N12,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,N13,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,N22,NX1,NY1,NZ1,TF2) 
      CALL FILTER_SPEC(DX,DY,N23,NX1,NY1,NZ1,TF2)
      CALL FILTER_SPEC(DX,DY,N33,NX1,NY1,NZ1,TF2)
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
	 
          N11(I,J,K)=CONST*(N11(I,J,K) - TF2_2*S_HAT(I,J,K)*S11_HAT(I,J,K))
          N12(I,J,K)=CONST*(N12(I,J,K) - TF2_2*S_HAT(I,J,K)*S12_HAT(I,J,K))
          N13(I,J,K)=CONST*(N13(I,J,K) - TF2_2*S_HAT(I,J,K)*S13_HAT(I,J,K))
          N22(I,J,K)=CONST*(N22(I,J,K) - TF2_2*S_HAT(I,J,K)*S22_HAT(I,J,K))
          N23(I,J,K)=CONST*(N23(I,J,K) - TF2_2*S_HAT(I,J,K)*S23_HAT(I,J,K))
          N33(I,J,K)=CONST*(N33(I,J,K) - TF2_2*S_HAT(I,J,K)*S33_HAT(I,J,K))
!       Deviatoric part of the stress
          TRACE=(N11(I,J,K)+N22(I,J,K)+N33(I,J,K))/3.0
          N11(I,J,K)=N11(I,J,K)-TRACE
          N22(I,J,K)=N22(I,J,K)-TRACE
          N33(I,J,K)=N33(I,J,K)-TRACE
! Contracting the tensors
          LM(I,J,K)=L11(I,J,K)*M11(I,J,K)+L22(I,J,K)*M22(I,J,K)+ &
                    L33(I,J,K)*M33(I,J,K)+                       &
                2.*(L12(I,J,K)*M12(I,J,K)+L13(I,J,K)*M13(I,J,K)+ &
                    L23(I,J,K)*M23(I,J,K)) 
          MM(I,J,K)=M11(I,J,K)**2+M22(I,J,K)**2+M33(I,J,K)**2+   &
                2.*(M12(I,J,K)**2+M13(I,J,K)**2+M23(I,J,K)**2)   

          QN(I,J,K)=Q11(I,J,K)*N11(I,J,K)+Q22(I,J,K)*N22(I,J,K)+ &
                    Q33(I,J,K)*N33(I,J,K)+                       &
                2.*(Q12(I,J,K)*N12(I,J,K)+Q13(I,J,K)*N13(I,J,K)+ &
                    Q23(I,J,K)*N23(I,J,K))
          NN(I,J,K)=N11(I,J,K)**2+N22(I,J,K)**2+N33(I,J,K)**2+   &
                2.*(N12(I,J,K)**2+N13(I,J,K)**2+N23(I,J,K)**2)
        END DO
      END DO
    END DO

    DEALLOCATE(U_BAR,V_BAR,W_BAR,U_HAT,V_HAT,W_HAT)
    DEALLOCATE(L11,L22,L33,L12,L13,L23,Q11,Q22,Q33,Q12,Q13,Q23)
    DEALLOCATE(S11_BAR,S22_BAR,S33_BAR,S12_BAR,S13_BAR,S23_BAR,S_BAR, &
               S11_HAT,S22_HAT,S33_HAT,S12_HAT,S13_HAT,S23_HAT,S_HAT)
    DEALLOCATE(M11,M22,M33,M12,M13,M23,N11,N22,N33,N12,N13,N23)
! for first test filter
! the lines below check if this is the first time step at which we compute P (P is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ratio = cs^2 =CS2
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
      CALL ASSEM_ALL(PLM,NX1,NY1,NZ1,VART1)
      CALL ASSEM_ALL(PMM,NX1,NY1,NZ1,VART2)
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
            XP=DMIN1(DMAX1(XP,XI(1)),XI(NXT))
            YP=DMIN1(DMAX1(YP,YI(1)),YI(NYT))
            ZP=DMIN1(DMAX1(ZP,ZI(1)),ZI(NZT))
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART1,PLMP)
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART2,PMMP)
            PLM(I,J,K)=DMAX1(EPSI*LM(I,J,K)+(1.0D0-EPSI)*PLMP,ZERO)
            PMM(I,J,K)=EPSI*MM(I,J,K)+(1.0e0-EPSI)*PMMP  
          END DO
        END DO
      END DO
      DEALLOCATE(VART1,VART2)
    END IF
    DEALLOCATE(LM,MM)
! for second test filter
! the lines below check if this is the first time step at which we compute F (F is equivalent to I in bou-Zeid et al 2005)
! if it is the first time step it initializes the denominator as the current MM and the numerator such as the ration = cs^2 =CS20
    IF(LAG_START.EQ.1)THEN
      DO J=1,NY
        DO K=1,NZ
          DO I=1,NX
            PNN(I,J,K)=DMAX1(ZERO,NN(I,J,K))
            PQN(I,J,K)=PNN(I,J,K)*CS20               
          END DO
        END DO
      END DO
    ELSE  ! CALCULATE PQN & PNN USING LAGRANGIAN AVERAGING    
      ALLOCATE(VART1(NXT,NYT,NZT),VART2(NXT,NYT,NZT))
      CALL ASSEM_ALL(PQN,NX1,NY1,NZ1,VART1)
      CALL ASSEM_ALL(PNN,NX1,NY1,NZ1,VART2)
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            DEL=SQRT(DX*DY)
            T=1.5d0*DEL*ABS(PQN(I,J,K)*PNN(I,J,K))**(-1.0e0/8.0e0)
            T=DMAX1(1D-24,T)	    ! clip to avoid numerical problem is tend to zero
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
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART1,PQNP)
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART2,PNNP)
            PQN(I,J,K)=DMAX1(EPSI*QN(I,J,K)+(1.0D0-EPSI)*PQNP,ZERO)
            PNN(I,J,K)=EPSI*NN(I,J,K)+(1.0e0-EPSI)*PNNP
          END DO
        END DO
      END DO
      DEALLOCATE(VART1,VART2)
    END IF
    DEALLOCATE(QN,NN)
!---CALCULATE SMAGORINSKY CONSTANT AND EDDY VISCOSITY------------------ 
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          DEL=(DX*DY)**(1.0/2.0)
          CS2_2=PLM(I,J,K)/PMM(I,J,K)
          CS2_2=DMAX1(ZERO,CS2_2)    ! clip to avoid numerical problem is tend to zero
          CS2_4=PQN(I,J,K)/PNN(I,J,K)
          CS2_4=DMAX1(ZERO,CS2_4)
! Line with logs is not needed if test filter =2delta and second test filter =4delta
          BETA(I,J,K)=(CS2_4/CS2_2)**(log(tf1)/(log(tf2)-log(tf1)))
          BETA(I,J,K)=DMAX1(BETA(I,J,K),1.0D0/(tf1*tf2))    ! clipping Beta at 1/8
          CS2(I,J,K)=CS2_2/BETA(I,J,K)
     
          NU(I,J,K)=CS2(I,J,K)*(DEL**2)*S(I,J,K)                ! SGS eddy viscosity
          DISSIP(I,J,K)=-CS2(I,J,K)*(DEL**2)*(S(I,J,K)**3)        ! SGS dissipation
        END DO
      END DO
    END DO
    CALL GET_BC(NX,NY,NZ,NU,NBX,NBY,NBZ,0,I_BC)

    LAG_START=0

!---OUTPUT SPATIAL AVERAGED NU,CS AND DISSIP
    IF(SCREEN_LEVEL.EQ.2)THEN
      CALL AVE_S_GLOBAL(NU,NX1,NY1,NZ1,NU_SA,1)
      CALL AVE_S_GLOBAL(CS2,NX1,NY1,NZ1,CS2_SA,1)
      CALL AVE_S_GLOBAL(BETA,NX1,NY1,NZ1,BETA_SA,2)
      CALL AVE_S_GLOBAL(DISSIP,NX1,NY1,NZ1,DISSIP_SA,1)
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
      CALL DEALIAS(U,NX1,NY1,NZ1) 
      CALL DEALIAS(V,NX1,NY1,NZ1)
      CALL DEALIAS(W,NX1,NY1,NZ1)
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
      CALL FILTER(DX,DY,DZ,U_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,V_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,W_BAR,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,U_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,V_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,W_BAR,NX1,NY1,NZ1,TF1)
    END IF
! computing Lij as defined in Bou-Zeid et al 2005 
    IF(IFILTER.NE.2)THEN            
      CALL FILTER(DX,DY,DZ,L11,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L12,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L13,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L22,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L23,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,L33,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,L11,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,L12,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L13,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L22,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,L23,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,L33,NX1,NY1,NZ1,TF1)
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
      CALL FILTER(DX,DY,DZ,S11_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S12_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S13_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S22_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S23_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,S33_BAR,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,S11_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,S12_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S13_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S22_BAR,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,S23_BAR,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,S33_BAR,NX1,NY1,NZ1,TF1)
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
      CALL FILTER(DX,DY,DZ,M11,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M12,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M13,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M22,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M23,NX1,NY1,NZ1,TF1)
      CALL FILTER(DX,DY,DZ,M33,NX1,NY1,NZ1,TF1)
    ELSE
      CALL FILTER_SPEC(DX,DY,M11,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,M12,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M13,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M22,NX1,NY1,NZ1,TF1) 
      CALL FILTER_SPEC(DX,DY,M23,NX1,NY1,NZ1,TF1)
      CALL FILTER_SPEC(DX,DY,M33,NX1,NY1,NZ1,TF1)
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
      CALL ASSEM_ALL(PLM,NX1,NY1,NZ1,VART1)
      CALL ASSEM_ALL(PMM,NX1,NY1,NZ1,VART2)
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
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART1,PLMP)
            CALL INTER_GLOBAL(XP,YP,ZP,NXT,NYT,NZT,XI,YI,ZI,1,VART2,PMMP)
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
      CALL AVE_S_GLOBAL(NU,NX1,NY1,NZ1,NU_SA,1)
      CALL AVE_S_GLOBAL(CS2,NX1,NY1,NZ1,CS2_SA,1)
      CALL AVE_S_GLOBAL(DISSIP,NX1,NY1,NZ1,DISSIP_SA,1)
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
      CALL AVE_S_GLOBAL(NU,NX1,NY1,NZ1,NU_SA,1)
      CALL AVE_S_GLOBAL(DISSIP,NX1,NY1,NZ1,DISSIP_SA,1)
      IF(MYID.EQ.0)THEN
        PRINT*,'SPATIAL AVERAGED NU=',NU_SA,' EPS=',DISSIP_SA
      END IF
    END IF

    END SUBROUTINE
  
  END MODULE
