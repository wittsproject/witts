! This module is for controling of the post process, including exporting //
! result files and statistics
!
  MODULE post_process

  USE mpi
  USE parameters
  USE field_shared   
  USE tools
  USE output

  CONTAINS
!=========================================================================!
!                  MAIN SUBROUTINE OF POST PROCESS                        !
!=========================================================================!
    SUBROUTINE POSTPROCESS()

    IMPLICIT NONE

    INTEGER :: I_3D_PRIME,I_3D_TAU,I_3D_SGS,N_1D,N_3D
    INTEGER :: I_RE_PRIME,I_RE_SGS,N_RE
    INTEGER :: STAT_FLAG,ORDER_STAT,N_SKIP
    REAL(KIND=DP) :: T_STAT_1,T_STAT_2 
      
    INTEGER :: I,J,K,N_AVE_1,N_AVE_2
    REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: VAR,VARA
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: UC,VC,WC
    REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: TAU11C,TAU22C,TAU33C, &
                                                  TAU12C,TAU13C,TAU23C
    CHARACTER*10,DIMENSION(:),ALLOCATABLE:: NAME_VAR
    CHARACTER*10 :: FNAME

    OPEN(1,FILE="post.in")
    READ(1,*) 
    READ(1,*) I_3D_PRIME   !(1: GENERATE 3D FIELD OF PRIME VARIABLES)
    READ(1,*) I_3D_TAU     !(1: GENERATE 3D FIELD OF STRESS TENSOR)
    READ(1,*) I_3D_SGS     !(1: GENERATE 3D FIELD OF SGS VARIABLES)
    READ(1,*) N_1D         !(TIME STEPS OF GENERATING AND SAVING 1D TEMPORAL DATA)
    READ(1,*) N_3D         !(TIME STEPS OF GENERATING AND SAVING 3D FIELD DATA)
    READ(1,*) 
    READ(1,*) I_RE_PRIME   !(1: GENERATE 3D FIELD OF PRIME VARIABLES)
    READ(1,*) I_RE_SGS     !(1: GENERATE 3D FIELD OF SGS VARIABLES)
    READ(1,*) N_RE         !(TIME STEPS OF GENERATING AND SAVING 3D FIELD DATA)
    READ(1,*)  
    READ(1,*) STAT_FLAG    !(FLAG OF TIME-AVERAGE BASED STATISTICS)
    READ(1,*) ORDER_STAT   !(ORDER OF STATISTICS)
    READ(1,*) N_SKIP       !(TIME STEPS SKIPPED)
    READ(1,*) T_STAT_1     !(START TIME OF 1ST-ORDER STATISTICS)
    READ(1,*) T_STAT_2     !(START TIME OF 2ND-ORDER STATISTICS)
    CLOSE(1) 
!---GET VELOCITY AND STRESS TENSOR AT THE CENTER OF CELL
    ALLOCATE(UC(NX1:NX2,NY1:NY2,NZ1:NZ2),VC(NX1:NX2,NY1:NY2,NZ1:NZ2),WC(NX1:NX2,NY1:NY2,NZ1:NZ2))
    ALLOCATE(TAU11C(NX1:NX2,NY1:NY2,NZ1:NZ2),TAU22C(NX1:NX2,NY1:NY2,NZ1:NZ2),TAU33C(NX1:NX2,NY1:NY2,NZ1:NZ2),  &
             TAU12C(NX1:NX2,NY1:NY2,NZ1:NZ2),TAU13C(NX1:NX2,NY1:NY2,NZ1:NZ2),TAU23C(NX1:NX2,NY1:NY2,NZ1:NZ2))
    IF(ICOLL.EQ.1)THEN
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            UC(I,J,K)=U(I,J,K)
            VC(I,J,K)=V(I,J,K)
            WC(I,J,K)=W(I,J,K)
          END DO
        END DO
      END DO
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            TAU11C(I,J,K)=TAU11(I,J,K)
            TAU22C(I,J,K)=TAU22(I,J,K)
            TAU33C(I,J,K)=TAU33(I,J,K)
            TAU12C(I,J,K)=TAU12(I,J,K)
            TAU13C(I,J,K)=TAU13(I,J,K)
            TAU23C(I,J,K)=TAU23(I,J,K)
          END DO
        END DO
      END DO
    ELSE
      DO I=1,NX
        DO J=1,NY
          DO K=1,NZ
            UC(I,J,K)=(U(I,J,K)+U(I+1,J,K))/2.0
            VC(I,J,K)=(V(I,J,K)+V(I,J+1,K))/2.0
            WC(I,J,K)=(W(I,J,K)+W(I,J,K+1))/2.0
            TAU11C(I,J,K)=TAU11(I,J,K)
            TAU22C(I,J,K)=TAU22(I,J,K)
            TAU33C(I,J,K)=TAU33(I,J,K)
            TAU12C(I,J,K)=(TAU12(I,J,K)+TAU12(I+1,J,K)+TAU12(I,J+1,K)+TAU12(I+1,J+1,K))/4.0
            TAU13C(I,J,K)=(TAU13(I,J,K)+TAU13(I+1,J,K)+TAU13(I,J,K+1)+TAU13(I+1,J,K+1))/4.0
            TAU23C(I,J,K)=(TAU23(I,J,K)+TAU23(I,J+1,K)+TAU23(I,J,K+1)+TAU23(I,J+1,K+1))/4.0
          END DO
        END DO
      END DO
    END IF
!---OBTAIN SOME CONTROL PARAMETERS
    IF(N.EQ.NT.OR.ISTOP.EQ.1.OR.TIME.GE.TIME_MAX)THEN
      I_END=1   ! FLAG IF THE RUN IS ABOUT TO END
    ELSE
      I_END=0
    END IF

    IF(MYID.EQ.0)THEN
      OPEN(1,FILE="post_process.echo")
      WRITE(1,*)'ISTOP =',ISTOP
      WRITE(1,*)'I_END =',I_END
      WRITE(1,*)
      WRITE(1,*)'N =',N
      WRITE(1,*)'NT =',NT
      WRITE(1,*)'N_1D =', N_1D
      WRITE(1,*)'N_3D =', N_3D
      WRITE(1,*)'N_RE =', N_RE  
      WRITE(1,*)
      WRITE(1,*)'I_3D_PRIME =', I_3D_PRIME 
      WRITE(1,*)'I_3D_TAU =', I_3D_TAU                     
      WRITE(1,*)'I_3D_SGS =', I_3D_SGS
      WRITE(1,*)'I_RE_PRIME =', I_RE_PRIME 
      WRITE(1,*)'I_RE_SGS =', I_RE_SGS 
      WRITE(1,*)
      WRITE(1,*)'STAT_FLGA =', STAT_FLAG
      WRITE(1,*)'ORDER_STAT =', ORDER_STAT
      WRITE(1,*)'N_SKIP =', N_SKIP 
      WRITE(1,*)'T_STAT_1 =', T_STAT_1 
      WRITE(1,*)'T_STAT_2 =', T_STAT_2
      CLOSE(1)    
    END IF
!---FOR 1D TEMPORAL OUTPUT
    IF(MOD(N,N_1D).EQ.0)THEN
      CALL OUTPUT_WALL()
    END IF   
!---FOR 3D FIELD OUTPUT-------------------------------------------------------
    IF(I_END.EQ.1.OR.MOD(N,N_3D).EQ.0)THEN
      IF(I_3D_PRIME.EQ.1)THEN    !  export prime variables
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for prime variables.'
        END IF
        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,6),NAME_VAR(6))
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              VAR(I,J,K,1)=UC(I,J,K)
              VAR(I,J,K,2)=VC(I,J,K)
              VAR(I,J,K,3)=WC(I,J,K) 
              VAR(I,J,K,4)=TE(I,J,K)
              VAR(I,J,K,5)=PD(I,J,K)
              VAR(I,J,K,6)=RHO(I,J,K)
            END DO
          END DO
        END DO

        NAME_VAR(1)="U"
        NAME_VAR(2)="V"
        NAME_VAR(3)="W"
        NAME_VAR(4)="TE"
        NAME_VAR(5)="PD"
        NAME_VAR(6)="RHO"
        FNAME="PRIM"        
        CALL OUTPUT_3D(6,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,0)
        DEALLOCATE(VAR,NAME_VAR)
      END IF

      IF(I_3D_TAU.EQ.1)THEN  ! export stress tensor
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for stress tensor.'
        END IF
        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,6),NAME_VAR(6))
        VAR(:,:,:,1)=TAU11C
        VAR(:,:,:,2)=TAU22C
        VAR(:,:,:,3)=TAU33C 
        VAR(:,:,:,4)=TAU12C
        VAR(:,:,:,5)=TAU13C
        VAR(:,:,:,6)=TAU23C
        NAME_VAR(1)="T11"
        NAME_VAR(2)="T22"
        NAME_VAR(3)="T33"
        NAME_VAR(4)="T12"
        NAME_VAR(5)="T13"
        NAME_VAR(6)="T23"
        FNAME="STES"        
        CALL OUTPUT_3D(6,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,0)
        DEALLOCATE(VAR,NAME_VAR)
      END IF

      IF(I_3D_SGS.EQ.1.AND.ITYPE.EQ.1)THEN  ! export sgs variables
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for SGS variables.'
        END IF
        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,3),NAME_VAR(3))
        VAR(:,:,:,1)=CS2
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
             VAR(I,J,K,2)=NU(I,J,K)
            END DO
          END DO
        END DO
        VAR(:,:,:,3)=DISSIP
        NAME_VAR(1)="CS2"
        NAME_VAR(2)="NU"
        NAME_VAR(3)="DISSIP"
        FNAME="SGS"        
        CALL OUTPUT_3D(3,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,0)
        DEALLOCATE(VAR,NAME_VAR)
      END IF
    END IF
!---FOR RESTART FILES-----------------------------------------------
    IF(I_END.EQ.1.OR.MOD(N,N_RE).EQ.0)THEN
      IF(I_RE_PRIME.EQ.1)THEN    !  export prime variables
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Write restart file for prime variables.'
        END IF
        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,5),NAME_VAR(5))
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              VAR(I,J,K,1)=UC(I,J,K)
              VAR(I,J,K,2)=VC(I,J,K)
              VAR(I,J,K,3)=WC(I,J,K)
              VAR(I,J,K,4)=TE(I,J,K)
              VAR(I,J,K,5)=RHO(I,J,K)
            END DO
          END DO
        END DO
        NAME_VAR(1)="U"
        NAME_VAR(2)="V"
        NAME_VAR(3)="W"
        NAME_VAR(4)="TE"
        NAME_VAR(5)="RHO"
        FNAME="PRIM"        
        CALL OUTPUT_3D(5,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,1)
        DEALLOCATE(VAR,NAME_VAR)
      END IF

      IF(I_RE_SGS.EQ.1.AND.ITYPE.EQ.1.AND.(ISGS.EQ.2.OR.ISGS.EQ.3))THEN  ! export SGS information for Lagrangian-type models
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Write restart file for SGS variables.'
        END IF

        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,4),NAME_VAR(4))
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              VAR(I,J,K,1)=PLM(I,J,K)
              VAR(I,J,K,2)=PMM(I,J,K)
              VAR(I,J,K,3)=PQN(I,J,K)
              VAR(I,J,K,4)=PNN(I,J,K)
            END DO
          END DO
        END DO
        NAME_VAR(1)="PLM"
        NAME_VAR(2)="PMM"
        NAME_VAR(3)="PQN"
        NAME_VAR(4)="PNN"
        FNAME="SGS"        
        CALL OUTPUT_3D(4,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,1)
        DEALLOCATE(VAR,NAME_VAR)
      END IF
    END IF
!---FOR STATISTICS-----------------------------------------------
    IF(STAT_FLAG.EQ.1)THEN
       IF(TIME.GE.T_STAT_1.AND.MOD(N,N_SKIP).EQ.0)THEN
        OPEN(1,FILE="TIME_AVRAGE_STEP.DAT")
        READ(1,*)N_AVE_1
        CLOSE(1)

        ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,6),VARA(NX1:NX2,NY1:NY2,NZ1:NZ2,6),NAME_VAR(6))
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              VAR(I,J,K,1)=UC(I,J,K)
              VAR(I,J,K,2)=VC(I,J,K)
              VAR(I,J,K,3)=WC(I,J,K) 
              VAR(I,J,K,4)=TE(I,J,K)
              VAR(I,J,K,5)=PD(I,J,K)
              VAR(I,J,K,6)=RHO(I,J,K)
              VARA(I,J,K,1)=UA(I,J,K)
              VARA(I,J,K,2)=VA(I,J,K)
              VARA(I,J,K,3)=WA(I,J,K)
              VARA(I,J,K,4)=TEA(I,J,K)
              VARA(I,I,K,5)=PDA(I,J,K)
              VARA(I,J,K,6)=RHOA(I,J,K) 
            END DO
          END DO
        END DO       
        NAME_VAR(1)="U"
        NAME_VAR(2)="V"
        NAME_VAR(3)="W"
        NAME_VAR(4)="TE"
        NAME_VAR(5)="PD"
        NAME_VAR(6)="RHO"
        FNAME="PRIME"    
       
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Perform time average for prime variables.'
        END IF

        CALL TIME_AVERAGE(6,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,VARA,N_AVE_1,T_STAT_1)

        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              UA(I,J,K)  = VARA(I,J,K,1)
              VA(I,J,K)  = VARA(I,J,K,2)
              WA(I,J,K)  = VARA(I,J,K,3) 
              TEA(I,J,K) = VARA(I,J,K,4)
              PDA(I,J,K) = VARA(I,J,K,5)
              RHOA(I,J,K)= VARA(I,J,K,6)
            END DO
          END DO
        END DO

        DEALLOCATE(VAR,VARA,NAME_VAR)
       
        N_AVE_1=N_AVE_1+1

        OPEN(1,FILE="TIME_AVRAGE_STEP.DAT")
        WRITE(1,*)N_AVE_1
        WRITE(1,*)N_AVE_2
        CLOSE(1)
      END IF

      IF(ORDER_STAT.EQ.2)THEN
        IF(TIME.GE.T_STAT_2.AND.MOD(N,N_SKIP).EQ.0)THEN   
          OPEN(1,FILE="TIME_AVRAGE_STEP.DAT")
          READ(1,*)N_AVE_1
          READ(1,*)N_AVE_2
          CLOSE(1)
          ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,12),VARA(NX1:NX2,NY1:NY2,NZ1:NZ2,12),NAME_VAR(12))
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                VAR(I,J,K,1)=(UC(I,J,K)-UA(I,J,K))**2
                VAR(I,J,K,2)=(VC(I,J,K)-VA(I,J,K))**2
                VAR(I,J,K,3)=(WC(I,J,K)-WA(I,J,K))**2 
                VAR(I,J,K,4)=(UC(I,J,K)-UA(I,J,K))*(VC(I,J,K)-VA(I,J,K))
                VAR(I,J,K,5)=(UC(I,J,K)-UA(I,J,K))*(WC(I,J,K)-WA(I,J,K))
                VAR(I,J,K,6)=(VC(I,J,K)-VA(I,J,K))*(WC(I,J,K)-WA(I,J,K))
                VAR(I,J,K,7)=(UC(I,J,K)-UA(I,J,K))*(TE(I,J,K)-TEA(I,J,K))
                VAR(I,J,K,8)=(VC(I,J,K)-VA(I,J,K))*(TE(I,J,K)-TEA(I,J,K))
                VAR(I,J,K,9)=(WC(I,J,K)-WA(I,J,K))*(TE(I,J,K)-TEA(I,J,K))
                VAR(I,J,K,10)=(UC(I,J,K)-UA(I,J,K))*(PD(I,J,K)-PDA(I,J,K))
                VAR(I,J,K,11)=(VC(I,J,K)-VA(I,J,K))*(PD(I,J,K)-PDA(I,J,K))
                VAR(I,J,K,12)=(WC(I,J,K)-WA(I,J,K))*(PD(I,J,K)-PDA(I,J,K))

                VARA(I,J,K,1)=UUA(I,J,K)
                VARA(I,J,K,2)=VVA(I,J,K)
                VARA(I,J,K,3)=WWA(I,J,K)
                VARA(I,J,K,4)=UVA(I,J,K)
                VARA(I,J,K,5)=UWA(I,J,K)
                VARA(I,J,K,6)=VWA(I,J,K)
                VARA(I,J,K,7)=UTA(I,J,K)
                VARA(I,J,K,8)=VTA(I,J,K)
                VARA(I,J,K,9)=WTA(I,J,K)
                VARA(I,J,K,10)=UPA(I,J,K)
                VARA(I,J,K,11)=VPA(I,J,K)
                VARA(I,J,K,12)=WPA(I,J,K)
              END DO
            END DO
          END DO
        
          NAME_VAR(1)="UU"
          NAME_VAR(2)="VV"
          NAME_VAR(3)="WW"      
          NAME_VAR(4)="UV"
          NAME_VAR(5)="UW"
          NAME_VAR(6)="VW"
          NAME_VAR(7)="UT"
          NAME_VAR(8)="VT"
          NAME_VAR(9)="WT" 
          NAME_VAR(10)="UP"
          NAME_VAR(11)="VP"
          NAME_VAR(12)="WP"
          FNAME="VARI"   
       
          IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
            PRINT*,'Perform time average for variances.'
          END IF
 
          CALL TIME_AVERAGE(12,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,VARA,N_AVE_2,T_STAT_2)

          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                UUA(I,J,K) = VARA(I,J,K,1)
                VVA(I,J,K) = VARA(I,J,K,2)
                WWA(I,J,K) = VARA(I,J,K,3) 
                UVA(I,J,K) = VARA(I,J,K,4)
                UWA(I,J,K) = VARA(I,J,K,5)
                VWA(I,J,K) = VARA(I,J,K,6)
                UTA(I,J,K) = VARA(I,J,K,7)
                VTA(I,J,K) = VARA(I,J,K,8)
                WTA(I,J,K) = VARA(I,J,K,9) 
                UPA(I,J,K) = VARA(I,J,K,10)
                VPA(I,J,K) = VARA(I,J,K,11)
                WPA(I,J,K) = VARA(I,J,K,12)
              END DO
            END DO
          END DO

          DEALLOCATE(VAR,VARA,NAME_VAR)
 
          IF(ITYPE.EQ.1)THEN    ! FOR LES, GET SGS CONTRIBUTIONS
            ALLOCATE(VAR(NX1:NX2,NY1:NY2,NZ1:NZ2,8),VARA(NX1:NX2,NY1:NY2,NZ1:NZ2,8),NAME_VAR(8))
            DO I=1,NX
              DO J=1,NY
                DO K=1,NZ
                  VAR(I,J,K,1)=TAU11C(I,J,K)
                  VAR(I,J,K,2)=TAU22C(I,J,K)
                  VAR(I,J,K,3)=TAU33C(I,J,K)
                  VAR(I,J,K,4)=TAU12C(I,J,K)
                  VAR(I,J,K,5)=TAU13C(I,J,K)
                  VAR(I,J,K,6)=TAU23C(I,J,K)
                  VAR(I,J,K,7)=NU(I,J,K)
                  VAR(I,J,K,8)=DISSIP(I,J,K)

                  VARA(I,J,K,1)=TAU11A(I,J,K)
                  VARA(I,J,K,2)=TAU22A(I,J,K)
                  VARA(I,J,K,3)=TAU33A(I,J,K)
                  VARA(I,J,K,4)=TAU12A(I,J,K)
                  VARA(I,J,K,5)=TAU13A(I,J,K)
                  VARA(I,J,K,6)=TAU23A(I,J,K)
                  VARA(I,J,K,7)=NUA(I,J,K)
                  VARA(I,J,K,8)=DISA(I,J,K)
                END DO
              END DO
            END DO     

            NAME_VAR(1)="T11"
            NAME_VAR(2)="T22"
            NAME_VAR(3)="T33"
            NAME_VAR(4)="T12"
            NAME_VAR(5)="T13"
            NAME_VAR(6)="T23"
            NAME_VAR(7)="NU"
            NAME_VAR(8)="DISSIP"
            FNAME="SGS"   

            IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
              PRINT*,'Perform time average for SGS variables.'
            END IF

        
            CALL TIME_AVERAGE(8,FNAME,NAME_VAR,VAR,NX1,NY1,NZ1,VARA,N_AVE_2,T_STAT_2)
            DO I=1,NX
              DO J=1,NY
                DO K=1,NZ  
                  TAU11A(I,J,K) = VARA(I,J,K,1)
                  TAU22A(I,J,K) = VARA(I,J,K,2)
                  TAU33A(I,J,K) = VARA(I,J,K,3) 
                  TAU12A(I,J,K) = VARA(I,J,K,4)
                  TAU13A(I,J,K) = VARA(I,J,K,5)
                  TAU23A(I,J,K) = VARA(I,J,K,6)
                  NUA(I,J,K)    = VARA(I,J,K,7)
                  DISA(I,J,K)   = VARA(I,J,K,8)
                END DO
              END DO
            END DO

            DEALLOCATE(VAR,VARA,NAME_VAR)
          END IF
       
          N_AVE_2=N_AVE_2+1

          OPEN(1,FILE="TIME_AVRAGE_STEP.DAT")
          WRITE(1,*)N_AVE_1
          WRITE(1,*)N_AVE_2
          CLOSE(1)
        END IF
      END IF
    END IF
  
    DEALLOCATE(UC,VC,WC,TAU11C,TAU22C,TAU33C,TAU12C,TAU13C,TAU23C)

    END SUBROUTINE
!=========================================================================!
!                  TIME AVERAGE OF PRIME VARIABLES                        !
!=========================================================================!
!   NUM_VAR:  NUMBER OF VARIABLES
!   FNAME:    KEY STRING OF THE FILE NAME, E.G. "VELO", "SGS", ETC. (<= 4 CHARACTERS)
!   NAME_VAR: VARIABLE NAME (<= 4 CHARACTERS). IT IS A 1-D ARRAY.
!   VAR:      VARIABLE DATA THAT WILL BE EXPORTED. IT IS A 4-D ARRAY.
!   N_TAVE:   NUMBER OF TOTAL SAMPLES
!   T_START:  START TIME
    SUBROUTINE TIME_AVERAGE(NUM_VAR,FNAME,NAME_VAR,VAR,SI1,SI2,SI3,VARA,N_TAVE,T_START)

    IMPLICIT NONE

    INTEGER :: SI1,SI2,SI3
    CHARACTER*10,DIMENSION(:):: NAME_VAR
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:,:):: VAR
    REAL(KIND=DP),DIMENSION(:,:,:,:):: VARA

    REAL(KIND=DP),DIMENSION(:,:,:,:),ALLOCATABLE:: TR
    REAL(KIND=DP) :: T_START,DUM

    CHARACTER*10 :: STR1,FNAME
    CHARACTER*10 :: STR2
    CHARACTER*20 :: FILENAME
    INTEGER NUM_VAR,N_TAVE,I,J,K,M,SAVE_RESTART
  
    IF(N.EQ.NSTART.AND.TIME.NE.T_START)THEN
      ALLOCATE(TR(NXT,NYT,NZT,NUM_VAR))
      STR1='_TAVE.DAT'
      FILENAME=TRIM(FNAME) // TRIM(STR1)

      OPEN(UNIT=1,FILE=FILENAME)   
      READ(1,*)
      READ(1,*) 
      DO J=1,NYT
        DO K=1,NZT
          DO I=1,NXT
            READ(1,*)DUM,DUM,DUM,(TR(I,J,K,M),M=1,NUM_VAR)
          ENDDO
        ENDDO
      ENDDO  
      CLOSE(1)

      DO J=1,NY
        DO K=1,NZ
          DO I=1,NX
            DO M=1,NUM_VAR
              VARA(I,J,K,M) = TR(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ,M)
            END DO
          END DO
        ENDDO
      ENDDO
      DEALLOCATE(TR)
    END IF
    
    VARA=(VARA*N_TAVE+VAR)/(N_TAVE+1)

    IF(I_END.EQ.1)THEN
      CALL OUTPUT_3D(NUM_VAR,FNAME,NAME_VAR,VARA,SI1,SI2,SI3,0,  &
                     FILENAME)
    END IF

    END SUBROUTINE

  END MODULE
