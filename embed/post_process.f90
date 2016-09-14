! This module is for controling of the post process, including exporting //
! result files and statistics
!
  MODULE post_process

  USE mpi
  USE parameters
  USE class_shared   
  USE tools
  USE output

  CONTAINS
!=========================================================================!
!                  MAIN SUBROUTINE OF POST PROCESS                        !
!=========================================================================!
    SUBROUTINE POSTPROCESS()

    IMPLICIT NONE

    INTEGER :: I_3D_PRIME,I_3D_TAU,I_3D_HF,I_3D_SGS,N_1D,N_3D
    INTEGER :: I_RE_PRIME,I_RE_SGS,N_RE
    INTEGER :: STAT_FLAG,ORDER_STAT,N_SKIP
    REAL(KIND=DP) :: T_STAT_1,T_STAT_2 
      
    INTEGER :: I,J,K,N_AVE_1,N_AVE_2,M,MM
    REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE:: VAR,VARA
    CHARACTER*10,DIMENSION(:),ALLOCATABLE:: NAME_VAR
    CHARACTER*10 :: FNAME

    OPEN(1,FILE="post.in")
    READ(1,*) 
    READ(1,*) I_3D_PRIME   !(1: GENERATE 3D FIELD OF PRIME VARIABLES)
    READ(1,*) I_3D_TAU     !(1: GENERATE 3D FIELD OF STRESS TENSOR)
    READ(1,*) I_3D_HF      !(1: GENERATE 3D FIELD OF HEAT FLUX)
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
!-----Export prime variables       
      IF(I_3D_PRIME.EQ.1)THEN    
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for prime variables.'
        END IF
        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,9),NAME_VAR(9))

        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"
        NAME_VAR(4)="U"
        NAME_VAR(5)="V"
        NAME_VAR(6)="W"
        NAME_VAR(7)="TE"
        NAME_VAR(8)="PD"
        NAME_VAR(9)="PHI"
        FNAME="PRIM"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            VAR(MM,4)=CELL_FV(M)%CELL_VEL(1)
            VAR(MM,5)=CELL_FV(M)%CELL_VEL(2)
            VAR(MM,6)=CELL_FV(M)%CELL_VEL(3)
            VAR(MM,7)=CELL_FV(M)%CELL_TE
            VAR(MM,8)=CELL_FV(M)%CELL_PD
            VAR(MM,9)=CELL_FV(M)%CELL_PHI
          END IF
        END DO
       
        CALL OUTPUT_3D(9,FNAME,NAME_VAR,VAR,0,TOTAL_CELL_ACTIVE) 
        DEALLOCATE(VAR,NAME_VAR)
      END IF
     
!-----Export stress tensor
      IF(I_3D_TAU.EQ.1)THEN  
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for stress tensor.'
        END IF
        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,9),NAME_VAR(9))

        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"
        NAME_VAR(4)="T11"
        NAME_VAR(5)="T22"
        NAME_VAR(6)="T33"
        NAME_VAR(7)="T12"
        NAME_VAR(8)="T13"
        NAME_VAR(9)="T23"
        FNAME="STES"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            DO I=1,6
              VAR(MM,3+I)=CELL_FV(M)%CELL_TAU(I)
            END DO
          END IF
        END DO
   
        CALL OUTPUT_3D(9,FNAME,NAME_VAR,VAR,0,TOTAL_CELL_ACTIVE)        
        DEALLOCATE(VAR,NAME_VAR)
      END IF

!-----Export heat flux
      IF(I_3D_TAU.EQ.1)THEN  
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for heat flux.'
        END IF
        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,6),NAME_VAR(6))

        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"
        NAME_VAR(4)="Q1"
        NAME_VAR(5)="Q2"
        NAME_VAR(6)="Q3"
        FNAME="HFLUX"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            DO I=1,3
              VAR(MM,3+I)=CELL_FV(M)%CELL_HF(I)
            END DO
          END IF
        END DO

        CALL OUTPUT_3D(6,FNAME,NAME_VAR,VAR,0,TOTAL_CELL_ACTIVE)         
        DEALLOCATE(VAR,NAME_VAR)
      END IF
     
!-----Export sgs variables
      IF(I_3D_SGS.EQ.1.AND.ITYPE.EQ.1)THEN  
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Export 3d field data for SGS variables.'
        END IF
        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,5),NAME_VAR(5))
        
        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"
        NAME_VAR(4)="NU"
        NAME_VAR(5)="EPS"
        FNAME="SGS"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            VAR(MM,4)=CELL_FV(M)%CELL_NU
            VAR(MM,4)=CELL_FV(M)%CELL_EPS
          END IF
        END DO
       
        CALL OUTPUT_3D(5,FNAME,NAME_VAR,VAR,0,TOTAL_CELL_ACTIVE)
        DEALLOCATE(VAR,NAME_VAR)
      END IF
    END IF
!---FOR RESTART FILES-----------------------------------------------
    IF(I_END.EQ.1.OR.MOD(N,N_RE).EQ.0)THEN
      IF(I_RE_PRIME.EQ.1)THEN    !  export prime variables
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Write restart file for prime variables.'
        END IF
       
        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,8),NAME_VAR(8))

        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"        
        NAME_VAR(4)="U"
        NAME_VAR(5)="V"
        NAME_VAR(6)="W"
        NAME_VAR(7)="TE"
        NAME_VAR(8)="PHI"
        FNAME="PRIM"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            VAR(MM,4)=CELL_FV(M)%CELL_VEL(1)
            VAR(MM,5)=CELL_FV(M)%CELL_VEL(2)
            VAR(MM,6)=CELL_FV(M)%CELL_VEL(3)
            VAR(MM,7)=CELL_FV(M)%CELL_TE
            VAR(MM,8)=CELL_FV(M)%CELL_PHI
          END IF
        END DO
       
        CALL OUTPUT_3D(8,FNAME,NAME_VAR,VAR,1,TOTAL_CELL_ACTIVE)
        DEALLOCATE(VAR,NAME_VAR)
      END IF

      IF(I_RE_SGS.EQ.1.AND.ITYPE.EQ.1.AND.(ISGS.EQ.2.OR.ISGS.EQ.3))THEN  ! export SGS information for Lagrangian-type models
        IF(MYID.EQ.0.AND.SCREEN_LEVEL.GT.1)THEN
          PRINT*,'Write restart file for SGS variables.'
        END IF

        ALLOCATE(VAR(TOTAL_CELL_ACTIVE,7),NAME_VAR(7))
        
        NAME_VAR(1)="X"
        NAME_VAR(2)="Y"
        NAME_VAR(3)="Z"        
        NAME_VAR(4)="PLM"
        NAME_VAR(5)="PMM"
        NAME_VAR(6)="PQN"
        NAME_VAR(7)="PNN"
        FNAME="SGS"

        MM=0
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
            MM=MM+1
       
            VAR(MM,1)=CELL_FV(M)%CELL_X
            VAR(MM,2)=CELL_FV(M)%CELL_Y
            VAR(MM,3)=CELL_FV(M)%CELL_Z
            VAR(MM,4)=CELL_FV(M)%CELL_PLM
            VAR(MM,5)=CELL_FV(M)%CELL_PMM
            VAR(MM,6)=CELL_FV(M)%CELL_PQN
            VAR(MM,7)=CELL_FV(M)%CELL_PNN
          END IF
        END DO
       
        CALL OUTPUT_3D(7,FNAME,NAME_VAR,VAR,1,TOTAL_CELL_ACTIVE)
        DEALLOCATE(VAR,NAME_VAR)
      END IF
    END IF

    END SUBROUTINE

  END MODULE
