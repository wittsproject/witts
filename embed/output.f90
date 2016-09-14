! This module is for exporting output files
!
  MODULE output

  USE mpi  
  USE parameters
  USE class_shared
  USE tools
                          
  IMPLICIT NONE
  
  CONTAINS
!=========================================================================!
!                EXPORT FIELD DATA OF THE FLOW VARIABLES                  !
!=========================================================================!
!   NUM_VAR:  NUMBER OF VARIABLES
!   FNAME:    KEY STRING OF THE FILE NAME, E.G. "VELO", "SGS", ETC. (<= 4 CHARACTERS)
!   NAME_VAR: VARIABLE NAME (<= 4 CHARACTERS). IT IS A 1-D ARRAY.
!   VAR:      VARIABLE DATA THAT WILL BE EXPORTED. IT IS A 2-D ARRAY.    
!   SAVE_RESTART: 
!             =1: ONLY SAVE THE INPUT DATA FOR RESTART USE
!             =0: ONLY EXPORT THE INPUT DATA IN *.DAT FILE
!             =2: EXPORT THE INPUT DATA FOR BOTH RESTART AND *.DAT FILE
!   NUM_CELL: NUMBER OF CELL ON EACH RANK    
!   FNAME_IN: IF IT IS PASSED, THEN IT WILL BE USED AS OUTPUT FILENAME
!             OTHERWISE, THE FILENAME WILL BE GIVEN AS DEFAULT.
!
!   GLOBAL VARIABLES USED: N,TIME,NX,NY,NZ,NXT,NYT,NZT,MYID,XI,YI,ZI
!
!   THIS SUBROUTINE WILL USE THOSE INPUT TO GENERATE FOLLOWING FILES CONTAINING
!   3-D FIELD DATA:
!            N_FNAME.DAT (IN THE TECPLOT FORMAT)
!   OR/AND   RESTART_FNAME (FOR RESTART USE)
    SUBROUTINE OUTPUT_3D(NUM_VAR,FNAME,NAME_VAR,VAR,SAVE_RESTART,NUM_CELL,  &
                         FNAME_IN)   !  THIS LINE FOR OPTIONAL
 
    IMPLICIT NONE

    INTEGER :: SI1,SI2,SI3
    CHARACTER*10,DIMENSION(:):: NAME_VAR
    REAL(KIND=DP),DIMENSION(:,:) :: VAR
    REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE:: TR

    CHARACTER*10 :: STR1,FNAME
    CHARACTER*10 :: STR2,STR3
    CHARACTER*20 :: FILENAME
    CHARACTER*20, OPTIONAL :: FNAME_IN
    INTEGER:: NUM_VAR,NUM_CELL,NUM_CELL_TOTAL
    INTEGER:: I,J,K,M,SAVE_RESTART

    CALL MPI_ALLREDUCE(NUM_CELL,NUM_CELL_TOTAL,1,MPI_INTEGER, &
                       MPI_SUM,MPI_COMM_WORLD,IERR)
    
    ALLOCATE(TR(NUM_CELL_TOTAL,NUM_VAR))
!---SAVE VELOCITY AND SCALAR
    DO M=1,NUM_VAR
      CALL ASSEM(VAR(:,M),TR(:,M),NUM_CELL)
    END DO

    IF(PRESENT(FNAME_IN))THEN
      FILENAME=FNAME_IN
    ELSE          
      STR1='.DAT'
      WRITE(STR2,'(I6.6,A1)')N,'_'
      FILENAME=TRIM(STR2) // TRIM(FNAME) // TRIM(STR1)
      IF(SAVE_RESTART.EQ.1.OR.SAVE_RESTART.EQ.3)THEN
        STR2='RESTART_'
        WRITE(STR3,'(I6.6)')N
        FILENAME=TRIM(STR2) // TRIM(FNAME) // TRIM(STR3)
      END IF
    END IF

!   OUTPUT DOMAIN DATA
    IF(MYID.EQ.0)THEN
      IF(SAVE_RESTART.EQ.1.OR.SAVE_RESTART.EQ.3)THEN    !  SAVE THE DOMAIN DATA FOR RESTART
        OPEN(UNIT=1,FILE=FILENAME)  
        WRITE(1,*)N,TIME
        DO I=1,NUM_CELL_TOTAL
          WRITE(1,*)(TR(I,M),M=1,NUM_VAR)
        ENDDO
        CLOSE(1)
      ELSE IF(SAVE_RESTART.EQ.0.OR.SAVE_RESTART.EQ.3)THEN  !  EXPORT 3D FIELD
        OPEN(UNIT=1,FILE=FILENAME)  
        WRITE(1,*)'VARIABLES=',('"',TRIM(NAME_VAR(I)),'",',I=1,NUM_VAR-1), &
                  '"',TRIM(NAME_VAR(NUM_VAR)),'"'
        WRITE(1,*)'ZONE T="',TIME,'", F=POINT'
        DO I=1,NUM_CELL_TOTAL
          WRITE(1,*)(TR(I,M),M=1,NUM_VAR)
        ENDDO
      CLOSE(1)
      END IF        
    END IF

    DEALLOCATE(TR)

    END SUBROUTINE OUTPUT_3D

!=========================================================================!
!                              WALL MODEL                                 !
!=========================================================================!
    SUBROUTINE OUTPUT_WALL()

    IMPLICIT NONE
    REAL(KIND=DP) :: TAU_W,Q_W,TE_W,VEL_W,RHO_W
    REAL(KIND=DP) :: TAU_W_GLOBAL,Q_W_GLOBAL,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL
    INTEGER :: M
    LOGICAL :: EXIST
!---AT BOUNDARY #1----------------------------------------------------------    
    IF(IWALL(1).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
      
      IF(MYIDX.EQ.0)THEN    
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_X.LT.CELL_FV(M)%CELL_DX)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(4)**2+ &
                             CELL_FV(M)%CELL_TAU(5)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(1)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(2)**2+  &
                             CELL_FV(M)%CELL_VEL(3)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF
        END DO
        TAU_W=TAU_W/(NY+2)/(NZ+2)
        Q_W=Q_W/(NY+2)/(NZ+2)
        TE_W=TE_W/(NY+2)/(NZ+2)
        VEL_W=VEL_W/(NY+2)/(NZ+2)
        RHO_W=RHO_W/(NY+2)/(NZ+2)
      END IF
     
      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPY*NPZ)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPY*NPZ)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPY*NPZ)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPY*NPZ)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPY*NPZ)
      
      IF(MYID.EQ.0)THEN
        INQUIRE(FILE="WALL_1.DAT", EXIST=EXIST)
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_1.DAT',ACCESS="APPEND")
        ELSE
          OPEN(1,FILE="WALL_1.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'           
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)
      END IF         
    END IF
!---AT BOUNDARY #2----------------------------------------------------------       
    IF(IWALL(2).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
      
      IF(MYIDX.EQ.NPX-1)THEN
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_X.GT.LX-CELL_FV(M)%CELL_DX)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(4)**2+ &
                             CELL_FV(M)%CELL_TAU(5)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(1)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(2)**2+  &
                             CELL_FV(M)%CELL_VEL(3)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF
        END DO        
       
        TAU_W=TAU_W/(NY+2)/(NZ+2)
        Q_W=Q_W/(NY+2)/(NZ+2)
        TE_W=TE_W/(NY+2)/(NZ+2)
        VEL_W=VEL_W/(NY+2)/(NZ+2)
        RHO_W=RHO_W/(NY+2)/(NZ+2)       
      END IF

      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPY*NPZ)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPY*NPZ)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPY*NPZ)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPY*NPZ)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPY*NPZ)
      
      IF(MYID.EQ.0)THEN
        INQUIRE(FILE="WALL_2.DAT", EXIST=EXIST)
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_2.DAT',ACCESS="APPEND")
        ELSE
          OPEN(1,FILE="WALL_2.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'           
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)        
      END IF
    END IF     
!---AT BOUNDARY #3---------------------------------------------------------- 
    IF(IWALL(3).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
      
      IF(MYIDY.EQ.0)THEN
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_Y.LT.CELL_FV(M)%CELL_DY)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(4)**2+ &
                             CELL_FV(M)%CELL_TAU(6)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(2)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(1)**2+  &
                             CELL_FV(M)%CELL_VEL(3)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF   
        END DO 
        
        TAU_W=TAU_W/(NX+2)/(NZ+2)
        Q_W=Q_W/(NX+2)/(NZ+2)
        TE_W=TE_W/(NX+2)/(NZ+2)
        VEL_W=VEL_W/(NX+2)/(NZ+2)
        RHO_W=RHO_W/(NX+2)/(NZ+2)
      END IF
     
      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPX*NPZ)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPX*NPZ)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPX*NPZ)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPX*NPZ)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPX*NPZ)
      
      IF(MYID.EQ.0)THEN
        INQUIRE(FILE="WALL_3.DAT", EXIST=EXIST)
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_3.DAT',ACCESS="APPEND")
        ELSE
          OPEN(1,FILE="WALL_3.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'           
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)        
      END IF
    END IF     
!---AT BOUNDARY #4---------------------------------------------------------- 
    IF(IWALL(4).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
      
      IF(MYIDY.EQ.NPY-1)THEN
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_Y.GT.LY+CELL_FV(M)%CELL_DY)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(4)**2+ &
                             CELL_FV(M)%CELL_TAU(6)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(2)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(1)**2+  &
                             CELL_FV(M)%CELL_VEL(3)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF   
        END DO
       
        TAU_W=TAU_W/(NX+2)/(NZ+2)
        Q_W=Q_W/(NX+2)/(NZ+2)
        TE_W=TE_W/(NX+2)/(NZ+2)
        VEL_W=VEL_W/(NX+2)/(NZ+2)
        RHO_W=RHO_W/(NX+2)/(NZ+2)       
      END IF

      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPX*NPZ)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPX*NPZ)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPX*NPZ)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPX*NPZ)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPX*NPZ)
      
      IF(MYID.EQ.0)THEN
        INQUIRE(FILE="WALL_4.DAT", EXIST=EXIST)
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_4.DAT',ACCESS="APPEND")           
        ELSE
          OPEN(1,FILE="WALL_4.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)        
      END IF
    END IF     
!---AT BOUNDARY #5---------------------------------------------------------- 
    IF(IWALL(5).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
      
      IF(MYIDZ.EQ.0)THEN
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_Z.LT.CELL_FV(M)%CELL_DZ)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(5)**2+ &
                             CELL_FV(M)%CELL_TAU(6)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(3)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(1)**2+  &
                             CELL_FV(M)%CELL_VEL(2)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF   
        END DO
        
        TAU_W=TAU_W/(NX+2)/(NY+2)
        Q_W=Q_W/(NX+2)/(NY+2)
        TE_W=TE_W/(NX+2)/(NY+2)
        VEL_W=VEL_W/(NX+2)/(NY+2)
        RHO_W=RHO_W/(NX+2)/(NY+2)
      END IF

      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPX*NPY)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPX*NPY)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPX*NPY)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPX*NPY)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPX*NPY)
      
      IF(MYID.EQ.0)THEN
         INQUIRE(FILE="WALL_5.DAT", EXIST=EXIST)         
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_5.DAT',ACCESS="APPEND")
        ELSE
          OPEN(1,FILE="WALL_5.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'           
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)        
      END IF     
    END IF
!---AT BOUNDARY #6----------------------------------------------------------  
    IF(IWALL(6).EQ.1)THEN
      TAU_W=0.0
      Q_W=0.0
      TE_W=0.0
      VEL_W=0.0
      RHO_W=0.0
       
      IF(MYIDZ.EQ.NPZ-1)THEN
        DO M=1,TOTAL_CELL
          IF(CELL_FV(M)%CELL_ACTIVE.EQ.1.AND.  &
             CELL_FV(M)%CELL_Z.GT.LZ+CELL_FV(M)%CELL_DZ)THEN    
            TAU_W=TAU_W+SQRT(CELL_FV(M)%CELL_TAU(5)**2+ &
                             CELL_FV(M)%CELL_TAU(6)**2)
            Q_W=Q_W+CELL_FV(M)%CELL_HF(3)
            TE_W=TE_W+CELL_FV(M)%CELL_TE
            VEL_W=VEL_W+SQRT(CELL_FV(M)%CELL_VEL(1)**2+  &
                             CELL_FV(M)%CELL_VEL(2)**2)
            RHO_W=RHO_W+CELL_FV(M)%CELL_RHO
          END IF   
        END DO
        
        TAU_W=TAU_W/(NX+2)/(NY+2)
        Q_W=Q_W/(NX+2)/(NY+2)
        TE_W=TE_W/(NX+2)/(NY+2)
        VEL_W=VEL_W/(NX+2)/(NY+2)
        RHO_W=RHO_W/(NX+2)/(NY+2)
      END IF

      CALL MPI_ALLREDUCE(TAU_W,TAU_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(Q_W,Q_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(TE_W,TE_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(VEL_W,VEL_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(RHO_W,RHO_W_GLOBAL,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)       
      TAU_W_GLOBAL=TAU_W_GLOBAL/(NPX*NPY)
      Q_W_GLOBAL=Q_W_GLOBAL/(NPX*NPY)
      TE_W_GLOBAL=TE_W_GLOBAL/(NPX*NPY)
      VEL_W_GLOBAL=VEL_W_GLOBAL/(NPX*NPY)
      RHO_W_GLOBAL=RHO_W_GLOBAL/(NPX*NPY)
      
      IF(MYID.EQ.0)THEN
        INQUIRE(FILE="WALL_6.DAT", EXIST=EXIST)
        IF(EXIST)THEN
          OPEN(1,FILE='WALL_6.DAT',ACCESS="APPEND")           
        ELSE
          OPEN(1,FILE="WALL_6.DAT")
          WRITE(1,*)'VARIABLES="Time","Te_near","Vel_Tan_near","Density_near","Tau_wall","U*","Q_wall","Obukhov_Length"'
        END IF
        WRITE(1,*)TIME,TE_W_GLOBAL,VEL_W_GLOBAL,RHO_W_GLOBAL,TAU_W_GLOBAL,SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL),Q_W_GLOBAL, &
             -(SQRT(TAU_W_GLOBAL/RHO_W_GLOBAL)**3.0)*TE_W_GLOBAL/(KAPPA*G*Q_W_GLOBAL+1.0E-12)
        CLOSE(1)        
      END IF      
    END IF
   
    END SUBROUTINE OUTPUT_WALL    
  END MODULE
