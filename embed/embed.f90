! This module declares the sizes of dynamic global arrays.
! Also, it will deallocate those arrays at the end of simulation.
  MODULE EMBED

  USE parameters
  USE field_shared
  USE cell

  IMPLICIT NONE

  CONTAINS
!-----------------------------------------------------------------!
!                                                                 !
!-----------------------------------------------------------------!
    SUBROUTINE EMBED_INITIAL()
    IMPLICIT NONE

    TYPE(CELL),DIMENSION(:),ALLOCATABLE:: CELL_EM
    INTEGER:: I,J,K,II,JJ,KK,M
    INTEGER:: EM_LEVEL_TOTAL
    INTEGER:: NUM_CELL_LIMIT
    INTEGER,DIMENSION(:),ALLOCATABLE:: EM_NUM
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: LX_EM_1,LX_EM_2,&
                                         LY_EM_1,LY_EM_2,&
                                         LZ_EM_1,LZ_EM_2

    ALLOCATE(CELL_EM(NUM_CELL_LIMIT))

    OPEN(1,FILE='embed.in')
    READ(1,*)EM_LEVEL_TOTAL
    READ(1,*)NUM_CELL_LIMIT

    ALLOCATE(EM_NUM(EM_LEVEL_TOTAL))
    ALLOCATE(LX_EM_1(EM_LEVEL_TOTAL,10),LX_EM_2(EM_LEVEL_TOTAL,10), &
             LY_EM_1(EM_LEVEL_TOTAL,10),LY_EM_2(EM_LEVEL_TOTAL,10), &  
             LZ_EM_1(EM_LEVEL_TOTAL,10),LZ_EM_2(EM_LEVEL_TOTAL,10))
    
    DO I=1,EM_LEVEL_TOTAL
      READ(1,*)
      READ(1,*) 
      READ(1,*)EM_NUM(I)

      EM_NUM(I)=MIN(EM_NUM(I),10) ! CUTOFF AT 10
      
      DO J=1,EM_NUM(I)
        READ(1,*) 
        READ(1,*)LX_EM_1(I,J)
        READ(1,*)LX_EM_2(I,J)
        READ(1,*)LY_EM_1(I,J)
        READ(1,*)LY_EM_2(I,J)
        READ(1,*)LZ_EM_1(I,J)
        READ(1,*)LZ_EM_2(I,J)
      END DO  
    END DO  
!---SETUP THE BASIC CELLS
    TOTAL_CELL=0
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LE.NUM_CELL_LIMIT)THEN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
            CELL_EM(TOTAL_CELL)%CELL_EMID=0
            CELL_EM(TOTAL_CELL)%CELL_GHOST=0
            CELL_EM(TOTAL_CELL)%CELL_SPLIT=0
            CELL_EM(TOTAL_CELL)%CELL_NUM_VAR=NUM_VAR

            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=0    
            CELL_EM(TOTAL_CELL)%CELL_PID=MYID
            CELL_EM(TOTAL_CELL)%CELL_NEAR=0
            
            CELL_EM(TOTAL_CELL)%CELL_PARENT=0
            CELL_EM(TOTAL_CELL)%CELL_CHILD=0

            CELL_EM(TOTAL_CELL)%CELL_DX=DX0
            CELL_EM(TOTAL_CELL)%CELL_DY=DY0
            CELL_EM(TOTAL_CELL)%CELL_DZ=DZ0

            CELL_EM(TOTAL_CELL)%CELL_X=XI(I+MYIDX*NX)
            CELL_EM(TOTAL_CELL)%CELL_Y=YI(J+MYIDY*NY)
            CELL_EM(TOTAL_CELL)%CELL_Z=ZI(K+MYIDZ*NZ)

            CELL_EM(TOTAL_CELL)%CELL_VAR(1)=U(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(2)=V(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(3)=W(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(4)=TE(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(5)=PD(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(6)=NU(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(7)=MU(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(8)=RHO(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VAR(9)=PHI(I,J,K)
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO
      END DO
    END DO
!---MAKE THE EMBEDDING
    DO I=1,EM_LEVEL_TOTAL
      DO J=1,EM_NUM(I)

        TOTAL_CELL0=TOTAL_CELL

        DO M=1,TOTAL_CELL0
          IF(CELL_EM(M)%CELL_X.GE.LX_EM_1(I,J).AND.CELL_EM(M)%CELL_X.LE.LX_EM_2(I,J).AND. &
             CELL_EM(M)%CELL_Y.GE.LY_EM_1(I,J).AND.CELL_EM(M)%CELL_Y.LE.LY_EM_2(I,J).AND. &
             CELL_EM(M)%CELL_Z.GE.LZ_EM_1(I,J).AND.CELL_EM(M)%CELL_Z.LE.LZ_EM_2(I,J).AND. &
             CELL_EM(M)%CELL_EMID.EQ.I-1)THEN  ! CHECK IF THE CELL IS QUALIFIED TO BE FURTHER REFINED

            DX=DX0/2.0**I
            DY=DY0/2.0**I
            DZ=DZ0/2.0**I

            CELL_EM(M)%CELL_SPLIT=1

            IF(TOTAL_CELL+8.GT.NUM_CELL_LIMIT)THEN
              EXIT
            END IF

            DO II=1,2
              DO JJ=1,2
                DO KK=1,2
                  TOTAL_CELL=TAOTAL_CELL+1
                  
                  CELL_EM(TOTAL_CELL)%CEL_INDEX=TOTAL_CELL
                  CELL_EM(TOTAL_CELL)%CELL_EMID=I
                  CELL_EM(TOTAL_CELL)%CELL_GHOST=0

                  CELL_EM(TOTAL_CELL)%CELL_PARENT=M
                  CELL_EM(M)%CELL_CHILD(II,JJ,KK)=TOTAL_CELL
                  
                  CELL_EM(TOTAL_CELL)%CELL_DX=DX
                  CELL_EM(TOTAL_CELL)%CELL_DY=DY
                  CELL_EM(TOTAL_CELL)%CELL_DZ=DZ

                  CELL_EM(TOTAL_CELL)%CELL_X=CELL_EM(M)%CELL_X-DX/2.0+DX*(II-1)
                  CELL_EM(TOTAL_CELL)%CELL_Y=CELL_EM(M)%CELL_Y-DY/2.0+DY*(II-1)
                  CELL_EM(TOTAL_CELL)%CELL_Z=CELL_EM(M)%CELL_Z-DZ/2.0+DZ*(II-1)

                  DO MM=1,NUM_VAR
                    CELL_EM(TOTAL_CELL)%CELL_VAR(MM)=CELL_EM(M)%CELL_VAR(MM)
                 END DO
                END DO
              END DO
            END DO            
          END IF
        END DO

      END DO         
    END DO
!---GENERATE GHOST CELLS
    CALL GENERATE_GHOST_CELL()
!---GET THE NEIGHBORS FOR EACH CELL
    DO M=1,TOTAL_CELL 
      CALL NEIGHBOR_INDEX(M) 
    END DO
!---GET VALUES ON THE GHOST CELLS
    DO M=1,TOTAL_CELL
      IF(CELL_EM(M)%CELL_GHOST.EQ.1)THEN
        DO II=1,NUM_VAR
          CALL GHOST_BOUNDARY(M,II)
        END DO
      END IF
    END DO
   
    DEALLOCATE(EM_NUM)
    DEALLOCATE(LX_EM_1,LX_EM_2,LY_EM_1,LY_EM_2,LZ_EM_1,LZ_EM_2)
    
    END SUBROUTINE EMBED_INITIAL
  
!-------------------------------------------------------------------!
!                       GENERATE GHOST CELLS                        !
!-------------------------------------------------------------------!    
    SUBROUTINE GENERATE_GHOST_CELL()
    IMPLICIT NONE
    TYPE(CELL),DIMENSION(:),ALLOCATABLE:: BUF_SEND,BUF_RECE
    INTEGER:: IUP,IDOWN,TOTAL_CELL0
    REAL(KIND=DP):: DIS,DIS0

    TOTAL_CELL0=TOTAL_CELL
    
    DO M=1,TOTAL_CELL0
!---X DIRECTION--------------------
      IF(CELL_EM(M)%CELL_NEI_X(1)=0)THEN
        IUP=MYID+1
        IDOWN=MYID-1          

        IF(MYIDX.EQ.NPX-1)THEN
          IUP=MYID-(NPX-1)
        END IF
        IF(MYIDX.EQ.0)THEN
          IDOWN=MYID+(NPX-1)
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_X(2))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_X(2))%CELL_NEI_X(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=I
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_EM(M)%CELL_NEI_X(2)=0)THEN
        IUP=MYID+1
        IDOWN=MYID-1          

        IF(MYIDX.EQ.NPX-1)THEN
          IUP=MYID-(NPX-1)
        END IF
        IF(MYIDX.EQ.0)THEN
          IDOWN=MYID+(NPX-1)
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_X(1))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_X(1))%CELL_NEI_X(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)        
      END IF
!---Y DIRECTION--------------------
      IF(CELL_EM(M)%CELL_NEI_Y(1)=0)THEN          
        IUP=MYID+NPX*NPZ
        IDOWN=MYID-NPX*NPZ

        IF(MYIDY.EQ.NPY-1)THEN
          IUP=MYIDX+MYIDZ*NPX
        END IF
        IF(MYIDY.EQ.0)THEN
          IDOWN=MYIDX+(NPY-1)*NPX*NPZ+MYIDZ*NPX
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_Y(2))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_Y(2))%CELL_NEI_Y(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_EM(M)%CELL_NEI_Y(2)=0)THEN
        IUP=MYID+NPX*NPZ
        IDOWN=MYID-NPX*NPZ

        IF(MYIDY.EQ.NPY-1)THEN
          IUP=MYIDX+MYIDZ*NPX
        END IF
        IF(MYIDY.EQ.0)THEN
          IDOWN=MYIDX+(NPY-1)*NPX*NPZ+MYIDZ*NPX
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_Y(1))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_Y(1))%CELL_NEI_Y(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)        
      END IF
!---Z DIRECTION--------------------
      IF(CELL_EM(M)%CELL_NEI_Y(1)=0)THEN          
        IUP=MYID+NPX
        IDOWN=MYID-NPX        

        IF(MYIDZ.EQ.NPZ-1)THEN
          IUP=MYIDX+MYIDY*NPX*NPZ
        END IF
        IF(MYIDZ.EQ.0)THEN
          IDOWN=MYIDX+MYIDY*NPX*NPZ+(NPZ-1)*NPX
        END IF
               
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_Z(2))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_Z(2))%CELL_NEI_Z(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_EM(M)%CELL_NEI_Z(2)=0)THEN
        IUP=MYID+NPX
        IDOWN=MYID-NPX        

        IF(MYIDZ.EQ.NPZ-1)THEN
          IUP=MYIDX+MYIDY*NPX*NPZ
        END IF
        IF(MYIDZ.EQ.0)THEN
          IDOWN=MYIDX+MYIDY*NPX*NPZ+(NPZ-1)*NPX
        END IF

        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_EM(M)
        BUF_SEND(2)=CELL_EM(CELL_EM(M)%CELL_NEI_Z(1))
        BUF_SEND(3)=CELL_EM(CELL_EM(CELL_EM(M)%CELL_NEI_Z(1))%CELL_NEI_Z(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_EM(TOTAL_CELL)=BUF_RECE(I)
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            CELL_EM(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_EM(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)        
      END IF
    END DO
!---FOR EACH GHOST CELL, GET THE INDEX OF THE NEAREST NON-GHOST CELL
    DO M=1,TOTAL_CELL
      IF(CELL_EM(M)%CELL_GHOST.EQ.1)THEN
        DIS0=1.0D6
        DO MM=1,TOTAL_CELL
          IF(CELL_EM(MM)%CELL_GHOST_EQ.0)THEN
            DIS=SQRT((CELL_EM(MM)%CELL_X-CELL_EM(M)%CELL_X)**2+ &
                     (CELL_EM(MM)%CELL_Y-CELL_EM(M)%CELL_Y)**2+ &
                     (CELL_EM(MM)%CELL_Z-CELL_EM(M)%CELL_Z)**2)
            IF(DIS.LT.DIS0)THEN
              DIS0=DIS
              CELL_EM(M)%CELL_NEAR=MM
            END IF
          END IF
        END DO
      END IF
    END DO          
 
    END SUBROUTINE GENERATE_GHOST_CELL
  
!-------------------------------------------------------------------!
!                  GET BC FOR AN INNER GHOST CELL                   !
!-------------------------------------------------------------------!    
    SUBROUTINE GETBC_GHOST_SEND(INDEX,NUM) 
    IMPLICIT NONE
    INTEGER:: INDEX,NUM
    INTEGER:: INDEX_SEND,PID_SEND,WIN,IERR

    INDEX_SEND=CELL_EM(INDEX)%CELL_INDEX_ORIGIN
    PID_SEND=CELL_EM(INDEX)%CELL_PID

    CALL MPI_GET(CELL_EM(INDEX)%CELL_VAR(NUM),1,MPI_DOUBLE_PRECISION, &
                 PID_SEND,CELL_EM(INDEX_SEND)%CELL_VAR(NUM),          &
                 1,MPI_DOUBLE_PRECISION,WIN,IERR)   
 
    END SUBROUTINE GETBC_GHOST_SEND
  
!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE NEIGHBOR_INDEX(INDEX)
    IMPLICIT NONE

    INTEGER:: INDEX
    INTEGER:: M

    DO M=1,TOTAL_CELL
      IF(M.NE.INDEX)THEN 
        DISX=CELL_EM(M)%CELL_X-CELL_EM(INDEX)%CELL_X
        DISY=CELL_EM(M)%CELL_Y-CELL_EM(INDEX)%CELL_Y
        DISZ=CELL_EM(M)%CELL_Z-CELL_EM(INDEX)%CELL_Z        
        IF(CELL_EM(INDEX)%CELL_EMID.EQ.CELL_EM(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT THE SAME LEVEL
          !  X DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*1.01.AND. &
             ABS(DISY).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
            IF(DISX.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_X(1)=M
              CELL_EM(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_X(2)=M
              CELL_EM(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF  
          !  Y DIRECTION
          IF(ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*1.01.AND. &
             ABS(DISX).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
            IF(DISY.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Y(1)=M
              CELL_EM(M)%CELL_NEI_Y(2)=INDEX               
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Y(2)=M
              CELL_EM(M)%CELL_NEI_Y(1)=INDEX                
            END IF
          END IF  
          !  Z DIRECTION
          IF(ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*1.01.AND. &
             ABS(DISX).LT.ZERO.AND.ABS(DISY).LT.ZERO)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Z(1)=M
              CELL_EM(M)%CELL_NEI_Z(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Z(2)=M
              CELL_EM(M)%CELL_NEI_Z(1)=INDEX 
            END IF
          END IF
        ELSE IF(CELL_EM(INDEX)%CELL_EMID+1.EQ.CELL_EM(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT ONE LEVEL ABOVE
           !  X DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.75*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.25*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.25*1.01)THEN
            IF(DISX.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_X(1)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_X(2)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF
           !  Y DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.75*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.25*1.01)THEN
            IF(DISY.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Y(1)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_Y(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Y(2)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_Y(1)=INDEX 
            END IF
          END IF
          !  Z DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.25*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.75*1.01)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Z(1)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_Z(2)=INDEX                
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Z(2)=CELL_EM(M)%CELL_PARENT
              CELL_EM(M)%CELL_NEI_Z(1)=INDEX  
            END IF
          END IF
        ELSE IF(CELL_EM(INDEX)%CELL_EMID-1.EQ.CELL_EM(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT ONE LEVEL BELOW
        !  X DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*1.5*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.5*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.5*1.01)THEN
            IF(DISX.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_X(1)=M
              CELL_EM(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_X(2)=M
              CELL_EM(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF
          !  Y DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.5*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*1.5*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.5*1.01)THEN
            IF(DISY.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Y(1)=M
              CELL_EM(M)%CELL_NEI_Y(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Y(2)=M
              CELL_EM(M)%CELL_NEI_Y(1)=INDEX 
            END IF
          END IF
          !  Z DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.5*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.5*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*1.5*1.01)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Z(1)=M
              CELL_EM(M)%CELL_NEI_Z(2)=INDEX
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Z(2)=M
              CELL_EM(M)%CELL_NEI_Z(1)=INDEX
            END IF
          END IF           
        END IF

      END IF
    END DO
    END SUBROUTINE NEIGHBOR_INDEX
!-------------------------------------------------------------------!
!                     SET BC ON THE GHOST CELLS                     !
!-------------------------------------------------------------------!
    SUBROUTINE GHOST_BOUNDARY(INDEX,NUM)
    IMPLICIT NONE

    INTEGER:: INDEX,BC_TYPE
    INTEGER:: NUM,ID_INNER
    INTEGER:: I,ID_L,ID_LL,ID_R,ID_RR

 
    IF(CELL_EM(INDEX)%CELL_GHOST.NE.1)THEN
      RETURN
    END IF

    IF(NUM.LT.1)THEN
      NUM=1
    END IF

    IF(NUM.GT.NUM_VAR)THEN
      NUM=NUM_VAR
    END IF
!---JUDGE IF THE CELL IS AN INNER GHOST CELL OR AN OUTER GHOST CELL
    IF(CELL_EM(INDEX)%CELL_X.LT.0.0.OR.CELL_EM(INDEX)%CELL_X.GT.LX.OR. &
       CELL_EM(INDEX)%CELL_Y.LT.0.0.OR.CELL_EM(INDEX)%CELL_Y.GT.LY.OR. &
       CELL_EM(INDEX)%CELL_Z.LT.0.0.OR.CELL_EM(INDEX)%CELL_Z.GT.LZ)THEN
      ID_INNER=0
    ELSE
      ID_INNER=1
    END IF
!---FOR INNER GHOST CELLS
    IF(ID_INNER.EQ.1)THEN
      CALL GETBC_GHOST_SEND(INDEX,NUM)
    ELSE
!---FOR OUTER GHOST CELLS
      IF(CELL_EM(INDEX)%CELL_X.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,1,BC(NUM,1),BV(NUM,1),NUM)
      ELSE IF(CELL_EM(INDEX)%CELL_X.GT.LX)THEN
        CALL GETBC_CELL(INDEX,2,BC(NUM,2),BV(NUM,2),NUM)
      END IF
      
      IF(CELL_EM(INDEX)%CELL_Y.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,3,BC(NUM,3),BV(NUM,3),NUM)
      ELSE IF(CELL_EM(INDEX)%CELL_Y.GT.LY)THEN
        CALL GETBC_CELL(INDEX,4,BC(NUM,4),BV(NUM,4),NUM)
      END IF

      IF(CELL_EM(INDEX)%CELL_Z.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,5,BC(NUM,5),BV(NUM,5),NUM)
      ELSE IF(CELL_EM(INDEX)%CELL_Z.GT.LZ)THEN
        CALL GETBC_CELL(INDEX,6,BC(NUM,6),BV(NUM,6),NUM)
      END IF
    END IF

    END SUBROUTINE GHOST_BOUNDARY
!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE GETBC_CELL(INDEX,ID,BC_LOCAL,BV_LOCAL,NUM)

    IMPLICIT NONE
    INTEGER:: NUM
    INTEGER:: INDEX_NEAR_1,INDEX_NEAR_2,INDEX_NEAR_3
    INTEGER:: BC_LOCAL,BV_LOCAL
    REAL(KIND=DP):: V0,V1,D0,D1   

    IF(ID.EQ.1)THEN
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_X(2)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_X(2)
    ELSE IF(ID.EQ.2)THEN     
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_X(1)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_X(1)
    ELSE IF(ID.EQ.3)THEN
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_Y(2)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_Y(2)
    ELSE IF(ID.EQ.4)THEN
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_Y(1)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_Y(1)
    ELSE IF(ID.EQ.5)THEN
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_Z(2)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_Z(2)
    ELSE IF(ID.EQ.6)THEN
      INDEX_NEAR_1=CELL_EM(INDEX)%CELL_NEAR     ! FIND THE NEAREST CELLS
      INDEX_NEAR_2=CELL_EM(INDEX_NEAR_1)%CELL_NEI_Z(1)
      INDEX_NEAR_3=CELL_EM(INDEX_NEAR_2)%CELL_NEI_Z(1)
    END IF

    IF(NUM.LE.4)THEN   ! FOR U, V, W AND TE
      BC_TYPE=BC_LOCAL
    ELSE IF(NUM.EQ.5)THEN ! FOR PD
      BC_TYPE=2
    ELSE
      BC_TYPE=5
    END IF
!---DIRICHLET BC-------------------------------------------------------------------          
    IF(BC_TYPE.EQ.1.OR.BC_TYPE.EQ.10)THEN   
      IF(ID.EQ.1.OR.ID.EQ.2)THEN
        CELL_EM(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_EM(INDEX)%CELL_X)/ABS(CELL_EM(INDEX_NEAR_1)%CELL_X)* &
                               (BV_LOCAL-CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM))
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        CELL_EM(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_EM(INDEX)%CELL_Y)/ABS(CELL_EM(INDEX_NEAR_1)%CELL_Y)* &
                               (BV_LOCAL-CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM))
      ELSE
        CELL_EM(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_EM(INDEX)%CELL_Z)/ABS(CELL_EM(INDEX_NEAR_1)%CELL_Z)* &
                               (BV_LOCAL-CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM))
      END IF
!---NEUMANN BC---------------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.2)THEN                
      IF(ID.EQ.1.OR.ID.EQ.2)THEN
        IF(ABS(CELL_EM(INDEX)%CELL_X).LT.CELL_EM(INDEX)%CELL_DX)THEN                 ! 1ST CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_X).LT.CELL_EM(INDEX)%CELL_DX*1.5*1.01)THEN   ! 2ND CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_X).LT.CELL_EM(INDEX)%CELL_DX*2.5*1.01)THEN   ! 3RD CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF 
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        IF(ABS(CELL_EM(INDEX)%CELL_Y).LT.CELL_EM(INDEX)%CELL_DY)THEN                 ! 1ST CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_Y).LT.CELL_EM(INDEX)%CELL_DY*1.5*1.01)THEN   ! 2ND CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_Y).LT.CELL_EM(INDEX)%CELL_DY*2.5*1.01)THEN   ! 3RD CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF
      ELSE
        IF(ABS(CELL_EM(INDEX)%CELL_Z).LT.CELL_EM(INDEX)%CELL_DZ)THEN                 ! 1ST CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_Z).LT.CELL_EM(INDEX)%CELL_DZ*1.5*1.01)THEN   ! 2ND CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_EM(INDEX)%CELL_Z).LT.CELL_EM(INDEX)%CELL_DZ*2.5*1.01)THEN   ! 3RD CELL
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF
      END IF
!---PERIODIC BC---------------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.3)THEN                
      CALL GETBC_GHOST_SEND(INDEX,NUM) 
!---2D DIRICHLET BC-----------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.4)THEN                
 
!---LINEAR EXTRAPOLATION------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.5)THEN                
      V0=CELL_EM(INDEX_NEAR_1)%CELL_VAR(NUM)
      V1=CELL_EM(INDEX_NEAR_2)%CELL_VAR(NUM)
      IF(ID.EQ.1.OR.ID.EQ.2)THEN        
        D1=CELL_EM(INDEX_NEAR_1)%CELL_DX
        D2=ABS(CELL_EM(INDEX_NEAR_1)%CELL_X-CELL_EM(INDEX)%CELL_X)
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        D1=CELL_EM(INDEX_NEAR_1)%CELL_DY
        D2=ABS(CELL_EM(INDEX_NEAR_1)%CELL_Y-CELL_EM(INDEX)%CELL_Y)
      ELSE 
        D1=CELL_EM(INDEX_NEAR_1)%CELL_DZ
        D2=ABS(CELL_EM(INDEX_NEAR_1)%CELL_Z-CELL_EM(INDEX)%CELL_Z)
      END IF
      CELL_EM(INDEX)%CELL_VAR(NUM)=V0+D2/D1*(V0-V1)
    END IF  
  
    END SUBROUTINE GETBC_CELL
!-------------------------------------------------------------------!
!       GET VALUES ON PARENT CELL BY MERGING FROM CHILD CELLS       !
!-------------------------------------------------------------------!
    SUBROUTINE CELL_VALUE_MERGE(INDEX,NUM)
    IMPLICIT NONE
    INTEGER:: INDEX,NUM

    IF(CELL_EM(INDEX)%CELL_SPLIT.EQ.0.OR.CELL_EM(INDEX)%CELL_GHOST.EQ.1)THEN
      RETURN
    END IF

    CELL_EM(INDEX)%CELL_VAR(NUM)=0.0
    DO I=1,2
      DO J=1,2
        DO K=1,2
          CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX)%CELL_VAR(NUM)+ &
            CELL_EM(CELL_EM(INDEX)%CELL_CHILD(I,J,K))%CELL_VAR(NUM)
        END DO
      END DO
    END DO

    CELL_EM(INDEX)%CELL_VAR(NUM)=CELL_EM(INDEX)%CELL_VAR(NUM)/8.0
      
    END SUBROUTINE CELL_VALUE_SPLIT
     
  END MODULE
