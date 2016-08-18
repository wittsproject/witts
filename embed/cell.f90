! This module declares the sizes of dynamic global arrays.
! Also, it will deallocate those arrays at the end of simulation.
  MODULE CLASS_CELL

  USE parameter
  USE field_shared
  USE class_shared
  
  CONTAINS
!-------------------------------------------------------------------!
!                      INITIALIZE CELLS                             !
!-------------------------------------------------------------------!
    SUBROUTINE CELL_INITIAL()
    IMPLICIT NONE

    ALLOCATE(CELL_FV(NUM_CELL_LIMIT))   

!---SETUP THE BASIC CELLS
    TOTAL_CELL=0
    DO I=1,NX
      DO J=1,NY
        DO K=1,NZ
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LE.NUM_CELL_LIMIT)THEN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
            CELL_FV(TOTAL_CELL)%CELL_EMID=0
            CELL_FV(TOTAL_CELL)%CELL_GHOST=0
            CELL_FV(TOTAL_CELL)%CELL_SPLIT=0
            CELL_FV(TOTAL_CELL)%CELL_NUM_VAR=NUM_VAR

            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=0    
            CELL_FV(TOTAL_CELL)%CELL_PID=MYID
            CELL_FV(TOTAL_CELL)%CELL_NEAR=0
            
            CELL_FV(TOTAL_CELL)%CELL_PARENT=0
            CELL_FV(TOTAL_CELL)%CELL_CHILD=0

            CELL_FV(TOTAL_CELL)%CELL_DX=DX0
            CELL_FV(TOTAL_CELL)%CELL_DY=DY0
            CELL_FV(TOTAL_CELL)%CELL_DZ=DZ0

            CELL_FV(TOTAL_CELL)%CELL_X=XI(I+MYIDX*NX)
            CELL_FV(TOTAL_CELL)%CELL_Y=YI(J+MYIDY*NY)
            CELL_FV(TOTAL_CELL)%CELL_Z=ZI(K+MYIDZ*NZ)

            CELL_FV(TOTAL_CELL)%CELL_VAR(1)=U(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(2)=V(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(3)=W(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(4)=TE(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(5)=PD(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(6)=NU(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(7)=MU(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(8)=RHO(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VAR(9)=PHI(I,J,K)
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO
      END DO
    END DO

    END SUBROUTINE CELL_INITIAL
  
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
      IF(CELL_FV(M)%CELL_NEI_X(1)=0)THEN
        IUP=MYID+1
        IDOWN=MYID-1          

        IF(MYIDX.EQ.NPX-1)THEN
          IUP=MYID-(NPX-1)
        END IF
        IF(MYIDX.EQ.0)THEN
          IDOWN=MYID+(NPX-1)
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_X(2))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_X(2))%CELL_NEI_X(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=I
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_FV(M)%CELL_NEI_X(2)=0)THEN
        IUP=MYID+1
        IDOWN=MYID-1          

        IF(MYIDX.EQ.NPX-1)THEN
          IUP=MYID-(NPX-1)
        END IF
        IF(MYIDX.EQ.0)THEN
          IDOWN=MYID+(NPX-1)
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_X(1))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_X(1))%CELL_NEI_X(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=1
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)        
      END IF
!---Y DIRECTION--------------------
      IF(CELL_FV(M)%CELL_NEI_Y(1)=0)THEN          
        IUP=MYID+NPX*NPZ
        IDOWN=MYID-NPX*NPZ

        IF(MYIDY.EQ.NPY-1)THEN
          IUP=MYIDX+MYIDZ*NPX
        END IF
        IF(MYIDY.EQ.0)THEN
          IDOWN=MYIDX+(NPY-1)*NPX*NPZ+MYIDZ*NPX
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_Y(2))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_Y(2))%CELL_NEI_Y(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=1
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_FV(M)%CELL_NEI_Y(2)=0)THEN
        IUP=MYID+NPX*NPZ
        IDOWN=MYID-NPX*NPZ

        IF(MYIDY.EQ.NPY-1)THEN
          IUP=MYIDX+MYIDZ*NPX
        END IF
        IF(MYIDY.EQ.0)THEN
          IDOWN=MYIDX+(NPY-1)*NPX*NPZ+MYIDZ*NPX
        END IF
       
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_Y(1))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_Y(1))%CELL_NEI_Y(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=1
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)        
      END IF
!---Z DIRECTION--------------------
      IF(CELL_FV(M)%CELL_NEI_Y(1)=0)THEN          
        IUP=MYID+NPX
        IDOWN=MYID-NPX        

        IF(MYIDZ.EQ.NPZ-1)THEN
          IUP=MYIDX+MYIDY*NPX*NPZ
        END IF
        IF(MYIDZ.EQ.0)THEN
          IDOWN=MYIDX+MYIDY*NPX*NPZ+(NPZ-1)*NPX
        END IF
               
        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_Z(2))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_Z(2))%CELL_NEI_Z(2))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,idown,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,iup,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=1
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO

        DEALLOCATE(BUF_SEND,BUF_RECE)
        
      ELSE IF(CELL_FV(M)%CELL_NEI_Z(2)=0)THEN
        IUP=MYID+NPX
        IDOWN=MYID-NPX        

        IF(MYIDZ.EQ.NPZ-1)THEN
          IUP=MYIDX+MYIDY*NPX*NPZ
        END IF
        IF(MYIDZ.EQ.0)THEN
          IDOWN=MYIDX+MYIDY*NPX*NPZ+(NPZ-1)*NPX
        END IF

        ALLOCATE(BUF_SEND(3),BUF_RECE(3))
        BUF_SEND(1)=CELL_FV(M)
        BUF_SEND(2)=CELL_FV(CELL_FV(M)%CELL_NEI_Z(1))
        BUF_SEND(3)=CELL_FV(CELL_FV(CELL_FV(M)%CELL_NEI_Z(1))%CELL_NEI_Z(1))
       
        CALL mpi_isend(BUF_SEND,3,mpi_TYPE,iup,1,   &
                       mpi_comm_world,isd1,ierR)
        CALL mpi_irecv(BUF_RECE,3,mpi_TYPE,idown,1, &
                       mpi_comm_world,irv1,ierR)
        CALL mpi_wait(isd1,stat,ierR)
        CALL mpi_wait(irv1,stat,ierR)

        DO I=1,3
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LT.NUM_CELL_LIMIT)THEN 
            CELL_FV(TOTAL_CELL)=BUF_RECE(I)
            CELL_FV(TOTAL_CELL)%CELL_GHOST=1
            CELL_FV(TOTAL_CELL)%CELL_INDEX_ORIGIN=CELL_FV(TOTAL_CELL)%CELL_INDEX  ! STORE THE ORIGINAL INDEX IN THE CELL_INDEX_ORIGIN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
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
      IF(CELL_FV(M)%CELL_GHOST.EQ.1)THEN
        DIS0=1.0D6
        DO MM=1,TOTAL_CELL
          IF(CELL_FV(MM)%CELL_GHOST_EQ.0)THEN
            DIS=SQRT((CELL_FV(MM)%CELL_X-CELL_FV(M)%CELL_X)**2+ &
                     (CELL_FV(MM)%CELL_Y-CELL_FV(M)%CELL_Y)**2+ &
                     (CELL_FV(MM)%CELL_Z-CELL_FV(M)%CELL_Z)**2)
            IF(DIS.LT.DIS0)THEN
              DIS0=DIS
              CELL_FV(M)%CELL_NEAR=MM
            END IF
          END IF
        END DO
      END IF
    END DO          
 
    END SUBROUTINE GENERATE_GHOST_CELL
  
!----------------------------------------------------!
!         GET BC FOR AN INNER GHOST CELL             !
!----------------------------------------------------!    
    SUBROUTINE GETBC_GHOST_INNER(INDEX,NUM) 
    IMPLICIT NONE
    INTEGER:: INDEX,NUM
    INTEGER:: INDEX_SEND,PID_SEND,WIN,IERR

    INDEX_SEND=CELL_FV(INDEX)%CELL_INDEX_ORIGIN
    PID_SEND=CELL_FV(INDEX)%CELL_PID

    CALL MPI_GET(CELL_FV(INDEX)%CELL_VAR(NUM),1,MPI_DOUBLE_PRECISION, &
                 PID_SEND,CELL_FV(INDEX_SEND)%CELL_VAR(NUM),          &
                 1,MPI_DOUBLE_PRECISION,WIN,IERR)   
 
    END SUBROUTINE GETBC_GHOST_INNER

!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE NEIGHBOR_INDEX(INDEX)
    IMPLICIT NONE

    INTEGER:: INDEX
    INTEGER:: M

    DO M=1,TOTAL_CELL
      IF(M.NE.INDEX)THEN 
        DISX=CELL_FV(M)%CELL_X-CELL_FV(INDEX)%CELL_X
        DISY=CELL_FV(M)%CELL_Y-CELL_FV(INDEX)%CELL_Y
        DISZ=CELL_FV(M)%CELL_Z-CELL_FV(INDEX)%CELL_Z        
        IF(CELL_FV(INDEX)%CELL_EMID.EQ.CELL_FV(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT THE SAME LEVEL
          !  X DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*1.01.AND. &
             ABS(DISY).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
            IF(DISX.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_X(1)=M
              CELL_FV(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_X(2)=M
              CELL_FV(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF  
          !  Y DIRECTION
          IF(ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*1.01.AND. &
             ABS(DISX).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
            IF(DISY.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Y(1)=M
              CELL_FV(M)%CELL_NEI_Y(2)=INDEX               
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Y(2)=M
              CELL_FV(M)%CELL_NEI_Y(1)=INDEX                
            END IF
          END IF  
          !  Z DIRECTION
          IF(ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*1.01.AND. &
             ABS(DISX).LT.ZERO.AND.ABS(DISY).LT.ZERO)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Z(1)=M
              CELL_FV(M)%CELL_NEI_Z(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Z(2)=M
              CELL_FV(M)%CELL_NEI_Z(1)=INDEX 
            END IF
          END IF
        ELSE IF(CELL_FV(INDEX)%CELL_EMID+1.EQ.CELL_FV(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT ONE LEVEL ABOVE
           !  X DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*0.75*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*0.25*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*0.25*1.01)THEN
            IF(DISX.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_X(1)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_X(2)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF
           !  Y DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*0.75*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*0.25*1.01)THEN
            IF(DISY.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Y(1)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_Y(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Y(2)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_Y(1)=INDEX 
            END IF
          END IF
          !  Z DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*0.25*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*0.75*1.01)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Z(1)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_Z(2)=INDEX                
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Z(2)=CELL_FV(M)%CELL_PARENT
              CELL_FV(M)%CELL_NEI_Z(1)=INDEX  
            END IF
          END IF
        ELSE IF(CELL_FV(INDEX)%CELL_EMID-1.EQ.CELL_FV(M)%CELL_EMID)THEN ! LOOP OVER CELLS AT ONE LEVEL BELOW
        !  X DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*1.5*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*0.5*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*0.5*1.01)THEN
            IF(DISX.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_X(1)=M
              CELL_FV(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_X(2)=M
              CELL_FV(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF
          !  Y DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*0.5*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*1.5*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*0.5*1.01)THEN
            IF(DISY.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Y(1)=M
              CELL_FV(M)%CELL_NEI_Y(2)=INDEX 
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Y(2)=M
              CELL_FV(M)%CELL_NEI_Y(1)=INDEX 
            END IF
          END IF
          !  Z DIRECTION
          IF(ABS(DISX).LT.CELL_FV(INDEX)%CELL_DX*0.5*1.01.AND. &
             ABS(DISY).LT.CELL_FV(INDEX)%CELL_DY*0.5*1.01.AND. &
             ABS(DISZ).LT.CELL_FV(INDEX)%CELL_DZ*1.5*1.01)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_FV(INDEX)%CELL_NEI_Z(1)=M
              CELL_FV(M)%CELL_NEI_Z(2)=INDEX
            ELSE
              CELL_FV(INDEX)%CELL_NEI_Z(2)=M
              CELL_FV(M)%CELL_NEI_Z(1)=INDEX
            END IF
          END IF           
        END IF

      END IF
    END DO
    END SUBROUTINE NEIGHBOR_INDEX

!-------------------------------------------------------------------!
!                     SET BC ON THE GHOST CELLS                     !
!-------------------------------------------------------------------!
!
    SUBROUTINE GHOST_BOUNDARY(INDEX,NUM)
    IMPLICIT NONE

    INTEGER:: INDEX,BC_TYPE
    INTEGER:: NUM,ID_INNER
    INTEGER:: I,ID_L,ID_LL,ID_R,ID_RR

 
    IF(CELL_FV(INDEX)%CELL_GHOST.NE.1)THEN
      RETURN
    END IF

    IF(NUM.LT.1)THEN
      NUM=1
    END IF

    IF(NUM.GT.NUM_VAR)THEN
      NUM=NUM_VAR
    END IF
!---JUDGE IF THE CELL IS AN INNER GHOST CELL OR AN OUTER GHOST CELL
    IF(CELL_FV(INDEX)%CELL_X.LT.0.0.OR.CELL_FV(INDEX)%CELL_X.GT.LX.OR. &
       CELL_FV(INDEX)%CELL_Y.LT.0.0.OR.CELL_FV(INDEX)%CELL_Y.GT.LY.OR. &
       CELL_FV(INDEX)%CELL_Z.LT.0.0.OR.CELL_FV(INDEX)%CELL_Z.GT.LZ)THEN
      ID_INNER=0
    ELSE
      ID_INNER=1
    END IF
!---FOR INNER GHOST CELLS
    IF(ID_INNER.EQ.1)THEN
      CALL GETBC_GHOST_INNER(INDEX,NUM)
    ELSE
!---FOR OUTER GHOST CELLS
      IF(CELL_FV(INDEX)%CELL_X.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,1,BC(NUM,1),BV(NUM,1),NUM)
      ELSE IF(CELL_FV(INDEX)%CELL_X.GT.LX)THEN
        CALL GETBC_CELL(INDEX,2,BC(NUM,2),BV(NUM,2),NUM)
      END IF
      
      IF(CELL_FV(INDEX)%CELL_Y.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,3,BC(NUM,3),BV(NUM,3),NUM)
      ELSE IF(CELL_FV(INDEX)%CELL_Y.GT.LY)THEN
        CALL GETBC_CELL(INDEX,4,BC(NUM,4),BV(NUM,4),NUM)
      END IF

      IF(CELL_FV(INDEX)%CELL_Z.LT.0.0)THEN
        CALL GETBC_CELL(INDEX,5,BC(NUM,5),BV(NUM,5),NUM)
      ELSE IF(CELL_FV(INDEX)%CELL_Z.GT.LZ)THEN
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
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_X(2)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_X(2)
    ELSE IF(ID.EQ.2)THEN     
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_X(1)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_X(1)
    ELSE IF(ID.EQ.3)THEN
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_Y(2)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_Y(2)
    ELSE IF(ID.EQ.4)THEN
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_Y(1)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_Y(1)
    ELSE IF(ID.EQ.5)THEN
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_Z(2)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_Z(2)
    ELSE IF(ID.EQ.6)THEN
      INDEX_NEAR_1=CELL_FV(INDEX)%CELL_NEAR     ! FIND THE NEAREST 3 CELLS
      INDEX_NEAR_2=CELL_FV(INDEX_NEAR_1)%CELL_NEI_Z(1)
      INDEX_NEAR_3=CELL_FV(INDEX_NEAR_2)%CELL_NEI_Z(1)
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
        CELL_FV(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_FV(INDEX)%CELL_X)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_X)* &
                               (BV_LOCAL-CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM))
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        CELL_FV(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_FV(INDEX)%CELL_Y)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_Y)* &
                               (BV_LOCAL-CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM))
      ELSE
        CELL_FV(INDEX)%CELL_VAR(NUM)=BV_LOCAL+ &
                               ABS(CELL_FV(INDEX)%CELL_Z)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_Z)* &
                               (BV_LOCAL-CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM))
      END IF
!---NEUMANN BC---------------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.2)THEN                
      IF(ID.EQ.1.OR.ID.EQ.2)THEN
        IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX)THEN                 ! 1ST CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX*1.5*1.01)THEN   ! 2ND CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX*2.5*1.01)THEN   ! 3RD CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF 
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY)THEN                 ! 1ST CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY*1.5*1.01)THEN   ! 2ND CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY*2.5*1.01)THEN   ! 3RD CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF
      ELSE
        IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ)THEN                 ! 1ST CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ*1.5*1.01)THEN   ! 2ND CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_2)%CELL_VAR(NUM)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ*2.5*1.01)THEN   ! 3RD CELL
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX_NEAR_3)%CELL_VAR(NUM)
        END IF
      END IF
!---PERIODIC BC---------------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.3)THEN                
      CALL GETBC_GHOST_INNER(INDEX,NUM) 
!---2D DIRICHLET BC-----------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.4)THEN                
 
!---LINEAR EXTRAPOLATION------------------------------------------------------------
    ELSE IF(BC_TYPE.EQ.5)THEN                
      V0=CELL_FV(INDEX_NEAR_1)%CELL_VAR(NUM)
      V1=CELL_FV(INDEX_NEAR_2)%CELL_VAR(NUM)
      IF(ID.EQ.1.OR.ID.EQ.2)THEN        
        D1=CELL_FV(INDEX_NEAR_1)%CELL_DX
        D2=ABS(CELL_FV(INDEX_NEAR_1)%CELL_X-CELL_FV(INDEX)%CELL_X)
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        D1=CELL_FV(INDEX_NEAR_1)%CELL_DY
        D2=ABS(CELL_FV(INDEX_NEAR_1)%CELL_Y-CELL_FV(INDEX)%CELL_Y)
      ELSE 
        D1=CELL_FV(INDEX_NEAR_1)%CELL_DZ
        D2=ABS(CELL_FV(INDEX_NEAR_1)%CELL_Z-CELL_FV(INDEX)%CELL_Z)
      END IF
      CELL_FV(INDEX)%CELL_VAR(NUM)=V0+D2/D1*(V0-V1)
    END IF  
  
    END SUBROUTINE GETBC_CELL
   
!---------------------------------------------------!
!      OBTAIN THE DERIVATIVE OF CELL VARIABLES      !
!---------------------------------------------------!
! INDEX: INDEX OF THE CELL
! ID=1: X DERIVATIVE; =2: Y DERIVATIVE; =3: Z DERIVATIVE
! SCHEME=1: CENTRAL SCHEME; =2 UPWIND SCHEME
! ISTAG=1: GET DERIVATIVE AT A STAGGERED LOCATION
!      =0: GET DERIVATIVE AT THE EXACT LOCATION
! ORDER: ORDER OF NUMERICAL ACCURACY    
    REAL(KIND=DP) FUNCTION DERIV_CELL(INDEX,NUM,ID,SCHEME,ISTAG,ORDER)
    IMPLICIT NONE
    INTEGER:: INDEX,NUM,ID,SCHEME,ISTAG,ORDER
    INTEGER:: I,J
    INTEGER:: IND(-ORDER+1:ORDER-1)
    REAL(KIND=DP),DIMENSION(:,:,:) ALLOCATABLE:: VAR

    NB=ORDER-1
    
    ALLOCATE(VAR(-NB:NB,-NB:NB,-NB:NB)
    
    CALL CELL_TO_STRUCT(INDEX,NB,NUM,VAR)    
!---X DERIVATIVE---------------------------    
    IF(ID.EQ.1)THEN      
      IF(SCHEME.EQ.1)THEN             
        DERIV_CELL=DERIV_X(VAR,-NB,-NB,-NB,0,0,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DX)
      ELSE
        !  Upwind scheme should be added here  
      END IF
!---Y DERIVATIVE---------------------------
    ELSE IF(ID.EQ.2)THEN
      IF(SCHEME.EQ.1)THEN       
        DERIV_CELL=DERIV_Y(VAR,-NB,-NB,-NB,0,0,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DY)
      ELSE
        !  Upwind scheme should be added here  
      END IF
!---Z DERIVATIVE---------------------------
    ELSE      
      IF(SCHEME.EQ.1)THEN     
        DERIV_CELL=DERIV_Z(VAR,-NB,-NB,-NB,0,0,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DZ)
      ELSE
        !  Upwind scheme should be added here  
      END IF      
    END IF
   
    DEALLOCATE(VAR)
    
    END FUNCTION DERIV_CELL

!---------------------------------------------------------!
!    TRANSFORM THE CELL ARRAY TO A 3D STRUCTURED ARRAY    !
!---------------------------------------------------------!    
    SUBROUTINE CELL_TO_STRUCT(INDEX,NB,NUM,VAR)

    IMPLICIT NONE

    INTEGER:: INDEX,NB,NUM
    REAL(KIND=DP),DIMENSION(-NB:,-NB:,-NB:):: VAR
    INTEGER:: IND(-NB:NB,-NB:NB,-NB:NB)
    INTEGER:: I,J,K,M,LEVEL
    
!---GET THE INDEX FOR EACH NEIGHBORING CELL         
    IND(0,0,0)=INDEX
    LEVEL=CELL_FV(INDEX)%CELL_EMID
    
    DO I=1,NB
      IND(I,0,0)=CELL_FV(IND(I-1,0,0))%CELL_NEI_X(2)      
      IF(CELL_FV(IND(I,0,0))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
        DO M=I,NB
          IND(M,0,0)=IND(I+INT((M-I)/2),0,0)
        END DO
        EXIT
      END IF 
    END DO
   
    DO I=-1,-NB,-1
      IND(I,0,0)=CELL_FV(IND(I+1,0,0))%CELL_NEI_X(1)
      IF(CELL_FV(IND(I,0,0))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
        DO M=I,-NB,-1
          IND(M,0,0)=IND(I+INT((M-I)/2),0,0)
        END DO
        EXIT
      END IF       
    END DO

    ILOOP: DO I=-NB,NB
      JLOOP1: DO J=1,NB
        IND(I,J,0)=CELL_FV(IND(I,J-1,0))%CELL_NEI_Y(2)
        IF(CELL_FV(IND(I,J,0))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
          DO M=J,NB
            IND(I,M,0)=IND(I,J+INT((M-J)/2),0)
          END DO
          EXIT JLOOP1
        END IF       
      END DO JLOOP1
     
      JLOOP2: DO J=-1,-NB,-1
        IND(I,J,0)=CELL_FV(IND(I,J+1,0))%CELL_NEI_Y(1)
        IF(CELL_FV(IND(I,J,0))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
          DO M=J,-NB,-1
            IND(I,M,0)=IND(I,J+INT((M-J)/2),0)
          END DO
          EXIT JLOOP2
        END IF         
      END DO JLOOP2
    END DO ILOOP

    ILOOP: DO I=-NB,NB
      JLOOP: DO J=-NB,NB 
        KLOOP1: DO K=1,NB
          IND(I,J,K)=CELL_FV(IND(I,J,K-1))%CELL_NEI_Z(2)
          IF(CELL_FV(IND(I,J,K))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
            DO M=K,NB
              IND(I,J,M)=IND(I,J,K+INT((M-K)/2))
            END DO
            EXIT KLOOP1
          END IF           
        END DO KLOOP1

        KLOOP2: DO K=-1,-NB,-1
          IND(I,J,K)=CELL_FV(IND(I,J,K+1))%CELL_NEI_Z(1)
          IF(CELL_FV(IND(I,J,K))%CELL_EMID.LT.LEVEL)THEN  ! If encounter a lower-level cell
            DO M=K,-NB,-1
              IND(I,J,M)=IND(I,J,K+INT((M-K)/2))
            END DO
            EXIT JLOOP2
          END IF           
        END DO KLOOP2   
      END DO JLOOP
    END DO ILOOP 
!---GET THE CORRESPONDING VARIABLE VALUES
    DO I=-NB,NB
      DO J=-NB,NVB
        DO K=-NB,NB
           VAR(I,J,K)=CELL_FV(IND(I,J,K))%CELL_VAR(NUM)
        END DO
      END DO     
   END DO

   END SUBROUTINE CELL_TO_STRUCT
   
  END MODULE CLASS_CELL
