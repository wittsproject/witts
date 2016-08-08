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

    ALLOCATE(CELL_EM(NUM_CELL_LIMIT))


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

            CELL_EM(TOTAL_CELL)%CELL_MASTER_INDEX=0          

            CELL_EM(TOTAL_CELL)%CELL_DX=DX0
            CELL_EM(TOTAL_CELL)%CELL_DY=DY0
            CELL_EM(TOTAL_CELL)%CELL_DZ=DZ0

            CELL_EM(TOTAL_CELL)%CELL_X=XI(I+MYIDX*NX)
            CELL_EM(TOTAL_CELL)%CELL_Y=YI(J+MYIDY*NY)
            CELL_EM(TOTAL_CELL)%CELL_Z=ZI(K+MYIDZ*NZ)

            CELL_EM(TOTAL_CELL)%CELL_VEL(1)=U(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VEL(2)=V(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_VEL(3)=W(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_TEM=TE(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_PD=PD(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_NU=NU(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_MU=MU(I,J,K)
            CELL_EM(TOTAL_CELL)%CELL_RHO=RHO(I,J,K)
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO
      END DO
    END DO
!---MAKE THE EMBEDDING
    DO I=1,EM_LEVEL
      DO J=1,EM_NUM(I)

        TOTAL_CELL0=TOTAL_CELL

        DO M=1,TOTAL_CELL0
          IF(CELL_EM(M)%CELL_X.GE.LX_EM(I,J).AND.CELL_EM(M)%CELL_X.GE.LX_EM(I,J).AND. &
             CELL_EM(M)%CELL_Y.GE.LY_EM(I,J).AND.CELL_EM(M)%CELL_Y.GE.LY_EM(I,J).AND. &
             CELL_EM(M)%CELL_Z.GE.LZ_EM(I,J).AND.CELL_EM(M)%CELL_Z.GE.LZ_EM(I,J).AND. &
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

                  CELL_EM(TOTAL_CELL)%CELL_MASTER_INDEX=M         

                  CELL_EM(TOTAL_CELL)%CELL_DX=DX
                  CELL_EM(TOTAL_CELL)%CELL_DY=DY
                  CELL_EM(TOTAL_CELL)%CELL_DZ=DZ

                  CELL_EM(TOTAL_CELL)%CELL_X=CELL_EM(M)%CELL_X-DX/2.0+DX*(II-1)
                  CELL_EM(TOTAL_CELL)%CELL_Y=CELL_EM(M)%CELL_Y-DY/2.0+DY*(II-1)
                  CELL_EM(TOTAL_CELL)%CELL_Z=CELL_EM(M)%CELL_Z-DZ/2.0+DZ*(II-1)

                  CELL_EM(TOTAL_CELL)%CELL_VEL(1)=CELL_EM(M)%CELL_VEL(1)
                  CELL_EM(TOTAL_CELL)%CELL_VEL(2)=CELL_EM(M)%CELL_VEL(2)
                  CELL_EM(TOTAL_CELL)%CELL_VEL(3)=CELL_EM(M)%CELL_VEL(3)

                  CELL_EM(TOTAL_CELL)%CELL_TEM=CELL_EM(M)%CELL_TEM
                  CELL_EM(TOTAL_CELL)%CELL_PD=CELL_EM(M)%CELL_PD
                  CELL_EM(TOTAL_CELL)%CELL_NU=CELL_EM(M)%CELL_NU
                  CELL_EM(TOTAL_CELL)%CELL_MU=CELL_EM(M)%CELL_MU
                  CELL_EM(TOTAL_CELL)%CELL_RHO=CELL_EM(M)%CELL_RHO   
                END DO
              END DO
            END DO            
          END IF
        END DO

      END DO         
    END DO
!---GET THE NEIGHBORS FOR EACH CELL
    DO M=1,TOTAL_CELL 
      CALL NEIGHBOR_INDEX(CELL_EM,M,TOTAL_CELL) 
    END DO
!---GENERATE GHOST CELLS
    
!-------------------------------------------------------------------!
!                  GENERATE GHOST CELLS AND SET BC                  !
!-------------------------------------------------------------------!
    SUBROUTINE GETBC_CELL(CELL_EM,TOTAL_CELL)
    IMPLICIT NONE
      
    TYPE(CELL),DIMENSION(:):: CELL_EM
    TYPE(CELL),DIMENSION(:),ALLOCATABLE:: BUF_SEND,BUF_RECE
    INTEGER :: TOTAL_CELL

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
            CELL_EM(TOTAL_CELL)%CELL_GHOST=1
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
!---UPDATE NEIGHBOR INDEX
    DO M=1,TOTAL_CELL
      IF(CELL_EM(M)%CELL_GHOST.EQ.1)THEN ! LOOP OVER GHOST CELLS
        CALL NEIGHBOR_INDEX(CELL_EM,M,TOTAL_CELL)
      END IF
    END DO
!---APPLY BOUNDARY CONDITIONS

 
    END SUBROUTINE GETBC_CELL           
!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE NEIGHBOR_INDEX(CELL_EM,INDEX,TOTAL_CELL)
    IMPLICIT NONE

    TYPE(CELL),DIMENSION(:):: CELL_EM
    INTEGER:: INDEX,TOTAL_CELL
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
              CELL_EM(INDEX)%CELL_NEI_X(1)=CELL_EM(M)%CELL_MASTER_INDEX
              CELL_EM(M)%CELL_NEI_X(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_X(2)=CELL_EM(M)%CELL_MASTER_INDEX
              CELL_EM(M)%CELL_NEI_X(1)=INDEX 
            END IF
          END IF
           !  Y DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.75*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.25*1.01)THEN
            IF(DISY.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Y(1)=CELL_EM(M)%CELL_MASTER_INDEX
              CELL_EM(M)%CELL_NEI_Y(2)=INDEX 
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Y(2)=CELL_EM(M)%CELL_MASTER_INDEX
              CELL_EM(M)%CELL_NEI_Y(1)=INDEX 
            END IF
          END IF
          !  Z DIRECTION
          IF(ABS(DISX).LT.CELL_EM(INDEX)%CELL_DX*0.25*1.01.AND. &
             ABS(DISY).LT.CELL_EM(INDEX)%CELL_DY*0.25*1.01.AND. &
             ABS(DISZ).LT.CELL_EM(INDEX)%CELL_DZ*0.75*1.01)THEN
            IF(DISZ.LT.0.0)THEN
              CELL_EM(INDEX)%CELL_NEI_Z(1)=CELL_EM(M)%CELL_MASTER_INDEX
              CELL_EM(M)%CELL_NEI_Z(2)=INDEX                
            ELSE
              CELL_EM(INDEX)%CELL_NEI_Z(2)=CELL_EM(M)%CELL_MASTER_INDEX
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
  
  END MODULE
