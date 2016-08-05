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
    DO I=NX1,NX2
      DO J=NY1,NY2
        DO K=NZ1,NZ2
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LE.NUM_CELL_LIMIT)THEN
            CELL_EM(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
            CELL_EM(TOTAL_CELL)%CELL_EMID=0
            IF(I.GE.1.AND.I.LE.NX.AND.J.GE.1.AND.J.LE.NY.AND. &
               K.GE.1.AND.K.LE.NZ)THEN
              CELL_EM(TOTAL_CELL)%CELL_GHOST=0
            ELSE
              CELL_EM(TOTAL_CELL)%CELL_GHOST=1
            END IF

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
             CELL_EM(M)%CELL_EMID.EQ.I-1.AND.CELL_EM(M)%CELL_GHOST.EQ.0)THEN  ! CHECK IF THE CELL IS QUALIFIED TO BE FURTHER REFINED

            DX=DX0/2.0**I
            DY=DY0/2.0**I
            DZ=DZ0/2.0**I

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

      CELL_EM(M)%CELL_NEI_X=0
      CELL_EM(M)%CELL_NEI_Y=0
      CELL_EM(M)%CELL_NEI_Z=0

      IF(CELL_EM(M)%CELL_GHOST.EQ.0)THEN
        DO II=1,TOTAL_CELL
          IF(M.NE.II)THEN
            DISX=CELL_EM(II)%CELL_X-CELL_EM(M)%CELL_X
            DISY=CELL_EM(II)%CELL_Y-CELL_EM(M)%CELL_Y
            DISZ=CELL_EM(II)%CELL_Z-CELL_EM(M)%CELL_Z
            !  X DIRECTION
            IF(ABS(DISX).LT.CELL_EM(M)%CELL_DX*1.01.AND. &
               ABS(DISY).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
              IF(DISX.LT.0.0)THEN
                CELL_EM(M)%CELL_NEI_X(1)=II
              ELSE
                CELL_EM(M)%CELL_NEI_X(2)=II
              END IF
            END IF  
            !  Y DIRECTION
            IF(ABS(DISY).LT.CELL_EM(M)%CELL_DY*1.01.AND. &
               ABS(DISX).LT.ZERO.AND.ABS(DISZ).LT.ZERO)THEN
              IF(DISY.LT.0.0)THEN
                CELL_EM(M)%CELL_NEI_Y(1)=II
              ELSE
                CELL_EM(M)%CELL_NEI_Y(2)=II
              END IF
            END IF  
            !  Z DIRECTION
            IF(ABS(DISZ).LT.CELL_EM(M)%CELL_DZ*1.01.AND. &
               ABS(DISX).LT.ZERO.AND.ABS(DISY).LT.ZERO)THEN
              IF(DISZ.LT.0.0)THEN
                CELL_EM(M)%CELL_NEI_Y(1)=II
              ELSE
                CELL_EM(M)%CELL_NEI_Y(2)=II
              END IF
            END IF 
          END IF
        END DO
      END DO
    END IF
!---GENERATE GHOST CELLS FOR THE EMBEDDED CELLS
    
   
  END MODULE
