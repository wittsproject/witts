! This module declares the sizes of dynamic global arrays.
! Also, it will deallocate those arrays at the end of simulation.
  MODULE EMBED

  USE parameters
  USE field_shared
  USE class_shared
  USE class_cell

  IMPLICIT NONE

  CONTAINS
!-----------------------------------------------------------------!
!                                                                 !
!-----------------------------------------------------------------!
    SUBROUTINE EMBED_INITIAL()
    IMPLICIT NONE

    INTEGER:: I,J,K,II,JJ,KK,M
    INTEGER:: EM_LEVEL_TOTAL
    INTEGER:: NUM_CELL_LIMIT
    INTEGER,DIMENSION(:),ALLOCATABLE:: EM_NUM
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: LX_EM_1,LX_EM_2,&
                                         LY_EM_1,LY_EM_2,&
                                         LZ_EM_1,LZ_EM_2

    OPEN(1,FILE='embed.in')
    READ(1,*)EM_LEVEL_TOTAL

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

    ALLOCATE(CELL_FV(NUM_CELL_LIMIT))   
!---MAKE THE EMBEDDING
    DO I=1,EM_LEVEL_TOTAL
      DO J=1,EM_NUM(I)

        TOTAL_CELL0=TOTAL_CELL

        DO M=1,TOTAL_CELL0
          IF(CELL_FV(M)%CELL_X.GE.LX_EM_1(I,J).AND.CELL_FV(M)%CELL_X.LE.LX_EM_2(I,J).AND. &
             CELL_FV(M)%CELL_Y.GE.LY_EM_1(I,J).AND.CELL_FV(M)%CELL_Y.LE.LY_EM_2(I,J).AND. &
             CELL_FV(M)%CELL_Z.GE.LZ_EM_1(I,J).AND.CELL_FV(M)%CELL_Z.LE.LZ_EM_2(I,J).AND. &
             CELL_FV(M)%CELL_EMID.EQ.I-1)THEN  ! CHECK IF THE CELL IS QUALIFIED TO BE FURTHER REFINED

            DX=DX0/2.0**I
            DY=DY0/2.0**I
            DZ=DZ0/2.0**I

            CELL_FV(M)%CELL_SPLIT=1

            IF(TOTAL_CELL+8.GT.NUM_CELL_LIMIT)THEN
              EXIT
            END IF

            DO II=1,2
              DO JJ=1,2
                DO KK=1,2
                  TOTAL_CELL=TAOTAL_CELL+1
                  
                  CELL_FV(TOTAL_CELL)%CEL_INDEX=TOTAL_CELL
                  CELL_FV(TOTAL_CELL)%CELL_EMID=I
                  CELL_FV(TOTAL_CELL)%CELL_GHOST=0

                  CELL_FV(TOTAL_CELL)%CELL_PARENT=M
                  CELL_FV(M)%CELL_CHILD(II,JJ,KK)=TOTAL_CELL
                  
                  CELL_FV(TOTAL_CELL)%CELL_DX=DX
                  CELL_FV(TOTAL_CELL)%CELL_DY=DY
                  CELL_FV(TOTAL_CELL)%CELL_DZ=DZ

                  CELL_FV(TOTAL_CELL)%CELL_X=CELL_FV(M)%CELL_X-DX/2.0+DX*(II-1)
                  CELL_FV(TOTAL_CELL)%CELL_Y=CELL_FV(M)%CELL_Y-DY/2.0+DY*(II-1)
                  CELL_FV(TOTAL_CELL)%CELL_Z=CELL_FV(M)%CELL_Z-DZ/2.0+DZ*(II-1)

                  DO MM=1,NUM_VAR
                    CELL_FV(TOTAL_CELL)%CELL_VAR(MM)=CELL_FV(M)%CELL_VAR(MM)
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
      IF(CELL_FV(M)%CELL_GHOST.EQ.1)THEN
        DO II=1,NUM_VAR
          CALL GHOST_BOUNDARY(M,II)
        END DO
      END IF
    END DO
   
    DEALLOCATE(EM_NUM)
    DEALLOCATE(LX_EM_1,LX_EM_2,LY_EM_1,LY_EM_2,LZ_EM_1,LZ_EM_2)
    
    END SUBROUTINE EMBED_INITIAL

!-------------------------------------------------------------------!
!       GET VALUES ON PARENT CELL BY MERGING FROM CHILD CELLS       !
!-------------------------------------------------------------------!
    SUBROUTINE CELL_VALUE_MERGE(INDEX,NUM)
    IMPLICIT NONE
    INTEGER:: INDEX,NUM

    IF(CELL_FV(INDEX)%CELL_SPLIT.EQ.0.OR.CELL_FV(INDEX)%CELL_GHOST.EQ.1)THEN
      RETURN
    END IF

    CELL_FV(INDEX)%CELL_VAR(NUM)=0.0
    DO I=1,2
      DO J=1,2
        DO K=1,2
          CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX)%CELL_VAR(NUM)+ &
            CELL_FV(CELL_FV(INDEX)%CELL_CHILD(I,J,K))%CELL_VAR(NUM)
        END DO
      END DO
    END DO

    CELL_FV(INDEX)%CELL_VAR(NUM)=CELL_FV(INDEX)%CELL_VAR(NUM)/8.0
      
    END SUBROUTINE CELL_VALUE_MERGE    
  END MODULE
