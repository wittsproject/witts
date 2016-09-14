! This module declares the sizes of dynamic global arrays.
! Also, it will deallocate those arrays at the end of simulation.
  MODULE class_cell

  USE mpi  
  USE parameters
  USE field_shared
  USE class_shared
  
  CONTAINS
!-------------------------------------------------------------------!
!                      INITIALIZE CELLS                             !
!-------------------------------------------------------------------!
    SUBROUTINE CELL_INITIAL()
    IMPLICIT NONE
    INTEGER:: I,J,K,M
      

    ALLOCATE(CELL_FV(NUM_CELL_LIMIT))   

!---SETUP THE BASIC CELLS
    TOTAL_CELL=0
    DO I=NX1,NX2
      DO J=NY1,NY2
        DO K=NZ1,NZ2
          TOTAL_CELL=TOTAL_CELL+1
          IF(TOTAL_CELL.LE.NUM_CELL_LIMIT)THEN
            CELL_FV(TOTAL_CELL)%CELL_INDEX=TOTAL_CELL
            CELL_FV(TOTAL_CELL)%CELL_EMID=0
            
            IF(I.GE.1.AND.I.LE.NX.AND.J.GE.1.AND.J.LE.NY.AND. &
               K.GE.1.AND.K.LE.NZ)THEN 
              CELL_FV(TOTAL_CELL)%CELL_GHOST=0
            ELSE
              CELL_FV(TOTAL_CELL)%CELL_GHOST=1
           END IF
           
            CELL_FV(TOTAL_CELL)%CELL_SPLIT=0
    
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

            CELL_FV(TOTAL_CELL)%CELL_VEL(1)=U(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VEL(2)=V(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_VEL(3)=W(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_TE=TE(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_RHO=RHO(I,J,K)
            CELL_FV(TOTAL_CELL)%CELL_PHI=PHI(I,J,K)
          ELSE
            PRINT*,'ERROR: the total number of basic cells exceeds the limit.'
            CALL MPI_FINALIZE(IERR)
            STOP
          END IF
        END DO
      END DO
    END DO
!---INITIALIZE THE GRID EMBEDDING
    IF(IEMBED.EQ.1)THEN
      CALL EMBED_INITIAL()      
    END IF

    DO I=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.0.AND. CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
        CELL_FV(M)%CELL_ACTIVE=1
      ELSE
        CELL_FV(M)%CELL_ACTIVE=0
      END IF
    END DO 
!---GET THE NEIGHBORS FOR EACH CELL
    DO M=1,TOTAL_CELL 
      CALL NEIGHBOR_INDEX(M) 
    END DO
!---GET THE BOUNDARY VALUES ON GHOST CELLS
    DO M=1,9
      CALL GHOST_BOUNDARY(9)
    END DO
   
    END SUBROUTINE CELL_INITIAL 
!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE NEIGHBOR_INDEX(INDEX)
    IMPLICIT NONE

    INTEGER:: INDEX
    INTEGER:: M
    REAL(KIND=DP):: DISX,DISY,DISZ,ZERO

    ZERO=1E-12

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
    SUBROUTINE GHOST_BOUNDARY(VAR, &
                              BC0,BV0)   !  THIS LINE FOR OPTIONAL
    IMPLICIT NONE

    INTEGER:: M,BC_TYPE
    INTEGER:: NUM,ID_INNER
    INTEGER:: I,ID_L,ID_LL,ID_R,ID_RR
    REAL(KIND=DP),DIMENSION(:):: VAR
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: X_LOC,Y_LOC,Z_LOC, &
                                             X_CELL,Y_CELL,Z_CELL, &
                                             VAR_LOC,VAR_CELL, &  
                                             XT,YT,ZT,VART
    REAL(KIND=DP):: DX,DY,DZ,ZERO

    INTEGER:: BC(6)
    REAL(KIN=DP):: BV(6)
    INTEGER, OPTIONAL:: BC0(6)
    REAL(KIN=DP),OPTIONAL:: BV0(6)

    CALL UPDATE_INNER_GHOST(VAR)

    DATA ZERO /1.0E-6/

    IF(PRESENT(BC0))THEN
      BC=BC0
    ELSE
      BC=5
    END IF

    IF(PRESENT(BV0))THEN
      BV=BV0
    ELSE
      BV=0.0
    END IF
!---ASSEMBLE CELL ARRAY-----------------------------------------
    ALLOCATE(X_LOC(TOTAL_CELL_ACTIVE))
    ALLOCATE(Y_LOC(TOTAL_CELL_ACTIVE))   
    ALLOCATE(Z_LOC(TOTAL_CELL_ACTIVE)) 
    ALLOCATE(VAR_LOC(TOTAL_CELL_ACTIVE))

    ALLOCATE(XT(GLOBAL_CELL_ACTIVE))
    ALLOCATE(YT(GLOBAL_CELL_ACTIVE))   
    ALLOCATE(ZT(GLOBAL_CELL_ACTIVE)) 
    ALLOCATE(VART(GLOBAL_CELL_ACTIVE))
    
    MM=0
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_ACTIVE.EQ.1)THEN
        MM=MM+1          
        X_LOC(MM)=CELL_FV(M)%CELL_X
        Y_LOC(MM)=CELL_FV(M)%CELL_Y
        Z_LOC(MM)=CELL_FV(M)%CELL_Z
        VAR_LOC(MM)=VAR(M)        
      END IF
    END DO

    CALL ASSEM(X_LOC,XT,TOTAL_CELL_ACTIVE)
    CALL ASSEM(Y_LOC,YT,TOTAL_CELL_ACTIVE)
    CALL ASSEM(Z_LOC,ZT,TOTAL_CELL_ACTIVE)
    CALL ASSEM(VAR_LOC,VART,TOTAL_CELL_ACTIVE)
  
    DEALLOCATE(X_LOC,Y_LOC,Z_LOC,VAR_LOC)
   
!---UPDATE VALUES ON GHOST CELLS-----------------------------------    
    DO M=1,TOTAL_CELL
      IF(CELL_FV(M)%CELL_GHOST.EQ.1.AND.CELL_FV(M)%CELL_SPLIT.EQ.0)THEN  
        X0=CELL_FV(M)%CELL_X
        Y0=CELL_FV(M)%CELL_Y
        Z0=CELL_FV(M)%CELL_Z
!---FOR OUTRT GHOST CELL WITH PERIODIC BC, USE THE COORDINATE OF ITS INNER PAIR
        IF(X0.LT.0.0.AND.BC(1).EQ.3)THEN
          X0=X0+LX
        ELSE IF(X0.GT.LX.AND.BC(2).EQ.3)THEN
          X0=X0-LX
        END IF

        IF(Y0.LT.0.0.AND.BC(3).EQ.3)THEN
          Y0=Y0+LY
        ELSE IF(Y0.GT.LY.AND.BC(4).EQ.3)THEN
          Y0=Y0-LY
        END IF

        IF(Z0.LT.0.0.AND.BC(5).EQ.3)THEN
          Z0=Z0+LZ
        ELSE IF(Z0.GT.LZ.AND.BC(6).EQ.3)THEN
          Z0=Z0-LZ
        END IF
!---JUDGE IF THE CELL IS AN INNER GHOST CELL OR AN OUTER GHOST CELL
        IF(X0.LT.0.0.OR.X0.GT.LX.OR.Y0.LT.0.0.OR.Y0.GT.LY.OR. &
           Z0.LT.0.0.OR.Z0.GT.LZ)THEN
          ID_INNER=1
        ELSE
          ID_INNER=0
        END IF
!---FOR INNER GHOST CELLS, COPY VALUES FROM ITS PAIR
        IF(ID_INNER.EQ.1)THEN
 JLOOP:   DO MM=1,GLOBAL_CELL_ACTIVE
            DX=ABS(XT(MM)-CELL_FV(M)%CELL_X)
            DY=ABS(YT(MM)-CELL_FV(M)%CELL_Y)
            DZ=ABS(ZT(MM)-CELL_FV(M)%CELL_Z)
            IF(DX.LT.ZERO.AND.DY.LT.ZERO.AND.DZ.LT.ZERO)THEN
              VAR(M)=VART(MM)
              EXIT JLOOP
            END IF
          END DO JLOOP
!---FOR OUTER GHOST CELLS
        ELSE
          IF(CELL_FV(M)%CELL_X.LT.0.0)THEN
            IF(BC(1).EQ.3)THEN  ! PERIODIC BC
              
            ELSE
              CALL GETBC_CELL(VAR,M,1,BC(1),BV(1))
            END IF
          ELSE IF(CELL_FV(M)%CELL_X.GT.LX)THEN
            CALL GETBC_CELL(VAR,M,2,BC(2),BV(2))
          END IF
      
          IF(CELL_FV(M)%CELL_Y.LT.0.0)THEN
            CALL GETBC_CELL(VAR,M,3,BC(3),BV(3))
          ELSE IF(CELL_FV(M)%CELL_Y.GT.LY)THEN
            CALL GETBC_CELL(VAR,M,4,BC(4),BV(4))
          END IF

          IF(CELL_FV(M)%CELL_Z.LT.0.0)THEN
            CALL GETBC_CELL(VAR,M,5,BC(5),BV(5))
          ELSE IF(CELL_FV(M)%CELL_Z.GT.LZ)THEN
            CALL GETBC_CELL(VAR,M,6,BC(6),BV(6))
          END IF
        END IF
      END IF
    END DO

    DEALLOCATE(XT,YT,ZT,VART)

    END SUBROUTINE GHOST_BOUNDARY
!---------------------------------------------------!
!      OBTAIN THE INDEX OF NEIGHBORING CELLS        !
!---------------------------------------------------!
    SUBROUTINE GETBC_CELL(VAR,INDEX,ID,BC_LOCAL,BV_LOCAL)

    IMPLICIT NONE
    INTEGER:: INDEX,ID
    INTEGER:: INDEX_NEAR_1,INDEX_NEAR_2,INDEX_NEAR_3
    INTEGER:: BC_LOCAL
    REAL(KIND=DP),DIMENSION(:):: VAR
    REAL(KIND=DP):: BV_LOCAL
    REAL(KIND=DP):: V0,V1,D0,D1,D2   

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

!---DIRICHLET BC-------------------------------------------------------------------          
    IF(BC_LOCAL.EQ.1.OR.BC_LOCAL.EQ.10)THEN   
      IF(ID.EQ.1.OR.ID.EQ.2)THEN
        VAR(INDEX)=BV_LOCAL+ &
                   ABS(CELL_FV(INDEX)%CELL_X)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_X)* &
                   (BV_LOCAL-VAR(INDEX_NEAR_1))
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        VAR(INDEX)=BV_LOCAL+ &
                   ABS(CELL_FV(INDEX)%CELL_Y)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_Y)* &
                   (BV_LOCAL-VAR(INDEX_NEAR_1))
      ELSE
        VAR(INDEX)=BV_LOCAL+ &
                   ABS(CELL_FV(INDEX)%CELL_Z)/ABS(CELL_FV(INDEX_NEAR_1)%CELL_Z)* &
                   (BV_LOCAL-VAR(INDEX_NEAR_1))
      END IF
!---NEUMANN BC---------------------------------------------------------------------
    ELSE IF(BC_LOCAL.EQ.2)THEN                
      IF(ID.EQ.1.OR.ID.EQ.2)THEN
        IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX)THEN                 ! 1ST CELL
          VAR(INDEX)=VAR(INDEX_NEAR_1)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX*1.5*1.01)THEN   ! 2ND CELL
          VAR(INDEX)=VAR(INDEX_NEAR_2)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_X).LT.CELL_FV(INDEX)%CELL_DX*2.5*1.01)THEN   ! 3RD CELL
          VAR(INDEX)=VAR(INDEX_NEAR_3)
        END IF 
      ELSE IF(ID.EQ.3.OR.ID.EQ.4)THEN
        IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY)THEN                 ! 1ST CELL
          VAR(INDEX)=VAR(INDEX_NEAR_1)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY*1.5*1.01)THEN   ! 2ND CELL
          VAR(INDEX)=VAR(INDEX_NEAR_2)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Y).LT.CELL_FV(INDEX)%CELL_DY*2.5*1.01)THEN   ! 3RD CELL
          VAR(INDEX)=VAR(INDEX_NEAR_3)
        END IF
      ELSE
        IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ)THEN                 ! 1ST CELL
          VAR(INDEX)=VAR(INDEX_NEAR_1)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ*1.5*1.01)THEN   ! 2ND CELL
          VAR(INDEX)=VAR(INDEX_NEAR_2)
        ELSE IF(ABS(CELL_FV(INDEX)%CELL_Z).LT.CELL_FV(INDEX)%CELL_DZ*2.5*1.01)THEN   ! 3RD CELL
          VAR(INDEX)=VAR(INDEX_NEAR_3)
        END IF
      END IF
!---2D DIRICHLET BC-----------------------------------------------------------------
    ELSE IF(BC_LOCAL.EQ.4)THEN                
 
!---LINEAR EXTRAPOLATION------------------------------------------------------------
    ELSE IF(BC_LOCAL.EQ.5)THEN                
      V0=VAR(INDEX_NEAR_1)
      V1=VAR(INDEX_NEAR_2)
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
      VAR(INDEX)=V0+D2/D1*(V0-V1)
    END IF  
  
    END SUBROUTINE GETBC_CELL
!---------------------------------------------------------!
!    TRANSFORM THE CELL ARRAY TO A 3D STRUCTURED ARRAY    !
!---------------------------------------------------------!    
    SUBROUTINE CELL_TO_STRUCT(VAR_CELL,INDEX,NB,VAR)

    IMPLICIT NONE

    INTEGER:: INDEX,NB
    REAL(KIND=DP),DIMENSION(:):: VAR_CELL
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

    ILOOP2: DO I=-NB,NB
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
            EXIT KLOOP2
          END IF           
        END DO KLOOP2   
      END DO JLOOP
    END DO ILOOP2
!---GET THE CORRESPONDING VARIABLE VALUES
    DO I=-NB,NB
      DO J=-NB,NB
        DO K=-NB,NB  
          VAR(I,J,K)=VAR_CELL(IND(I,J,K))
        END DO
      END DO     
    END DO

    END SUBROUTINE CELL_TO_STRUCT
!-------------------------------------------------------------------!
!           GET THE CELL NUMBER OF NEIGHBORING CELLS                !
!-------------------------------------------------------------------!
!  This function is used to get the index of neighboring cells
!  INDEX: Index of the cell
!  ID: =1, x direction; =2, y direction; =3, z direction
!  NUM: Number of index skipped (positive or negative)  
   FUNCTION LOOKUP_NEI(INDEX,ID,NUM)
   IMPLICIT NONE
   INTEGER:: INDEX,ID,NUM,IND
   INTEGER:: LOOKUP_NEI,I
   
   IND=INDEX

   IF(NUM.EQ.0)THEN
     LOOKUP_NEI=INDEX
     RETURN
   END IF  
!--X DIRECTION
   IF(ID.EQ.1)THEN
     IF(NUM.GT.0)THEN 
       DO I=1,NUM
         IND=CELL_FV(IND)%CELL_NEI_X(2)       
       END DO
     ELSE
       DO I=-1,NUM,-1
         IND=CELL_FV(IND)%CELL_NEI_X(1)       
       END DO
     END IF
!--Y DIRECTION    
   ELSE IF(ID.EQ.2)THEN
     IF(NUM.GT.0)THEN 
       DO I=1,NUM
         IND=CELL_FV(IND)%CELL_NEI_Y(2)       
       END DO
     ELSE
       DO I=-1,NUM,-1
         IND=CELL_FV(IND)%CELL_NEI_Y(1)       
       END DO
     END IF
!--Z DIRECTION    
   ELSE
     IF(NUM.GT.0)THEN 
       DO I=1,NUM
         IND=CELL_FV(IND)%CELL_NEI_Z(2)       
       END DO
     ELSE
       DO I=-1,NUM,-1
         IND=CELL_FV(IND)%CELL_NEI_Z(1)       
       END DO
     END IF
   END IF

   LOOKUP_NEI=IND

   END FUNCTION LOOKUP_NEI
!-------------------------------------------------!
!               GET THE CELL SKIP                 !           
!-------------------------------------------------!
!   This function is used to obtain the number of
!   active cells that has to be skipped when read or export
!   a global array   
    INTEGER FUNCTION CELL_SKIP(COUNT)
    IMPLICIT NONE
    INTEGER:: NCPU,M,TRANI,COUNT
    INTEGER,DIMENSION(:),ALLOCATABLE:: RANK_CELL
   
!---OBTAIN THE CELL COUNT ON EACH RANK
    NCPU=NPX*NPY*NPZ-1
    ALLOCATE(RANK_CELL(0:NCPU))
      
    DO M=0,NCPU
      IF(MYID.EQ.M)THEN
        RANK_CELL(M)=COUNT
      ELSE
        RANK_CELL(M)=0 
      END IF
    END DO

    DO M=0,NCPU 
       CALL MPI_ALLREDUCE(RANK_CELL(M),TRANI,1,MPI_INTEGER,MPI_MAX, &
                          MPI_COMM_WORLD,IERR) 
       RANK_CELL(M)=TRANI
    END DO

    CELL_SKIP=0
    DO M=1,MYID
      CELL_SKIP=CELL_SKIP+RANK_CELL(M-1)
    END DO
   
    DEALLOCATE(RANK_CELL)
    
    END FUNCTION CELL_SKIP
 
  END MODULE CLASS_CELL
