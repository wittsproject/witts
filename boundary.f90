! This module is used for setting boundary conditions
!
  MODULE boundary

  USE mpi
  USE parameters, ONLY: DP,NBX,NBY,NBZ,ICOLL,MYID,MYIDX,MYIDY,MYIDZ,NPX,NPY,NPZ,BC,BV
  USE field_shared

  CONTAINS
!=========================================================================!
!                        VELOCITY FIELD BC                                !
!=========================================================================!
    SUBROUTINE BOUNDARY_VEL(NX,NY,NZ)

    IMPLICIT NONE
    INTEGER ::  NX,NY,NZ

    IF(ICOLL.EQ.0)THEN
      CALL GET_BC(NX,NY,NZ,U,NBX,NBY,NBZ,1, &
                  BC(1,1),BV(1,1),U0)  
      CALL GET_BC(NX,NY,NZ,V,NBX,NBY,NBZ,2, &
                  BC(1,2),BV(1,2),V0)
      CALL GET_BC(NX,NY,NZ,W,NBX,NBY,NBZ,3, &
                  BC(1,3),BV(1,3),W0)
    ELSE
      CALL GET_BC(NX,NY,NZ,U,NBX,NBY,NBZ,0, &
                  BC(1,1),BV(1,1),U0)
      CALL GET_BC(NX,NY,NZ,V,NBX,NBY,NBZ,0, &
                  BC(1,2),BV(1,2),V0)
      CALL GET_BC(NX,NY,NZ,W,NBX,NBY,NBZ,0, &
                  BC(1,3),BV(1,3),W0)
    END IF

    END SUBROUTINE
!*************************************************************************!
!                      GET BOUNDARY CONDITIONS                            !
!*************************************************************************!
!     ID=1: AT CELL'S X FACE
!       =2: AT CELL'S Y FACE
!       =3: AT CELL'S Z FACE
!       =0: AT CELL'S CENTER
!     BC=1: DIRICHLET BC
!       =2: NEUMANN BC
!       =3: PERIODIC BC
!       =4: USE EXTERNAL INPUT (A0)
!       =5: LINEAR EXTRAPOLATION
!     OPTIONAL ARGUMENTS:
!     IBC: IF NOT EXPLICITLY SPECIFIED, THEN USE DEFAULT VALUE
!     A0: EXTERNAL BC INPUT, ONLY WHEN BC=4
!     VB1, .., VB6: BOUNDARY VALUES, ONLY WHEN BC=1
      SUBROUTINE GET_BC(NX,NY,NZ,A,MX,MY,MZ,ID,  &
                        IBC,VB,A0)    ! THIS LINE IS FOR OPTIONAL ARGUMENTS

      IMPLICIT NONE
!      INCLUDE "mpif.h"
      INTEGER :: MX,MY,MZ,ID,NX,NY,NZ
      REAL(KIND=DP), DIMENSION(1-MX:,1-MY:,1-MZ:):: A
      REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE:: BUFS0X,BUFS1X,BUFR0X,BUFR1X, &
                                            BUFS0Y,BUFS1Y,BUFR0Y,BUFR1Y, &
                                            BUFS0Z,BUFS1Z,BUFR0Z,BUFR1Z
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE,OPTIONAL:: A0      
      REAL(KIND=DP),OPTIONAL :: VB(6)
      INTEGER,OPTIONAL :: IBC(6)
      INTEGER :: BCI(6),I,J,K,IL,IERR,ISD1,IRV1,ISD2,IRV2,IYIKSIZE
      INTEGER :: IUP,IDOWN,STAT(MPI_STATUS_SIZE)

      BCI=5  ! THE DEFAULT BC TYPE IS LINEAR EXTRAPOLATION
      IF(PRESENT(IBC))THEN
      DO I=1,6
          BCI(I)=IBC(I)
        END DO
      END IF

      ALLOCATE(BUFS0X(1-MY:NY+MY,1-MZ:NZ+MZ,MX),BUFS1X(1-MY:NY+MY,1-MZ:NZ+MZ,MX),&
               BUFR0X(1-MY:NY+MY,1-MZ:NZ+MZ,MX),BUFR1X(1-MY:NY+MY,1-MZ:NZ+MZ,MX))
      ALLOCATE(BUFS0Y(1-MX:NX+MX,1-MZ:NZ+MZ,MY),BUFS1Y(1-MX:NX+MX,1-MZ:NZ+MZ,MY),&
               BUFR0Y(1-MX:NX+MX,1-MZ:NZ+MZ,MY),BUFR1Y(1-MX:NX+MX,1-MZ:NZ+MZ,MY))
      ALLOCATE(BUFS0Z(1-MX:NX+MX,1-MY:NY+MY,MZ),BUFS1Z(1-MX:NX+MX,1-MY:NY+MY,MZ),&
               BUFR0Z(1-MX:NX+MX,1-MY:NY+MY,MZ),BUFR1Z(1-MX:NX+MX,1-MY:NY+MY,MZ))
!-----FOR X DIRECTION----------------------------------------------------------
      IUP=MYID+1
      IDOWN=MYID-1

      IYIKSIZE =(NY+2*MY)*(NZ+2*MZ)*MX

      DO J=1-MY,NY+MY
        DO K=1-MZ,NZ+MZ
          DO IL=1,MX
            BUFS0X(J,K,IL)=A(NX+1-IL,J,K)
            BUFS1X(J,K,IL)=A(IL,J,K)
          END DO
        END DO
      END DO
!-----PERIODIC BC
      IF(MYIDX.EQ.NPX-1)THEN
        IUP=MYID-(NPX-1)
      END IF
      IF(MYIDX.EQ.0)THEN
        IDOWN=MYID+(NPX-1)
      END IF

      CALL mpi_isend(bufs0x,iyiksize,mpi_DOUBLE_PRECISION,iup,1,   &
                     mpi_comm_world,isd1,ierR)
      CALL mpi_irecv(bufr0x,iyiksize,mpi_double_precision,idown,1, &
                     mpi_comm_world,irv1,ierR)
      CALL mpi_wait(isd1,stat,ierR)
      CALL mpi_wait(irv1,stat,ierR)

      CALL mpi_isend(bufs1x,iyiksize,mpi_double_precision,idown,1, &
                     mpi_comm_world,isd2,ierR)
      CALL mpi_irecv(bufr1x,iyiksize,mpi_double_precision,iup,1,   &
                     mpi_comm_world,irv2,ierR)
      CALL mpi_wait(isd2,stat,ierR)
      CALL mpi_wait(irv2,stat,ierR)
      DO K=1,NZ
        DO J=1,NY
          DO IL=1,MX
            A(1-IL,J,K)=BUFR0X(J,K,IL)
            A(NX+IL,J,K)=BUFR1X(J,K,IL)
          END DO
        END DO
      END DO
!-----CORRECTION OF GLOBAL BC
      IF(MYIDX.EQ.0)THEN
        DO K=1,NZ
          DO J=1,NY
            DO IL=1,MX
              IF(BCI(1).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.1)THEN         !  FOR X FACE
                  A(1,J,K)=VB(1)
                  A(1-IL,J,K)=VB(1)*2.0-A(1+IL,J,K)
                ELSE                     
                  A(1-IL,J,K)=VB(1)*2.0-A(IL,J,K)
                END IF
              ELSE IF(BCI(1).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.1)THEN         !  FOR X FACE
                  A(1-IL,J,K)=A(1+IL,J,K)
                ELSE                              
                  A(1-IL,J,K)=A(IL,J,K)
                END IF 
              ELSE IF(BCI(1).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(1-IL,J,K)=A0(IL,J,K)  
              ELSE IF(BCI(1).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(1-IL,J,K)=A(1-IL+1,J,K)*2-A(1-IL+2,J,K)         
              END IF
            END DO
          END DO
        END DO
      END IF
      
      IF(MYIDX.EQ.NPX-1)THEN
        DO K=1,NZ
          DO J=1,NY
            DO IL=1,MX
              IF(BCI(2).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.1)THEN         !  FOR X FACE
                  IF(IL.EQ.1)THEN
                    A(NX+1,J,K)=VB(2)
                  ELSE
                    A(NX+IL,J,K)=VB(2)*2.0-A(NX+2-IL,J,K)
                  END IF
                ELSE                              
                  A(NX+IL,J,K)=VB(2)*2.0-A(NX+1-IL,J,K)
                END IF
              ELSE IF(BCI(2).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.1)THEN         !  FOR X FACE
                  IF(IL.EQ.1)THEN
                    A(NX+1,J,K)=A(NX,J,K)
                  ELSE
                    A(NX+IL,J,K)=A(NX+2-IL,J,K)
                  END IF
                ELSE                              
                  A(NX+IL,J,K)=A(NX+1-IL,J,K)
                END IF 
              ELSE IF(BCI(2).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(NX+IL,J,K)=A0(IL,J,K)   
              ELSE IF(BCI(2).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(NX+IL,J,K)=A(NX+IL-1,J,K)*2-A(NX+IL-2,J,K)                  
              END IF
            END DO
          END DO
        END DO
      END IF
!-----Y DIRECTION-----------------------------------------------
      IUP=MYID+NPX*NPZ
      IDOWN=MYID-NPX*NPZ

      IYIKSIZE=(NX+2*MX)*(NZ+2*MZ)*MY

      DO I=1-MX,NX+MX
        DO K=1-MZ,NZ+MZ
          DO IL=1,MY
            BUFS0Y(I,K,IL)=A(I,NY+1-IL,K)
            BUFS1Y(I,K,IL)=A(I,IL,K)
          END DO
        END DO
      END DO
!-----PERIODIC BC
      IF(MYIDY.EQ.NPY-1)THEN
        IUP=MYIDX+MYIDZ*NPX
      END IF
      IF(MYIDY.EQ.0)THEN
        IDOWN=MYIDX+(NPY-1)*NPX*NPZ+MYIDZ*NPX
      END IF

      CALL mpi_isend(bufs0y,iyiksize,mpi_double_precision,iup,1,  &
                     mpi_comm_world,isd1,ierR)
      CALL mpi_irecv(bufr0y,iyiksize,mpi_double_precision,idown,1,&
                     mpi_comm_world,irv1,ierR)
      CALL mpi_wait(isd1,stat,ierR)
      CALL mpi_wait(irv1,stat,ierR)

      CALL mpi_isend(bufs1y,iyiksize,mpi_double_precision,idown,1,&
                     mpi_comm_world,isd2,ierR)
      CALL mpi_irecv(bufr1y,iyiksize,mpi_double_precision,iup,1,  &
                     mpi_comm_world,irv2,ierR)
      CALL mpi_wait(isd2,stat,ierR)
      CALL mpi_wait(irv2,stat,ierR)
      DO K=1,NZ
        DO I=1-MX,NX+MX
          DO IL=1,MY
            A(I,1-IL,K)=BUFR0Y(I,K,IL)
            A(I,NY+IL,K)=BUFR1Y(I,K,IL)
          ENDDO
        END DO
      END DO
!-----CORRECTION OF GLOBAL BC
      IF(MYIDY.EQ.0)THEN
        DO K=1,NZ
          DO I=1-MX,NX+MX
            DO IL=1,MY
              IF(BCI(3).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.2)THEN         !  FOR Y FACE
                  A(I,1,K)=VB(3)
                  A(I,1-IL,K)=VB(3)*2.0-A(I,1+IL,K)
                ELSE                            
                  A(I,1-IL,K)=VB(3)*2.0-A(I,IL,K)
                END IF
              ELSE IF(BCI(3).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.2)THEN         !  FOR Y FACE
                  A(I,1-IL,K)=A(I,1+IL,K)
                ELSE                              
                  A(I,1-IL,K)=A(I,IL,K)
                END IF 
              ELSE IF(BCI(3).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(I,1-IL,K)=A0(I,IL,K)  
              ELSE IF(BCI(3).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(I,1-IL,K)=A(I,1-IL+1,K)*2-A(I,1-IL+2,K)                 
              END IF
            END DO
          END DO
        END DO
      END IF

      IF(MYIDY.EQ.NPY-1)THEN
        DO K=1,NZ
          DO I=1-MX,NX+MX
            DO IL=1,MY
              IF(BCI(4).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.2)THEN         !  FOR Y FACE
                  IF(IL.EQ.1)THEN
                    A(I,NY+1,K)=VB(4)
                  ELSE
                    A(I,NY+IL,K)=VB(4)*2.0-A(I,NY+2-IL,K)
                  END IF
                ELSE                              
                  A(I,NY+IL,K)=VB(4)*2.0-A(I,NY+1-IL,K)
                END IF
              ELSE IF(BCI(4).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.2)THEN         !  FOR V ON STAGGERED GRID
                  IF(IL.EQ.1)THEN
                    A(I,NY+1,K)=A(I,NY,K)
                  ELSE
                    A(I,NY+IL,K)=A(I,NY+2-IL,K)
                  END IF
                ELSE                        
                  A(I,NY+IL,K)=A(I,NY+1-IL,K)
                END IF   
              ELSE IF(BCI(4).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(I,NY+IL,K)=A0(I,IL,K)  
              ELSE IF(BCI(4).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(I,NY+IL,K)=A(I,NY+IL-1,K)*2-A(I,NY+IL-2,K)       
              END IF
            END DO
          END DO
        END DO
      END IF
!-----Z DIRECTION------------------------------------------------
      IUP=MYID+NPX
      IDOWN=MYID-NPX

      IYIKSIZE=(NY+2*MY)*(NX+2*MX)*MZ

      DO I=1-MX,NX+MX
        DO J=1-MY,NY+MY
          DO IL=1,MZ
            BUFS0Z(I,J,IL)=A(I,J,NZ+1-IL)
            BUFS1Z(I,J,IL)=A(I,J,IL)
          END DO
        END DO
      END DO
!-----PERIODIC BC
      IF(MYIDZ.EQ.NPZ-1)THEN
        IUP=MYIDX+MYIDY*NPX*NPZ
      END IF
      IF(MYIDZ.EQ.0)THEN
        IDOWN=MYIDX+MYIDY*NPX*NPZ+(NPZ-1)*NPX
      END IF

      CALL mpi_isend(bufs0z,iyiksize,mpi_double_precision,iup,1,   &
                     mpi_comm_world,isd1,ierR)
      CALL mpi_irecv(bufr0z,iyiksize,mpi_double_precision,idown,1, &
                     mpi_comm_world,irv1,ierR)
      CALL mpi_wait(isd1,stat,ierR)
      CALL mpi_wait(irv1,stat,ierR)

      CALL mpi_isend(bufs1z,iyiksize,mpi_double_precision,idown,1, &
                     mpi_comm_world,isd2,ierR)
      CALL mpi_irecv(bufr1z,iyiksize,mpi_double_precision,iup,1,   &
                     mpi_comm_world,irv2,ierR)
      CALL mpi_wait(isd2,stat,ierR)
      CALL mpi_wait(irv2,stat,ierR)

      DO I=1-MX,NX+MX
        DO J=1-MY,NY+MY
          DO IL=1,MZ
            A(I,J,1-IL)=BUFR0Z(I,J,IL)
            A(I,J,NZ+IL)=BUFR1Z(I,J,IL)
          END DO
        END DO
      END DO
!-----CORRECTION OF GLOBAL BC
      IF(MYIDZ.EQ.0)THEN
        DO J=1-MY,NY+MY
          DO I=1-MX,NX+MX
            DO IL=1,MZ
              IF(BCI(5).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.3)THEN         !  FOR Z FACE
                  A(I,J,1)=VB(5)
                  A(I,J,1-IL)=VB(5)*2.0-A(I,J,1+IL)
                ELSE                             
                  A(I,J,1-IL)=VB(5)*2.0-A(I,J,IL)
                END IF
              ELSE IF(BCI(5).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.3)THEN         ! FOR W ON STAGGERED GRID
                  A(I,J,1-IL)=A(I,J,1+IL)
                ELSE                     
                  A(I,J,1-IL)=A(I,J,IL)
                END IF 
              ELSE IF(BCI(5).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(I,J,1-IL)=A0(I,J,IL)                   
              ELSE IF(BCI(5).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(I,J,1-IL)=A(I,J,1-IL+1)*2-A(I,J,1-IL+2)              
              END IF
            END DO
          END DO
        END DO
      END IF
      
      IF(MYIDZ.EQ.NPZ-1)THEN
        DO J=1-MY,NY+MY
          DO I=1-MX,NX+MX
            DO IL=1,MZ
              IF(BCI(6).EQ.1)THEN        !  DIRICHLET BC WITH FIXED VALUE
                IF(ID.EQ.3)THEN         !  FOR Z FACE
                  IF(IL.EQ.1)THEN
                    A(I,J,NZ+1)=VB(6)
                  ELSE
                    A(I,J,NZ+IL)=VB(6)*2.0-A(I,J,NZ+2-IL)
                  END IF
                ELSE                             
                  A(I,J,NZ+IL)=VB(6)*2.0-A(I,J,NZ+1-IL)
                END IF
              ELSE IF(BCI(6).EQ.2)THEN   !  NEUMANN BC 
                IF(ID.EQ.3)THEN         !  FOR Z FACE
                  IF(IL.EQ.1)THEN
                    A(I,J,NZ+1)=A(I,J,NZ)
                  ELSE
                    A(I,J,NZ+IL)=A(I,J,NZ+2-IL)
                  END IF
                ELSE                            
                  A(I,J,NZ+IL)=A(I,J,NZ+1-IL)
                END IF       
              ELSE IF(BCI(6).EQ.4)THEN   !  USE EXTERNAL INPUT 
                A(I,J,NZ+IL)=A0(I,J,IL)   
              ELSE IF(BCI(6).EQ.5)THEN   !  LINEAR EXTRAPOLATION
                A(I,J,NZ+IL)=A(I,J,NZ+IL-1)*2-A(I,J,NZ+IL-2)              
              END IF      
            END DO
          END DO
        END DO
      END IF

      DEALLOCATE(BUFS0X,BUFS1X,BUFR0X,BUFR1X, &
                 BUFS0Y,BUFS1Y,BUFR0Y,BUFR1Y, &
                 BUFS0Z,BUFS1Z,BUFR0Z,BUFR1Z)

      END SUBROUTINE 

  END MODULE
