! This is the main prcess of a simulation.
!
     PROGRAM MAIN

     USE mpi
     USE parameters
     USE field_shared
     USE allocation
     USE initial
     USE inflow
     USE dynamic_dt
     USE boundary
     USE momentum
     USE scalar
     USE post_process

     IMPLICIT NONE

!      INTERFACE 
!        SUBROUTINE READ_INPUT()
!        END SUBROUTINE 
!      END INTERFACE

      INTEGER :: NCPU,NPROC,I,J,K

      CALL MPI_INIT(IERR)     
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)

      
!-----READ INPUT FILES (*.IN)
      CALL READ_INPUT()
!-----INITIAL DIAGNOSE------------------------------------
      NCPU=NPX*NPY*NPZ
      IF(NPROC.NE.NCPU) THEN
        IF(MYID.EQ.0)THEN
          PRINT *,"WRONG NUMBER OF CPUS:",NPROC
        END IF
        STOP
      ENDIF
      IF(MOD(NXT,NPX).NE.0 .OR. MOD(NYT,NPY).NE.0.OR. & 
         MOD(NZT,NPZ).NE.0) THEN
        PRINT *,"NX,NY,NZ SHOULD BE DIVIDED BY NPX,NPY,NPZ EXACTLY"
        STOP
      ENDIF
!-----INITIALIZATION--------------------------------------
      NX1=1-NBX
      NX2=NX+NBX
      NY1=1-NBY
      NY2=NY+NBY
      NZ1=1-NBZ
      NZ2=NZ+NBZ
!-----ALLOCATE GLOBAL ARRAYS
      CALL ALLOCATE_GLOBAL()

      CALL INITIALIZE()
!-----CHECK IF TIME MAX OR NT IS REACHED----------------
      IF(TIME.GT.TIME_MAX.OR.NSTART.GE.NT)THEN
        IF(MYID.EQ.0)THEN
          PRINT*,'STOP: time max or total time step is reached'
        END IF
        CALL MPI_FINALIZE(IERR)
        STOP
      END IF
!-----MAIN LOOP-------------------------------------------
      DO N=NSTART,NT
        IF(N.EQ.NSTART.OR.MOD(N,5).EQ.0)THEN
          CALL READ_INPUT()  ! READ *.in FILES
        END IF   

        IF(MYID.EQ.0)THEN
          PRINT*,'-----------------------------------------'
        END IF

        CALL GENDT(DX0,DY0,DZ0)  ! DYNAMICALLY SET TIME STEP

        IF(MYID.EQ.0)THEN
          PRINT*,'TIME STEP=',N,' TIME=',TIME,' DT=',DT
       END IF
!------IF IINFLOW=1, READ INFLOW FILES
       IF(IINFLOW_READ.EQ.1)THEN
         CALL INFLOW_READ()
         IF(ISTOP.EQ.1)THEN
           EXIT
         END IF
       END IF
!------SOLVE THE MOMENTUM EQUATION    
        IF(ITDER.EQ.1)THEN        !  A-B SCHEME
          CALL MOMENTUM_AB(DX0,DY0,DZ0)
        ELSE IF(ITDER.EQ.2)THEN   !  R-K SCHEME
          CALL MOMENTUM_RK(DX0,DY0,DZ0)
        ELSE IF(ITDER.EQ.20)THEN  !  SEMI R-K SCHEME (SEMI: PROJECTION IS OUTSIDE THE LOOP)
          CALL MOMENTUM_SEMI_RK(DX0,DY0,DZ0)
       END IF      
!------SOLVE THE TEMPERATURE EQUATION
        IF(ITEMP.EQ.1)THEN
          CALL SCALAR_WRAP(TE,NX1,NY1,NZ1,DX0,DY0,DZ0)
        END IF               
!------POST PROCESS AND STATISTICS
        CALL POSTPROCESS()

        TIME=TIME+DT      

        IF(ISTOP.EQ.1.OR.TIME.GT.TIME_MAX)THEN
          EXIT
        END IF
      ENDDO    

      CALL DEALLOCATE_GLOBAL()

      CALL MPI_FINALIZE(IERR)

      END PROGRAM
!***********************************************************************!
!                SUBROUTINE OF READING INPUT FILES                      !
!***********************************************************************!
      SUBROUTINE READ_INPUT()

      USE parameters
      USE boundary
      USE inflow
      
      IMPLICIT NONE
      INTEGER :: I,J

!      IF(MYID.EQ.0)THEN
!        PRINT*,'Read "parameters.in"'
!      END IF
      OPEN(1,FILE="parameter.in")
      READ(1,*)
      READ(1,*) LX
      READ(1,*) LY
      READ(1,*) LZ
      READ(1,*)
      READ(1,*) RHO0
      READ(1,*) MU0
      READ(1,*) LATITUDE
      READ(1,*) TR
      READ(1,*) USTAR0
      READ(1,*) UG0
      READ(1,*) VG0
      READ(1,*) WG0
      READ(1,*)
      READ(1,*) G
      READ(1,*) KAPPA
      READ(1,*) PI
      READ(1,*) OMEGA_EARTH
      CLOSE(1)

!      IF(MYID.EQ.0)THEN
!        PRINT*,'Read "solver.in"'
!      END IF
      OPEN(1, FILE="solver.in")
      READ(1,*) IINCOM  
      READ(1,*) ITYPE     
      READ(1,*) ISCHEME   
      READ(1,*) ICOLL     
      READ(1,*) 
      READ(1,*) ITDER
      READ(1,*) ORDER_TIM
      READ(1,*) CFL
      READ(1,*) DT_MAX
      READ(1,*)
      READ(1,*) IFLUX
      READ(1,*) ILIMIT
      READ(1,*) ORDER_CON 
      READ(1,*) BLEND_CON 
      READ(1,*)
      READ(1,*) IIBM      
      READ(1,*) ITURBINE  
      READ(1,*) ICORI     
      READ(1,*) IGEOB
      READ(1,*) IBOUS       
      READ(1,*)
      READ(1,*) ORDER_POI 
      READ(1,*) ISCHE_POI 
      READ(1,*) NMUL_POI  
      READ(1,*) NITE_POI  
      READ(1,*) TOLE_POI  
      READ(1,*)
      READ(1,*) ITEMP
      READ(1,*) ISCALAR        
      CLOSE(1)

!      IF(MYID.EQ.0)THEN
!        PRINT*,'Read "control.in"'
!      END IF
      OPEN(1,FILE="control.in")
      READ(1,*) ISTOP  
      READ(1,*)
      READ(1,*) IRESTART
      READ(1,*)
      READ(1,*) IPROLX
      READ(1,*) IPROLY
      READ(1,*) IPROLZ
      READ(1,*)
      READ(1,*) IINFLOW_READ
      READ(1,*) IINFLOW_WRITE
      READ(1,*)      
      READ(1,*) NXT
      READ(1,*) NYT 
      READ(1,*) NZT 
      READ(1,*) NPX 
      READ(1,*) NPY 
      READ(1,*) NPZ
      READ(1,*)
      READ(1,*) TIME_MAX
      READ(1,*) NT
      READ(1,*)
      READ(1,*) SCREEN_LEVEL
      CLOSE(1)

!      IF(MYID.EQ.0)THEN
!        PRINT*,'Read "boundary.in"'
!      END IF
      OPEN(1,FILE="boundary.in")
      READ(1,*) NBX
      READ(1,*) NBY
      READ(1,*) NBZ       
      DO I=1,6
        READ(1,*)
        DO J=1,5
          READ(1,*) BC(I,J)
        END DO
        DO J=1,5
          READ(1,*) BV(I,J)
        END DO
        READ(1,*) IWALL(I)
        READ(1,*) QS(I)
        READ(1,*) Z0(I)
      END DO
      CLOSE(1)

      CALL UPDATE_BV()  !  Update boundary values when bc=10, from module boundary

      OPEN(1,FILE="sgs.in")
      READ(1,*) ISGS
      READ(1,*) CS0
      READ(1,*) IFILTER
      READ(1,*) N_SGS_SKIP
      CLOSE(1)

      IF(IINFLOW_READ.EQ.1)THEN  ! Read inflow control coefficients
        CALL INFLOW_COE()
      END IF
      
      NX=NXT/NPX
      NY=NYT/NPY
      NZ=NZT/NPZ

      END SUBROUTINE           
                                                              


