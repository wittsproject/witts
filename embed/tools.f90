! This module contains numerical tools for general use purpose.
!
  MODULE tools

  USE mpi 
  USE parameters, ONLY: DP,NX,NY,NZ,NPX,NPY,NPZ,NXT,NYT,NZT, &
                       IERR,MYID,MYIDX,MYIDY,MYIDZ,LX,LY,LZ,PI, &
                       TOTAL_CELL

  USE class_shared
  USE class_cell

  IMPLICIT NONE
!  INCLUDE "fftw.f"
  
  CONTAINS

!============================================================!
!        SUBROUTINE OF GETTING THE INDEX OF PROCESSOR        !
!============================================================!
!     MYID=MYIDX+MYIDY*NPX*NPZ+MYIDZ*NPX
      SUBROUTINE GETMYID()
      IMPLICIT NONE

      MYIDX=MOD(MYID,NPX)
      MYIDY=INT(MYID/(NPX*NPZ))
      MYIDZ=INT((MYID-(NPX*NPZ)*MYIDY)/NPX)

      END SUBROUTINE GETMYID
!============================================================!
!        ASSEMBLE THE GLOBAL ARRAY WITH ACTIVE CELLS         !
!============================================================!
      SUBROUTINE ASSEM(F,F_GLOBAL,COUNT)
      IMPLICIT NONE
      INTEGER:: COUNT,COUNT_TOTAL,M  
      REAL(KIND=DP),DIMENSION(:):: F,F_GLOBAL
      REAL(KIND=DP),DIMENSION(:),ALLOCATABLE::VAR
      
      CALL MPI_ALLREDUCE(COUNT,COUNT_TOTAL,1,MPI_INTEGER,MPI_SUM, &
                         MPI_COMM_WORLD,IERR) 
      ALLOCATE(VAR(COUNT_TOTAL))

      VAR=0.0

      DO M=1,COUNT
        VAR(CELL_SKIP(COUNT)+M)=F(M)
      END DO

      CALL MPI_ALLREDUCE(VAR,F_GLOBAL,COUNT_TOTAL,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERR) 

      DEALLOCATE(VAR)
      
      END SUBROUTINE      
!============================================================!
!        SPATIAL AVERAGE OVER THE GLOBAL DOMAIN              !
!============================================================!
!     ID=1: ARITHMETIC AVERAGE
!     ID=2: GEOMETRIC AVERAGE
      SUBROUTINE AVE_S_GLOBAL(F,SI1,SI2,SI3,FA,ID)
     
      IMPLICIT NONE

      INTEGER :: I,J,K,SI1,SI2,SI3,ID,NUM
      REAL(KIND=DP), DIMENSION(SI1:,SI2:,SI3:) :: F
      REAL(KIND=DP) :: FAI,FA

      IF(ID.EQ.1)THEN
        FAI=0.0
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              FAI=FAI+F(I,J,K)
            END DO
          END DO
        END DO

        CALL MPI_ALLREDUCE(FAI,FA,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                           MPI_COMM_WORLD,IERR)
        FA=FA/(NXT*NYT*NZT)
      ELSE IF(ID.EQ.2)THEN
        FAI=1.0
        NUM=0
        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              IF(ABS(F(I,J,K)).LT.10.0d0.AND.ABS(F(I,J,K)).GT.1.0E-5)THEN
                FAI=DMIN1(FAI*F(I,J,K),1.0D100)
                NUM=NUM+1
              ELSE
                FAI=FAI*1.0d0
              END IF
            END DO
          END DO
        END DO

        IF(NUM.EQ.0)THEN
          FAI=1.0
        ELSE
          FAI=ABS(FAI)**(1.0d0/DFLOAT(NUM))
        END IF
        
        CALL MPI_ALLREDUCE(FAI,FA,1,MPI_DOUBLE_PRECISION,MPI_PROD, &
                           MPI_COMM_WORLD,IERR)
        FA=FA**(1.0d0/DFLOAT(NPX*NPY*NPZ))
      END IF

      END SUBROUTINE
!============================================================!
!         SUBROUTINE OF PLANAR AVERAGE                       !
!============================================================!
!     ID=1: ARITHMETIC AVERAGE
!     ID=2: GEOMETRIC AVERAGE
!     IDIR=1: ALONG X DIRECTION
!         =2: ALONG Y DIRECTION
!         =3: ALONG Z DIRECTION
      SUBROUTINE AVE_P(F,SI1,SI2,SI3,FA,ID,IDIR)

      IMPLICIT NONE
      INTEGER :: I,J,K,M,II,JJ,KK,SI1,SI2,SI3
      INTEGER :: ICNT,ID,IDIR,NUM
      REAL(KIND=DP):: FX,FY,FZ,SUMF,SUMA,FAT
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: F
      REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: FAI,FA
!-----ALONG X DIRECTION-------------------------------
      IF(IDIR.EQ.1)THEN
        ALLOCATE(FAI(NX))

        IF(ID.EQ.1)THEN
          DO I=1,NX
            SUMF=0.0
            SUMA=0.0
            DO K=1,NZ
              DO J=1,NY
                SUMF=SUMF+F(I,J,K)
              END DO
            END DO
            FAI(I)=SUMF/(NZ*NY)
            FA(I)=0.0
          END DO
        ELSE
          DO I=1,NX
            SUMF=0.0
            DO K=1,NZ
              DO J=1,NY
                SUMF=SUMF+LOG(DMAX1(F(I,J,K),1.0D-12))
              END DO
            END DO
            FAI(I)=SUMF/(NZ*NY)
            FA(I)=0.0
          END DO
        END IF
!       AVERAGE OVER PROCESSORS
        ICNT=0
    
 1      IF(ID.EQ.1)THEN
          DO I=1,NX
            IF(MYIDX.EQ.ICNT)THEN
              FX=FAI(I)
            ELSE
              FX=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FX,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR) 

            IF(MYIDX.EQ.ICNT)THEN            
              FA(I)=FAT/(NPZ*NPY)
            END IF
          END DO
        ELSE
          DO I=1,NX
            IF(MYIDX.EQ.ICNT)THEN
              FX=FAI(I)
            ELSE
              FX=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FX,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR) 
           
            IF(MYIDX.EQ.ICNT)THEN            
              FA(I)=EXP(FAT/(NPZ*NPY))
            END IF
          END DO 
        END IF

        ICNT=ICNT+1
        IF(ICNT.LT.NPX)THEN
          GOTO 1
        END IF     
        DEALLOCATE(FAI)
!-----ALONG Y DIRECTION-------------------------------
      ELSE IF(IDIR.EQ.2)THEN
        ALLOCATE(FAI(NY))

        IF(ID.EQ.1)THEN
          DO J=1,NY
            SUMF=0.0
            SUMA=0.0
            DO K=1,NZ
              DO I=1,NX
                SUMF=SUMF+F(I,J,K)
              END DO
            END DO
            FAI(J)=SUMF/(NZ*NX)
            FA(J)=0.0
          END DO
        ELSE
          DO J=1,NY
            SUMF=0.0
            DO K=1,NZ
              DO I=1,NX
                SUMF=SUMF+LOG(DMAX1(F(I,J,K),1.0D-12))
              END DO
            END DO
            FAI(J)=SUMF/(NZ*NX)
            FA(J)=0.0
          END DO
        END IF
!       AVERAGE OVER PROCESSORS
        ICNT=0
    
 2      IF(ID.EQ.1)THEN
          DO J=1,NY
            IF(MYIDY.EQ.ICNT)THEN
              FY=FAI(J)
            ELSE
              FY=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FY,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR)          

            IF(MYIDY.EQ.ICNT)THEN            
              FA(J)=FAT/(NPZ*NPX)
            END IF
          END DO
        ELSE
          DO J=1,NY
            IF(MYIDY.EQ.ICNT)THEN
              FY=FAI(J)
            ELSE
              FY=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FY,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR) 
          
            IF(MYIDY.EQ.ICNT)THEN            
              FA(J)=EXP(FAT/(NPZ*NPX))
            END IF
          END DO 
        END IF

        ICNT=ICNT+1
        IF(ICNT.LT.NPY)THEN
          GOTO 2
        END IF     
        DEALLOCATE(FAI)
!-----ALONG Z DIRECTION-------------------------------
      ELSE
!     GET PLANE AVERAGE ON EACH PROCESSOR
        ALLOCATE(FAI(NZ))

        IF(ID.EQ.1)THEN
          DO K=1,NZ
            SUMF=0.0
            SUMA=0.0
            DO I=1,NX
              DO J=1,NY
                SUMF=SUMF+F(I,J,K)
              END DO
            END DO
            FAI(K)=SUMF/(NX*NY)
            FA(K)=0.0
          END DO
        ELSE
          DO K=1,NZ
            SUMF=0.0
            DO I=1,NX
              DO J=1,NY
                SUMF=SUMF+LOG(DMAX1(F(I,J,K),1.0D-12))
              END DO
            END DO
            FAI(K)=SUMF/(NX*NY)
            FA(K)=0.0
          END DO
        END IF
!       AVERAGE OVER PROCESSORS
        ICNT=0
    
 3      IF(ID.EQ.1)THEN
          DO K=1,NZ
            IF(MYIDZ.EQ.ICNT)THEN
              FZ=FAI(K)
            ELSE
              FZ=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FZ,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR) 

            IF(MYIDZ.EQ.ICNT)THEN            
              FA(K)=FAT/(NPX*NPY)
            END IF
          END DO
        ELSE
          DO K=1,NZ
            IF(MYIDZ.EQ.ICNT)THEN
              FZ=FAI(K)
            ELSE
              FZ=0.0
            ENDIF
            CALL MPI_ALLREDUCE(FZ,FAT,1,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERR) 

            IF(MYIDZ.EQ.ICNT)THEN            
              FA(K)=EXP(FAT/(NPX*NPY))
            END IF
          END DO 
        END IF

        ICNT=ICNT+1
        IF(ICNT.LT.NPZ)THEN
          GOTO 3
        END IF     
        DEALLOCATE(FAI)
      END IF

      END SUBROUTINE 
!============================================================!    
!            PLANAR AVERAGE OVER THE TOTAL LENGTH            !
!============================================================!
      SUBROUTINE AVE_P_GLOBAL(F,SI1,SI2,SI3,FA,ID,IDIR)

      IMPLICIT NONE

      INTEGER :: ID,IDIR,SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: F
      REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: FAI,FA
!-----GET SPTIAL AVERAGE ON EACH PROCESSOR
      IF(IDIR.EQ.1)THEN
        ALLOCATE(FAI(NX))
        CALL AVE_P(F,SI1,SI2,SI3,FAI,ID,IDIR) 
        CALL ASSEM_1D_X(FAI,SI1,FA)    
        DEALLOCATE(FAI) 
      ELSE IF(IDIR.EQ.2)THEN
        ALLOCATE(FAI(NY))
        CALL AVE_P(F,SI1,SI2,SI3,FAI,ID,IDIR) 
        CALL ASSEM_1D_Y(FAI,SI2,FA)    
        DEALLOCATE(FAI)
      ELSE
        ALLOCATE(FAI(NZ))
        CALL AVE_P(F,SI1,SI2,SI3,FAI,ID,IDIR) 
        CALL ASSEM_1D_Z(FAI,SI3,FA)    
        DEALLOCATE(FAI)
      END IF

      END SUBROUTINE 
!============================================================!
!    SUBROUTINE OF DEALIASING USING TRUNCATING (2/3 RULE)    !
!============================================================!
      SUBROUTINE DEALIAS(U,SI1,SI2,SI3)
      IMPLICIT NONE
!      INCLUDE "fftw3.f"

      INTEGER :: I,J,K,SI1,SI2,SI3
      REAL(KIND=DP):: U(SI1:,SI2:,SI3:)
!-----FOR FFTW
      INTEGER*8 :: plan_f, plan_b
      REAL(KIND=DP):: KX,KY
      REAL(KIND=DP):: tmpv(nx,ny)
      COMPLEX*16 :: TMPVH(NX/2+1,NY),TMPVX(NX/2+1,NY),TMPVZ(NX/2+1,NY)

 !     CALL dfftw_plan_dft_r2c_2d(plan_f,nx,ny,tmpv,tmpvh,FFTW_ESTIMATE)
 !     CALL dfftw_plan_dft_c2r_2d(plan_b,nx,ny,tmpvh,tmpv,FFTW_ESTIMATE)

      DO K=1,NZ
!-----FORWARD 2D FFT IN THE HORIZONTAL PALNE-------------------------
        DO I=1,NX
          DO J=1,NY
            TMPV(I,J)=U(I,J,K)
          END DO
        END DO
!        call dfftw_execute_dft_r2c(plan_f,tmpv,tmpvh)     ! forward fft
!-----CALCULATE FIRST DERIVATIVES IN THE SPECTRAL SPACE
        DO J=1,NY
          DO I=1,NX/2+1
  	    KX=DBLE(PI*2.0/LX*(I-1))
            IF(J.LE.NY/2)THEN
  	      KY=DBLE(PI*2.0/LY*(J-1))
            ELSE
              KY=-DBLE(PI*2.0/LY*(NY-(J-1)))
            END IF
            IF(ABS(KX).GT.PI*2.0/LX*(NX/3).OR.ABS(KY).GT.PI*2.0/LY*(NY/3))THEN
              TMPVH(I,J)=0.0
            END IF
          ENDDO
        ENDDO
!-----BACKWARD 2D FFT TO GET DERIVATIVCES IN THE PHYSICAL SPACE-------
        DO I=1,NX
          DO J=1,NY
            TMPV(I,J)=0.0
          END DO
        END DO
!        call dfftw_execute_dft_c2r(plan_b,tmpvh,tmpv)     
        DO I=1,NX
          DO J=1,NY
            U(I,J,K)=TMPV(I,J)/(NX*NY)                   ! rescale the results
          END DO
        END DO  
      END DO

!      call dfftw_destroy_plan(plan_f)
!      call dfftw_destroy_plan(plan_b)

      END SUBROUTINE 
!*****************************************************************************!
!                     LAGRANGIAN POLYNIMIAL INTERPOLATION (3D)                !
!*****************************************************************************!
      SUBROUTINE INTER_GLOBAL(XI,YI,ZI,NX0,NY0,NZ0,X,Y,Z,N,F,SI1,SI2,SI3,FI)

        IMPLICIT NONE
        INTEGER :: I,J,K,NX0,NY0,NZ0,I1,I2,J1,J2,K1,K2,N,SI1,SI2,SI3
	REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: F
        REAL(KIND=DP),DIMENSION(SI1:) :: X,Y,Z
!        REAL(KIND=DP) :: X(:),Y(:),Z(:)
	REAL(KIND=DP) :: LX(N+1),LY(N+1),LZ(N+1)
	REAL(KIND=DP) :: XI,YI,ZI,FI,XT,YT,ZT
	REAL(KIND=DP) :: ZERO

	ZERO = 1.0E-12

!-----BOUND THE POINT INSIDE THE COMPUTATIONAL DOMAIN
        XT=DMIN1(DMAX1(XI,X(1)),X(NX0+1))
        YT=DMIN1(DMAX1(YI,Y(1)),Y(NY0+1))
        ZT=DMIN1(DMAX1(ZI,Z(1)),Z(NZ0+1))
!-----DETERMINE THE NEIGHBORING POINTS
	DO I=1,NX0
	  IF((X(I)-XT).LE.SQRT(ZERO).AND.X(I+1).GT.XT)THEN
            IF((I-(N+1)/2+1).LE.1)THEN
              I1=1
              I2=1+N
            ELSE IF((I+(N+1)/2).GE.NX0)THEN
              I2=NX0
              I1=NX0-N
            ELSE
	      I1=I-(N+1)/2+1
	      I2=I+(N+1)/2
            END IF
	    EXIT
	  END IF
	END DO

        DO J=1,NY0
	  IF((Y(J)-YT).LE.SQRT(ZERO).AND.Y(J+1).GT.YT)THEN 
            IF((J-(N+1)/2+1).LE.1)THEN
              J1=1
              J2=1+N
            ELSE IF((J+(N+1)/2).GE.NY0)THEN
              J2=NY0
              J1=NY0-N
            ELSE
	      J1=J-(N+1)/2+1
	      J2=J+(N+1)/2
            END IF
	    EXIT
	  END IF
	END DO

        DO K=1,NZ0
	  IF((Z(K)-ZT).LE.SQRT(ZERO).AND.Z(K+1).GT.ZT)THEN
            IF((K-(N+1)/2+1).LE.1)THEN
              K1=1
              K2=1+N
            ELSE IF((K+(N+1)/2).GE.NZ0)THEN
              K2=NZ0
              K1=NZ0-N
            ELSE
	      K1=K-(N+1)/2+1
	      K2=K+(N+1)/2
            END IF
	    EXIT
	  END IF
	END DO
!-----USE LAGRANGIAN SCHEME TO DO THE INTERPOLATION
        FI=0. 

        I1=MAX(I1,1)
        I2=MIN(I2,NX0)
        J1=MAX(J1,1)
        J2=MIN(J2,NY0)
        K1=MAX(K1,1)
        K2=MIN(K2,NZ0)

        DO I=I1,I2
	  K=I-I1+1
          LX(K)=1.
	  DO J=I1,I2
	    IF(I.NE.J)THEN
	      LX(K)=LX(K)*(XT-X(J))/(X(I)-X(J))
	    END IF
	  END DO
	END DO

	DO I=J1,J2
	  K=I-J1+1
          LY(K)=1.
	  DO J=J1,J2
	    IF(I.NE.J)THEN
	      LY(K)=LY(K)*(YT-Y(J))/(Y(I)-Y(J))
	    END IF
	  END DO
	END DO

	DO I=K1,K2
	  K=I-K1+1
          LZ(K)=1.
	  DO J=K1,K2
	    IF(I.NE.J)THEN
	      LZ(K)=LZ(K)*(ZT-Z(J))/(Z(I)-Z(J))
	    END IF
	  END DO
	END DO

        FI=0.0
	DO I=I1,I2
	  DO J=J1,J2
	    DO K=K1,K2
	      FI=FI+LX(I-I1+1)*LY(J-J1+1)*LZ(K-K1+1)*F(I,J,K)
	    END DO
	  END DO
	END DO
        
	END SUBROUTINE 
!*****************************************************************************!
!                     LAGRANGIAN POLYNIMIAL INTERPOLATION (3D)                !
!*****************************************************************************!
      SUBROUTINE INTER_CELL(XI,YI,ZI,VAR,FI)

      IMPLICIT NONE
      
      INTEGER :: M
      REAL(KIND=DP),DIMENSION(:):: VAR
      REAL(KIND=DP):: XI,YI,ZI,BX1,BX2,BY1,BY2,BZ1,BZ2
      REAL(KIND=DP):: FI0,FI

      FI0=0.0
      DO M=1,TOTAL_CELL
        IF(CELL_FV(M)%CELL_SPLIT.EQ.0)THEN
          BX1=CELL_FV(M)%CELL_X-CELL_FV(M)%CELL_DX/2.0
          BX2=CELL_FV(M)%CELL_X+CELL_FV(M)%CELL_DX/2.0
          BY1=CELL_FV(M)%CELL_Y-CELL_FV(M)%CELL_DY/2.0
          BY2=CELL_FV(M)%CELL_Y+CELL_FV(M)%CELL_DY/2.0
          BZ1=CELL_FV(M)%CELL_Z-CELL_FV(M)%CELL_DZ/2.0
          BZ2=CELL_FV(M)%CELL_Z+CELL_FV(M)%CELL_DZ/2.0
          IF(XI.GE.BX1.AND.XI.LT.BX2.AND. &
             YI.GE.BY1.AND.YI.LT.BY2.AND. &
             ZI.GE.BZ1.AND.ZI.LT.BZ2)THEN
            FI0=VAR(M)
            EXIT
          END IF
        END IF
      END DO

      CALL MPI_ALLREDUCE(FI0,FI,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,IERR)
       
      END SUBROUTINE
!*************************************************************************!
!                  SUBROUTINE OF GENERATING RANDOM NUMBERS                !
!*************************************************************************!
        SUBROUTINE RANDOM(NX,NY,NZ,U,SI1,SI2,SI3)
        IMPLICIT NONE
!        INCLUDE "mpif.h"
        INTEGER :: NX,NY,NZ,NPX,NPY,NPZ,I,J,K,SI1,SI2,SI3
        INTEGER :: MYIDX,MYIDY,MYIDZ,NXT,NYT,NZT
        REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: U
        REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: UF
        REAL(KIND=DP):: MU,SUM,DU,R

        MU=0.0   !MEAN VALUE IN SPANWISE DIRECTION

        R=0.5

        ALLOCATE(UF(NX*NPX,NY*NPY,NZ*NPZ))
        DO K=1,NZT
          DO J=1,NYT
            DO I=1,NXT
              CALL RAND_NUM(R)
              UF(I,J,K)=R
1             IF(ABS(UF(I,J,K)).GT.0.1)THEN
                UF(I,J,K)=UF(I,J,K)/1.5
                GOTO 1
              END IF
            END DO
          END DO
        END DO

        DO K=1,NZT
	  SUM=0.0
          DO I=1,NXT
            DO J=1,NYT
              SUM=SUM+UF(I,J,K)
            END DO
          END DO
          DU=MU-SUM
          DO I=1,NXT
            DO J=1,NYT
              UF(I,J,K)=UF(I,J,K)+DU/(NXT*NYT)
            END DO
          END DO  
!---------TEST
	  SUM=0.0
          DO I=1,NXT
            DO J=1,NYT
              SUM=SUM+UF(I,J,K)
            END DO
          END DO     
        END DO

        DO I=1,NX
          DO J=1,NY
            DO K=1,NZ
              U(I,J,K)=UF(I+MYIDX*NX,J+MYIDY*NY,K+MYIDZ*NZ)
            END DO
          END DO
        END DO

        DEALLOCATE(UF)
      END SUBROUTINE
!*************************************************************************!
!                  FUNCTION OF GENERATING RANDOM NUMBERS                  !
!*************************************************************************!
        SUBROUTINE RAND_NUM(R)

        IMPLICIT NONE
	REAL(KIND=DP):: S,U,V,R,M,NRND1
	S=65536.0
	U=2053.0
	V=13849.0
	M=R/S
	R=R-M*S
	R=U*R+V
	M=R/S
	R=R-M*S
	NRND1=R/S
      END SUBROUTINE 
!=========================================================================!
!                       SUBROUTINE OF LES FILTERING                       !
!=========================================================================!
!     IFILTER=1: GAUSSIAN FUNCTION
!            =2: SHARP SPECTRAL FUNCTION (IN THE PHYSICAL SPACE)
!            =3: BOX FILTER (DEFAULT)
!     NDEL: FILTER WIDTH
      SUBROUTINE FILTER(DX,DY,DZ,U,SI1,SI2,SI3,NDEL,  &
                        IFILTER)  ! THIS LINE FOR OPTIONAL ARGUMENTS

      IMPLICIT NONE
      INTEGER :: I,J,K,ID,BC_ID,SI1,SI2,SI3
      INTEGER,OPTIONAL :: IFILTER
      INTEGER :: I1,I2,J1,J2,K1,K2,NB,M
      REAL(KIND=DP):: DX,DY,DZ
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:) :: U
      REAL(KIND=DP):: UT,R,DELTA,NDEL,GS
      REAL(KIND=DP):: RX,RY,RZ,WX,WY,WZ,W
      REAL(KIND=DP):: EPS

      EPS=1.E-12

      NB=NDEL/2

      IF(PRESENT(IFILTER))THEN
        ID=IFILTER
      ELSE
        ID=3
      END IF
      
      I1=-NB
      I2= NB
      J1=-NB
      J2= NB
      K1=-NB
      K2= NB
!-----GET FILTER WIDTH
      DELTA=NDEL*(DX*DY*DZ)**(1.0/3.0)
!-----IMPLEMENT THE FILTERING
!  NUMERICAL INTEGRATION (1d): Integral=D/2*[F(1)+2*F(2)+...+2*F(N-1)+F(N)]
      UT=0.0
      GS=0.0
      DO I=I1,I2
        RX=ABS(I)*DX
        IF(I.EQ.I1.OR.I.EQ.I2)THEN
          WX=1.0
        ELSE
          WX=2.0
        END IF
       
        DO J=J1,J2
          RY=ABS(J)*DY
          IF(J.EQ.J1.OR.J.EQ.J2)THEN
            WY=1.0
          ELSE
            WY=2.0
          END IF
         
          DO K=K1,K2
            RZ=ABS(K)*DZ
            IF(K.EQ.K1.OR.K.EQ.K2)THEN
              WZ=1.0
            ELSE
              WZ=2.0
            END IF
           
            W=WX*WY*WZ
            
            UT=UT+U(I,J,K)* &
                  FILTER_G(DELTA,RX,RY,RZ,3,ID)*W*(DY/2.0)*(DX/2.0)*(DZ/2.0)
            GS=GS+FILTER_G(DELTA,RX,RY,RZ,3,ID)*W*(DY/2.0)*(DX/2.0)*(DZ/2.0)
          END DO
        END DO
      END DO
      U(0,0,0)=UT/GS

      END SUBROUTINE FILTER
!=========================================================================!
!                      FUNCTION OF FILTER FUNCTIONS                       !
!=========================================================================!
!       ID=1: GAUSSIAN FUNCTION
!       ID=2: SHARP SPECTRAL FUNCTION
!       ID=3: BOX FILTER
!       N: DIMENSION
        REAL(KIND=DP) FUNCTION FILTER_G(DELTA,RX,RY,RZ,N,ID)
        IMPLICIT NONE
        INTEGER :: ID,N
        REAL(KIND=DP):: DELTA,RX,RY,RZ,Q
        REAL(KIND=DP):: FILTER_GX,FILTER_GY,FILTER_GZ
        REAL(KIND=DP):: PI,SIG2,KC

        DATA PI /3.1415926/

        IF(ID.EQ.1)THEN
          FILTER_G=SQRT(6.0/(PI*DELTA**2.0))**N*  &
                   EXP(-6.0*RX**2.0/DELTA**2.0)*  &
                   EXP(-6.0*RY**2.0/DELTA**2.0)*  &
                   EXP(-6.0*RZ**2.0/DELTA**2.0)
        ELSE IF(ID.EQ.2)THEN
          KC=PI/DELTA
          FILTER_G=SIN(KC*RX)/(KC*MAX(REAL(RX),REAL(1.E-16)))*  &
                   SIN(KC*RY)/(KC*MAX(REAL(RY),REAL(1.E-16)))*  &
                   SIN(KC*RZ)/(KC*MAX(REAL(RZ),REAL(1.E-16)))
        ELSE IF(ID.EQ.3)THEN      ! BOX FILTER
          IF(RX.LE.DELTA/2.0)THEN
            FILTER_GX=1.0/DELTA
          ELSE
            FILTER_GX=0.0
          END IF
         
          IF(RY.LE.DELTA/2.0)THEN
            FILTER_GY=1/DELTA
          ELSE
            FILTER_GY=0.0
          END IF
         
          IF(RZ.LE.DELTA/2.0)THEN
            FILTER_GZ=1/DELTA
          ELSE
            FILTER_GZ=0.0
          END IF
         
          FILTER_G=FILTER_GX*FILTER_GY*FILTER_GZ
        END IF
      END FUNCTION
!=========================================================================!
!              SUBROUTINE OF FILTERING IN THE WAVENUMBER DOMAIN           !
!=========================================================================!
      SUBROUTINE FILTER_SPEC(DX,DY,U,SI1,SI2,SI3,NDEL)

      IMPLICIT NONE
!      INCLUDE "fftw3.f"
      INTEGER :: I,J,K,SI1,SI2,SI3                                   
      REAL(KIND=DP):: U(SI1:,SI2:,SI3:) 
      REAL(KIND=DP):: DEL,NDEL
      REAL(KIND=DP):: DX,DY,GK,KC
      REAL(KIND=DP):: KCX,KCY,GKX,GKY
!-----FOR FFTW
      INTEGER*8 :: plan_f, plan_b
      REAL(KIND=DP):: KX,KY
      REAL(KIND=DP):: tmpv(nx,ny)
      COMPLEX*16 :: TMPVH(NX/2+1,NY),TMPVX(NX/2+1,NY),TMPVZ(NX/2+1,NY)     

!      CALL dfftw_plan_dft_r2c_2d(plan_f,nx,ny,tmpv,tmpvh,FFTW_ESTIMATE)
!      CALL dfftw_plan_dft_c2r_2d(plan_b,nx,ny,tmpvh,tmpv,FFTW_ESTIMATE)

      DEL=SQRT(DX*DY)
      KC=PI/(NDEL*DEL)

      DO K=1,NZ
!-----FORWARD 2D FFT IN THE HORIZONTAL PALNE-------------------------
        DO I=1,NX
          DO J=1,NY
            TMPV(I,J)=U(I,J,K)
          END DO
        END DO
!        call dfftw_execute_dft_r2c(plan_f,tmpv,tmpvh)     ! forward fft
!-----DO THE FILTERING USING SHARP SPECTRAL FILTER
        DO J=1,NY
          DO I=1,NX/2+1
  	    KX=DBLE(PI*2.0/LX*(I-1))
            IF(J.LE.NY/2)THEN
  	      KY=DBLE(PI*2.0/LY*(J-1))
            ELSE
              KY=-DBLE(PI*2.0/LY*(NY-(J-1)))
            END IF

            IF(SQRT(KX**2+KY**2)-KC.LT.0.0)THEN
              GK=1.0
            ELSE
              GK=0.0
            END IF
            TMPVH(I,J)=TMPVH(I,J)*GK
          END DO
        END DO
!-----BACKWARD 2D FFT TO GET FILTERED FUNCTION IN THE PHYSICAL SPACE-------
        DO I=1,NX
          DO J=1,NY
            TMPV(I,J)=0.0
          END DO
        END DO
!        call dfftw_execute_dft_c2r(plan_b,tmpvh,tmpv)
        DO I=1,NX
          DO J=1,NY
            U(I,J,K)=TMPV(I,J)/(NX*NY)                   ! rescale the results
          END DO
        END DO    
      END DO

!      call dfftw_destroy_plan(plan_f)
!      call dfftw_destroy_plan(plan_b)

      END SUBROUTINE 
!=========================================================================!
!	        SUBROUTINE FOR SOLVING TRIDIAGONAL SYSTEM  		  !
!=========================================================================!
!       THOMAS ALGORITHM IS USED
	SUBROUTINE THOMAS(A0,B0,C0,D0,X,N)
	IMPLICIT NONE
	INTEGER :: N,K
	REAL(KIND=DP):: A0(N),B0(N),C0(N),D0(N)
	REAL(KIND=DP):: A(N),B(N),C(N),D(N),X(N)
	REAL(KIND=DP):: M

	DO K=1,N
	  A(K)=A0(K)
	  B(K)=B0(K)
	  C(K)=C0(K)
	  D(K)=D0(K)
	END DO
!       FORWARD ELIMINATION PHASE
        DO K=2,N
	  M=A(K)/B(K-1)
	  B(K)=B(K)-M*C(K-1)
	  D(K)=D(K)-M*D(K-1)
	END DO
!       BACKWARD SUBSTITUTION PHASE
        X(N)=D(N)/B(N)
	DO K=N-1,1,-1
	  X(K)=(D(K)-C(K)*X(K+1))/B(K)
	END DO
        END SUBROUTINE 
!==========================================================!
!              FUNCTION OF NUMERICAL DAMPING               !
!==========================================================!
      REAL(KIND=DP) FUNCTION SITA(X,X0,LS,ID)
      USE parameters, ONLY: PI
      IMPLICIT NONE
      REAL(KIND=DP):: X,X0,LS
      REAL(KIND=DP):: A,Z,RLS_TIME
      INTEGER :: ID
      
      IF(ID.EQ.1)THEN  ! GAUSSIAN FUNCTION
        A=1.0
        Z=(X-X0)/LS*3.5
        SITA=A*EXP(-Z**2/2.0)
      ELSE IF(ID.EQ.2)THEN ! COSINE FUNCTION (PORTE-AGEL)
        RLS_TIME=10.0   ! RELAXATION TIME SCALE, CONSTANT 60s IS SET
        SITA=(1.0-COS(PI*(X-X0)/LS))/RLS_TIME
      END IF
      END FUNCTION
!**************************DERIVATIVES********************************!
!---------------------------------------------------!
!      OBTAIN THE DERIVATIVE OF CELL VARIABLES      !
!---------------------------------------------------!
! INDEX: INDEX OF THE CELL
! ID=1: X DERIVATIVE; =2: Y DERIVATIVE; =3: Z DERIVATIVE
! SCHEME=1: CENTRAL SCHEME; =2 UPWIND SCHEME
! ISTAG=1: GET DERIVATIVE AT A STAGGERED LOCATION
!      =0: GET DERIVATIVE AT THE EXACT LOCATION
!      =10: STAGGERED, BUT WITH OFFSET=1      
! ORDER: ORDER OF NUMERICAL ACCURACY    
    REAL(KIND=DP) FUNCTION DERIV_CELL(VAR_CELL,INDEX,ID,SCHEME,ISTAG,ORDER)
    IMPLICIT NONE
    INTEGER:: INDEX,ID,SCHEME,ISTAG,ORDER
    INTEGER:: I,J,NB,OFFSET
    INTEGER:: IND(-ORDER+1:ORDER-1)
    REAL(KIND=DP),DIMENSION(:):: VAR_CELL
    REAL(KIND=DP),DIMENSION(:,:,:), ALLOCATABLE:: VAR

    NB=ORDER-1
    
    ALLOCATE(VAR(-NB:NB,-NB:NB,-NB:NB))

    OFFSET=0

    IF(ISTAG.EQ.10)THEN
      OFFSET=1
    END IF

    CALL CELL_TO_STRUCT(VAR_CELL,INDEX,NB,VAR)    
!---X DERIVATIVE---------------------------    
    IF(ID.EQ.1)THEN      
      IF(SCHEME.EQ.1)THEN             
        DERIV_CELL=DERIV_X(VAR,-NB,-NB,-NB,OFFSET,0,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DX)
      ELSE
        !  Upwind scheme should be added here  
      END IF
!---Y DERIVATIVE---------------------------
    ELSE IF(ID.EQ.2)THEN
      IF(SCHEME.EQ.1)THEN       
        DERIV_CELL=DERIV_Y(VAR,-NB,-NB,-NB,0,OFFSET,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DY)
      ELSE
        !  Upwind scheme should be added here  
      END IF
!---Z DERIVATIVE---------------------------
    ELSE      
      IF(SCHEME.EQ.1)THEN     
        DERIV_CELL=DERIV_Z(VAR,-NB,-NB,-NB,0,OFFSET,0,ISTAG,ORDER,CELL_FV(INDEX)%CELL_DZ)
      ELSE
        !  Upwind scheme should be added here  
      END IF      
    END IF
   
    DEALLOCATE(VAR)
    
    END FUNCTION DERIV_CELL     
!==========================================================!
!           DERIVATIVE STENCIL IN X DIRECTION              !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,ND,DX)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DX

    IF(ND.EQ.1)THEN
      STENCIL_DX=(VAR(I,J,K)-VAR(I-1,J,K))/DX
    ELSE IF(ND.EQ.2)THEN
      STENCIL_DX=(VAR(I+1,J,K)-VAR(I-1,J,K))/(DX*2.0)
    ELSE IF(ND.EQ.3)THEN
      STENCIL_DX=(VAR(I+1,J,K)-VAR(I-2,J,K))/(DX*3.0)
    ELSE IF(ND.EQ.4)THEN
      STENCIL_DX=(VAR(I+2,J,K)-VAR(I-2,J,K))/(DX*4.0)
    END IF

    END FUNCTION
!==========================================================!
!           DERIVATIVE STENCIL IN Y DIRECTION              !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,ND,DY)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DY

    IF(ND.EQ.1)THEN
      STENCIL_DY=(VAR(I,J,K)-VAR(I,J-1,K))/DY
    ELSE IF(ND.EQ.2)THEN
      STENCIL_DY=(VAR(I,J+1,K)-VAR(I,J-1,K))/(DY*2.0)
    ELSE IF(ND.EQ.3)THEN
      STENCIL_DY=(VAR(I,J+1,K)-VAR(I,J-2,K))/(DY*3.0)
    ELSE IF(ND.EQ.4)THEN
      STENCIL_DY=(VAR(I,J+2,K)-VAR(I,J-2,K))/(DY*4.0)
    END IF

    END FUNCTION
!==========================================================!
!           DERIVATIVE STENCIL IN Z DIRECTION              !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,ND,DZ)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DZ

    IF(ND.EQ.1)THEN
      STENCIL_DZ=(VAR(I,J,K)-VAR(I,J,K-1))/DZ
    ELSE IF(ND.EQ.2)THEN
      STENCIL_DZ=(VAR(I,J,K+1)-VAR(I,J,K-1))/(DZ*2.0)
    ELSE IF(ND.EQ.3)THEN
      STENCIL_DZ=(VAR(I,J,K+1)-VAR(I,J,K-2))/(DZ*3.0)
    ELSE IF(ND.EQ.4)THEN
      STENCIL_DZ=(VAR(I,J,K+2)-VAR(I,J,K-2))/(DZ*4.0)
   END IF

   END FUNCTION STENCIL_DZ
!==========================================================!
!        DERIVATIVE IN X DIRECTION (CENTRAL SCHEME)        !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A STAGGERED LOCATION
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION DERIV_X(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER,DX)
    IMPLICIT NONE
    INTEGER :: I,J,K,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DX 

    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        DERIV_X=STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,1,DX)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_X=STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,1,DX)*9.0/8.0- &
                STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,3,DX)*1.0/8.0 
      END IF
    ELSE 
      IF(ORDER.EQ.2)THEN
        DERIV_X=STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,2,DX)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_X=STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,2,DX)*4.0/3.0- &
                STENCIL_DX(VAR,SI1,SI2,SI3,I,J,K,4,DX)*1.0/3.0 
      END IF 
    END IF
    END FUNCTION
!==========================================================!
!        DERIVATIVE IN Y DIRECTION (CENTRAL SCHEME)        !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A STAGGERED LOCATION
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION DERIV_Y(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER,DY)
    IMPLICIT NONE
    INTEGER :: I,J,K,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DY      

    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        DERIV_Y=STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,1,DY)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_Y=STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,1,DY)*9.0/8.0- &
                STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,3,DY)*1.0/8.0 
      END IF
    ELSE 
      IF(ORDER.EQ.2)THEN
        DERIV_Y=STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,2,DY)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_Y=STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,2,DY)*4.0/3.0- &
                STENCIL_DY(VAR,SI1,SI2,SI3,I,J,K,4,DY)*1.0/3.0 
      END IF 
    END IF
    END FUNCTION
!==========================================================!
!        DERIVATIVE IN Z DIRECTION (CENTRAL SCHEME)        !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A STAGGERED LOCATION
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION DERIV_Z(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER,DZ)
    IMPLICIT NONE
    INTEGER :: I,J,K,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: DZ

    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        DERIV_Z=STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,1,DZ)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_Z=STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,1,DZ)*9.0/8.0- &
                STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,3,DZ)*1.0/8.0 
      END IF
    ELSE 
      IF(ORDER.EQ.2)THEN
        DERIV_Z=STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,2,DZ)
      ELSE IF(ORDER.EQ.4)THEN
        DERIV_Z=STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,2,DZ)*4.0/3.0- &
                STENCIL_DZ(VAR,SI1,SI2,SI3,I,J,K,4,DZ)*1.0/3.0 
      END IF 
    END IF
    END FUNCTION
!**************************INTERPOLATION********************************!
!==========================================================!
!           INTERPOLATION STENCIL IN X DIRECTION           !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_IX(VAR,SI1,SI2,SI3,I,J,K,ND)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR

    IF(ND.EQ.1)THEN
      STENCIL_IX=(VAR(I,J,K)+VAR(I-1,J,K))/2.0
    ELSE IF(ND.EQ.2)THEN
      STENCIL_IX=(VAR(I+1,J,K)+VAR(I-1,J,K))/2.0
    ELSE IF(ND.EQ.3)THEN
      STENCIL_IX=(VAR(I+1,J,K)+VAR(I-2,J,K))/2.0
    ELSE IF(ND.EQ.4)THEN
      STENCIL_IX=(VAR(I+2,J,K)+VAR(I-2,J,K))/2.0
    END IF

    END FUNCTION
!==========================================================!
!           DERIVATIVE STENCIL IN Y DIRECTION              !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_IY(VAR,SI1,SI2,SI3,I,J,K,ND)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR

    IF(ND.EQ.1)THEN
      STENCIL_IY=(VAR(I,J,K)+VAR(I,J-1,K))/2.0
    ELSE IF(ND.EQ.2)THEN
      STENCIL_IY=(VAR(I,J+1,K)+VAR(I,J-1,K))/2.0
    ELSE IF(ND.EQ.3)THEN
      STENCIL_IY=(VAR(I,J+1,K)+VAR(I,J-2,K))/2.0
    ELSE IF(ND.EQ.4)THEN
      STENCIL_IY=(VAR(I,J+2,K)+VAR(I,J-2,K))/2.0
    END IF

    END FUNCTION
!==========================================================!
!           DERIVATIVE STENCIL IN Z DIRECTION              !
!==========================================================!
    REAL(KIND=DP) FUNCTION STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,ND)
    IMPLICIT NONE
    INTEGER :: I,J,K,ISTAG,ND,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
 
    IF(ND.EQ.1)THEN
      STENCIL_IZ=(VAR(I,J,K)+VAR(I,J,K-1))/2.0
    ELSE IF(ND.EQ.2)THEN
      STENCIL_IZ=(VAR(I,J,K+1)+VAR(I,J,K-1))/2.0
    ELSE IF(ND.EQ.3)THEN
      STENCIL_IZ=(VAR(I,J,K+1)+VAR(I,J,K-2))/2.0
    ELSE IF(ND.EQ.4)THEN
      STENCIL_IZ=(VAR(I,J,K+2)+VAR(I,J,K-2))/2.0
    END IF

    END FUNCTION
!==========================================================!
!         INTERPOLATION IN X DIR (CENTRAL SCHEME)          !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A DISTANCE D/2 BEFORE
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION VAR_INTER_X(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER)
    IMPLICIT NONE
    INTEGER :: I,J,K,IDIRE,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR

    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        VAR_INTER_X=STENCIL_IX(VAR,SI1,SI2,SI3,I,J,K,1)
      ELSE IF(ORDER.EQ.4)THEN
        VAR_INTER_X=STENCIL_IX(VAR,SI1,SI2,SI3,I,J,K,1)*9.0/8.0- &
                    STENCIL_IX(VAR,SI1,SI2,SI3,I,J,K,3)*1.0/8.0 
      END IF
    ELSE 
      VAR_INTER_X=VAR(I,J,K) 
    END IF   
    END FUNCTION
!==========================================================!
!         INTERPOLATION IN Y DIR (CENTRAL SCHEME)          !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A DISTANCE D/2 BEFORE
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION VAR_INTER_Y(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER)
    IMPLICIT NONE
    INTEGER :: I,J,K,IDIRE,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
   
    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        VAR_INTER_Y=STENCIL_IY(VAR,SI1,SI2,SI3,I,J,K,1)
      ELSE IF(ORDER.EQ.4)THEN
        VAR_INTER_Y=STENCIL_IY(VAR,SI1,SI2,SI3,I,J,K,1)*9.0/8.0- &
                    STENCIL_IY(VAR,SI1,SI2,SI3,I,J,K,3)*1.0/8.0 
      END IF
    ELSE 
      VAR_INTER_Y=VAR(I,J,K) 
    END IF   
    END FUNCTION
!==========================================================!
!         INTERPOLATION IN Z DIR (CENTRAL SCHEME)          !
!==========================================================!
!   ISTAG=1: GET DERIVATIVE AT A DISTANCE D/2 BEFORE
!        =0: GET DERIVATIVE AT THE EXACT LOCATION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION VAR_INTER_Z(VAR,SI1,SI2,SI3,I,J,K,ISTAG,ORDER)
    IMPLICIT NONE
    INTEGER :: I,J,K,IDIRE,ORDER,ISTAG,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    
    IF(ISTAG.EQ.1)THEN
      IF(ORDER.EQ.2)THEN
        VAR_INTER_Z=STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,1)
      ELSE IF(ORDER.EQ.4)THEN
        VAR_INTER_Z=STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,1)*9.0/8.0- &
                    STENCIL_IZ(VAR,SI1,SI2,SI3,I,J,K,3)*1.0/8.0 
      END IF
    ELSE 
      VAR_INTER_Z=VAR(I,J,K) 
    END IF   
    END FUNCTION
!--------------------------------------------------------------!
!            GET THE X INDEX NEAREST TO A POINT                !
!--------------------------------------------------------------!
    INTEGER FUNCTION INDEX_X(X0,XA)
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:):: XA
    REAL(KIND=DP):: X0
    INTEGER:: I
      
    INDEX_X=0   
    DO I=1,NX
      IF(XA(I+MYIDX*NX).LE.X0.AND.XA(I+MYIDX*NX+1).GT.X0)THEN
        INDEX_X=I
        EXIT
      END IF
    END DO     
    END FUNCTION
!--------------------------------------------------------------!
!            GET THE Y INDEX NEAREST TO A POINT                !
!--------------------------------------------------------------!
    INTEGER  FUNCTION INDEX_Y(Y0,YA)
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: YA
    REAL(KIND=DP):: Y0
    INTEGER:: I
      
    INDEX_Y=0  
    DO I=1,NY
      IF(YA(I+MYIDY*NY).LE.Y0.AND.YA(I+MYIDY*NY+1).GT.Y0)THEN
        INDEX_Y=I
        EXIT
      END IF
    END DO    
    END FUNCTION
!--------------------------------------------------------------!
!            GET THE Z INDEX NEAREST TO A POINT                !
!--------------------------------------------------------------!
    INTEGER  FUNCTION INDEX_Z(Z0,ZA)
    IMPLICIT NONE
    REAL(KIND=DP),DIMENSION(:),ALLOCATABLE:: ZA
    REAL(KIND=DP):: Z0
    INTEGER:: I
      
    INDEX_Z=0  
    DO I=1,NZ
      IF(ZA(I+MYIDZ*NZ).LE.Z0.AND.ZA(I+MYIDZ*NZ+1).GT.Z0)THEN
        INDEX_Z=I
        EXIT
      END IF
    END DO     
    END FUNCTION
!--------------------------------------------------------------!
!		        FUNCTION OF RMS                        !
!--------------------------------------------------------------!
      REAL(KIND=DP) FUNCTION RMS(F,NUM)
      IMPLICIT NONE      
      INTEGER:: NUM
      REAL(KIND=DP),DIMENSION(NUM)::F
      INTEGER:: M,NSUMT
      REAL(KIND=DP):: RMSF

      RMSF=0.0
      DO M=1,NUM        
        RMSF=RMSF+F(M)**2
      END DO

      CALL MPI_ALLREDUCE(RMSF,RMS,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(NUM,NSUMT,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,IERR)

      RMS=SQRT(RMS/NSUMT)

      END FUNCTION
!--------------------------------------------------------------!
!		          FUNCTION OF MAX                      !
!--------------------------------------------------------------!
      REAL(KIND=DP) FUNCTION MAXF(F,SI1,SI2,SI3,IM,JM,KM)
      IMPLICIT NONE      
      INTEGER:: SI1,SI2,SI3
      REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:)::F
      INTEGER:: IM,JM,KM,I,J,K
      REAL(KIND=DP):: MAXFF

      MAXFF=0.0
      DO I=1,IM
        DO J=1,JM
          DO K=1,KM
            MAXFF=MAX(ABS(F(I,J,K)),MAXFF)
          END DO
        END DO
      END DO

      CALL MPI_ALLREDUCE(MAXFF,MAXF,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                         MPI_COMM_WORLD,IERR)
     END FUNCTION MAXF
!---------------------------------------------------!
!          FUNCTION OF MATRIX MULTIPLICATION        !
!---------------------------------------------------!
    FUNCTION MATMUL_2D_1D(A,B,RANK1,RANK2)
    IMPLICIT NONE
    INTEGER:: RANK1,RANK2  
    REAL(KIND=DP):: A(1:RANK1,1:RANK2),B(1:RANK2)
    REAL(KIND=DP),DIMENSION(1:RANK1):: MATMUL_2D_1D
    INTEGER:: I,J
    
    DO I=1,RANK1
      MATMUL_2D_1D(I)=0.0
      DO J=1,RANK2
        MATMUL_2D_1D(I)=MATMUL_2D_1D(I)+A(I,J)*B(J)
      END DO
    END DO

    END FUNCTION MATMUL_2D_1D
!---------------------------------------------------!
!          FUNCTION OF MATRIX MULTIPLICATION        !
!---------------------------------------------------!
    FUNCTION MATMUL_1D_1D(A,B,RANK)
    IMPLICIT NONE
    INTEGER:: RANK
    REAL(KIND=DP):: A(1:RANK),B(1:RANK)
    REAL(KIND=DP):: MATMUL_1D_1D
    INTEGER:: J
    
    MATMUL_1D_1D=0.0
    DO J=1,RANK
      MATMUL_1D_1D=MATMUL_1D_1D+A(J)*B(J)
    END DO
    
    END FUNCTION MATMUL_1D_1D
  
  END MODULE
