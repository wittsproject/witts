! This module is used to calculate fluxes on each cell.
!
  MODULE flux
  
  USE parameters
  USE field_shared
  USE tools
  
  CONTAINS
!==========================================================!
!                GET VELOCITY FLUX AT FACES                !
!==========================================================!
!   IFLUX=1: CENTRAL SCHEME
!   IFLUX=2: MUSCL SCHEME
!   THIS SUBROUTINE IS ONLY USED WHEN THE COLLOCATED GRID IS USED (ICOLL=1)
    SUBROUTINE FLUX_VELOCITY(DX,DY,DZ,ORDER)
    IMPLICIT NONE
    INTEGER :: I,J,K,ORDER
    REAL(KIND=DP) :: DX,DY,DZ

    IF(ICOLL.EQ.1)THEN
      IF(IFLUX.EQ.1)THEN
        DO I=NX1+1,NX2-1
          DO J=NY1+1,NY2-1
            DO K=NZ1+1,NZ2-1
              UF(I,J,K)=VAR_INTER_X(U,NX1,NY1,NZ1,I,J,K,1,ORDER)
              VF(I,J,K)=VAR_INTER_Y(V,NX1,NY1,NZ1,I,J,K,1,ORDER)
              WF(I,J,K)=VAR_INTER_Z(W,NX1,NY1,NZ1,I,J,K,1,ORDER)
            END DO
          END DO
        END DO
      ELSE IF(IFLUX.EQ.2)THEN
        DO I=NX1+1,NX2-1
          DO J=NX1+1,NY2-1
            DO K=NZ1+1,NZ2-1
              UF(I,J,K)=FLUX_MUSCL(U,NX1,NY1,NZ1,U(I,J,K),I,J,K,1,ORDER,DX)
              VF(I,J,K)=FLUX_MUSCL(V,NX1,NY1,NZ1,V(I,J,K),I,J,K,2,ORDER,DY)
              WF(I,J,K)=FLUX_MUSCL(W,NX1,NY1,NZ1,W(I,J,K),I,J,K,3,ORDER,DZ)
            END DO
          END DO
        END DO       
      END IF
    ELSE   ! FOR THE STAGGERED GRID, VELOCITY FLUX HAS ALREADY BEEN DEFINED
      UF=U
      VF=V
      WF=W
    END IF

    END SUBROUTINE
!==========================================================!
!                      MUSCL SCHEME                        !
!==========================================================!
!   IDIRE=1: X DIRECTION; =2: Y DIRECTION; =3: Z DIRECTION
!   ORDER: ORDER OF ACCURACY OF THE CENTRAL SCHEME
    REAL(KIND=DP) FUNCTION FLUX_MUSCL(VAR,SI1,SI2,SI3,ADV,I,J,K,IDIRE,ORDER,D)
    IMPLICIT NONE
    INTEGER :: I,J,K,IDIRE,ORDER,SI1,SI2,SI3
    REAL(KIND=DP),DIMENSION(SI1:,SI2:,SI3:):: VAR
    REAL(KIND=DP) :: ADV,ZERO,BETA,R,D
    REAL(KIND=DP) :: VAR_F_L,VAR_F_R
    REAL(KIND=DP) :: VAR_F_UP,VAR_F_CE,VAR_F_HO
    
    DATA ZERO /1.E-8/

    IF(IDIRE.EQ.1)THEN          ! X DIRECTION
      VAR_F_L=VAR(I-1,J,K)+DERIV_X(VAR,SI1,SI2,SI3,I-1,J,K,0,ORDER,D)*D/2.0
      VAR_F_R=VAR(I,  J,K)-DERIV_X(VAR,SI1,SI2,SI3,I,J,K,0,ORDER,D)*D/2.0
      IF(ABS(VAR(I+1,J,K)-VAR(I,J,K)).GT.ZERO)THEN
        R=(VAR(I,J,K)-VAR(I-1,J,K))/(VAR(I+1,J,K)-VAR(I,J,K))
      ELSE
        R=(VAR(I,J,K)-VAR(I-1,J,K))/ZERO
      END IF
    ELSE IF(IDIRE.EQ.2)THEN     ! Y DIRECTION
      VAR_F_L=VAR(I,J-1,K)+DERIV_Y(VAR,SI1,SI2,SI3,I,J-1,K,0,ORDER,D)*D/2.0
      VAR_F_R=VAR(I,J,  K)-DERIV_Y(VAR,SI1,SI2,SI3,I,J,K,0,ORDER,D)*D/2.0
      IF(ABS(VAR(I,J+1,K)-VAR(I,J,K)).GT.ZERO)THEN
        R=(VAR(I,J,K)-VAR(I,J-1,K))/(VAR(I,J+1,K)-VAR(I,J,K))
      ELSE
        R=(VAR(I,J,K)-VAR(I,J-1,K))/ZERO
      END IF
    ELSE                        ! Z DIRECTION
      VAR_F_L=VAR(I,J,K-1)+DERIV_Z(VAR,SI1,SI2,SI3,I,J,K-1,0,ORDER,D)*D/2.0
      VAR_F_R=VAR(I,J,K  )-DERIV_Z(VAR,SI1,SI2,SI3,I,J,K,0,ORDER,D)*D/2.0
      IF(ABS(VAR(I,J,K+1)-VAR(I,J,K)).GT.ZERO)THEN
        R=(VAR(I,J,K)-VAR(I,J,K-1))/(VAR(I,J,K+1)-VAR(I,J,K))
      ELSE
        R=(VAR(I,J,K)-VAR(I,J,K-1))/ZERO
      END IF
    END IF
!   GET RECONSTRUCTED UPWIND TERM
    IF(ADV.GT.ZERO)THEN
      VAR_F_UP=VAR_F_L
    ELSE IF(ADV.LT.-ZERO)THEN
      VAR_F_UP=VAR_F_R
    ELSE
      VAR_F_UP=(VAR_F_L+VAR_F_R)/2.0
    END IF
!   GET THE CENTRAL DIFFERENCE TERM
    VAR_F_CE=(VAR_F_L+VAR_F_R)/2.0
!   GET THE BLENDED TERM
    VAR_F_HO=BETA*VAR_F_UP+(1.0-BETA)*VAR_F_CE
!   IMPLEMENT FLUX LIMITER 
    FLUX_MUSCL=VAR_F_UP+FL(R,ILIMIT)*(VAR_F_HO-VAR_F_UP)

    END FUNCTION
!==========================================================!
!                       FLUX LIMITER                       !
!==========================================================!
    REAL(KIND=DP) FUNCTION FL(R,ILIMIT)
    IMPLICIT NONE
    INTEGER,OPTIONAL :: ILIMIT 
    INTEGER :: ID
    REAL(KIND=DP) :: R,ZERO

    DATA ZERO /1E-8/

    ID=1
    IF(PRESENT(ILIMIT))THEN
      ID=ILIMIT
    END IF

    IF(ID.EQ.1)THEN   ! CHARM SCHEME (Zhou, 1995)  (DEFAULT)
      IF(R.GT.0)THEN
        FL=R*(R*3.0+1.0)/((R+1.0)**2+ZERO)
      ELSE
        FL=0.0
      END IF
    ELSE IF(ID.EQ.2)THEN ! HCUS SCHEME 
      FL=1.5*(R+ABS(R))/(R+2.0)
    ELSE IF(ID.EQ.3)THEN ! HQUICK SCHEME 
      FL=2.0*(R+ABS(R))/(R+3.0)
    ELSE IF(ID.EQ.4)THEN ! KOREN SCHEME
      FL=DMAX1(0.0D0,DMIN1(2.0D0*R,DMIN1((2.0+R)/3.0,2.0D0)))
    ELSE IF(ID.EQ.5)THEN ! MINMOD SCHEME
      FL=DMAX1(0.0D0,DMIN1(1.0D0,R))
    ELSE IF(ID.EQ.6)THEN ! MC SCHEME
      FL=DMAX1(0.0D0,DMIN1(2.0D0*R,DMIN1((1.0+R)/2.0,2.0D0)))
    ELSE IF(ID.EQ.7)THEN ! OSPRE SCHEME
      FL=1.5*(R**2+R)/(R**2+R+1)
    END IF

    END FUNCTION

    END MODULE
