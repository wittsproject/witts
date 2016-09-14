! This module is used to declare global arrays
!
  MODULE class_shared

  USE parameters, ONLY: DP
    
  IMPLICIT NONE

  TYPE :: CELL
    INTEGER:: CELL_INDEX                     ! INDEX OF THE CELL
    INTEGER:: CELL_EMID                      ! THE EMBEDDING LEVEL OF THE CELL
    INTEGER:: CELL_GHOST                     ! =1: GHOST CELL, =0: ACTIVE CELL
    INTEGER:: CELL_SPLIT                     ! =1: THIS CELL IS SPLITTED
    INTEGER:: CELL_ACTIVE                    ! =1: ACTIVE CELL, I.E., CELL_GHOST=0 & CELL_SPLIT=0
!---FOR GHOST CELL-------------------------------------------------------------------
    INTEGER:: CELL_PID                       ! RANK (PROCESSOR) ID
    INTEGER:: CELL_NEAR                      ! THE INDEX OF NEAREST NON-GHOST CELL
!---FOR EMBEDDED CELL----------------------------------------------------------------
    INTEGER:: CELL_PARENT                    ! INDEX OF THE PARENT CELL
    INTEGER:: CELL_CHILD(2,2,2)              ! INDEX OF THE CHILD CELL
!------------------------------------------------------------------------------------
    INTEGER:: CELL_NEI_X(2)                  ! THE INDEX OF NEIGHBORING CELLS
    INTEGER:: CELL_NEI_Y(2)
    INTEGER:: CELL_NEI_Z(2)

    INTEGER:: CELL_BOU_FLAG                  ! =0: INNER CELL, =1: BOUNDARY CELL
    INTEGER:: CELL_WALL                      ! =1: THE CELL ADJACENT TO WALL

    REAL(KIND=DP):: CELL_X,CELL_Y,CELL_Z     ! COORDINATES AT THE CELL CENTER

    REAL(KIND=DP):: CELL_DX,CELL_DY,CELL_DZ  ! SPACING OF THE CELL

    REAL(KIND=DP):: CELL_VEL(3)              ! VELOCITY COMPONENTS (1: U, 2: V, 3: W)
    REAL(KIND=DP):: CELL_TE                  ! TEMPERATURE
    REAL(KIND=DP):: CELL_PD                  ! DYNAMIC PRESSURE
    REAL(KIND=DP):: CELL_NU                  ! EDDY VISCOSITY
    REAL(KIND=DP):: CELL_EPS                 ! SGS DISSIPATION                             
    REAL(KIND=DP):: CELL_RHO                 ! DENSITY
    REAL(KIND=DP):: CELL_PHI                 ! LEVEL-SET FUNCTION
    REAL(KIND=DP):: CELL_FX                  ! FORCING TERMS
    REAL(KIND=DP):: CELL_FY 
    REAL(KIND=DP):: CELL_FZ 
    REAL(KIND=DP):: CELL_S(6)                ! STRAIN RATE TENSOR (1: S11, 2: S22, 3: S33, 4: S12, 5: S13, 6: S23)
    REAL(KIND=DP):: CELL_SS                  ! MODULUS OF THE STRAIN RATE TENSOR (S)
    REAL(KIND=DP):: CELL_TAU(6)              ! STRESS TENSOR (1: TAU11, 2: TAU22, 3: TAU33, 4: TAU12, 5: TAU13, 6: TAU23)
    REAL(KIND=DP):: CELL_HF(3)               ! HEAT FLUX (1: Q1, 2: Q2, 3: Q3)

    REAL(KIND=DP):: CELL_PLM                 ! LASD RELATED VARIABLES (PLM,PMM,PQN,PNN)
    REAL(KIND=DP):: CELL_PMM
    REAL(KIND=DP):: CELL_PQN
    REAL(KIND=DP):: CELL_PNN

    REAL(KIND=DP):: CELL_UF(3,2)             ! VELOCITY FLUX AT THE CELL FACES 
  END TYPE CELL

  TYPE(CELL),DIMENSION(:),ALLOCATABLE:: CELL_FV

 
  END MODULE class_shared

