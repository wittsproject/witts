! This is the user defined function of boundary condition
!
   MODULE UDF_BOUNDARY

   USE parameters

   CONTAINS
!---------------------------------------------------!
!        Function of calculating boundary value	    !
!---------------------------------------------------!
     REAL(KIND=DP) FUNCTION BOUNDARY_VALUE()

     IMPLICIT NONE

     REAL(KIND=DP) :: VB0,VB_RATE

     VB0=265.0	 
     VB_RATE=-0.0000694444

     BOUNDARY_VALUE=VB0+TIME*VB_RATE

     END FUNCTION
 
   END MODULE


     
