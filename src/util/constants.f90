MODULE constants
  USE kinds, ONLY: DP
  IMPLICIT NONE
  REAL(DP), PARAMETER::pi = 3.14159265358979323846_DP
  REAL(DP), PARAMETER::tpi = 2.0_DP*pi !< 2*pi
  COMPLEX(DP), PARAMETER::zi = CMPLX(0.0_DP, 1.0_DP, DP) !< imaginary unit

  INTEGER, PARAMETER::vec_0(3) = (/0, 0, 0/)
  !< (0, 0, 0) vector
END MODULE constants
