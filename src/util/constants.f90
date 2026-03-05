MODULE constants
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !--- Mathematical constants
  REAL(DP), PARAMETER::pi = ACOS(-1.0_DP)
  REAL(DP), PARAMETER::tpi = 2.0_DP*pi !< 2*pi
  REAL(dp), PARAMETER :: sqtpi = SQRT(tpi) !< sqrt(2pi)
  COMPLEX(DP), PARAMETER::zi = CMPLX(0.0_DP, 1.0_DP, DP) !< imaginary unit
  COMPLEX(dp), PARAMETER :: pi2zi = (0.0_DP, 2.0_DP)*pi !< 2pi*(zi)

  INTEGER, PARAMETER::vec_0(3) = (/0, 0, 0/)
  !< (0, 0, 0) vector

  ! --- Physical Constants & Unit Conversions ---
  REAL(DP), PARAMETER :: bohr2ang = 0.529177249_DP
  REAL(DP), PARAMETER :: ang2bohr = 1.0_DP/bohr2ang
  REAL(DP), PARAMETER :: ry2ev = 13.605698066_DP
  REAL(DP), PARAMETER :: ev2ry = 1.0_DP/ry2ev
  REAL(DP), PARAMETER :: ev2j = 1.602176634E-19_DP  ! 1 eV -> J
  REAL(DP), PARAMETER :: ang2m = 1.0E-10_DP         ! 1 Ang -> m

  REAL(DP), PARAMETER :: m_e = 9.10938356E-31_DP    ! electron mass [kg]
  REAL(DP), PARAMETER :: hbar = 1.054571817E-34_DP  ! [J s]
  REAL(DP), PARAMETER :: e_chg = 1.602176634E-19_DP ! [C]
END MODULE constants
