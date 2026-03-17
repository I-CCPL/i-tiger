MODULE constants
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !--- Mathematical constants
  REAL(DP), PARAMETER::pi = ACOS(-1.0_DP)
  REAL(DP), PARAMETER::tpi = 2.0_DP*pi !< 2*pi
  REAL(DP), PARAMETER :: sqtpi = SQRT(tpi) !< sqrt(2pi)
  COMPLEX(DP), PARAMETER::zero = CMPLX(0.0_DP, 0.0_DP, DP)
  COMPLEX(DP), PARAMETER::zi = CMPLX(0.0_DP, 1.0_DP, DP) !< imaginary unit
  COMPLEX(DP), PARAMETER :: pi2zi = (0.0_DP, 2.0_DP)*pi !< 2pi*(zi)

  INTEGER, PARAMETER::vec_0(3) = (/0, 0, 0/)
  !< (0, 0, 0) vector

  ! --- Physical Constants & Unit Conversions ---
  ! Unit Usage
  ! - Length: Angstrom (Ang)
  ! - Time  : femto second (fs)
  ! - Charge: elementary charge (e)
  ! - Energy: electron volt (eV)
  REAL(DP), PARAMETER :: bohr2ang = 0.529177249_DP
  REAL(DP), PARAMETER :: ang2bohr = 1.0_DP/bohr2ang
  REAL(DP), PARAMETER :: ry2ev = 13.605698066_DP
  REAL(DP), PARAMETER :: ev2ry = 1.0_DP/ry2ev
  REAL(DP), PARAMETER :: ev2j = 1.602176634E-19_DP
  REAL(DP), PARAMETER :: j2ev = 1.0_DP/ev2j
  REAL(DP), PARAMETER :: Ang2m = 1.0E-10_DP
  REAL(DP), PARAMETER :: m2Ang = 1.0E10_DP

  REAL(DP), PARAMETER :: h_SI = 6.62607015E-34_DP     !  [J s]
  REAL(DP), PARAMETER :: hbar_SI = h_SI/tpi           ! [J s]
  REAL(DP), PARAMETER :: h_evfs = 4.135667696_DP     ! Planck constant [eV fs]
  REAL(DP), PARAMETER :: hbar_evfs = 0.6582119514_DP  ! [eV fs]

  REAL(DP), PARAMETER :: SEC2FS = 1.0E15_DP           ! second to femto second
  REAL(DP), PARAMETER :: FS2SEC = 1.0E-15_DP          ! femto second to second

  REAL(DP), PARAMETER :: m_e_kg = 9.10938356E-31_DP        ! electron mass [kg]
  REAL(DP), PARAMETER :: m_e = m_e_kg*j2ev*(SEC2FS**2)/(m2Ang**2) ! electron mass [eV fs^2 / Ang^2]

  REAL(DP), PARAMETER :: e_chg = 1.602176634E-19_DP     ! [C]
END MODULE constants
