MODULE constants
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !--- Mathematical constants
  REAL(DP), PARAMETER::pi = ACOS(-1.0_DP)
  REAL(DP), PARAMETER::tpi = 2.0_DP*pi !< 2*pi
  REAL(DP), PARAMETER::sq_pi = SQRT(pi) !< sqrt(pi)
  REAL(DP), PARAMETER :: sq_tpi = SQRT(tpi) !< sqrt(2pi)

  COMPLEX(DP), PARAMETER::cmplx_0 = CMPLX(0.0_DP, 0.0_DP, DP)
  COMPLEX(DP), PARAMETER::cmplx_1 = CMPLX(1.0_DP, 0.0_DP, DP)
  COMPLEX(DP), PARAMETER::cmplx_i = CMPLX(0.0_DP, 1.0_DP, DP)

  INTEGER, PARAMETER::vec_0(3) = (/0, 0, 0/)
  !< (0, 0, 0) vector

  ! --- Physical Constants & Unit Conversions ---
  ! ARU: Atomic Rydberg Unit
  ! AU : Atomic Unit (Hartree)
  ! Unit Usage
  ! - Charge: elementary charge (e)
  REAL(DP), PARAMETER :: e_chg_au = 1.0_DP
  REAL(DP), PARAMETER :: e_chg_si = 1.602176634E-19_DP     ! [C]
  REAL(DP), PARAMETER :: e_chg_au2si = e_chg_si/e_chg_au  ! [C/e]
  REAL(DP), PARAMETER :: e_chg_si2au = e_chg_au/e_chg_si  ! [e/C]

  ! - Length: Angstrom (Ang)
  REAL(DP), PARAMETER :: bohr2ang = 0.529177249_DP
  REAL(DP), PARAMETER :: ang2bohr = 1.0_DP/bohr2ang
  REAL(DP), PARAMETER :: Ang2m = 1.0E-10_DP
  REAL(DP), PARAMETER :: m2Ang = 1.0E10_DP

  ! - Time  : femto second (fs)
  REAL(DP), PARAMETER :: SEC2FS = 1.0E15_DP           ! second to femto second
  REAL(DP), PARAMETER :: FS2SEC = 1.0E-15_DP          ! femto second to second
  REAL(DP), PARAMETER :: t_aru2fs = 0.048377687_DP

  ! - Energy: electron volt (eV)
  REAL(DP), PARAMETER :: ry2ev = 13.605698066_DP
  REAL(DP), PARAMETER :: ev2ry = 1.0_DP/ry2ev
  REAL(DP), PARAMETER :: ev2j = 1.602176634E-19_DP
  REAL(DP), PARAMETER :: j2ev = 1.0_DP/ev2j

  REAL(DP), PARAMETER :: h_SI = 6.62607015E-34_DP     !  [J s]
  REAL(DP), PARAMETER :: hbar_SI = h_SI/tpi           ! [J s]
  REAL(DP), PARAMETER :: h_evfs = 4.135667696_DP      ! Planck constant [eV fs]
  REAL(DP), PARAMETER :: hbar_eVfs = 0.6582119514_DP  ! [eV fs]

  REAL(DP), PARAMETER :: m_e_kg = 9.10938356E-31_DP               ! electron mass [kg]
  REAL(DP), PARAMETER :: m_e = m_e_kg*j2ev*(SEC2FS**2)/(m2Ang**2) ! electron mass [eV fs^2 / Ang^2]

  !> Vacuum permittivity
  REAL(DP), PARAMETER :: epsilon_0_SI = 8.8541878128E-12_DP             ! SI units [C/V/m]
  REAL(DP), PARAMETER :: epsilon_0_au = 1.0_DP/(4*pi)                   ! atomic units
  REAL(DP), PARAMETER :: epsilon_0 = epsilon_0_SI*e_chg_si2au/(m2Ang)   ! [e/V/Ang]
END MODULE constants
