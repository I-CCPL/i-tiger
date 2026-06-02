
MODULE f_params
  USE kinds, ONLY: DP, inf_DP
  !... b for boolean flags
  !... m for matrices
  IMPLICIT NONE
  !> User input parameters
  !... itg
  CHARACTER(LEN=8)::FFT_conv = 'simple'
  !< 'simple' 'TB' or 'Test' FFT convention for R2k
  INTEGER::convention = 0
  !< 0 for consider degenerate R0s and nearest R0+T (build_ws)
  !< 2 for consider nearest R0+T+shift_vectors (build_R)
  CHARACTER(LEN=20)::formula = 'gauge'
  !< 'gauge' or 'projection'

  LOGICAL::lBand = .FALSE.
  ! LOGICAL::lDOS = .FALSE.
  ! LOGICAL::lPDOS = .FALSE.
  ! INTEGER::dos_dE
  ! INTEGER::dos_Emin
  ! INTEGER::dos_Emax
  LOGICAL::lOAM = .FALSE.
  LOGICAL::lOAM_g = .FALSE.

  LOGICAL::lBerry = .FALSE.
  LOGICAL::lBerry_p = .FALSE.
  LOGICAL::lBerry_g = .FALSE.

  LOGICAL::lBCD = .FALSE.
  LOGICAL::lBCD_p = .FALSE.
  !< Berry Curvature Dipole
  REAL(DP)::dE_thr = 1D-8
  !< threshold for identifying degenerate states in eV
  REAL(DP)::dE_eta = 0.04
  !< broadening parameter for 1/dE
  REAL(DP)::E_fermi = 0.0_DP
  !< Now use for occupation number.

  REAL(DP)::Ef_min = inf_DP
  REAL(DP)::Ef_max = -inf_DP
  REAL(DP)::Ef_step = -1.0_DP
  INTEGER::Ef_nE
  !< Fermi energy range for BCD calculation.

  !... Nonlinear optics calculation parameters
  !... Frequency: eV units from input
  LOGICAL::lNLO = .FALSE.
  LOGICAL::lNLO_g = .FALSE.
  LOGICAL::lshift_g = .FALSE.
  REAL(DP)::NLO_Emin = 0.0_DP
  REAL(DP)::NLO_Emax = 0.0_DP
  REAL(DP)::NLO_dE = 0.0_DP
  INTEGER::NLO_nE
  REAL(DP)::NLO_eta = 0.01_DP
  REAL(DP)::NLO_w_thr = 5.0_DP
  !< speeding up frequency integration

  LOGICAL::lshift_E_g = .FALSE.
  REAL(DP)::shift_hw = -1.0_DP
  LOGICAL::ldielec_E_g = .FALSE.

END MODULE f_params
