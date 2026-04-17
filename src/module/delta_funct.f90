MODULE delta_func
  USE kinds, ONLY: DP
CONTAINS
  !> 1/x with eta for numerical stability
  PURE REAL(DP) FUNCTION dE_inv(dE, eta)
    REAL(DP), INTENT(IN) :: dE, eta
    dE_inv = dE/(dE**2 + eta**2)
  END FUNCTION dE_inv
  ! ==================================================
  !> Lorentzian broadening d(x)= eta/(pi*(x^2 + eta^2))
  PURE FUNCTION delta_lorentz(x, eta)
    USE constants, ONLY: pi
    REAL(DP), INTENT(IN) :: x(:), eta
    REAL(DP) :: delta_lorentz(SIZE(x))

    delta_lorentz(:) = eta/(pi*(x(:)**2 + eta**2))
  END FUNCTION delta_lorentz
  !> Gaussian broadening d(x)= exp(-(x/eta)^2)/(sqrt(pi)*eta)
  PURE FUNCTION delta_gaussian(x, eta)
    USE constants, ONLY: sq_pi
    REAL(DP), INTENT(IN) :: x(:), eta
    REAL(DP)::delta_gaussian(SIZE(x))
    REAL(DP)::inv_eta, pref
    inv_eta = 1.0_DP/eta
    pref = inv_eta/sq_pi
    delta_gaussian(:) = pref*EXP(-((x(:)*inv_eta)**2))
  END FUNCTION delta_gaussian
END MODULE delta_func
