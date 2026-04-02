MODULE delta_func
  USE kinds, ONLY: DP
CONTAINS
  PURE REAL(DP) FUNCTION delta_lorentz(x, eta)
    USE constants, ONLY: pi
    REAL(DP), INTENT(IN) :: x, eta
    delta_lorentz = eta/(pi*(x*x + eta*eta))
  END FUNCTION delta_lorentz
!
!> Gaussian broadening d(x)= exp(-(x/eta)^2)/(sqrt(pi)*eta)
  PURE REAL(DP) FUNCTION delta_gaussian(x, eta)
    USE constants, ONLY: sq_pi
    REAL(DP), INTENT(IN) :: x, eta
    delta_gaussian = EXP(-(x/eta)**2)/(sq_pi*eta)
  END FUNCTION delta_gaussian
END MODULE delta_func
