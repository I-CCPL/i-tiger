MODULE delta_func
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTERFACE w1gauss
    MODULE PROCEDURE w1gauss_scalar, w1gauss_vec
  END INTERFACE w1gauss
CONTAINS
  !> 1/x with eta for numerical stability
  PURE REAL(DP) FUNCTION dE_inv(dE, eta)
    REAL(DP), INTENT(IN) :: dE, eta
    dE_inv = dE/(dE**2 + eta**2)
  END FUNCTION dE_inv
  ! ==================================================
  PURE FUNCTION w0gauss(ldX, x, eta, n) RESULT(w0)
    !< Returns the occupation factor for a given energy difference x and broadening eta. \
    !< (n=-99): Fermi-Dirac distribution f(x) = 1/(exp(x/eta) + 1) \
    !< (n=0): Gaussian broadening f(x) = 0.5*erfc(x/eta) \
    INTEGER, INTENT(in) :: ldX
    REAL(dp), INTENT(in) :: x(ldX), eta
    INTEGER, INTENT(in) :: n
    REAL(dp) :: w0(ldX)
    IF (n == -99) THEN
      w0(:) = 1.0_DP/(EXP(x(:)/eta) + 1.0_DP)
    ELSE IF (n == 0) THEN
      w0(:) = 0.5_DP*ERFC(x(:)/eta)
    END IF
  END FUNCTION w0gauss
  ! ==================================================
  PURE FUNCTION w1gauss_scalar(x, eta, n) RESULT(w1)
    USE constants, ONLY: pi, sq_pi
    !< Returns the Delta function for a given energy difference x and broadening eta.
    !< (Derivative of -w0gauss) \
    !< (n=-99): Derivative of Fermi-Dirac distribution
    !<     d(x)= (1/4*eta)*sech^2(x/(2*eta)) \
    !< (n=-2): Lorentzian broadening d(x)= eta/(pi*(x^2 + eta^2)) \
    !< (n=0): Gaussian broadening d(x)= exp(-(x/eta)^2)/(sqrt(pi)*eta) \
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(in) :: x, eta
    REAL(dp) :: w1
    REAL(dp) :: inv_eta, pref
    IF (n == -99) THEN
      inv_eta = 1.0_DP/eta
      pref = inv_eta/4.0_DP
      w1 = pref/(COSH(x*inv_eta/2.0_DP)**2)
    ELSE IF (n == -2) THEN
      w1 = eta/(pi*(x**2 + eta**2))
    ELSE IF (n == 0) THEN
      inv_eta = 1.0_DP/eta
      pref = inv_eta/sq_pi
      w1 = pref*EXP(-((x*inv_eta)**2))
    END IF
  END FUNCTION w1gauss_scalar
  PURE FUNCTION w1gauss_vec(ldX, x, eta, n) RESULT(w1)
    USE constants, ONLY: pi, sq_pi
    !< Returns the Delta function for a given energy difference x and broadening eta.
    !< (Derivative of -w0gauss) \
    !< (n=-99): Derivative of Fermi-Dirac distribution
    !<     d(x)= (1/4*eta)*sech^2(x/(2*eta)) \
    !< (n=-2): Lorentzian broadening d(x)= eta/(pi*(x^2 + eta^2)) \
    !< (n=0): Gaussian broadening d(x)= exp(-(x/eta)^2)/(sqrt(pi)*eta) \
    INTEGER, INTENT(in) :: ldX, n
    REAL(dp), INTENT(in) :: x(ldX), eta
    REAL(dp) :: w1(ldX)
    REAL(dp) :: inv_eta, pref
    IF (n == -99) THEN
      inv_eta = 1.0_DP/eta
      pref = inv_eta/4.0_DP
      w1(:) = pref/(COSH(x(:)*inv_eta/2.0_DP)**2)
    ELSE IF (n == -2) THEN
      w1(:) = eta/(pi*(x(:)**2 + eta**2))
    ELSE IF (n == 0) THEN
      inv_eta = 1.0_DP/eta
      pref = inv_eta/sq_pi
      w1(:) = pref*EXP(-((x(:)*inv_eta)**2))
    END IF
  END FUNCTION w1gauss_vec
END MODULE delta_func
