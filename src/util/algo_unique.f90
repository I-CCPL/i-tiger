MODULE algo_unique
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  LOGICAL FUNCTION lex_less(a, b)
    !< Lexicographical comparison of two 3D vectors
    REAL(DP), INTENT(IN) :: a(3), b(3)
    IF (a(1) < b(1)) THEN
      lex_less = .TRUE.
    ELSEIF (a(1) > b(1)) THEN
      lex_less = .FALSE.
    ELSEIF (a(2) < b(2)) THEN
      lex_less = .TRUE.
    ELSEIF (a(2) > b(2)) THEN
      lex_less = .FALSE.
    ELSE
      lex_less = a(3) < b(3)
    END IF
  END FUNCTION

  !==================================================

  SUBROUTINE unique_vec3(v, tol, v_unique, nu)
    USE kinds, ONLY: eq_vec_real
    REAL(DP), INTENT(IN) :: v(:, :)
    !< v: list of 3D vectors (3, Nvec)
    REAL(DP), INTENT(IN) :: tol
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: v_unique(:, :)
    INTEGER, INTENT(OUT) :: nu
    REAL(DP), ALLOCATABLE :: tmp(:, :)
    INTEGER :: n, i

    n = SIZE(v, 2)

    ALLOCATE (tmp(3, n))
    tmp = v

    ! 1. sort lexicographically
    CALL qsort_vec(tmp, 1, n)

    ! 2. Number of unique
    nu = 1
    DO i = 2, n
      IF (.NOT. eq_vec_real(tmp(:, i), tmp(:, i - 1), tol)) THEN
        nu = nu + 1
      END IF
    END DO

    ALLOCATE (v_unique(3, nu))

    ! 3. extract unique
    v_unique(:, 1) = tmp(:, 1)
    nu = 1
    DO i = 2, n
      IF (.NOT. eq_vec_real(tmp(:, i), tmp(:, i - 1), tol)) THEN
        nu = nu + 1
        v_unique(:, nu) = tmp(:, i)
      END IF
    END DO

  END SUBROUTINE
  !--------------------------------------------------
  RECURSIVE SUBROUTINE qsort_vec(v, left, right)
    REAL(DP), INTENT(INOUT)::v(:, :)
    !< vector list (3, Nvec)
    INTEGER, INTENT(IN) :: left, right
    INTEGER :: i, j
    REAL(DP) :: pivot(3), temp(3)

    IF (left >= right) RETURN

    pivot = v(:, (left + right)/2)
    i = left
    j = right

    DO
      DO WHILE (lex_less(v(:, i), pivot))
        i = i + 1
      END DO
      DO WHILE (lex_less(pivot, v(:, j)))
        j = j - 1
      END DO
      IF (i <= j) THEN
        temp = v(:, i)
        v(:, i) = v(:, j)
        v(:, j) = temp
        i = i + 1
        j = j - 1
      END IF
      IF (i > j) EXIT
    END DO

    IF (left < j) CALL qsort_vec(v, left, j)
    IF (i < right) CALL qsort_vec(v, i, right)
  END SUBROUTINE qsort_vec

  !==================================================

  SUBROUTINE unique_vec3_inv(v, tol, u, nu, inv)
    USE kinds, ONLY: eq_vec_real
    REAL(DP), INTENT(IN) :: v(:, :)
    !< v: list of 3D vectors (3, Nvec)
    REAL(DP), INTENT(IN) :: tol
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: u(:, :)
    !< unique of v (3, nu)
    INTEGER, INTENT(OUT) :: nu
    INTEGER, INTENT(OUT) :: inv(:)
    !< mapping from v to u (Nvec)
    INTEGER, ALLOCATABLE::perm(:)
    INTEGER :: n, i

    n = SIZE(v, 2)

    ALLOCATE (perm(n))
    DO i = 1, n
      perm(i) = i
    END DO

    ! 1. sort permutation only
    CALL qsort_vec_perm(v, 1, n, perm)

    ! 2. Number of unique
    nu = 1
    DO i = 2, n
      IF (.NOT. eq_vec_real(v(:, perm(i)), v(:, perm(i - 1)), tol)) THEN
        nu = nu + 1
      END IF
    END DO

    ! extract unique and mapping
    ALLOCATE (u(3, nu))
    nu = 1
    inv(perm(1)) = nu
    u(:, nu) = v(:, perm(1))
    DO i = 2, n
      IF (.NOT. eq_vec_real(v(:, perm(i)), v(:, perm(i - 1)), tol)) THEN
        nu = nu + 1
        u(:, nu) = v(:, perm(i))
      END IF
      inv(perm(i)) = nu
    END DO
  END SUBROUTINE unique_vec3_inv
  !--------------------------------------------------
  RECURSIVE SUBROUTINE qsort_vec_perm(v, left, right, perm)
    REAL(DP), INTENT(IN)::v(:, :)
    INTEGER, INTENT(IN) :: left, right
    INTEGER, INTENT(INOUT) :: perm(:)
    INTEGER::i, j, piv_idx, tmp
    REAL(DP):: pivot(3)

    IF (left >= right) RETURN

    piv_idx = perm((left + right)/2)
    pivot = v(:, piv_idx)
    i = left
    j = right

    DO
      DO WHILE (lex_less(v(:, perm(i)), pivot))
        i = i + 1
      END DO
      DO WHILE (lex_less(pivot, v(:, perm(j))))
        j = j - 1
      END DO

      IF (i <= j) THEN
        tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
        i = i + 1
        j = j - 1
      END IF
      IF (i > j) EXIT
    END DO

    IF (left < j) CALL qsort_vec_perm(v, left, j, perm)
    IF (i < right) CALL qsort_vec_perm(v, i, right, perm)
  END SUBROUTINE qsort_vec_perm
END MODULE algo_unique
