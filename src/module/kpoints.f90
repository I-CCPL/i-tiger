MODULE kpoints
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC::kpoint_type
    INTEGER::nkpt
    !< Number of k-points in this node
    INTEGER::nktot
    !< Total number of k-points
    REAL(DP)::wk
    !< Weight of the k-point (1/nkpt)
    REAL(DP), ALLOCATABLE:: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nkpt)
    REAL(DP), ALLOCATABLE:: k_red(:, :)
    !< k-points in reduced coordinates (3, nkpt)
  CONTAINS
    PROCEDURE::build_path => build_kpath
    PROCEDURE::build_mesh => build_kmesh
  END TYPE kpoint_type

  TYPE(kpoint_type), PUBLIC::kpts
CONTAINS
  SUBROUTINE build_kpath(self, npath, skp, nkpps)
    USE cell, ONLY: red2cart
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER, INTENT(IN)::npath
    REAL(DP), INTENT(IN)::skp(3, npath)
    !< Edges for each path segment (3, npath)
    INTEGER, INTENT(IN)::nkpps(npath)
    !< Number of k-points for each path segment (npath)
    INTEGER::ikpt, ipath, ikpp, nkpp
    INTEGER::nkpt

    nkpt = SUM(nkpps(1:npath - 1)) + 1
    !< Total number of k-points along the path (including the last point)
    ALLOCATE (self%k_cart(3, nkpt))
    ALLOCATE (self%k_red(3, nkpt))

    ikpt = 1
    DO ipath = 1, npath - 1
      nkpp = nkpps(ipath)
      DO ikpp = 1, nkpp
        self%k_red(:, ikpt) = skp(1:3, ipath) + &
                              (skp(1:3, ipath + 1) - skp(1:3, ipath))*REAL(ikpp - 1, DP)/REAL(nkpp, DP)
        ikpt = ikpt + 1
      END DO
    END DO
    self%k_red(:, ikpt) = skp(:, npath)
    CALL red2cart(self%k_red, self%k_cart, nkpt)
    self%nkpt = nkpt
    self%nktot = nkpt
    self%wk = 1.0_DP/REAL(nkpt, DP)
  END SUBROUTINE build_kpath
  !
  SUBROUTINE build_kmesh(self, nk1, nk2, nk3, sk1, sk2, sk3)
    USE cell, ONLY: red2cart
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER, INTENT(IN)::nk1, nk2, nk3
    !< number of k-points align axis 1,2,3
    INTEGER, INTENT(IN)::sk1, sk2, sk3
    !< shift of k-points align axis 1,2,3 (0 or 1)
    REAL(DP)::kx, ky, kz
    INTEGER::ik1, ik2, ik3, ikpt
    INTEGER::nktot
    IF (nk1 <= 0 .OR. nk2 <= 0 .OR. nk3 <= 0) THEN
      CALL errore(1, 'build_kmesh', 'nk must be positive')
    ELSE IF (sk1 < 0 .OR. sk1 > 1 &
             .OR. sk2 < 0 .OR. sk2 > 1 &
             .OR. sk3 < 0 .OR. sk3 > 1) THEN
      CALL errore(1, 'build_kmesh', 'sk must be 0 or 1')
    END IF

    nktot = nk1*nk2*nk3
    ALLOCATE (self%k_cart(3, nktot))
    ALLOCATE (self%k_red(3, nktot))
    ikpt = 0
    DO ik1 = 1, nk1
      kx = (DBLE(ik1 - 1) + DBLE(sk1)/2)/nk1
      kx = kx - NINT(kx)
      DO ik2 = 1, nk2
        ky = (DBLE(ik2 - 1) + DBLE(sk2)/2)/nk2
        ky = ky - NINT(ky)
        DO ik3 = 1, nk3
          kz = (DBLE(ik3 - 1) + DBLE(sk3)/2)/nk3
          kz = kz - NINT(kz)

          ikpt = ikpt + 1
          self%k_red(:, ikpt) = (/kx, ky, kz/)
        END DO
      END DO
    END DO
    CALL red2cart(self%k_red, self%k_cart, nktot)
    self%nkpt = nktot
    self%nktot = nktot
    self%wk = 1.0_DP/REAL(nktot, DP)
  END SUBROUTINE build_kmesh
END MODULE kpoints
