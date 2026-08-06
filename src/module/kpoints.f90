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
    !< Weight of the k-point (1/nktot)
    REAL(DP), ALLOCATABLE:: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nktot)
    REAL(DP), ALLOCATABLE:: k_red(:, :)
    !< k-points in reduced coordinates (3, nktot)
    COMPLEX(DP), ALLOCATABLE::H_k(:, :)
    !< Hamiltonian in k-space (Nw, Nw)
    REAL(DP), ALLOCATABLE::eigval(:, :)
    !< Eigenvalues (Nw, nktot)
    COMPLEX(DP), ALLOCATABLE::eigvec(:, :)
    !< Eigenvectors (Nw, Nw)
  CONTAINS
    PROCEDURE::divide_k => divide_k_idx
    PROCEDURE::global_k => global_k_idx
    PROCEDURE::gather_r => gather_r_data
    PROCEDURE::gather_c => gather_c_data
    PROCEDURE::receive => rec_r_data
    PROCEDURE::build_path => build_kpath
    PROCEDURE::build_mesh => build_kmesh
    PROCEDURE::rotate => rotate_k
  END TYPE kpoint_type

  TYPE(kpoint_type), PUBLIC::t_kpt
  !< Interpolated k-point list.
  INTEGER, PUBLIC::t_iks
  !< Interpolated k-point index. ikpt for coarse, iks for dense.
  INTERFACE
    MODULE SUBROUTINE divide_k_idx(self)
      CLASS(kpoint_type), INTENT(INOUT)::self
    END SUBROUTINE divide_k_idx
    MODULE FUNCTION global_k_idx(self, k_local) RESULT(k_global)
      CLASS(kpoint_type), INTENT(IN)::self
      INTEGER, INTENT(IN)::k_local
      INTEGER :: k_global
    END FUNCTION global_k_idx
    MODULE SUBROUTINE gather_r_data(self, length, f_in, f_out)
      CLASS(kpoint_type), INTENT(INOUT) :: self
      INTEGER, INTENT(IN) :: length
      REAL(DP), INTENT(IN) :: f_in(length, self%nkpt)
      REAL(DP), INTENT(OUT) :: f_out(length, self%nktot)
    END SUBROUTINE gather_r_data
    MODULE SUBROUTINE gather_c_data(self, length, f_in, f_out)
      CLASS(kpoint_type), INTENT(INOUT) :: self
      INTEGER, INTENT(IN) :: length
      COMPLEX(DP), INTENT(IN) :: f_in(length, self%nkpt)
      COMPLEX(DP), INTENT(OUT) :: f_out(length, self%nktot)
    END SUBROUTINE gather_c_data
    MODULE SUBROUTINE rec_r_data(self, length, vec)
      CLASS(kpoint_type), INTENT(INOUT) :: self
      INTEGER, INTENT(IN) :: length
      REAL(DP), INTENT(INOUT) :: vec(length, self%nktot)
    END SUBROUTINE rec_r_data
  END INTERFACE
  INTEGER, PUBLIC, ALLOCATABLE::recvcounts(:)
  INTEGER, PUBLIC, ALLOCATABLE::displs(:)
CONTAINS
  ! ==================================================
  !... K point generation
  ! ==================================================
  SUBROUTINE build_kpath(self, npath, skp, nkpps)
    USE system, ONLY: red2cart_recip
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER, INTENT(IN)::npath
    REAL(DP), INTENT(IN)::skp(3, npath)
    !< Edges for each path segment (3, npath)
    INTEGER, INTENT(IN)::nkpps(npath)
    !< Number of k-points for each path segment (npath)
    INTEGER::ikpt, ipath, ikpp, nkpp
    INTEGER::nktot

    nktot = SUM(nkpps(1:npath - 1)) + 1
    !< Total number of k-points along the path (including the last point)
    IF (ALLOCATED(self%k_cart)) THEN
      CALL errore(1, 'build_kpath', 'k-points are already allocated.')
    END IF
    ALLOCATE (self%k_cart(3, nktot))
    ALLOCATE (self%k_red(3, nktot))

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
    CALL red2cart_recip(nktot, self%k_red, self%k_cart)
    self%nkpt = nktot
    self%nktot = nktot
    self%wk = 1.0_DP/REAL(nktot, DP)
  END SUBROUTINE build_kpath
  !
  SUBROUTINE build_kmesh(self, nk1, nk2, nk3, sk1, sk2, sk3)
    USE system, ONLY: red2cart_recip
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
    IF (ALLOCATED(self%k_cart)) THEN
      CALL errore(1, 'build_kmesh', 'k-points are already allocated.')
    END IF
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
    CALL red2cart_recip(nktot, self%k_red, self%k_cart)
    self%nkpt = nktot
    self%nktot = nktot
    self%wk = 1.0_DP/REAL(nktot, DP)
  END SUBROUTINE build_kmesh

  SUBROUTINE rotate_k(self, mat_in, mat_out)
    USE system, ONLY: Nw
    CLASS(kpoint_type), INTENT(INOUT)::self
    COMPLEX(DP), INTENT(IN)::mat_in(..)
    COMPLEX(DP), INTENT(OUT)::mat_out(..)
    INTEGER::ldX, ldY
    IF (.NOT. ALLOCATED(self%eigvec)) THEN
      CALL errore(1, 'rotate_k', 'Eigenvectors are not allocated.')
    END IF
    ldX = SIZE(mat_in)
    ldY = SIZE(mat_out)
    IF (ldX /= ldY) THEN
      CALL errore(1, 'rotate_k', 'Invalid size.')
    END IF
    ldX = ldX/Nw/Nw

    CALL rotate_3d(ldX, self%eigvec(:, :), mat_in, mat_out)
  END SUBROUTINE rotate_k
END MODULE kpoints

SUBROUTINE rotate_3d(ldX, eigvec, mat_in, mat_out)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0
  USE system, ONLY: Nw
  IMPLICIT NONE
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN)::eigvec(Nw, Nw)
  COMPLEX(DP), INTENT(IN)::mat_in(Nw, Nw, ldX)
  COMPLEX(DP), INTENT(OUT)::mat_out(Nw, Nw, ldX)
  COMPLEX(DP)::U_dag(Nw, Nw)
  INTEGER::idx
  CALL start_clock('rotate_3d')
  U_dag = TRANSPOSE(CONJG(eigvec))
  DO idx = 1, ldx
    mat_out(:, :, idx) = MATMUL( &
                         MATMUL(U_dag, mat_in(:, :, idx)) &
                         , eigvec)
  END DO
  CALL stop_clock('rotate_3d')
END SUBROUTINE rotate_3d
