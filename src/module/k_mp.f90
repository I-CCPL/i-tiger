SUBMODULE(kpoints) k_mp
  USE mp_global, ONLY: mp_rank, mp_size
  IMPLICIT NONE
CONTAINS
  MODULE SUBROUTINE divide_k_idx(self)
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER::ik_start, ik_end, i, q, r, m
    q = self%nktot/mp_size
    r = MOD(self%nktot, mp_size)
    self%nkpt = q
    IF (mp_rank < r) THEN
      self%nkpt = self%nkpt + 1
    END IF

    ik_start = global_k_idx(self, 1)
    ik_end = ik_start + self%nkpt - 1
    self%k_cart(:, 1:self%nkpt) = self%k_cart(:, ik_start:ik_end)
    self%k_red(:, 1:self%nkpt) = self%k_red(:, ik_start:ik_end)

    IF (mp_size > 1 .AND. mp_rank == 0) THEN
      ALLOCATE (recvcounts(mp_size))
      ALLOCATE (displs(mp_size))
      DO i = 1, mp_size
        IF (i <= r) THEN
          recvcounts(i) = q + 1
        ELSE
          recvcounts(i) = q
        END IF
        m = MIN(i - 1, r)
        displs(i) = q*(i - 1) + m
      END DO
    END IF
  END SUBROUTINE divide_k_idx

  MODULE FUNCTION global_k_idx(self, k_local) RESULT(k_global)
    CLASS(kpoint_type), INTENT(IN)::self
    INTEGER, INTENT(IN)::k_local
    INTEGER :: k_global
    INTEGER::q, r
    q = self%nktot/mp_size
    r = MOD(self%nktot, mp_size)
    k_global = mp_rank*q + MIN(mp_rank, r) + k_local
  END FUNCTION global_k_idx

  MODULE SUBROUTINE gather_l_data(self, length, f_in, f_out)
    USE kinds, ONLY: DP
    USE mp_global
    USE kpoints, ONLY: kpoint_type, &
                       recvcounts, displs
    IMPLICIT NONE
    CLASS(kpoint_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN)::length
    LOGICAL, INTENT(IN) :: f_in(length, self%nkpt)
    LOGICAL, INTENT(OUT) :: f_out(length, self%nktot)
#if defined (__MPI)
    INTEGER, ALLOCATABLE::f_recvcounts(:)
    INTEGER, ALLOCATABLE::f_displs(:)
    INTEGER::ikpt, info
    IF (mp_size == 1) THEN
      f_out(:, 1:self%nktot) = f_in(:, 1:self%nkpt)
      RETURN
    END IF

    IF (mp_rank == mp_root) THEN
      ALLOCATE (f_recvcounts(mp_size))
      ALLOCATE (f_displs(mp_size))
      f_recvcounts = recvcounts*length
      f_displs = displs*length
    ELSE
      ALLOCATE (f_recvcounts(0))
      ALLOCATE (f_displs(0))
    END IF

    CALL MPI_GATHERV(f_in, length*self%nkpt, MPI_LOGICAL, &
                     f_out, f_recvcounts, f_displs, MPI_LOGICAL, &
                     mp_root, mp_comm, info)
    IF (info /= 0) THEN
      CALL errore(info, 'gather_l_data', 'info<>0 in MPI_GATHERV')
    END IF

    IF (ALLOCATED(f_recvcounts)) DEALLOCATE (f_recvcounts)
    IF (ALLOCATED(f_displs)) DEALLOCATE (f_displs)
#endif
    RETURN
  END SUBROUTINE gather_l_data
  MODULE SUBROUTINE gather_r_data(self, length, f_in, f_out)
    USE kinds, ONLY: DP
    USE mp_global
    USE kpoints, ONLY: kpoint_type, &
                       recvcounts, displs
    IMPLICIT NONE
    CLASS(kpoint_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN)::length
    REAL(DP), INTENT(IN) :: f_in(length, self%nkpt)
    REAL(DP), INTENT(OUT) :: f_out(length, self%nktot)
#if defined (__MPI)
    INTEGER, ALLOCATABLE::f_recvcounts(:)
    INTEGER, ALLOCATABLE::f_displs(:)
    INTEGER::ikpt, info
    IF (mp_size == 1) THEN
      f_out(:, 1:self%nktot) = f_in(:, 1:self%nkpt)
      RETURN
    END IF

    IF (mp_rank == mp_root) THEN
      ALLOCATE (f_recvcounts(mp_size))
      ALLOCATE (f_displs(mp_size))
      f_recvcounts = recvcounts*length
      f_displs = displs*length
    ELSE
      ALLOCATE (f_recvcounts(0))
      ALLOCATE (f_displs(0))
    END IF

    CALL MPI_GATHERV(f_in, length*self%nkpt, MPI_DOUBLE_PRECISION, &
                     f_out, f_recvcounts, f_displs, MPI_DOUBLE_PRECISION, &
                     mp_root, mp_comm, info)
    IF (info /= 0) THEN
      CALL errore(info, 'gather_r_data', 'info<>0 in MPI_GATHERV')
    END IF

    IF (ALLOCATED(f_recvcounts)) DEALLOCATE (f_recvcounts)
    IF (ALLOCATED(f_displs)) DEALLOCATE (f_displs)
#endif
    RETURN
  END SUBROUTINE gather_r_data
  MODULE SUBROUTINE gather_c_data(self, length, f_in, f_out)
    USE kinds, ONLY: DP
    USE mp_global
    USE kpoints, ONLY: kpoint_type, &
                       recvcounts, displs
    IMPLICIT NONE

    CLASS(kpoint_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: length
    COMPLEX(DP), INTENT(IN) :: f_in(length, self%nkpt)
    COMPLEX(DP), INTENT(OUT) :: f_out(length, self%nktot)

#if defined (__MPI)
    INTEGER, ALLOCATABLE :: f_recvcounts(:)
    INTEGER, ALLOCATABLE :: f_displs(:)
    INTEGER :: info

    IF (mp_size == 1) THEN
      f_out(:, 1:self%nktot) = f_in(:, 1:self%nkpt)
      RETURN
    END IF

    IF (mp_rank == mp_root) THEN
      ALLOCATE (f_recvcounts(mp_size))
      ALLOCATE (f_displs(mp_size))

      f_recvcounts = recvcounts*length
      f_displs = displs*length
    ELSE
      ALLOCATE (f_recvcounts(0))
      ALLOCATE (f_displs(0))
    END IF

    CALL MPI_GATHERV( &
      f_in, length*self%nkpt, MPI_DOUBLE_COMPLEX, &
      f_out, f_recvcounts, f_displs, MPI_DOUBLE_COMPLEX, &
      mp_root, mp_comm, info)

    IF (info /= 0) THEN
      CALL errore(info, 'gather_c_data', &
                  'info<>0 in MPI_GATHERV')
    END IF

    IF (ALLOCATED(f_recvcounts)) DEALLOCATE (f_recvcounts)
    IF (ALLOCATED(f_displs)) DEALLOCATE (f_displs)
#else
    f_out(:, 1:self%nktot) = f_in(:, 1:self%nkpt)
#endif

    RETURN
  END SUBROUTINE gather_c_data
  MODULE SUBROUTINE rec_r_data(self, length, vec)
    USE kinds, ONLY: DP
    USE mp_global
    IMPLICIT NONE
    CLASS(kpoint_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: length
    REAL(DP), INTENT(INOUT) :: vec(length, self%nktot)
#if defined (__MPI)
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: i, q, r, info
    !
    IF (mp_size <= 1) RETURN
    CALL mp_barrier()
    !
    IF (mp_root /= mp_rank) THEN
      CALL MPI_SEND(vec, (length*self%nkpt), MPI_DOUBLE_PRECISION, &
                    0, 17, mp_comm, info)
      CALL errore(info, 'rec_r_data', 'info<>0 in send')
    END IF
    !
    DO i = 2, mp_size
      IF (mp_root == mp_rank) THEN
        CALL MPI_RECV(vec(1, displs(i) + 1), &
                      (length*recvcounts(i)), MPI_DOUBLE_PRECISION, &
                      (i - 1), 17, mp_comm, status, info)
        CALL errore(info, 'rec_r_data', 'info<>0 in recv')
      END IF
    END DO
#endif
    RETURN
  END SUBROUTINE rec_r_data
END SUBMODULE
