MODULE itg_k
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE itg_R, ONLY: R_vec
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  !> intermediate data in k space
  !> X_bar = U^+ X U
  !> X_k_H = X_bar only for Gauge-covariant X
  TYPE::k_data_type
    !... Hamiltonian and its derivatives
    COMPLEX(DP), ALLOCATABLE::mH_k_W(:, :)
    !< Hamiltonian (Nw, Nw)
    LOGICAL::bdH_k_W
    COMPLEX(DP), ALLOCATABLE::mdH_k_W(:, :, :)
    !< Hamiltonian derivatives (Nw, Nw, 3)
    LOGICAL::bdH_bar
    COMPLEX(DP), ALLOCATABLE::mdH_bar(:, :, :)
    !< Hamiltonian derivatives (Nw, Nw, 3)
    LOGICAL::bd2H_k_W
    COMPLEX(DP), ALLOCATABLE::md2H_k_W(:, :, :, :)
    !< Hamiltonian second derivatives (Nw, Nw, 3, 3)
    LOGICAL::bd2H_bar
    COMPLEX(DP), ALLOCATABLE::md2H_bar(:, :, :, :)
    !< Hamiltonian second derivatives (Nw, Nw, 3, 3)
    LOGICAL::bD_bar
    COMPLEX(DP), ALLOCATABLE::mD_bar(:, :, :)
    !< D_bar = U^+ dU = -dH_k_H/dE (Nw, Nw, 3)

    !... Berry connection and its derivatives
    LOGICAL::bA_k_W
    COMPLEX(DP), ALLOCATABLE::mA_k_W(:, :, :)
    !< Berry connection (Nw, Nw, 3)
    LOGICAL::bA_bar
    COMPLEX(DP), ALLOCATABLE::mA_bar(:, :, :)
    !< Berry connection (Nw, Nw, 3)
    LOGICAL::bdA_k_W
    COMPLEX(DP), ALLOCATABLE::mdA_k_W(:, :, :, :)
    !< Berry connection derivatives (Nw, Nw, 3, 3)
    LOGICAL::bdA_bar
    COMPLEX(DP), ALLOCATABLE::mdA_bar(:, :, :, :)
    !< Berry connection derivatives (Nw, Nw, 3, 3)
    LOGICAL::bd2A_k_W
    COMPLEX(DP), ALLOCATABLE::md2A_k_W(:, :, :, :, :)
    !< Berry connection second derivatives (Nw, Nw, 3, 3, 3)
    LOGICAL::bd2A_bar
    COMPLEX(DP), ALLOCATABLE::md2A_bar(:, :, :, :, :)
    !< Berry connection second derivatives (Nw, Nw, 3, 3, 3)
    LOGICAL::bO_k_W
    COMPLEX(DP), ALLOCATABLE::mO_k_W(:, :, :)
    !< Curl of Berry connection (Nw, Nw, 3)
    LOGICAL::bO_bar
    COMPLEX(DP), ALLOCATABLE::mO_bar(:, :, :)
    !< Curl of Berry connection (Nw, Nw, 3)
    LOGICAL::bdO_k_W
    COMPLEX(DP), ALLOCATABLE::mdO_k_W(:, :, :, :)
    !< Curl of Berry connection derivatives (Nw, Nw, 3, 3)
    LOGICAL::bdO_bar
    COMPLEX(DP), ALLOCATABLE::mdO_bar(:, :, :, :)
    !< Curl of Berry connection derivatives (Nw, Nw, 3, 3)
  END TYPE k_data_type
  TYPE(k_data_type)::k_data

CONTAINS
  SUBROUTINE set_k_flag()
    USE itg_R, ONLY: R_data
    IF (k_data%bD_bar) k_data%bdH_bar = .TRUE.
    IF (k_data%bdH_bar) k_data%bdH_k_W = .TRUE.
    IF (k_data%bd2H_bar) k_data%bd2H_k_W = .TRUE.

    IF (k_data%bA_bar) k_data%bA_k_W = .TRUE.
    IF (k_data%bdA_bar) k_data%bdA_k_W = .TRUE.
    IF (k_data%bd2A_bar) k_data%bd2A_k_W = .TRUE.
    IF (k_data%bO_bar) k_data%bO_k_W = .TRUE.
    IF (k_data%bdO_bar) k_data%bdO_k_W = .TRUE.
    !
    IF (k_data%bA_k_W) R_data%bA_R = .TRUE.
    IF (k_data%bdA_k_W) R_data%bA_R = .TRUE.
    IF (k_data%bd2A_k_W) R_data%bA_R = .TRUE.
    IF (k_data%bO_k_W) R_data%bA_R = .TRUE.
    IF (k_data%bdO_k_W) R_data%bA_R = .TRUE.
  END SUBROUTINE set_k_flag
  !
  SUBROUTINE allocate_k()
    ALLOCATE (t_kpt%H_k(Nw, Nw))
    ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
    ALLOCATE (t_kpt%eigvec(Nw, Nw))
    t_kpt%eigval = 0.0_DP

    IF (k_data%bdH_k_W) ALLOCATE (k_data%mdH_k_W(Nw, Nw, 3))
    IF (k_data%bdH_bar) ALLOCATE (k_data%mdH_bar(Nw, Nw, 3))
    IF (k_data%bd2H_k_W) ALLOCATE (k_data%md2H_k_W(Nw, Nw, 3, 3))
    IF (k_data%bd2H_bar) ALLOCATE (k_data%md2H_bar(Nw, Nw, 3, 3))
    IF (k_data%bD_bar) ALLOCATE (k_data%mD_bar(Nw, Nw, 3))

    IF (k_data%bA_k_W) ALLOCATE (k_data%mA_k_W(Nw, Nw, 3))
    IF (k_data%bA_bar) ALLOCATE (k_data%mA_bar(Nw, Nw, 3))
    IF (k_data%bdA_k_W) ALLOCATE (k_data%mdA_k_W(Nw, Nw, 3, 3))
    IF (k_data%bdA_bar) ALLOCATE (k_data%mdA_bar(Nw, Nw, 3, 3))
    IF (k_data%bd2A_k_W) ALLOCATE (k_data%md2A_k_W(Nw, Nw, 3, 3, 3))
    IF (k_data%bd2A_bar) ALLOCATE (k_data%md2A_bar(Nw, Nw, 3, 3, 3))
    IF (k_data%bO_k_W) ALLOCATE (k_data%mO_k_W(Nw, Nw, 3))
    IF (k_data%bO_bar) ALLOCATE (k_data%mO_bar(Nw, Nw, 3))
    IF (k_data%bdO_k_W) ALLOCATE (k_data%mdO_k_W(Nw, Nw, 3, 3))
    IF (k_data%bdO_bar) ALLOCATE (k_data%mdO_bar(Nw, Nw, 3, 3))
  END SUBROUTINE allocate_k
  !
  SUBROUTINE clear_k()
    IF (ALLOCATED(t_kpt%H_k)) DEALLOCATE (t_kpt%H_k)
    ! Used in write_f
    ! IF (ALLOCATED(t_kpt%eigval)) DEALLOCATE (t_kpt%eigval)
    IF (ALLOCATED(t_kpt%eigvec)) DEALLOCATE (t_kpt%eigvec)
    IF (ALLOCATED(k_data%mdH_k_W)) DEALLOCATE (k_data%mdH_k_W)
    IF (ALLOCATED(k_data%mdH_bar)) DEALLOCATE (k_data%mdH_bar)
    IF (ALLOCATED(k_data%md2H_k_W)) DEALLOCATE (k_data%md2H_k_W)
    IF (ALLOCATED(k_data%md2H_bar)) DEALLOCATE (k_data%md2H_bar)
    IF (ALLOCATED(k_data%mD_bar)) DEALLOCATE (k_data%mD_bar)

    IF (ALLOCATED(k_data%mA_k_W)) DEALLOCATE (k_data%mA_k_W)
    IF (ALLOCATED(k_data%mA_bar)) DEALLOCATE (k_data%mA_bar)
    IF (ALLOCATED(k_data%mdA_k_W)) DEALLOCATE (k_data%mdA_k_W)
    IF (ALLOCATED(k_data%mdA_bar)) DEALLOCATE (k_data%mdA_bar)
    IF (ALLOCATED(k_data%md2A_k_W)) DEALLOCATE (k_data%md2A_k_W)
    IF (ALLOCATED(k_data%md2A_bar)) DEALLOCATE (k_data%md2A_bar)
    IF (ALLOCATED(k_data%mO_k_W)) DEALLOCATE (k_data%mO_k_W)
    IF (ALLOCATED(k_data%mO_bar)) DEALLOCATE (k_data%mO_bar)
    IF (ALLOCATED(k_data%mdO_k_W)) DEALLOCATE (k_data%mdO_k_W)
    IF (ALLOCATED(k_data%mdO_bar)) DEALLOCATE (k_data%mdO_bar)
  END SUBROUTINE clear_k
  !
  SUBROUTINE make_k()
    USE fft_base, ONLY: fft_R2k, fft_R2k_vec
    USE itg_R, ONLY: R_data
    USE lin_eig_H, ONLY: eig_H
    USE kpoints, ONLY: t_iks, t_kpt
    INTEGER::iw
    CALL start_clock('make_k')

    !... Important: nonallocatable dummy is not present
    !...            if the actual argument is unallocated allocatable.
    !...            (F2008 12.5.2.12 / F2018 15.5.2.12 / F2023 15.5.2.13)
    !...            (No rule before F2008, errore would occur.)

    CALL fft_R2k(R_vec, R_data%mH_R, &
                 X_k=t_kpt%H_k, &
                 dX_k=k_data%mdH_k_W, &
                 d2X_k=k_data%md2H_k_W)
    DO iw = 1, Nw
      t_kpt%H_k(iw, iw) = REAL(t_kpt%H_k(iw, iw), DP)
    END DO
    CALL eig_H(Nw, t_kpt%H_k, t_kpt%eigval(:, t_iks), t_kpt%eigvec(:, :))

    IF (k_data%bdH_bar) CALL t_kpt%rotate(k_data%mdH_k_W, k_data%mdH_bar)
    IF (k_data%bd2H_bar) CALL t_kpt%rotate(k_data%md2H_k_W, k_data%md2H_bar)
    IF (k_data%bD_bar) CALL compute_D_bar(k_data%mdH_bar, t_kpt%eigval(:, t_iks), k_data%mD_bar)

    IF (R_data%bA_R) THEN
      CALL fft_R2k_vec(R_vec, R_data%mA_R, &
                       X_k=k_data%mA_k_W, &
                       dX_k=k_data%mdA_k_W, &
                       d2X_k=k_data%md2A_k_W, &
                       curl_X_k=k_data%mO_k_W, &
                       curl_dX_k=k_data%mdO_k_W)
      IF (k_data%bA_bar) CALL t_kpt%rotate(k_data%mA_k_W, k_data%mA_bar)
      IF (k_data%bdA_bar) CALL t_kpt%rotate(k_data%mdA_k_W, k_data%mdA_bar)
      IF (k_data%bd2A_bar) CALL t_kpt%rotate(k_data%md2A_k_W, k_data%md2A_bar)
      IF (k_data%bO_bar) CALL t_kpt%rotate(k_data%mO_k_W, k_data%mO_bar)
      IF (k_data%bdO_bar) CALL t_kpt%rotate(k_data%mdO_k_W, k_data%mdO_bar)
    END IF
    CALL stop_clock('make_k')
  END SUBROUTINE make_k
  !
  SUBROUTINE rotate_W2H(mat)
    USE lin_mat, ONLY: mat_mul
    COMPLEX(DP), INTENT(INOUT)::mat(Nw, Nw)
    COMPLEX(DP)::tmp(Nw, Nw)
    ! U^+ mat U
    ! CALL mat_mul(t_kpt%eigvec, 'C', mat, 'N', tmp)
    ! CALL mat_mul(tmp, 'N', t_kpt%eigvec, 'N', mat)
    !< 동치 확인 필요.
    CALL errore(1, 'rotate_W2H', 'Not implemented yet')
  END SUBROUTINE rotate_W2H
  SUBROUTINE rotate_H2W(mat)
    USE lin_mat, ONLY: mat_mul
    COMPLEX(DP), INTENT(INOUT)::mat(Nw, Nw)
    COMPLEX(DP)::tmp(Nw, Nw)
    ! U mat U^+
    ! CALL mat_mul(t_kpt%eigvec, 'N', mat, 'N', tmp)
    ! CALL mat_mul(tmp, 'N', t_kpt%eigvec, 'C', mat)
    CALL mat_mul(t_kpt%eigvec, 'N', mat, 'C', tmp)
    CALL mat_mul(t_kpt%eigvec, 'N', tmp, 'C', mat)
  END SUBROUTINE rotate_H2W
END MODULE itg_k
