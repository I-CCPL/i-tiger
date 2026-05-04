SUBROUTINE debug_k()
  USE io_global, ONLY: get_free_unit
  USE io_input, ONLY: debug => debug_k
  USE mp_global, ONLY: mp_rank
  USE itg_k
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  INTEGER::io_unit
  CHARACTER(len=30)::fname
  INTEGER::ikpt, m, n, a, b
  IF (.NOT. debug) RETURN
  io_unit = get_free_unit()
  WRITE (fname, '(A,I0,A)') "debug_k.H_k.", mp_rank, ".dat"
  OPEN (io_unit, file=fname)
  DO ikpt = 1, t_kpt%nkpt
    DO m = 1, Nw
      DO n = 1, Nw
        WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, t_kpt%H_k(m, n)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  WRITE (fname, '(A,I0,A)') "debug_k.eigval.", mp_rank, ".dat"
  OPEN (io_unit, file=fname)
  DO ikpt = 1, t_kpt%nkpt
    DO m = 1, Nw
      WRITE (io_unit, '(2I4, ES23.14E3)') ikpt, m, t_kpt%eigval(m, ikpt)
    END DO
  END DO
  CLOSE (io_unit)

  WRITE (fname, '(A,I0,A)') "debug_k.eigvec.", mp_rank, ".dat"
  OPEN (io_unit, file=fname)
  DO ikpt = 1, t_kpt%nkpt
    DO m = 1, Nw
      DO n = 1, Nw
        WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, t_kpt%eigvec(m, n)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  IF (ALLOCATED(A_k_W)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.A_k_W.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, A_k_W(m, n, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(A_bar)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.A_bar.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, A_bar(m, n, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(dH_k_W)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dH_k_W.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, dH_k_W(m, n, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(dH_bar)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dH_bar.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, dH_bar(m, n, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(v_k_H)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.v_k_H.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            WRITE (io_unit, '(3I4, 2ES23.14E3)') ikpt, m, n, v_k_H(m, n, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(dA_k_W)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dA_k_W.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            DO b = 1, 3
              WRITE (io_unit, '(5I4, 2ES23.14E3)') ikpt, m, n, a, b, dA_k_W(m, n, a, b)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(dA_bar)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dA_bar.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            DO b = 1, 3
              WRITE (io_unit, '(5I4, 2ES23.14E3)') ikpt, m, n, a, b, dA_bar(m, n, b, a)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(d2H_k_W)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.d2H_k_W.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            DO b = 1, 3
              WRITE (io_unit, '(5I4, 2ES23.14E3)') ikpt, m, n, a, b, d2H_k_W(m, n, b, a)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(d2H_bar)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.d2H_bar.", mp_rank, ".dat"
    OPEN (io_unit, file=fname)
    DO ikpt = 1, t_kpt%nkpt
      DO m = 1, Nw
        DO n = 1, Nw
          DO a = 1, 3
            DO b = 1, 3
              WRITE (io_unit, '(5I4, 2ES23.14E3)') ikpt, m, n, a, b, d2H_bar(m, n, b, a)
            END DO
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

END SUBROUTINE debug_k
