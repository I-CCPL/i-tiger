SUBROUTINE debug_k()
  USE io_global, ONLY: get_free_unit
  USE io_input, ONLY: debug => debug_k
  USE itg_k
  USE itg_f
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  INTEGER::io_unit
  CHARACTER(len=30)::fname
  INTEGER::ikpt, m, n, a, b
  IF (.NOT. debug) RETURN
  io_unit = get_free_unit()
  ikpt = t_kpt%global_k(t_iks)

  WRITE (fname, '(A,I0,A)') "debug_k.H_k.", ikpt, ".dat"
  OPEN (io_unit, file=fname)
  DO m = 1, Nw
    DO n = 1, Nw
      WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, t_kpt%H_k(m, n)
    END DO
  END DO
  CLOSE (io_unit)

  WRITE (fname, '(A,I0,A)') "debug_k.eigval.", ikpt, ".dat"
  OPEN (io_unit, file=fname)
  DO m = 1, Nw
    WRITE (io_unit, '(I4, ES23.14E3)') m, t_kpt%eigval(m, t_iks)
  END DO
  CLOSE (io_unit)

  WRITE (fname, '(A,I0,A)') "debug_k.eigvec.", ikpt, ".dat"
  OPEN (io_unit, file=fname)
  DO m = 1, Nw
    DO n = 1, Nw
      WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, t_kpt%eigvec(m, n)
    END DO
  END DO
  CLOSE (io_unit)

  IF (k_data%bA_k_W) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.A_k_W.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mA_k_W(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bA_bar) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.A_bar.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mA_bar(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bdH_k_W) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dH_k_W.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mdH_k_W(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bdH_bar) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dH_bar.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mdH_bar(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (ALLOCATED(v_k_H)) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.v_k_H.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, v_k_H(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bdA_k_W) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dA_k_W.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          DO b = 1, 3
            WRITE (io_unit, '(4I4, 2ES23.14E3)') m, n, a, b, k_data%mdA_k_W(m, n, a, b)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bdA_bar) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.dA_bar.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          DO b = 1, 3
            WRITE (io_unit, '(4I4, 2ES23.14E3)') m, n, a, b, k_data%mdA_bar(m, n, b, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bd2H_k_W) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.d2H_k_W.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          DO b = 1, 3
            WRITE (io_unit, '(4I4, 2ES23.14E3)') m, n, a, b, k_data%md2H_k_W(m, n, b, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bd2H_bar) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.d2H_bar.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          DO b = 1, 3
            WRITE (io_unit, '(4I4, 2ES23.14E3)') m, n, a, b, k_data%md2H_bar(m, n, b, a)
          END DO
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bO_k_W) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.O_k_W.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mO_k_W(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF

  IF (k_data%bO_bar) THEN
    WRITE (fname, '(A,I0,A)') "debug_k.O_bar.", ikpt, ".dat"
    OPEN (io_unit, file=fname)
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, '(2I4, 2ES23.14E3)') m, n, k_data%mO_bar(m, n, a)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF
CONTAINS
END SUBROUTINE debug_k
