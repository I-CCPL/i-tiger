SUBROUTINE debug_R()
  USE io_global, ONLY: get_free_unit
  USE io_input, ONLY: debug => debug_R
  USE itg_R
  IMPLICIT NONE
  INTEGER::io_unit
  INTEGER::irpt, iw, jw, a, b
  INTEGER::ir0pt, ideg, iuw
  IF (.NOT. debug) RETURN
  io_unit = get_free_unit()

  OPEN (io_unit, file="debug_R.R_red.dat")
  DO irpt = 1, R_vec%nRpt
    WRITE (io_unit, *) INT(R_vec%R_red(:, irpt))
  END DO
  CLOSE (io_unit)

  OPEN (io_unit, file="debug_R.R_cart.dat")
  DO irpt = 1, R_vec%nRpt
    WRITE (io_unit, *) R_vec%R_cart(:, irpt)
  END DO
  CLOSE (io_unit)

  OPEN (io_unit, file="debug_R.R0.dat")
  DO irpt = 1, R_vec%nRpt
    DO iw = 1, Nw
      DO jw = 1, Nw
        iuw = R_vec%shift_map_inv(iw, jw)
        ir0pt = R_vec%map_r2r0(1, iuw, irpt)
        ideg = R_vec%map_r2r0(2, iuw, irpt)
        WRITE (io_unit, '(5I3)') irpt, iw, jw, ir0pt, ideg
      END DO
    END DO
  END DO

  OPEN (io_unit, file="debug_R.shift_cart.dat")
  DO irpt = 1, R_vec%nRpt
    DO iw = 1, Nw
      DO jw = 1, Nw
        iuw = R_vec%shift_map_inv(iw, jw)
        WRITE (io_unit, *) R_vec%shift_cart_u(:, iuw)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (io_unit, file="debug_R.H_R.dat")
  DO irpt = 1, R_vec%nRpt
    DO iw = 1, Nw
      DO jw = 1, Nw
        WRITE (io_unit, *) H_R(iw, jw, irpt)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  IF (ALLOCATED(A_R)) THEN
    OPEN (io_unit, file="debug_R.A_R.dat")
    DO irpt = 1, R_vec%nRpt
      DO iw = 1, Nw
        DO jw = 1, Nw
          WRITE (io_unit, *) A_R(iw, jw, irpt, :)
        END DO
      END DO
    END DO
    CLOSE (io_unit)
  END IF
END SUBROUTINE debug_R
