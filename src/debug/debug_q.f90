SUBROUTINE debug_q()
  USE kinds, ONLY: DP
  USE io_global, ONLY: ionode, get_free_unit
  USE io_input, ONLY: debug => debug_q
  USE wannier90, ONLY: w90data
  USE itg_q
  USE system, ONLY: Nw, red2cart_recip
  IMPLICIT NONE
  INTEGER::io_unit
  INTEGER::ikpt, m, n, a, inb, jnb
  REAL(DP)::vec(3)
  IF (.NOT. debug) RETURN
  IF (.NOT. ionode) RETURN
  io_unit = get_free_unit()

  OPEN (unit=io_unit, file='debug_q.wcc.dat')
  DO m = 1, Nw
    WRITE (io_unit, *) w90data%wannier_center_cart(:, m)
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.win.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    WRITE (io_unit, '(2I6)') w90data%win_min(ikpt), w90data%ndimwin(ikpt)
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.eigval.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO n = 1, Nw
      WRITE (io_unit, *) w90data%eigval(n, ikpt)
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.v_mat.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO m = 1, w90data%nbnd
      DO n = 1, Nw
        WRITE (io_unit, *) w90data%v_matrix(m, n, ikpt)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.overlap.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO inb = 1, w90data%nnb
      DO m = 1, w90data%nbnd
        DO n = 1, w90data%nbnd
          WRITE (io_unit, *) w90data%overlap(m, n, inb, ikpt)
        END DO
      END DO
    END DO
  END DO

  OPEN (unit=io_unit, file='debug_q.wb.dat')
  DO inb = 1, w90data%nnb
    WRITE (io_unit, *) w90data%wb(inb)
  END DO

  OPEN (unit=io_unit, file='debug_q.bvec_red.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO inb = 1, w90data%nnb
      jnb = w90data%bvec_index(inb, ikpt)
      DO a = 1, w90data%nnb
        WRITE (io_unit, '(3F12.6)') w90data%bvec_red(:, jnb)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.bvec_cart.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO inb = 1, w90data%nnb
      jnb = w90data%bvec_index(inb, ikpt)
      CALL red2cart_recip(w90data%bvec_red(:, jnb), vec)
      DO a = 1, 3
        WRITE (io_unit, *) vec(a)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.Hq.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO m = 1, Nw
      DO n = 1, Nw
        WRITE (io_unit, *) w90data%Hq(m, n, ikpt)
      END DO
    END DO
  END DO
  CLOSE (io_unit)

  OPEN (unit=io_unit, file='debug_q.Aq.dat')
  DO ikpt = 1, w90data%kpts%nkpt
    DO m = 1, Nw
      DO n = 1, Nw
        DO a = 1, 3
          WRITE (io_unit, *) w90data%Aq(a, m, n, ikpt)
        END DO
      END DO
    END DO
  END DO
  CLOSE (io_unit)
END SUBROUTINE debug_q
