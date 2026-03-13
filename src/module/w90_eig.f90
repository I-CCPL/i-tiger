SUBMODULE(wannier90) w90_eig
CONTAINS
  MODULE SUBROUTINE read_w90_eig(self)
    !< Ref. wannier90/src/readwrite.F90
    USE io_global, ONLY: check_file
    USE mp_base, ONLY: mp_bcast
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::io_unit, ikpt, ibnd, jkpt, jbnd, ios
    CHARACTER(LEN=256)::msg
    !
    WRITE (stdout, '(2X, A)') 'Reading .eig file...'
    CALL check_file(TRIM(self%prefix)//'.eig')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(self%prefix)//'.eig', form='formatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_eig', 'Failed to open '//TRIM(self%prefix)//'.eig')
      ALLOCATE (self%eigval(self%nbnd, self%kpts%nkpt))

      DO ikpt = 1, self%kpts%nkpt
        DO ibnd = 1, self%nbnd
          READ (io_unit, *, iostat=ios) jbnd, jkpt, self%eigval(ibnd, ikpt)
          IF (ios /= 0) THEN
            WRITE (msg, '(A,1X,I0,1X,A,1X,I0)') 'Failed reading .eig at ibnd=', ibnd, 'ikpt=', ikpt
            CALL errore(ios, 'read_w90_eig', TRIM(msg))
          END IF
          IF (jbnd /= ibnd .OR. jkpt /= ikpt) THEN
            WRITE (msg, '(A,1X,I0,A,I0,A,I0,A,I0,A)') &
              'Unexpected (.eig) indices. Expected (', ibnd, ',', ikpt, '), got (', jbnd, ',', jkpt, ').'
            CALL errore(1, 'read_w90_eig', TRIM(msg))
          END IF
        END DO
      END DO
      CLOSE (io_unit)
    ELSE
      ALLOCATE (self%eigval(self%nbnd, self%kpts%nkpt))
    END IF
    CALL mp_bcast(self%eigval)
  END SUBROUTINE read_w90_eig
END SUBMODULE
