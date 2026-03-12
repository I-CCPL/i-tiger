MODULE env
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  CHARACTER(LEN=10)::itg_version = 'v0.0.0.19'
CONTAINS
  SUBROUTINE env_start(date, time)
    USE mp_global, ONLY: mp_start, mp_rank, mp_root, mp_size, mp_barrier
    USE io_global, ONLY: ionode, stdout, write_bold_line
    CHARACTER(LEN=*), INTENT(IN)::date, time
    CHARACTER(len=10)::cdate, ctime
    !
    CALL mp_start()
    ionode = (mp_rank == mp_root)
    IF (.NOT. ionode) THEN
      OPEN (unit=stdout, file='/dev/null', status='unknown')
    END IF
    !
    CALL current_date_time(cdate, ctime)
    WRITE (stdout, '(A)') 'Starting i-TIGER on '//TRIM(cdate)//' '//TRIM(ctime)
    WRITE (stdout, '(2X, 4A)') "= Compiled on ", date, ' ', time
    WRITE (stdout, '(2X,A,I0,A)') '= Running on ', mp_size, ' processors.'
    WRITE (stdout, *)
    CALL write_bold_line()
    WRITE (stdout, '(2X,A)') 'Incheon Tight-binding Induced Generalized Electronic Response'
    WRITE (stdout, '(2X,A,A)') 'i-TIGER ', itg_version
    CALL write_bold_line()
    WRITE (stdout, *)
    CALL mp_barrier()
  END SUBROUTINE env_start
  !
  SUBROUTINE env_end()
    USE mp_global, ONLY: mp_end
    CHARACTER(len=10)::cdate, ctime
    CALL mp_end()
    CALL current_date_time(cdate, ctime)
    WRITE (stdout, *)
    WRITE (stdout, '(A)') 'Finished i-TIGER on '//TRIM(cdate)//' '//TRIM(ctime)
    FLUSH (stdout)
  END SUBROUTINE env_end
  !
  SUBROUTINE current_date_time(cdate, ctime)
    CHARACTER(len=10), INTENT(OUT)::cdate
    CHARACTER(len=10), INTENT(OUT)::ctime
    CHARACTER(len=1), PARAMETER::delim1 = '-', delim2 = ':'
    INTEGER::date_time(8)
    CALL DATE_AND_TIME(values=date_time)
    WRITE (cdate, '(I4.4,A1,I2.2,A1,I2.2)') date_time(1), delim1, date_time(2), delim1, date_time(3)
    WRITE (ctime, '(I2.2,A1,I2.2,A1,I2.2)') date_time(5), delim2, date_time(6), delim2, date_time(7)
  END SUBROUTINE current_date_time
END MODULE
