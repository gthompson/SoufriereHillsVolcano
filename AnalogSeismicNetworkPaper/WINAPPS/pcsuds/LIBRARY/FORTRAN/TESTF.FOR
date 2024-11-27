!---  Test FORTRAN calls to SUDS library
!---  30-Nov-1993, RB

!----------------------------------------------------------------------

      INCLUDE 'suds.fi'   ! FORTRAN interface to SUDS library

!----------------------------------------------------------------------

      PROGRAM TEST

      IMPLICIT NONE

      ! SUDS structure definitions
      INCLUDE 'suds.inc'

      ! SUDS library function return types
      INTEGER*2 suds_open, suds_close
      CHARACTER*256 suds_get_err

      ! Locals
      INTEGER*2 i, test_write, test_read
      CHARACTER buf*256   ! Working string buffer
      INTEGER*2 fd(3)     ! File descriptors

      !--- Executable ---

      WRITE( *,'(A,A,F4.2/)' ) ' FORTRAN language SUDS library test',
     +  ', SUDS version ', SUDS_VERSION

      ! Test time functions ---
      WRITE( *, '(A)' ) ' Testing time functions...'
      CALL TIMETEST( )

      ! Test I/O functions ---
      WRITE( *, '(/A)' ) ' Testing I/O functions...'

      ! Open 3 files with read/write access
      DO i=1, 3
         buf = 'TEST'//CHAR(i+48)//'.SUD'//CHAR(0)
         fd(i) = SUDS_OPEN( buf, SUDS_CREATE )
         IF( fd(i) .EQ. FALSE ) THEN
            buf = SUDS_GET_ERR( )
            WRITE( *, '(A,A)' ) ' ', buf(:INDEX( buf, CHAR(0) ))
            STOP' '
         ENDIF
      ENDDO
      WRITE( *, '(A)' ) ' Files open'

      ! Write test
      WRITE( *, '(/A)' ) ' Press return for WRITE test'
      READ( *, '(A)' ) buf
      IF( TEST_WRITE( fd ) .EQ. FALSE ) THEN
         buf = SUDS_GET_ERR( )
         WRITE( *, '(A,A)' ) ' ', buf(:INDEX( buf, CHAR(0) ))
         STOP' '
      ENDIF

      ! Read, rewind and seek test
      WRITE( *, '(//A)' ) ' Press return for READ test'
      READ( *, '(A)' ) buf
      IF( TEST_READ( fd ) .EQ. FALSE ) THEN
         buf = SUDS_GET_ERR( )
         WRITE( *, '(A,A)' ) ' ', buf(:INDEX( buf, CHAR(0) ))
         STOP' '
      ENDIF

      ! Close the files
      DO i=1, 3
         IF( SUDS_CLOSE( fd(i) ) .EQ. FALSE ) THEN
            buf = SUDS_GET_ERR( )
            WRITE( *, '(A,A)' ) ' ', buf(:INDEX( buf, CHAR(0) ))
            STOP' '
         ENDIF
      ENDDO
      WRITE( *, '(/A)' ) ' Files closed'

      STOP' '

      END

!----------------------------------------------------------------------

      SUBROUTINE TIMETEST

      IMPLICIT NONE

      ! Library function return types
      REAL*8 MAKE_MSTIME
      CHARACTER*128 LIST_MSTIME
      INTEGER*2 DECODE_MSTIME, YRDAY, MNDAY, ISLEAP

      ! Locals
      INTEGER*2 yr, mn, dy, hr, mi, doy, dummy
      REAL*8 ms_time, sc
      CHARACTER*128 buf

      !--- Executable ---

      ! Make some time
      ms_time = MAKE_MSTIME( 1992, 8, 31, 13, 05, 12.345 )

      ! Note the the string is NULL (CHAR(0)) terminated
      buf = LIST_MSTIME( ms_time, 4 )
      WRITE( *, '(A,A)' ) ' ', buf(:INDEX(buf,CHAR(0)))
      buf = LIST_MSTIME( ms_time, 10 )
      WRITE( *, '(A,A)' ) ' ', buf(:INDEX(buf,CHAR(0)))

      ! Decode time back into is component parts
      dummy = DECODE_MSTIME( ms_time, yr, mn, dy, hr, mi, sc )

      ! Get the Day-of-year form month and day
      doy = YRDAY( mn, dy, ISLEAP( yr, 0 ) )

      ! Get month, day from doy
      dummy = MNDAY( doy, ISLEAP( yr, 0 ), mn, dy )

      WRITE( *, '(A,I4)' ) ' Year:        ', yr
      WRITE( *, '(A,I3.3)' ) ' Day of year: ', doy
      WRITE( *, '(A,I2.2)' ) ' Month:       ', mn
      WRITE( *, '(A,I2.2)' ) ' Day:         ', dy
      WRITE( *, '(A,I2.2)' ) ' Hour:        ', hr
      WRITE( *, '(A,I2.2)' ) ' Minute:      ', mi
      WRITE( *, '(A,F6.3)' ) ' Second:      ', sc

      RETURN

      END
