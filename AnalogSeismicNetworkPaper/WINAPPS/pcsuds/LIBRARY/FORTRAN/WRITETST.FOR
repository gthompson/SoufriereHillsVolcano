!---  Test Write functions
!---  09-Dec-1993, RB

!----------------------------------------------------------------------

      INCLUDE 'suds.fi'   ! FORTRAN interface to SUDS library

!----------------------------------------------------------------------

      INTEGER*2 FUNCTION TEST_WRITE( fd )

      IMPLICIT NONE

      ! Params
      INTEGER*2 fd(*)

      ! SUDS structure definitions
      INCLUDE 'suds.inc'

      ! SUDS library function return types
      INTEGER*2 suds_write, suds_write_data, suds_flush
      REAL*8 get_mstime

      ! SUDS structure used with library functions
      RECORD /suds_structs/ suds

      ! Locals
      INTEGER*2 i, j
      INTEGER*4 l, bytes, samps, error
      INTEGER*2 data [ALLOCATABLE, HUGE] (:)
      REAL*4 rate, freq, secs, PI, TWOPI, ist

      !--- Executable ---

      WRITE( *, '(A)' ) ' Testing write functions...'

      test_write = FALSE

      PI = 4.0*ATAN(1.0)
      TWOPI = 2.0*PI
      rate = 200.0
      freq = 1.0
      secs = 10

      samps = INT4( secs*rate )
      bytes = samps*2

      ist = GET_MSTIME( )

      ! Allocate the data buffer
      ALLOCATE( data(samps), STAT=error )
      IF( error .NE. 0 ) THEN
         WRITE( *, '(A)' ) ' ERROR: Not enough memory!'
         RETURN
      ENDIF

      ! Build a sine wave
      DO l=1, samps
         data(l) = INT(32767.0*SIN((l-1)*((TWOPI/rate)*freq)))
      ENDDO

      ! Write to all three files
      DO i=1, 3
         WRITE( *, '(/AI1\)' ) ' File #', i
         ! 3 channels per file
         DO j=1, 3
            WRITE( *, '(AI1\)' ) ', Chan #', j

            ! Each channel will consist of a STATIONCOMP, INSTRUMENT and
            ! DESCRIPTRACE structure with a waveform following DESCRIPTRACE
            suds.type = STATIONCOMP
            CALL SUDS_INIT( LOC(suds) )

            suds.sc.sc_name.st_name = 'STN'//CHAR(j+48)//CHAR(0)
            suds.sc.sc_name.component = 'v'
            suds.sc.channel = j
            suds.sc.st_lat = 37.08
            suds.sc.st_long = -111.67
            suds.sc.elev = 1286.25
            suds.sc.sensor_type = 'v'
            suds.sc.data_units = 'd'
            suds.sc.data_type = 'i'
            suds.sc.atod_gain = 32

            IF( SUDS_WRITE( fd(i), LOC(suds) ) .EQ. FALSE )
     +         RETURN

            suds.type = INSTRUMENT
            CALL SUDS_INIT( LOC(suds) )

            suds.in.in_name.st_name = 'STN'//CHAR(j+48)//CHAR(0)
            suds.in.in_name.component = 'v'
            suds.in.channel = j
            suds.in.in_serial = 1234
            suds.in.comps = 3
            suds.in.sens_type = 'v'
            suds.in.datatype = 'i'

            IF( SUDS_WRITE( fd(i), LOC(suds) ) .EQ. FALSE )
     +         RETURN

            suds.type = DESCRIPTRACE
            CALL SUDS_INIT( LOC(suds) )

            suds.dt.dt_name.st_name = 'STN'//CHAR(j+48)//CHAR(0)
            suds.dt.dt_name.component = 'v'
            suds.dt.begintime = ist
            suds.dt.datatype = 'i'
            suds.dt.rate = rate
            suds.dt.length = samps
            suds.dt.mindata = -32767.0
            suds.dt.maxdata = 32767.0
            suds.dt.avenoise = 0.0

            suds.data_len = bytes   ! Length of data to follow
            IF( SUDS_WRITE( fd(i), LOC(suds) ) .EQ. FALSE )
     +         RETURN

            IF( SUDS_WRITE_DATA( fd(i), LOC(data), bytes ) .EQ. FALSE )
     +         RETURN

            IF( SUDS_FLUSH( fd(i) ) .EQ. FALSE )
     +         RETURN

         ENDDO
      ENDDO

      ! Release memory
      DEALLOCATE( data )

      ! Success
      test_write = TRUE

      RETURN

      END
