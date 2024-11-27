!---  Test Read functions
!---  09-Dec-1993, RB

!----------------------------------------------------------------------

      INCLUDE 'suds.fi'   ! FORTRAN interface to SUDS library

!----------------------------------------------------------------------

      INTEGER*2 FUNCTION TEST_READ( fd )

      IMPLICIT NONE

      ! Params
      INTEGER*2 fd(*)

      ! SUDS structure definitions
      INCLUDE 'suds.inc'

      ! SUDS library function return types
      INTEGER*2 suds_read, suds_read_data, suds_rewind, suds_seek
      INTEGER*4 suds_pos

      ! SUDS structure used with library functions
      RECORD /suds_structs/ suds

      ! Locals
      INTEGER*2 i
      INTEGER*4 pos, l, error
      INTEGER*2 data [ALLOCATABLE, HUGE] (:)
      CHARACTER buf*256   ! Working string buffer

      !--- Executable ---

      WRITE( *, '(A)' ) ' Test read functions...'

      test_read = FALSE

      ! Read from all three files, structs only, no data
      DO i=1, 3

         ! Rewind the file
         IF( SUDS_REWIND( fd(i) ) .NE. TRUE )
     +      RETURN

         ! Read and display all structures
         DO WHILE( SUDS_READ( fd(i), LOC(suds) ) .NE. SUDS_EOF )
            ! Save position of 3rd channel for seek test
            IF( suds.type .EQ. STATIONCOMP .AND.
     +          suds.sc.channel .EQ. 3 ) THEN
               pos = SUDS_POS( fd(i) )
               IF( pos .EQ. SUDS_EOF )
     +            RETURN
            ENDIF
            CALL DISPLAY_STRUCT( suds )
         ENDDO
      ENDDO

      WRITE( *, '(/A)' ) ' Press return for SEEK test'
      READ( *, '(A)' ) buf

      IF( SUDS_SEEK( fd(3), pos ) .EQ. SUDS_EOF )
     +   RETURN

      DO WHILE( SUDS_READ( fd(3), LOC(suds) ) .NE. SUDS_EOF )
         CALL DISPLAY_STRUCT( suds )

         IF( suds.type .EQ. DESCRIPTRACE ) THEN
            ALLOCATE( data(suds.data_len), STAT=error )
            IF( error .NE. 0 ) THEN
               WRITE( *, '(A)' ) ' ERROR: Not enough memory!'
               RETURN
            ENDIF

            WRITE( *, '(A)' ) ' Waveform data:'

            IF( SUDS_READ_DATA( fd(3), LOC(data), suds.data_len )
     +         .EQ. FALSE ) RETURN

            DO l=1, suds.dt.length
               WRITE( *, '(A,I7\)' ) ' ', data(l)
               IF( MOD( i, 10 ) .EQ. 0 )
     +            WRITE( *, '(A)' )
            ENDDO

            DEALLOCATE( data )

         ENDIF

      ENDDO
      ! Success
      test_read = TRUE

      RETURN

      END

!----------------------------------------------------------------------

      SUBROUTINE DISPLAY_STRUCT( suds )

      IMPLICIT NONE

      ! SUDS structure definitions
      INCLUDE 'suds.inc'

      ! Params
      RECORD /suds_structs/ suds

      ! Return types
      CHARACTER*128 list_mstime

      SELECT CASE( suds.type )
         CASE( STATIONCOMP )
            WRITE(*,'(/A)')
     +         ' Stucture type:        SUDS_STATIONCOMP'
            WRITE(*,'(A,I7)')
     +         ' Data length (bytes):  ', suds.data_len
            WRITE(*,'(A,A)')
     +         ' Station:              ', suds.sc.sc_name.st_name
            WRITE(*,'(A,A)')
     +         ' Component:            ', suds.sc.sc_name.component
            WRITE(*,'(A,I7)')
     +         ' Incidence:            ', suds.sc.incid
            WRITE(*,'(A,I7)')
     +         ' Azimuth:              ', suds.sc.azim
            WRITE(*,'(A,F10.2)')
     +         ' Latitude:             ', suds.sc.st_lat
            WRITE(*,'(A,F10.2)')
     +         ' Longitude:            ', suds.sc.st_long
            WRITE(*,'(A,F10.2)')
     +         ' Elevation:            ', suds.sc.elev
            WRITE(*,'(A,A)')
     +         ' Sensor type:          ', suds.sc.sensor_type
            WRITE(*,'(A,A)')
     +         ' Data type:            ', suds.sc.data_type
            WRITE(*,'(A,A)')
     +         ' Data unit:            ', suds.sc.data_units
            WRITE(*,'(A,I7)')
     +         ' Gain:                 ', suds.sc.atod_gain

         CASE( INSTRUMENT )
            WRITE(*,'(/A)')
     +         ' Stucture type:        SUDS_INSTRUMENT'
            WRITE(*,'(A,I7)')
     +         ' Data length (bytes):  ', suds.data_len
            WRITE(*,'(A,A)')
     +         ' Station:              ', suds.in.in_name.st_name
            WRITE(*,'(A,A)')
     +         ' Component:            ', suds.in.in_name.component
            WRITE(*,'(A,I7)')
     +         ' Instrument serial #:  ', suds.in.in_serial
            WRITE(*,'(A,I7)')
     +         ' Channel:              ', suds.in.channel
            WRITE(*,'(A,I7)')
     +         ' # of components:      ', suds.in.comps
            WRITE(*,'(A,A)')
     +         ' Sensor type:          ', suds.in.sens_type
            WRITE(*,'(A,A)')
     +         ' Data type:            ', suds.in.datatype

         CASE( DESCRIPTRACE )
            WRITE(*,'(/A)')
     +         ' Stucture type:        SUDS_DESCRIPTRACE'
            WRITE(*,'(A,I7)')
     +         ' Data length (bytes):  ', suds.data_len
            WRITE(*,'(A,A)')
     +         ' Station:              ', suds.dt.dt_name.st_name
            WRITE(*,'(A,A)')
     +         ' Component:            ', suds.dt.dt_name.component
            WRITE(*,'(A,A20)')
     +         ' Initial sample time:  ',
     +         LIST_MSTIME( suds.dt.begintime, 10 )
            WRITE(*,'(A,I7)')
     +         ' Local time = GMT +    ', suds.dt.localtime
            WRITE(*,'(A,A)')
     +         ' Data type:            ', suds.dt.datatype
            WRITE(*,'(A,F10.2)')
     +         ' Sample rate:          ', suds.dt.rate
            WRITE(*,'(A,I7)')
     +         ' # of samples          ', suds.dt.length
            WRITE(*,'(A,F10.2)')
     +         ' Max data value:       ', suds.dt.maxdata
            WRITE(*,'(A,F10.2)')
     +         ' Min data value:       ', suds.dt.mindata
            WRITE(*,'(A,F10.2)')
     +         ' Zero level:           ', suds.dt.avenoise
            WRITE(*,'(A,F10.2)')
     +         ' IST correction:       ', suds.dt.time_correct
            WRITE(*,'(A,F10.2)')
     +         ' Rate correction:      ', suds.dt.rate_correct
         CASE DEFAULT
            WRITE(*,'(/A,I2)') ' Struct type: ', suds.type
      END SELECT

      RETURN

      END
