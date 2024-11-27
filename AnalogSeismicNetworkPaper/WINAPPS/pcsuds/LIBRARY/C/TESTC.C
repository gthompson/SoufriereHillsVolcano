// Test the SUDS library

// 09-Dec-1993, RLB

// Includes -----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "suds.h"

// Constants ----------------------------------------------------------
#define N_FILES     3            // Number of test files
#define N_WAVEFORMS 3            // Number of traces per file

#define SECONDS     20.0         // Length of traces (seconds)
#define FREQ        1.0          // Frequency (Hz)
#define RATE        200.0        // Sampling rate (Hz)
#define PEAK        32767        // Peak amplitude

#define TWO_PI      6.283185     // 2 * Pi
#define T_VOID      -2147472000  // Undefined or "Void" time

// Prototypes ---------------------------------------------------------
SH_INT write_files( VOID );
SH_INT copy_file( STRING *source, STRING *dest );
SH_INT read_files( VOID );
VOID display_struct( SUDS suds );

// Globals ------------------------------------------------------------
MS_TIME session_time;

//---------------------------------------------------------------------
main( int argc, char *argv[] ) {

   printf( "SUDS Library test, SUDS version %s\n\n", _SUDS_VERSION );

   session_time = get_mstime( );

   // Write out files
   printf( "Writing files...\n" );
   if( ! write_files( ) ) {
      printf( "\n%s\n", suds_get_err( ) );
      exit( 1 );
   }

   // Copy a file
   printf( "Copy test..." );
   if( ! copy_file( "test1.sud", "copy.sud" ) ) {
      printf( "\n%s\n", suds_get_err( ) );
      exit( 1 );
   }

   // Read data back in
   printf( "\nPress any key to read file(s) back in" );
   getch( );
   printf( "\n" );

   if( ! read_files( ) ) {
      printf( "\n%s\n", suds_get_err( ) );
      exit( 1 );
   }

   exit( 0 );
}

//---------------------------------------------------------------------
SH_INT copy_file( STRING *source, STRING *dest ) {
   SH_INT sfd, dfd, ret;
   LG_INT spos, cpos;
   SUDS suds;
   VOID _huge *ptr;

   if( ( sfd = suds_open( source, SUDS_READONLY ) ) == FALSE )
      return( FALSE );

   if( ( dfd = suds_open( dest, SUDS_CREATE ) ) == FALSE )
      return( FALSE );

   while( ( ret = suds_read( sfd, &suds ) ) != SUDS_EOF ) {
      if( ! ret )
         return( FALSE );

      if( ! suds_write( dfd, &suds ) )
         return( FALSE );

      if( suds.data_len > 0 ) {
         if( ( ptr = halloc( suds.data_len, 1 ) ) == NULL ) {
            printf( "ERROR: not enough memory!\n" );
            return( FALSE );
         }

         if( ! suds_read_data( sfd, ptr, suds.data_len ) )
            return( FALSE );

         if( ! suds_write_data( dfd, ptr, suds.data_len ) )
            return( FALSE );

         // Test suds_update function
         strcpy( suds.st.network, "TST" );
         if( ! suds_update( dfd, &suds ) )
            return( FALSE );
      }
   }

   if( ! suds_flush( dfd ) )
      return( FALSE );

   hfree( ptr );

   return( TRUE );
}

//---------------------------------------------------------------------
// Write out N_FILES SUDS files containing N_WAVEFORMS of sine waves.
// Files are written while open concurrently.

SH_INT write_files( VOID ) {
   SH_INT i, fd[N_FILES], chan;
   SH_INT _huge *ptr;
   LG_INT samps, l, bytes;
   CHAR filespec[_MAX_PATH], comp;
   CHANSETENTRY set[N_WAVEFORMS];
   SUDS suds;

   // Open the file(s)
   for( i=0; i<N_FILES; i++ ) {
      sprintf( filespec, "TEST%d.SUD", i+1 );
      if( ( fd[i] = suds_open( filespec, SUDS_CREATE ) ) == FALSE )
         return( FALSE );
   }

   // Compute storage
   samps = (LG_INT)(RATE*SECONDS);
   bytes = samps*(LG_INT)sizeof(SH_INT);

   // Allocate trace buffer
   if( ( ptr = (SH_INT _huge *)halloc( samps, sizeof(SH_INT) ) ) == NULL ) {
      fprintf( stderr, "\nERROR: Not enough memory!\n" );
      return( FALSE );
   }

   // Put a sine wave in the buffer
   printf( "Building trace...\n" );
   for( l=0; l<samps; l++ )
      *(ptr+l) = (SH_INT)(sin( (FLOAT)l * ((TWO_PI/RATE) * FREQ) ) * PEAK);

   // Output traces into each file
   for( i=0; i<N_FILES; i++ ) {

      printf( "TEST%d.SUD\n", i );

      // Output each channel
      for( chan=1; chan<=N_WAVEFORMS; chan++ ) {

         switch( chan-1 % 3 ) {
            case 0:
               comp = 'v';
               break;
            case 1:
               comp = 'n';
               break;
            case 2:
               comp = 'e';
               break;
         }

         // Stationcomp struct
         suds.type = STATIONCOMP;
         suds_init( &suds );

         sprintf( suds.sc.sc_name.st_name, "STN%1d", chan );
         suds.sc.channel = chan;
         suds.sc.sc_name.component = comp;
         suds.sc.st_lat = 37.08;
         suds.sc.st_long = -111.67;
         suds.sc.elev = 1286.25;
         suds.sc.sensor_type = 'v';
         suds.sc.data_type = 'i';
         suds.sc.data_units = 'd';
         suds.sc.atod_gain = 32;
         suds.sc.effective = (ST_TIME)session_time;

         if( ! suds_write( fd[i], &suds ) )
            return( FALSE );

         // Instrument struct
         suds.type = INSTRUMENT;
         suds_init( &suds );

         sprintf( suds.in.in_name.st_name, "STN%1d", chan );
         suds.in.channel = chan;
         suds.in.in_name.component = comp;
         suds.in.in_serial = 1234;
         suds.in.comps = N_WAVEFORMS;
         suds.in.sens_type = 'v';
         suds.in.datatype = 'i';
         suds.in.effective = (ST_TIME)session_time;

         if( ! suds_write( fd[i], &suds ) )
            return( FALSE );

         // Descriptrace struct
         suds.type = DESCRIPTRACE;
         suds_init( &suds );

         sprintf( suds.dt.dt_name.st_name, "STN%1d", chan );
         suds.dt.dt_name.component = comp;
         suds.dt.begintime = session_time;
         suds.dt.datatype = 'i';
         suds.dt.rate = RATE;
         suds.dt.length = samps;
         suds.dt.mindata = -32767.0;
         suds.dt.maxdata = 32767.0;
         suds.dt.avenoise = 0.0;

         suds.data_len = bytes;
         if( ! suds_write( fd[i], &suds ) )
            return( FALSE );

         // Write out the data
         if( ! suds_write_data( fd[i], ptr, bytes ) )
            return( FALSE );
         if( ! suds_flush( fd[i] ) )
            return( FALSE );

         // Channel set entry
         set[chan-1].inst_num = 1234;
         set[chan-1].stream_num = i+1;
         set[chan-1].chan_num = chan;
         set[chan-1].st.network[0] = '\0';
         sprintf( set[chan-1].st.st_name, "STN%1d", chan );
         set[chan-1].st.component = comp;
         set[chan-1].st.inst_type = 0;
      }

      // ChannelSet struct
      suds.type = CHANSET;
      suds_init( &suds );

      suds.cs.type = 1;
      suds.cs.entries = N_WAVEFORMS;
      strcpy( suds.cs.network, "BW" );
      strcpy( suds.cs.name, "STN" );
      suds.cs.active = (ST_TIME)session_time - 60;

      suds.data_len = sizeof(CHANSETENTRY) * N_WAVEFORMS;
      if( ! suds_write( fd[i], &suds ) )
         return( FALSE );

      // Write out the data
      if( ! suds_write_data( fd[i], set, suds.data_len ) )
         return( FALSE );
      if( ! suds_flush( fd[i] ) )
         return( FALSE );

   }

   // Free memory
   hfree( ptr );

   // Close the file(s)
   for( i=0; i<N_FILES; i++ ) {
      if( ! suds_close( fd[i] ) )
         return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
SH_INT read_files( VOID ) {
   SH_INT i, fd[N_WAVEFORMS], ret;
   int _huge *ptr;
   long samps, l, pos = 1;
   char filespec[_MAX_PATH];
   SUDS suds;
   CHANSETENTRY *set = NULL;

   // Open the file(s)
   for( i=0; i<N_FILES; i++ ) {
      sprintf( filespec, "TEST%d.SUD", i+1 );
      if( ( fd[i] = suds_open( filespec, SUDS_READONLY ) ) == FALSE )
         return( FALSE );
   }

   // Read structures only
   printf( "\nReading stuctures only, skipping all data\n" );
   for( i=0; i<N_FILES; i++ ) {
      printf( "\nFile #%d\n", i+1 );
      while( ( ret = suds_read( fd[i], &suds ) ) != SUDS_EOF ) {
         if( ! ret )
            return( FALSE );
         display_struct( suds );
      }
   }

   // Rewind the first file
   if( ! suds_rewind( fd[0] ) )
      return( FALSE );

   printf( "Press any key for complete dump of first file\n\n" );
   getch( );

   again:

   while( ( ret = suds_read( fd[0], &suds ) ) != SUDS_EOF ) {
      if( ! ret )
         return( FALSE );

      display_struct( suds );

      if( suds.type == STATIONCOMP && pos ) {
         if( suds.sc.channel == 3 ) {
            pos = suds_pos( fd[0] );
            printf( "Position of previous struct saved for seek test...\n\n" );
         }
      }

      if( suds.data_len > 0 ) {
         switch( suds.type ) {
            case DESCRIPTRACE:
               samps = suds.data_len / sizeof(SH_INT);
               // Allocate trace buffer
               if( ( ptr = (SH_INT _huge *)halloc( samps, sizeof(SH_INT) ) ) == NULL ) {
                  fprintf( stderr, "\nERROR: Not enough memory!\n" );
                  return( FALSE );
               }
               if( ( ret = suds_read_data( fd[0], ptr, suds.data_len ) ) == FALSE )
                  return( FALSE );
               if( ret == SUDS_EOF )
                  break;
               printf( "Waveform data:\n" );
               for( l=0; l<samps; l++ ) {
                  printf( " %6d", *(ptr+l) );
                  if( (l+1) % 10 == 0 )
                     printf( "\n" );
               }
               printf( "\n" );
               hfree( ptr );
               break;

            case CHANSET:
               if( (set = (CHANSETENTRY *)malloc( (size_t)suds.data_len )) == NULL ) {
                  fprintf( stderr, "\nERROR: Not enough memory!\n" );
                  return( FALSE );
               }
               if( ( ret = suds_read_data( fd[0], set, suds.data_len ) ) == FALSE )
                  return( FALSE );
               if( ret == SUDS_EOF )
                  break;
               for( i=0; i<suds.cs.entries; i++ ) {
                  printf( "   Channel entry:    %d\n", i+1 );
                  printf( "      Instrument#:      %d\n", (set+i)->inst_num );
                  printf( "      Stream#:          %d\n", (set+i)->stream_num );
                  printf( "      Channel#:         %d\n", (set+i)->chan_num );
                  printf( "      Station ID:       %s\n", (set+i)->st.st_name );
                  printf( "      Component ID:     %c\n", (set+i)->st.component );
               }
               printf( "\n" );
               free( set );
               break;
         }
      }
   }

   if( pos ) {
      printf( "\nPress any key to seek and re-read channel 3\n\n" );
      getch( );
      if( ! suds_seek( fd[0], pos ) )
         return( FALSE );
      pos = 0;
      goto again;
   }

   // Close the file(s)
   for( i=0; i<N_FILES; i++ ) {
      if( ! suds_close( fd[i] ) )
         return( FALSE );
   }

   return( TRUE );
}


//---------------------------------------------------------------------
VOID display_struct( SUDS suds ) {

   switch( suds.type ) {
      case DESCRIPTRACE:
         printf( "   Type:             SUDS_DESCRIPTRACE\n" );
         printf( "   Data length:      %ld bytes\n", suds.data_len );
         printf( "   Station:          %s\n", suds.dt.dt_name.st_name );
         printf( "   Component:        %c\n", suds.dt.dt_name.component );
         printf( "   IST:              %s\n", list_mstime( suds.dt.begintime, 10 ) );
         printf( "   Localtime:        %d\n", suds.dt.localtime );
         printf( "   Data type:        %c\n", suds.dt.datatype );
         printf( "   Sample rate:      %f\n", suds.dt.rate );
         printf( "   # of samples:     %ld\n", suds.dt.length );
         printf( "   Min data value:   %f\n", suds.dt.mindata );
         printf( "   Max data value:   %f\n", suds.dt.maxdata );
         printf( "   Zero level:       %f\n", suds.dt.avenoise );
         printf( "   Time correction:  %lf\n", suds.dt.time_correct );
         printf( "   Rate correction:  %lf\n\n", suds.dt.rate_correct );
         break;
      case INSTRUMENT:
         printf( "   Type:             SUDS_INSTRUMENT\n" );
         printf( "   Data length:      %ld bytes\n", suds.data_len );
         printf( "   Station:          %s\n", suds.in.in_name.st_name );
         printf( "   Component:        %c\n", suds.in.in_name.component );
         printf( "   Instrument #:     %d\n", suds.in.in_serial );
         printf( "   Channel:          %d\n", suds.in.channel );
         printf( "   # of components:  %d\n", suds.in.comps );
         printf( "   Sensor type:      %c\n", suds.in.sens_type );
         printf( "   Data type:        %c\n", suds.in.datatype );
         printf( "   Effective time:   %s\n\n", list_mstime( suds.in.effective, 10 ) );
         break;
      case STATIONCOMP:
         printf( "   Type:             SUDS_STATIONCOMP\n" );
         printf( "   Data length:      %ld bytes\n", suds.data_len );
         printf( "   Station:          %s\n", suds.sc.sc_name.st_name );
         printf( "   Component:        %c\n", suds.sc.sc_name.component );
         printf( "   Incidence:        %d\n", suds.sc.incid );
         printf( "   Azimuth:          %d\n", suds.sc.azim );
         printf( "   Latitude:         %lf\n", suds.sc.st_lat );
         printf( "   Longitude:        %lf\n", suds.sc.st_long );
         printf( "   Elevation:        %f\n", suds.sc.elev );
         printf( "   Sensor type:      %c\n", suds.sc.sensor_type );
         printf( "   Data type:        %c\n", suds.sc.data_type );
         printf( "   Data unit:        %c\n", suds.sc.data_units );
         printf( "   A/D gain:         %dX\n", suds.sc.atod_gain );
         printf( "   Effective time:   %s\n\n", list_mstime( suds.sc.effective, 10 ) );
         break;
      case CHANSET:
         printf( "   Type:             SUDS_CHANSET\n" );
         printf( "   Set type:         %d\n", suds.cs.type );
         printf( "   Set entries:      %d\n", suds.cs.entries );
         printf( "   Network name:     %s\n", suds.cs.network );
         printf( "   Set name:         %s\n", suds.cs.name );
         printf( "   Active at:        %s\n", list_mstime( suds.cs.active, 10 ) );
         if( suds.cs.inactive != T_VOID )
            printf( "   Inactive at:      %s\n", list_mstime( suds.cs.inactive, 10 ) );
         else
            printf( "   Inactive at:      undefined\n" );
         break;
      default:
         printf( "   Structure type    %d\n\n", suds.type );
         break;
   }

   return;
}
