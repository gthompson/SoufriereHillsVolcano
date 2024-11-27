// SUDS I/O routines, callable for C, FORTRAN, PASCAL, etc...

// This module is a component of the SUDS library.

// Written in Microsoft C/C++ 7.00, Large memory model by;

//    Robert Banfill,
//    Small Systems Support
//    2 Boston Harbor Place
//    P.O. Box 410205
//    Big Water,  UT  84741-0205
//    (801) 675-5827 Voice
//    (801) 675-3780 FAX

// Revision History:

// Version 1.40 Sep-1992, RLB
//    Complete rebuild, started over.

// Version 1.41 Dec-1992, RLB
//    Fixed access mode problems, minor cleanup

// Version 1.42 14-Feb-1993, RLB
//    Fixed problem when first struct in file has data following it.
//    Added suds_abs_pos( ) function.

// Version 1.43 02-Jun-1993, RLB
//    Fixed seek problems when working with multiple files
//    added suds_update( ) function;

// Includes -----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <share.h>

#define _SUDS_IO_
#include "suds.h"

// Prototypes ---------------------------------------------------------

// Open and Close file functions
int suds_open( char *filespec, int mode );
int suds_close( int fd );

// Read functions
int suds_read( int fd, SUDS *suds );
int suds_read_data( int fd, void _huge *ptr, long len );

// Write functions
int suds_write( int fd, SUDS *suds );
int suds_write_data( int fd, void _huge *ptr, long len );
int suds_update( int fd, SUDS *suds );
int suds_flush( int fd );

// Navigation functions
long suds_pos( int fd );
long suds_abs_pos( int fd );
int suds_seek( int fd, long pos );
int suds_rewind( int fd );

// Retrieve error messages
char *suds_get_err( void );

// Helper functions
int _valid_fd( int fd );
void _s_upper( char *buffer );

// Constants ----------------------------------------------------------
#define SUDS_MAX_FILES 10       // Number of concurrently open SUDS files allowed
#define SUDS_IO_BLK    32767L   // Size of block used for huge I/O operations

// User defined types -------------------------------------------------
typedef struct _FILE_TAG {
   FILE *stream;              // Stream pointer to file
   char filespec[_MAX_PATH];  // Filespec of file
   SUDS_STRUCTTAG tag;        // Struct tag for current struct
   long offset;               // Offset to structure in bytes
   int  access;               // 0=not open, 1=read only, 2=read/write, 3=append
   int  seek_flag;            // Last operation on file was a seek
} FILE_TAG;

// Globals ------------------------------------------------------------
char _suds_err_buf[256];               // Error message buffer
FILE_TAG _suds_files[SUDS_MAX_FILES];  // File tags
int _suds_nfiles = 0;                  // Number of open files

extern int _suds_size_struct[];        // Structure sizes (defined in SUDSFORM.H)

//---------------------------------------------------------------------
int suds_open( char *filespec, int mode ) {
   register i;
   int fd;
   char access[4];

   // Initialize file tags
   if( _suds_nfiles == 0 ) {
      for( i=0; i<SUDS_MAX_FILES; i++ ) {
         _suds_files[i].stream = NULL;
         _suds_files[i].filespec[0] = '\0';
         _suds_files[i].offset = 0L;
         _suds_files[i].access = 0;
         _suds_files[i].seek_flag = FALSE;
      }
   }

   // Make sure we have an available file descriptor
   if( _suds_nfiles >= SUDS_MAX_FILES ) {
      sprintf( _suds_err_buf, "ERROR: To many open files!" );
      return( FALSE );
   }

   // Find an available file descriptor
   for( i=1; i<=SUDS_MAX_FILES; i++ ) {
      if( _suds_files[i-1].stream == NULL ) {
         fd = i;
         break;
      }
   }

   // Setup access mode
   switch( mode ) {
      case SUDS_READWRITE:
      case SUDS_APPEND:
         strcpy( access, "r+b" );
         break;
      case SUDS_CREATE:
         strcpy( access, "w+b" );
         break;
      default:
         mode = SUDS_READONLY;
         strcpy( access, "rb" );
         break;
   }

   // Open the file
   if( ( _suds_files[fd-1].stream = _fsopen( filespec, access, _SH_COMPAT ) ) == NULL ) {
      sprintf( _suds_err_buf, "ERROR: Unable to open file: %s\n       %s",
         filespec, _strerror( NULL ) );
      return( FALSE );
   }

   // Move to end of file in append mode
   if( mode == SUDS_APPEND ) {
      if( fseek( _suds_files[fd-1].stream, 0L, SEEK_END ) != 0 ) {
         sprintf( _suds_err_buf, "ERROR: Seek error in: %s\n       %s",
            filespec, _strerror( NULL ) );
         return( FALSE );
      }
   }

   // Increment number of open files
   _suds_nfiles++;

   // Fill in the file tag and return the file descriptor
   _fullpath( _suds_files[fd-1].filespec, filespec, _MAX_PATH );
   _s_upper( _suds_files[fd-1].filespec );
   _suds_files[fd-1].offset = ftell( _suds_files[fd-1].stream );
   _suds_files[fd-1].access = mode;
   _suds_files[fd-1].tag.id_struct = 0;
   _suds_files[fd-1].tag.len_struct = 0;
   _suds_files[fd-1].tag.len_data = 0;
   _suds_files[fd-1].seek_flag = FALSE;

   return( fd );
}

//---------------------------------------------------------------------
int suds_close( int fd ) {
   int ret;

   if( ! _valid_fd( fd ) )
      return( FALSE );

   // Close the file
   if( ( ret = fclose( _suds_files[fd-1].stream ) ) == EOF )
      sprintf( _suds_err_buf, "ERROR: Unable to close file: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );

   // Initialize the file tag for future use
   _suds_files[fd-1].stream = NULL;
   _suds_files[fd-1].filespec[0] = '\0';
   _suds_files[fd-1].offset = 0L;
   _suds_files[fd-1].access = 0;
   _suds_files[fd-1].seek_flag = FALSE;

   // Decrement number of open files
   _suds_nfiles--;

   // Check for error on close
   if( ret == EOF )
      return( FALSE );
   else
      return( TRUE );
}

//---------------------------------------------------------------------
int suds_read( int fd, SUDS *suds ) {
   void *ptr;

   suds->type = 0;
   suds->data_len = 0L;
   ptr = &suds->st;

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   // If last operation was a seek, there is no need to align
   if( _suds_files[fd-1].seek_flag ) {
      _suds_files[fd-1].seek_flag = FALSE;
   }
   else {
      // Align the file, in case any data following previous struct was not read
      if( _suds_files[fd-1].offset != 0L ||
         (_suds_files[fd-1].offset == 0L &&
          _suds_files[fd-1].tag.len_data != 0) ) {
         if( fseek( _suds_files[fd-1].stream,
                    _suds_files[fd-1].offset + (long)sizeof(SUDS_STRUCTTAG) +
                    _suds_files[fd-1].tag.len_struct + _suds_files[fd-1].tag.len_data,
                     SEEK_SET ) != 0 ) {
            if( feof( _suds_files[fd-1].stream ) )
               return( SUDS_EOF );
            sprintf( _suds_err_buf, "ERROR: Seek error: %s\n       %s",
               _suds_files[fd-1].filespec, _strerror( NULL ) );
            return( FALSE );
         }
      }
      // Save new file position
      if( ( _suds_files[fd-1].offset = ftell( _suds_files[fd-1].stream ) ) == -1L ) {
         sprintf( _suds_err_buf, "ERROR: Unable to determine position in: %s", _suds_files[fd-1].filespec );
         return( FALSE );
      }
   }

   // Read in the tag
   if( fread( &_suds_files[fd-1].tag, 1, sizeof(SUDS_STRUCTTAG), _suds_files[fd-1].stream  ) != sizeof(SUDS_STRUCTTAG) ) {
      if( feof( _suds_files[fd-1].stream ) )
         return( SUDS_EOF );
      sprintf( _suds_err_buf, "ERROR: Read error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   // Check file sync
   if( _suds_files[fd-1].tag.sync != 'S' ) {
      sprintf( _suds_err_buf, "ERROR: File out of sync: %s\n       May not be a SUDS file.", _suds_files[fd-1].filespec );
      return( FALSE );
   }

   // Check hardware data format
   if( _suds_files[fd-1].tag.machine != MACHINE ) {
      sprintf( _suds_err_buf, "ERROR: Unrecognized machine type in: %s", _suds_files[fd-1].filespec );
      return( FALSE );
   }

   // Read in the structure
   if( fread( ptr, 1, (int)_suds_files[fd-1].tag.len_struct, _suds_files[fd-1].stream  ) != (int)_suds_files[fd-1].tag.len_struct ) {
      if( feof( _suds_files[fd-1].stream ) )
         return( SUDS_EOF );
      sprintf( _suds_err_buf, "ERROR: Read error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   suds->type = _suds_files[fd-1].tag.id_struct;
   suds->data_len = _suds_files[fd-1].tag.len_data;

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_read_data( int fd, void _huge *ptr, long len ) {
   char _huge *p;
   long i, bytes_read;

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   // Read in data
   p = (char *)ptr;
   bytes_read = 0L;
   for( i=0; i<len/SUDS_IO_BLK; i++ ) {
      bytes_read += fread( p, 1, (size_t)SUDS_IO_BLK, _suds_files[fd-1].stream );
      if( bytes_read != (i+1)*SUDS_IO_BLK ) {
         sprintf( _suds_err_buf, "ERROR: Read error: %s\n       %s",
            _suds_files[fd-1].filespec, _strerror( NULL ) );
         return( FALSE );
      }
      p += SUDS_IO_BLK;
   }
   bytes_read += fread( p, 1, (size_t)(len-bytes_read), _suds_files[fd-1].stream );
   if( bytes_read != len ) {
      sprintf( _suds_err_buf, "ERROR: Read error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_write( int fd, SUDS *suds ) {
   void *ptr;

   ptr = &suds->st;

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   if( _suds_files[fd-1].access == 1 ) {
      sprintf( _suds_err_buf, "ERROR: File open in read only mode: %s", _suds_files[fd-1].filespec );
      return( FALSE );
   }

   if( suds->type < 1 || suds->type > TOTAL_STRUCTS ) {
      sprintf( _suds_err_buf, "ERROR: Unknown structure type: %d", suds->type );
      return( FALSE );
   }

   // Save new file position
   if( ( _suds_files[fd-1].offset = ftell( _suds_files[fd-1].stream ) ) == -1L ) {
      sprintf( _suds_err_buf, "ERROR: Unable to determine position in: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   _suds_files[fd-1].tag.sync = 'S';
   _suds_files[fd-1].tag.machine = MACHINE;
   _suds_files[fd-1].tag.id_struct = suds->type;
   _suds_files[fd-1].tag.len_struct = (long)_suds_size_struct[suds->type];
   _suds_files[fd-1].tag.len_data = suds->data_len;

   if( fwrite( &_suds_files[fd-1].tag, 1, sizeof(SUDS_STRUCTTAG), _suds_files[fd-1].stream  ) != sizeof(SUDS_STRUCTTAG) ) {
      sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   if( fwrite( ptr, 1, (int)_suds_size_struct[suds->type], _suds_files[fd-1].stream  ) != (int)_suds_size_struct[suds->type] ) {
      sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_write_data( int fd, void _huge *ptr, long len ) {
   char _huge *p;
   long i, bytes_written;

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   // Write out data
   p = (char *)ptr;
   bytes_written = 0L;
   for( i=0; i<len/SUDS_IO_BLK; i++ ) {
      bytes_written += fwrite( p, 1, (size_t)SUDS_IO_BLK, _suds_files[fd-1].stream );
      if( bytes_written != (i+1)*SUDS_IO_BLK ) {
         sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
            _suds_files[fd-1].filespec, _strerror( NULL ) );
         return( FALSE );
      }
      p += SUDS_IO_BLK;
   }
   bytes_written += fwrite( p, 1, (size_t)(len-bytes_written), _suds_files[fd-1].stream );
   if( bytes_written != len ) {
      sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_update( int fd, SUDS *suds ) {
   void *ptr;
   long pos;

   ptr = &suds->st;

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   if( _suds_files[fd-1].access == 1 ) {
      sprintf( _suds_err_buf, "ERROR: File open in read only mode: %s", _suds_files[fd-1].filespec );
      return( FALSE );
   }

   if( suds->type < 1 || suds->type > TOTAL_STRUCTS ) {
      sprintf( _suds_err_buf, "ERROR: Unknown structure type: %d", suds->type );
      return( FALSE );
   }

   if( ( pos = ftell( _suds_files[fd-1].stream ) ) == -1L ) {
      sprintf( _suds_err_buf, "ERROR: Unable to determine position in: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   if( fseek( _suds_files[fd-1].stream, _suds_files[fd-1].offset, SEEK_SET ) != 0 ) {
      if( feof( _suds_files[fd-1].stream ) )
         return( SUDS_EOF );
      sprintf( _suds_err_buf, "ERROR: Seek error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   _suds_files[fd-1].tag.sync = 'S';
   _suds_files[fd-1].tag.machine = MACHINE;
   _suds_files[fd-1].tag.id_struct = suds->type;
   _suds_files[fd-1].tag.len_struct = (long)_suds_size_struct[suds->type];
   _suds_files[fd-1].tag.len_data = suds->data_len;

   if( fwrite( &_suds_files[fd-1].tag, 1, sizeof(SUDS_STRUCTTAG), _suds_files[fd-1].stream  ) != sizeof(SUDS_STRUCTTAG) ) {
      sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   if( fwrite( ptr, 1, (int)_suds_size_struct[suds->type], _suds_files[fd-1].stream  ) != (int)_suds_size_struct[suds->type] ) {
      sprintf( _suds_err_buf, "ERROR: Write error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   if( fseek( _suds_files[fd-1].stream, pos, SEEK_SET ) != 0 ) {
      if( feof( _suds_files[fd-1].stream ) )
         return( SUDS_EOF );
      sprintf( _suds_err_buf, "ERROR: Seek error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_flush( int fd ) {

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   if( fflush( _suds_files[fd-1].stream ) != 0 ) {
      sprintf( _suds_err_buf, "ERROR: Unable to flush: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
long suds_pos( int fd ) {

   _suds_err_buf[0] = '\0';

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( (long)SUDS_EOF );

   // Return position of last structure read
   return( _suds_files[fd-1].offset );
}

//---------------------------------------------------------------------
long suds_abs_pos( int fd ) {

   _suds_err_buf[0] = '\0';

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( (long)SUDS_EOF );

   // Return absolute position
   return( ftell( _suds_files[fd-1].stream ) );
}

//---------------------------------------------------------------------
int suds_seek( int fd, long pos ) {

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   // Seek
   if( fseek( _suds_files[fd-1].stream, pos, SEEK_SET ) != 0 ) {
      if( feof( _suds_files[fd-1].stream ) )
         return( SUDS_EOF );
      sprintf( _suds_err_buf, "ERROR: Seek error: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   // Save position in file tag
   if( ( _suds_files[fd-1].offset = ftell( _suds_files[fd-1].stream ) ) == -1L ) {
      sprintf( _suds_err_buf, "ERROR: Unable to determine position in: %s\n       %s",
         _suds_files[fd-1].filespec, _strerror( NULL ) );
      return( FALSE );
   }

   if( pos == 0 ) {
      _suds_files[fd-1].tag.id_struct = 0;
      _suds_files[fd-1].tag.len_struct = 0;
      _suds_files[fd-1].tag.len_data = 0;
   }

   // Set seek flag
   _suds_files[fd-1].seek_flag = TRUE;

   return( TRUE );
}

//---------------------------------------------------------------------
int suds_rewind( int fd ) {

   // Check out file descriptor
   if( ! _valid_fd( fd ) )
      return( FALSE );

   rewind( _suds_files[fd-1].stream );
   _suds_files[fd-1].offset = 0L;
   _suds_files[fd-1].tag.id_struct = 0;
   _suds_files[fd-1].tag.len_struct = 0;
   _suds_files[fd-1].tag.len_data = 0;
   _suds_files[fd-1].seek_flag = TRUE;

   return( TRUE );
}

//---------------------------------------------------------------------
char *suds_get_err( void ) {
   return( _suds_err_buf );
}

//---------------------------------------------------------------------
// Helper functions follow
//---------------------------------------------------------------------
int _valid_fd( int fd ) {

   if( fd < 1 || fd > SUDS_MAX_FILES ||
       _suds_files[fd-1].stream == NULL ||
       _suds_files[fd-1].access == 0 ) {
      sprintf( _suds_err_buf, "ERROR: Invalid file descriptor!" );
      return( FALSE );
   }

   return( TRUE );
}

//---------------------------------------------------------------------
void _s_upper( char *buffer ) {
   char *p;
   for( p = buffer; *p; p++ )
      *p = toupper( *p );
}
