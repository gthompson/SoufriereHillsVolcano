// SUDS time routines

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
// Version 1.45 Apr-1999,RLB
//    Y2K Compliant, added parse_mstime( ).

// Version 1.40 Sep-1992, RLB
//    Cleaned up original source provided by Bruce Julian and PLWard,
//    U. S. Geological Survey, Menlo Park, Ca. 94025

// Includes -----------------------------------------------------------
#include <sys\types.h>
#include <sys\timeb.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define _SUDS_TIME_
#include "suds.h"

// Constants ----------------------------------------------------------
#define BASEJDN         2440588       // 1 January 1970

// Calendar codes
#define GREGORIAN       1             // Modern calendar
#define JULIAN          0             // Old-style calendar

// Time intervals
#define SECOND          1L
#define MINUTE          (60L*SECOND)
#define HOUR            (60L*MINUTE)
#define DAY             (24L*HOUR)

#define MAXFORM 12

// Prototypes, for external use ---------------------------------------
double get_mstime( void );
double make_mstime( int, int, int, int, int, double );
int    decode_mstime( double, int *, int *, int *, int *, int *, double * );
char   *list_mstime( double, int );
double parse_mstime( char * buf, int mndy );
int    yrday( int, int, int );
void   mnday( int, int, int *, int * );
int    isleap( int, int );

// Prototypes, for internal use ---------------------------------------
long   jdn( int, int, int, int );
void   date( long, int *, int *, int *, int );
char *_next_field( char * ptr );

// Macros -------------------------------------------------------------
// floor(x/y), where x, y>0 are integers, using integer arithmetic
#define qfloor( x, y ) ( x>0 ? (x)/y : -((y-1-(x))/y) )

// Globals ------------------------------------------------------------
int  eom[2][15] = {
   { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365, 396, 424 },
   { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366, 397, 425 },
};
char *mon[]={ "", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

char outstr[34];              // Char buffer for list_mstime( )
extern char *_suds_err_buf;   // Error message buffer

// Recognize leap years -----------------------------------------------
int isleap( int yr, int cal) {
   int l;

   if( yr < 0 )
      yr++;
   l = ( yr%4 == 0 );
   if( cal == GREGORIAN )
      l = l && ( yr%100 != 0 || yr%400 == 0 );

   return( l );
}

// Compute day of year ------------------------------------------------
int yrday( int mo, int day, int lp ) {

   return( eom[lp][mo-1] + day );
}

// Calculate month, day from day of year ------------------------------
void mnday( int d, int lp, int *pm, int *pd ) {
   int  i;
   int  *et = eom[lp];

   for( i=1; d>et[i]; i++ );
   *pm = i;
   *pd = d - et[i-1];

   return;
}

// Compute Julian Day number from calendar date -----------------------
long jdn( int yr, int mo, int day, int cal ) {
   long ret;

   if( yr < 0 )
      yr++;

   // Move Jan. & Feb. to end of previous year
   if( mo <= 2 ) {
      --yr;
      mo += 12;
   }

   ret = qfloor( (long)(4*365+1) * (yr+4712), 4 ) + eom[0][mo-1] + day +
      ( cal==GREGORIAN ? -qfloor(yr, 100) + qfloor(yr, 400) + 2 : 0 );

   return( ret );
}

// Compute calendar date from Julian Day Number -----------------------
void date( long n, int *py, int *pm, int *pd, int cal ) {
   // n              Julian day number
   // *py, *pm, *pd  Year, month, day
   // cal            Calendar type

   long d, t;
   int y;

// Find position within cycles that are nd days long
#define CYCLE( n, nd ) { t=qfloor(d-1,nd); y+=t*n; d-=t*nd; }

// The same, with bound on cycle number
#define LCYCLE( n, nd, l ) { t=qfloor(d-1,nd); if (t>l) t=l; y+=t*n; d-=t*nd; }

   y = -4799;
   if( cal == GREGORIAN ) {
      d = n + 31739;            // JD -31739 = 31 Dec 4801 B.C.
      CYCLE( 400, 146097 )      // Four-century cycle
      LCYCLE( 100, 36524, 3 )   // 100-year cycle
   }
   else
      d = n + 31777;            // JD -31777 = 31 Dec 4801 B.C.

   CYCLE( 4, 1461 )             // Four-year cycle
   LCYCLE( 1, 365, 3 )          // Yearly cycle

   if( y <= 0 )
      --y;
   *py = y;

   mnday( (int)d, isleap(y, cal), pm, pd );

   return;
}

// Return system time -------------------------------------------------
double get_mstime( void ) {
   double systime;      
   
#ifdef _MSC_VER
   struct _timeb tp;
      
   _ftime (&tp);
   
   if( tp.dstflag )
      systime = (double)tp.time + ((double)tp.millitm / 1000.0) - (((double)tp.timezone - 60.0) * 60.0);
   else
      systime = (double)tp.time + ((double)tp.millitm / 1000.0) - ((double)tp.timezone * 60.0);
#else
   struct timeb tp;
      
   ftime (&tp);

   if( tp.dstflag )
      systime = (double)tp.time + ((double)tp.millitm / 1000.0) - (((double)tp._timezone - 60.0) * 60.0);
   else
      systime = (double)tp.time + ((double)tp.millitm / 1000.0) - ((double)tp._timezone * 60.0);
#endif    

   return( systime );
}

// Build ms_time from components --------------------------------------
double make_mstime( int year, int month, int day, int hour,
                    int min, double second) {
   int err;
   int *et = eom[isleap(year,GREGORIAN)];

   err = 0;
   if( year == 0 ) {
      sprintf( _suds_err_buf,"ERROR: Bad year value: %d", year );
      err = 1;
   }

   if( month == 0 )
      mnday( day, isleap( year, GREGORIAN ), &month, &day );

   if( month < 1 || month > 12 ) {
      sprintf( _suds_err_buf,"ERROR: Bad month value: %d", month );
      err = 1;
   }

   if( day < 1 || day > et[month]-et[month-1] ) {
      sprintf( _suds_err_buf,"ERROR: Bad day value: %d", day );
      err = 1;
   }

   if( err )
      return( 0.0 );

   return( (double) DAY*(jdn(year,month,day,GREGORIAN)-BASEJDN) +
      hour*HOUR + min*MINUTE + second*SECOND );
}

// Decode mstime -------------------------------------------------------
int decode_mstime( double time, int *year, int *month, int *day,
                   int *hour, int *min, double *second) {
   long d;

   d = (long)floor( (long)time/DAY );
   date( d+BASEJDN, year, month, day, GREGORIAN );
   time -= (double)(d*DAY);

   *hour = (int)((long)time/HOUR);
   time -= (double)((long)(*hour)*HOUR);

   *min = (int)((long)time/MINUTE);

   *second = time - (double)((long)(*min)*MINUTE);

   return( 1 );
}

// Display mstime ------------------------------------------------------
char *list_mstime( double time, int form ) {
   int year, month, day, hour, min;
   double second;

   if( form < 0 || form > MAXFORM )
      form=0;

   if( decode_mstime( time, &year, &month, &day, &hour, &min, &second ) == 0 )
      return( 0 );

   switch( form ) {
      case 0:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d%02d%02d%02d%02d%06.3f",
                  year, month, day, hour, min, second );
         break;
      case 1:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d%02d%02d%02d%02d%02d",
                  year, month, day, hour, min, (int)second );
         break;
      case 2:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d %02d %02d% 02d %02d %06.3f",
                  year, month, day, hour, min, second );
         break;
      case 3:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d %02d %02d% 02d %02d %02d",
                  year, month, day, hour, min, (int)second );
         break;
      case 4:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d/%02d/%02d %02d:%02d:%06.3f",
                  month, day, year, hour, min, second );
         break;
      case 5:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d/%02d/%02d %02d:%02d:%02d",
                  month, day, year, hour, min, (int)second );
         break;
      case 6:
         sprintf( outstr, "%s %d, %d %02d:%02d %06.3f GMT",
                  mon[month], day, year, hour, min, second );
         break;
      case 7:
         sprintf( outstr, "%s %d, %d %02d:%02d %02d GMT",
                  mon[month], day, year, hour, min, (int)second );
         break;
      case 8:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d %03d %02d:%02d:%06.03f",
                 year, yrday( month, day, isleap( year, 1 ) ), hour, min, second );
         break;
      case 9:
         sprintf( outstr, "%08lX", (long)time );
         break;
      case 10:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d*%03d+%02d:%02d:%06.03f",
                  year, yrday( month, day, isleap( year, 1 ) ), hour, min, second );
         break;
      case 11:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d/%02d/%02d %02d:%02d:%09.6f",
                  month, day, year, hour, min, second );
         break;
      case 12:
         year = (year < 2000 ? year - 1900 : year - 2000 );
         sprintf( outstr, "%02d*%03d+%02d:%02d:%09.06f",
                  year, yrday( month, day, isleap( year, 1 ) ), hour, min, second );
         break;
   }

   return( outstr );
}
                                                                 
/*--------------------------------------------------------------------- */
double parse_mstime( char * buf, int mndy )
   {
   char *ptr;
   int year, mon, day, doy, hour, min;
   double sec, time;

   ptr = buf;

   // Default year, month, and day to current time, hour, min, and sec to 0
   time = get_mstime(  );
   decode_mstime( time, &year, &mon, &day, &hour, &min, &sec );
   doy = yrday( mon, day, isleap( year, 1 ) );
   hour = 0;
   min = 0;
   sec = 0.0;

   // First field
   if ( *ptr )
      {
      if ( mndy )
         mon = atoi( ptr );            // Month
      else
         {
         year = atoi( ptr );           // Year  
         if ( year < 70 )
         	year += 2000;
         else if ( year < 100 )
            year += 1900;
         }
      }

   // Second field
   ptr = _next_field( ptr );
   if ( *ptr )
      {                                                           
      if ( mndy )
         day = atoi( ptr );            // Day
      else
         doy = atoi( ptr );            // Day of year
      }

   // Third field
   ptr = _next_field( ptr );
   if ( *ptr )
      {
      if ( mndy )
         {
         year = atoi( ptr );           // Year
         if ( year < 70 )
         	year += 2000;
         else if ( year < 100 )
            year += 1900;
         }
      else
         hour = atoi( ptr );           // Hour
      }

   // Fourth field
   ptr = _next_field( ptr );
   if ( *ptr )
      {
      if ( mndy )
         hour = atoi( ptr );           // Hour
      else
         min = atoi( ptr );            // Minute
      }

   // Fifth field
   ptr = _next_field( ptr );
   if ( *ptr )
      {
      if ( mndy )
         min = atoi( ptr );            // Minute
      else
         sec = atof( ptr );            // Second
      }

   if ( mndy )
      {
      // Sixth field
      ptr = _next_field( ptr );
      if ( *ptr )
         {
         sec = atof( ptr );            // Second
         }
      }
   else
      mnday( doy, isleap( year, 1 ), &mon, &day );

   time = make_mstime( year, mon, day, hour, min, sec );

   return ( time );
   }

/*--------------------------------------------------------------------- */
char *_next_field( char * ptr )
   {

   if ( *ptr )
      {
      while ( *ptr && isdigit( *ptr ) )
         ptr++;

      while ( *ptr && !isdigit( *ptr ) )
         ptr++;
      }

   return ( ptr );
   }

                                                                 