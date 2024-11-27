// SUDS Structure initialization routines

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
//    Started fresh

// Includes -----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#define _SUDS_INIT_
#include "suds.h"
#include "sudsform.h"

// Local prototypes ---------------------------------------------------
void suds_init( SUDS *suds );
void _suds_init_struct( void *ptr, int type );
void _suds_init_value( int i, void *ptr );

//---------------------------------------------------------------------
void suds_init( SUDS *suds ) {

   suds->data_len = 0L;

   _suds_init_struct( &suds->st, suds->type );

   return;
}

//---------------------------------------------------------------------
void _suds_init_struct( void *ptr, int type ) {
   register i;
   void *ptr1, *ptr2;

   ptr1 = ptr;

   for( i=_suds_beg_struct[type]; _suds_form[i].fstype==type; i++ ) {

      ptr2 = (char *)ptr1+_suds_form[i].offset;

      if( _suds_form[i].nextfstype != 0 )
         _suds_init_struct( ptr2, (int)_suds_form[i].nextfstype );
      else
         _suds_init_value( i, ptr2 );
   }
}

//---------------------------------------------------------------------
void _suds_init_value( int i, void *ptr ) {

   switch( _suds_form[i].ftype ) {
      case CHR:
         *(char *)ptr = _suds_form[i].initval[0];
         break;
      case BTS:
      case MIN:
         *(char *)ptr = atoi( _suds_form[i].initval );
          break;
      case STR:
         strcpy( ptr, _suds_form[i].initval );
         break;
      case SHT:
         *(short *)ptr = atoi( _suds_form[i].initval );
         break;
      case LNG:
      case STT:
         *(long *)ptr = atol( _suds_form[i].initval );
         break;
      case FLT:
         *(float *)ptr = (float)atof( _suds_form[i].initval );
         break;
      case LLT:
      case DBL:
      case MST:
         *(double *)ptr = atof( _suds_form[i].initval );
         break;
      case BTW:
         *(unsigned short *)ptr = atoi( _suds_form[i].initval );
         break;
      case CAL:
      case CPX:
      case STI:
      default:
         break;
   }
}
