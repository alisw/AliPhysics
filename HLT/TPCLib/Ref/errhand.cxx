/************************* Start of ERRHAND.C ************************
 *
 * This is a general purpose error handler used with every program in
 * the book.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "errhand.h"
/*
  #ifdef __STDC__
  
  
  #else
  #ifdef __UNIX__
  void fatal_error( fmt, va_alist )
  char *fmt;
  va_dcl
  #else
  void fatal_error( fmt )
  char *fmt;
  #endif
  #endif
*/
void fatal_error( char *fmt, ... )
{
    va_list argptr;

    va_start( argptr, fmt );
    printf( "Fatal error: " );
    vprintf( fmt, argptr );
    va_end( argptr );
    exit( -1 );
}

/************************** End of ERRHAND.C *************************/


