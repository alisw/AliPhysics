/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  1997/02/04 17:35:05  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE INTRAC
  CERN PROGLIB# Z044    INTRAC          .VERSION KERNFOR  4.38  931108
*/
#ifdef WIN32
#include <io.h>
#endif

#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
int type_of_call intrac_()
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call intrac()
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call INTRAC()
#endif
{
    return (((int) isatty(0)!=0) ? 1 : 0) ;
}
/*> END <----------------------------------------------------------*/
