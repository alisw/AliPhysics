/*
 * $Id$
 *
 * $Log$
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
#if defined(CERNLIB_QX_SC)
int intrac_()
#endif
#if defined(CERNLIB_QXNO_SC)
int intrac()
#endif
#if defined(CERNLIB_QXCAPT)
int INTRAC()
#endif
{
#ifdef WIN32
    return (((int) isatty(0)!=0) ? 1 : 0) ;
#else
    return ((int) isatty(0));
#endif
}
/*> END <----------------------------------------------------------*/
