/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:21  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:29  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:23  mclareni
 * Kernlib
 * 
 */
#include "kerngen/pilot.h"
#if defined(CERNLIB_WINNT)
#include "wntgs/intrac.c"
#elif defined(CERNLIB_QMDOS)
#include "dosgs/intrac.c"
#else
/*>    ROUTINE INTRAC
  CERN PROGLIB# Z044    INTRAC          .VERSION KERNFOR  4.39  940228
*/
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
    return ((int) isatty(0));
}
/*> END <----------------------------------------------------------*/
#endif
