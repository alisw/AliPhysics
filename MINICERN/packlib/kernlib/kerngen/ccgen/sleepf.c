/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:43  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:43  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#if defined(CERNLIB_WINNT)
#include "wntgs/sleepf.c"
#elif defined(CERNLIB_QMDOS)
#include "dosgs/sleepf.c"
#else
/*>    ROUTINE SLEEPF (NSECS)
  CERN PROGLIB# Z265    SLEEPF          .VERSION KERNFOR  4.26  910313
  ORIG. 22/02/91, JZ
  Fortran interface routine to sleep
*/
#if defined(CERNLIB_QX_SC)
void sleepf_(seconds)
#endif
#if defined(CERNLIB_QXNO_SC)
void sleepf(seconds)
#endif
#if defined(CERNLIB_QXCAPT)
void SLEEPF(seconds)
#endif
      int  *seconds;
{
      void sleep();
      int  secu;

      secu = *seconds;
      sleep(secu);
      return;
}
/*> END <----------------------------------------------------------*/
#endif
