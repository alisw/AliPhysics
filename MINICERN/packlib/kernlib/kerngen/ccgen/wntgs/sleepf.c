/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  1997/02/04 17:35:06  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SLEEPF (NSECS)
  CERN PROGLIB# Z265    SLEEPF          .VERSION KERNFOR  4.26  910313
  ORIG. 22/02/91, JZ
  Fortran interface routine to sleep
*/
#include <windows.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call sleepf_(seconds)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call sleepf(seconds)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call SLEEPF(seconds)
#endif
      int  *seconds;
{
      void sleep();
      int  secu;

      secu = *seconds;
#ifdef WIN32
      Sleep(secu*1000);
#else
      sleep(secu);
#endif
      return;
}
/*> END <----------------------------------------------------------*/
