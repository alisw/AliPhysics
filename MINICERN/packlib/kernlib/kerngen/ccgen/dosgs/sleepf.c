/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SLEEPF (NSECS)
  CERN PROGLIB# Z265    SLEEPF          .VERSION KERNFOR  4.26  910313
  ORIG. 22/02/91, JZ
  Fortran interface routine to sleep
*/
#ifdef WIN32
#include <windows.h>
#endif
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
#ifdef WIN32
      Sleep(secu*1000));
#else
      sleep(secu);
#endif
      return;
}
/*> END <----------------------------------------------------------*/
