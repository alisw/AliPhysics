/*
 * $Id$
 *
 * $Log$
 * Revision 1.6  1997/12/19 16:36:09  mclareni
 * After 2000, the date in ID, ND will have the old format with the year as 2 digits
 *
 * Revision 1.5  1997/11/05 10:35:32  mclareni
 * Remove the last WNT mod
 *
 * Revision 1.3  1997/09/02 14:26:35  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:15  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:26  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
 * Kernlib
 *  
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE DATIMH
  CERN PROGLIB# Z007    DATIMH          .VERSION KERNFOR  4.40  940929
*/
#if !defined(CERNLIB_QMOS9)
#include <sys/types.h>
#endif
#include <time.h>

#if defined(CERNLIB_QX_SC)
void type_of_call datimh_(dh, th)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call datimh(dh, th)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call DATIMH(dh, th)
#endif
   char dh[7], th[7];
{
      struct tm *tp;

   time_t tloc = time(0);
   tp = localtime(&tloc);
   dh[0] = tp->tm_mday / 10 + '0';
   dh[1] = tp->tm_mday % 10 + '0';
   dh[2] = '/';
   dh[3] = (tp->tm_mon + 1) / 10 + '0';
   dh[4] = (tp->tm_mon + 1) % 10 + '0';
   dh[5] = '/';
   dh[6] = (tp->tm_year % 100) / 10 + '0';
   dh[7] = (tp->tm_year % 100) % 10 + '0';
   th[0] = tp->tm_hour / 10 + '0';
   th[1] = tp->tm_hour % 10 + '0';
   th[2] = '.';
   th[3] = tp->tm_min  / 10 + '0';
   th[4] = tp->tm_min  % 10 + '0';
   th[5] = '.';
   th[6] = tp->tm_sec  / 10 + '0';
   th[7] = tp->tm_sec  % 10 + '0';
   return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_DATIMH
#undef CERNLIB_TCGEN_DATIMH
#endif
