/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:30  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE GMTIMEF (CLOCK, TARR)
  CERN PROGLIB# Z265    GMTIMEF         .VERSION KERNIRT  1.06  930811
  ORIG. 14/03/91, RDM
  Fortran interface routine to gmtime

     CLOCK  encoded time (returned by, e.g. STATF)
     TARR   decoded time (INTEGER TARR(9))
*/
#include <stdio.h>
#include <time.h>

#if defined(CERNLIB_QX_SC)
void gmtimef_(clock, tarr)
#endif
#if defined(CERNLIB_QXNO_SC)
void gmtimef(clock, tarr)
#endif
#if defined(CERNLIB_QXCAPT)
void GMTIMEF(clock, tarr)
#endif
      time_t *clock;
      int    *tarr;
{
    struct tm *gmtime(), *tm;

    tm = gmtime(clock);
    tarr[0] = tm->tm_sec;
    tarr[2] = tm->tm_min;
    tarr[4] = tm->tm_hour;
    tarr[6] = tm->tm_mday;
    tarr[8] = tm->tm_mon;
    tarr[10] = tm->tm_year;
    tarr[12] = tm->tm_wday;
    tarr[14] = tm->tm_yday;
    tarr[16] = tm->tm_isdst;
}
/*> END <----------------------------------------------------------*/
