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
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMIRTD)
#include "irtdgs/gmtimef.c"
#else
/*>    ROUTINE GMTIMEF (CLOCK, TARR)
  CERN PROGLIB# Z265    GMTIMEF         .VERSION KERNFOR  4.32  920229
  ORIG. 14/03/91, RDM
  Fortran interface routine to gmtime

     CLOCK  encoded time (returned by, e.g. STATF)
     TARR   decoded time (INTEGER TARR(9))
*/
#include <stdio.h>
#include <time.h>
#if defined(CERNLIB_QMAPO)||defined(CERNLIB_QMALT)
#include <sys/types.h>
#endif

#if defined(CERNLIB_QX_SC)
void type_of_call gmtimef_(clock, tarr)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call gmtimef(clock, tarr)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call GMTIMEF(clock, tarr)
#endif
      time_t *clock;
      int    *tarr;
{
    struct tm *gmtime(), *tm;

    tm = gmtime(clock);
    tarr[0] = tm->tm_sec;
    tarr[1] = tm->tm_min;
    tarr[2] = tm->tm_hour;
    tarr[3] = tm->tm_mday;
    tarr[4] = tm->tm_mon;
    tarr[5] = tm->tm_year;
    tarr[6] = tm->tm_wday;
    tarr[7] = tm->tm_yday;
    tarr[8] = tm->tm_isdst;
}
/*> END <----------------------------------------------------------*/
#endif
