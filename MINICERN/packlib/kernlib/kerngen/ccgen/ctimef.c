/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:14  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:24  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CTIMEF (CLOCK, STIME)
  CERN PROGLIB# Z265    CTIMEF          .VERSION KERNFOR  4.36  930602
  ORIG. 14/03/91, RDM
  Fortran interface routine to ctime

     CLOCK  encoded time (returned by, e.g. STATF)
     STIME  decoded time string of length 24 (CHARACTER*24 STIME)
*/
#if defined(CERNLIB_QX_SC)
void type_of_call ctimef_(clock, stime)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ctimef(clock, stime)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CTIMEF(clock, stime)
#endif
int  *clock;
char *stime;
{
    char *ctime();

    strncpy(stime,ctime(clock),24);
}
/*> END <----------------------------------------------------------*/
