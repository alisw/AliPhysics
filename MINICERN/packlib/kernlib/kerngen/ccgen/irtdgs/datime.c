/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:30  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE DATIME
  CERN PROGLIB# Z007    DATIME          .VERSION KERNIRT  1.06  930811
*/
#include <sys/types.h>
#include <time.h>

#define slate slate_
struct { int  inum[79]; } slate_;
void datime_(id, it)
#if defined(CERNLIB_QXCAPT)
#define slate SLATE
struct { int  inum[79]; } slate;
void DATIME(id, it)
#endif
   int  *id, *it;
{
      struct tm *tp;

   time_t tloc = time(0);
   tp = localtime(&tloc);
   slate.inum[0] = tp->tm_year + 1900;
   slate.inum[2] = tp->tm_mon + 1;
   slate.inum[4] = tp->tm_mday;
   slate.inum[6] = tp->tm_hour;
   slate.inum[8] = tp->tm_min;
   slate.inum[10] = tp->tm_sec;
   *id  = tp->tm_year * 10000;
   *id += (tp->tm_mon + 1) * 100;
   *id += tp->tm_mday;
   *it  = tp->tm_hour * 100;
   *it += tp->tm_min;
   return;
}
/*> END <----------------------------------------------------------*/
