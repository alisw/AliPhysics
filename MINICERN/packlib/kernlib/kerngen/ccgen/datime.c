/*
 * $Id$
 *
 * $Log$
 * Revision 1.5  1997/12/19 16:36:06  mclareni
 * After 2000, the date in ID, ND will have the old format with the year as 2 digits
 *
 * Revision 1.4  1997/12/15 16:52:27  mclareni
 * Make length of structure slate 40, same as common slate in all the other routines
 *
 * Revision 1.3  1997/09/02 14:26:35  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:15  mclareni
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

#if defined(CERNLIB_QMIRTD)
#include "irtdgs/datime.c"
#else
/*>    ROUTINE DATIME
  CERN PROGLIB# Z007    DATIME          .VERSION KERNFOR  4.40  940929
*/
#if !defined(CERNLIB_QMOS9)
#include <sys/types.h>
#endif
#include <time.h>

#if defined(CERNLIB_QX_SC)
#define slate slate_
struct { int  inum[40]; } slate_;
void type_of_call datime_(id, it)
#endif
#if defined(CERNLIB_QXNO_SC)
struct { int  inum[40]; } slate;
void type_of_call datime(id, it)
#endif
#if defined(CERNLIB_QXCAPT)
#define slate SLATE
struct { int  inum[40]; } SLATE;
void type_of_call DATIME(id, it)
#endif
   int  *id, *it;
{
      struct tm *tp;
#if defined(CERNLIB_QMAPO)
      int   nsl;
#endif

#if (defined(CERNLIB_QMAPO))&&(defined(CERNLIB_QX_SC))
      void  type_of_call toslat_();
#endif
#if (defined(CERNLIB_QMAPO))&&(defined(CERNLIB_QXNO_SC))
      void  type_of_call toslat();
#endif

#if defined(CERNLIB_QXCAPT)
      void  type_of_call TOSLAT();
#endif

   time_t tloc = time(0);
   tp = localtime(&tloc);
   slate.inum[0] = tp->tm_year + 1900;
   slate.inum[1] = tp->tm_mon + 1;
   slate.inum[2] = tp->tm_mday;
   slate.inum[3] = tp->tm_hour;
   slate.inum[4] = tp->tm_min;
   slate.inum[5] = tp->tm_sec;
#if defined(CERNLIB_QMAPO)
      nsl = 6;
#endif
#if (defined(CERNLIB_QMAPO))&&(defined(CERNLIB_QX_SC))
      toslat_ (slate.inum, &nsl);
#endif
#if (defined(CERNLIB_QMAPO))&&(defined(CERNLIB_QXNO_SC))
      toslat (slate.inum, &nsl);
#endif
   *id  = (tp->tm_year % 100 ) * 10000;
   *id += (tp->tm_mon + 1) * 100;
   *id += tp->tm_mday;
   *it  = tp->tm_hour * 100;
   *it += tp->tm_min;
   return;
}
/*> END <----------------------------------------------------------*/
#endif
