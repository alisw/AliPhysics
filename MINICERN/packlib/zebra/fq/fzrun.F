*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZRUN (LUNPAR,JRUNP,NUHP,IUHEAD)

C-    Write start/end of run
C-        JRUN > 0  start this run
C-             = 0  bump current run number by 1
C-             < 0  end of run

#include "zebra/mqsysh.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNPAR(9), JRUNP(9), NUHP(9), IUHEAD(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZRU, 4HN    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZRUN  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZRUN   ')
#endif


*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      JRUN  = JRUNP(1)

      CALL FZOUT (LUNPAR,0,JRUN,-1,'FZRUN',1,NUHP,IUHEAD)

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
