*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:16  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZROOT

C-    RESET TO ROOT LEVEL FOR QNEXT

#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZRO, 4HOT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZROOT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZROOT  ')
#endif

#include "zebra/q_jbit.inc"

      IF (JQLEV.LT.0)              RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"

      IF (JQLEV.EQ.0)              GO TO 31
   21 CALL JZOUT ('??? ')
      IF (JQLEV.GT.0)        GO TO 21

      IQ(KQS+LQSV+2) = IQ(KQS+LQSV+2) + 1

#if defined(CERNLIB_QDEBUG)
      IF (JBIT(JQLLEV,13).NE.0)  WRITE (IQLOG,9024)
 9024 FORMAT    (/' =======  JZROOT done.')
#endif

   31 CONTINUE
#if defined(CERNLIB_JZTIME)
#include "zebra/jztimout.inc"
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
