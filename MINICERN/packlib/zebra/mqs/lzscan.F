*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION LZSCAN (IXDIVP,LGOP)

C-    Find first bank in division IXDIV after bank at LGO
C-    user called

#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXDIVP(9), LGOP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HLZSC, 4HAN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HLZSCAN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'LZSCAN  ')
#endif

#include "zebra/q_jbit.inc"

      IXIN = IXDIVP(1)
      LGO  = LGOP(1)

      CALL MZSDIV (IXIN,4)
      LSTA = LQSTA(KQT+JQDIVI)
      IQNX = LSTA
      LEND = LQEND(KQT+JQDIVI)

      IF (LGO.EQ.0)                GO TO 31

      IF (LGO.LT.LSTA)             GO TO 91
      IF (LGO.GE.LEND)             GO TO 91

      CALL MZCHLS (-7,LGO)
      IF (IQFOUL.NE.0)             GO TO 92

   31 LN   = IQNX
      IF (LN.GE.LEND)              GO TO 39
      CALL MZCHLN (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 93
      IF (IQND.LT.0)               GO TO 31
      IQUEST(1) = JBIT(IQ(KQS+IQLS),IQDROP)
      LZSCAN = IQLS
      RETURN

   39 LZSCAN = 0
      RETURN

C------            ERROR CONDITIONS

   93 NQCASE = 1
      NQFATA = 1
      IQUEST(17) = LN
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 6
      IQUEST(11) = IXIN
      IQUEST(12) = 0
      IQUEST(13) = 0
      IQUEST(14) = LGO
      IQUEST(15) = LSTA
      IQUEST(16) = LEND
#include "zebra/qtrace.inc"
#include "zebra/qtofatal.inc"
      LZSCAN = 0
      END
*      ==================================================
#include "zebra/qcardl.inc"
