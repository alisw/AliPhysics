*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZMEMORY)
      SUBROUTINE FZMEMO (LUNP,MBUF,NWBUF)

C-    Connect user buffer for M mode, User called

#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcf.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), MBUF(99)

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZME, 4HMO   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZMEMO /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZMEMO  ')
#endif

#include "zebra/q_locf.inc"

*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      LUN   = LUNP(1)

      CALL FZLOC (LUN,0)
      IF (LUNF.EQ.0)    LOGLVF=1
      IF (IFIFOF.NE.3)  LOGLVF=1

      IQUEST(2) = LOCF (MBUF(1))
      IQUEST(3) = IQUEST(2) - LQASTO
      IQUEST(7) = NWBUF
      IQUEST(8) = IQUEST(7)
#if !defined(CERNLIB_B32)
      IF (IUPAKF.EQ.0)  THEN
#endif
#if !defined(CERNLIB_B32)
          IQUEST(8) = (IQBITW*IQUEST(8)) / 32
#endif
#if !defined(CERNLIB_B32)
        ENDIF
#endif
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.1)   WRITE (IQLOG,9024) LUN
 9024 FORMAT (1X/' FZMEMO.  Connect user buffer for LUN=',I3)

#endif
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVF.GE.1)  THEN
          IQUEST(4) = LOCF (LQ(IQUEST(3)))
          WRITE (IQLOG,9824) IQUEST(2),IQUEST(3),IQUEST(4)
          WRITE (IQLOG,9825) IQUEST(7),IQUEST(8)
        ENDIF
 9824 FORMAT (1X/' DEVZE FZMEMO,  LOCF(BUF) / LBUF /'
     F,' LOCF(LQ(LBUF)) = ',3I9)
 9825 FORMAT (16X,'NW machine / 32-bit =',2I6)

#endif
      IF (LUNF.EQ.0)               GO TO 91
      IF (IFIFOF.NE.3)             GO TO 92

      IQ(KQSP+LQFF+1) = IQUEST(3)
      IQ(KQSP+LQFF+8) = IQUEST(3)
      IQ(KQSP+LQFF+9) = IQUEST(8)
      GO TO 999

C-------------------------------------------------
C-                 Error handling
C-------------------------------------------------

   91 IQUEST(1) = 1
      CALL ZFATAM ('FZMEMO - Stream not existing.')

   92 IQUEST(1) = 2
      CALL ZFATAM ('FZMEMO - Stream not ready for Memory mode.')

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
