*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZCHANNEL)
      SUBROUTINE FZHOOK (LUNP,TARGET,BUF)

C-    Connect user routine TARGET for Channel mode, User called

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcf.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), BUF(9)
      EXTERNAL     TARGET

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZHO, 4HOK   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZHOOK /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZHOOK  ')
#endif

*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      LUN   = LUNP(1)

      CALL FZLOC (LUN,0)
      IF (LUNF.EQ.0)    LOGLVF=1
      IF (IACMOF.NE.3)  LOGLVF=1

#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.1)   WRITE (IQLOG,9024) LUN
 9024 FORMAT (1X/' FZHOOK.  Connect user routine for LUN=',I3)

#endif
      IF (LUNF.EQ.0)               GO TO 91
      IF (IACMOF.NE.3)             GO TO 92

      IQUEST(2) = JUMPAD (TARGET)
      IQ(KQSP+LQFF+1) = IQUEST(2)
      IF (LUNF.EQ.LUNI)  IADOPI = IQUEST(2)
      IF (LUNF.EQ.LUNX)  IADOPX = IQUEST(2)
      GO TO 999

C-------------------------------------------------
C-                 Error handling
C-------------------------------------------------

   91 IQUEST(1) = 1
      CALL ZFATAM ('FZHOOK - Stream not existing.')

   92 IQUEST(1) = 2
      CALL ZFATAM ('FZHOOK - Stream not ready for Channel mode.')

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
