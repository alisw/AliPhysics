*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:26:33  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/04/18 16:12:44  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABC

C-    Construct table of banks with value IQTVAL in status bit IQTBIT
C-
C-    LQMTC1   start adr of first bank
C-    LQMTC2   end+1 adr of  last bank
C-    LQTC1    first table word available
C-    LQTC2    last  table word available

C-    Return       IQPART   if table incomplete
C-                 NQFRTC   number of words occupied by the holes
C-                 NQLIVE   number of live banks

#include "zebra/mzcn.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LS,IQLS), (LNX,IQNX)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZTA, 4HBC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZTABC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZTABC  ')
#endif

#include "zebra/q_jbit.inc"
#if !defined(CERNLIB_QDEBUG)
#include "zebra/q_jbyt.inc"
#endif


#include "zebra/qtrace.inc"

      NQFRTC = 0
      NQLIVE = 0
      N      = 0

      LNX  = LQMTC1
      LQTE = LQTC1
      NEW  = IQTVAL

      LQ(LQTE)   = LNX
      LQ(LQTE+2) = 0

C----              Next bank

   21 MODE = NEW
   22 LN   = LNX
      IF (LN.GE.LQMTC2)            GO TO 41
      N = N + 1
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLN (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 91
#endif
#if !defined(CERNLIB_QDEBUG)
      NST = JBYT (LQ(KQS+LN),1,16) - 11
      IF (NST.GE.0)  THEN
          LS  = LN + NST
          LNX = LS + IQ(KQS+LS-1) + 9
        ELSE
          LS  = LN - 8
          LNX = LN + NST + 11
        ENDIF
#endif
      NEW = JBIT (IQ(KQS+LS),IQTBIT)
      IF (NEW.EQ.MODE)             GO TO 22
      IF (NEW.EQ.IQTVAL)           GO TO 36

C--                Start of a dead group

      NQLIVE = NQLIVE + N - 1
      LQ(LQTE+1) = LN
      LQ(LQTE+3) = 1
      LQTE = LQTE + 4
      GO TO 21

C--                Start of a live group

   36 NQFRTC = NQFRTC + (LN - LQ(LQTE-3))
      LQ(LQTE)   = LN
      LQ(LQTE+2) = 0
      N = 1
      IF (LQTE.LT.LQTC2)           GO TO 21
      CALL MZTABH
      IF (IQPART.EQ.0)             GO TO 21
      IQPART = 1
      LN     = LQMTC2

C----              Close live end group

   41 IF (NEW.NE.IQTVAL)           GO TO 43
      NQLIVE = NQLIVE + N
      LQ(LQTE+1) = LN
      LQ(LQTE+3) = 0
      GO TO 45

C--                Close dead end group

   43 NQFRTC = NQFRTC + (LN-LQ(LQTE-3))
      LQ(LQTE)   = LN
      LQ(LQTE+1) = LN
      LQ(LQTE+2) = 0
      LQ(LQTE+3) = 0

   45 LQTE = LQTE + 4
#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

   91 NQCASE = 1
      NQFATA = 3
      IQUEST(11) = LN
      IQUEST(12) = LQMTC1
      IQUEST(13) = LQMTC2
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
