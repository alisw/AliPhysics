*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZMOVE

C-    Garbage collector - memory move

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZMO, 4HVE   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZMOVE /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZMOVE  ')
#endif

#include "zebra/q_jbit.inc"


#include "zebra/qtrace.inc"

C---------------   UPDATE DIVISION TABLES
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9821) LQMTA,LQMTE
 9821 FORMAT (1X/' DEVZE MZMOVE,  Update division tables, LQMTA/LQMTE='
     F,2I12)
#endif

      LMT  = LQMTA

   23 JDIV = LQ(LMT)
      IACT = LQ(LMT+1)
      NSHF = LQ(LMT+2)
      IF (IACT.EQ.4)               GO TO 26
      IF (IACT.NE.3)               GO TO 31

C----              Division with garbage collection
C--                New start adr

      L  = LQ(LMT+3)
      LT = LQ(LMT+5) + LQRTA
      N  = LQ(LT+2)
      LQSTA(KQT+JDIV) = L + N

C--                New end adr

      L  = LQ(LMT+4)
      LT = LQ(LMT+6) + LQRTA - 4
      N  = LQ(LT+2)
      LQEND(KQT+JDIV) = L + N
      GO TO 36

C----              Wiped division

   26 MODE = JBIT (IQMODE(KQT+JDIV),1)
      IF (MODE.NE.0)               GO TO 28

C--                Division mode forward

      LQSTA(KQT+JDIV) = LQSTA(KQT+JDIV) + NSHF
      GO TO 29

C--                Division mode reverse

   28 LQSTA(KQT+JDIV) = LQEND(KQT+JDIV) + NSHF
   29 LQEND(KQT+JDIV) = LQSTA(KQT+JDIV)
      GO TO 36

C----              Inactive division

   31 IF (NSHF.EQ.0)               GO TO 37
      LQSTA(KQT+JDIV) = LQSTA(KQT+JDIV) + NSHF
      LQEND(KQT+JDIV) = LQEND(KQT+JDIV) + NSHF
   36 CONTINUE
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9836) JDIV,LQSTA(KQT+JDIV),LQEND(KQT+JDIV)
 9836 FORMAT (16X,'JDIV / LQSTA / LQEND =',I4,2I8)
#endif
   37 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 23
      IF (NQNOOP.NE.0)             GO TO 999

C---------------   MEMORY MOVE
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9861) LQTA,LQTE
 9861 FORMAT (1X/' DEVZE MZMOVE,  Execute memory move, LQTA/LQTE=',2I12)
#endif
      IF (LQTE.LE.LQTA)            GO TO 999

      LTF = LQTA

   61 NREL = LQ(LTF+2)
      IF (NREL)              64, 68, 71

C--                Left squeeze

   64 LOLD = LQ(LTF)
      LNEW = LOLD + NREL
      NW   = LQ(LTF+1) - LOLD
      IF (NW.EQ.0)                 GO TO 68
      CALL UCOPY (LQ(KQS+LOLD),LQ(KQS+LNEW),NW)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9865) LOLD, LNEW, NW
 9865 FORMAT (16X,'Copy  From / To / N=',3I8)
#endif
   68 LTF = LTF + 4
      IF (LTF.NE.LQTE)             GO TO 61
      GO TO  999

C----              Right squeeze

   71 LTFN = LTF
   72 LTFN = LTFN + 4
      IF (LTFN.EQ.LQTE)            GO TO 76
      IF (LQ(LTFN+2).GT.0)         GO TO 72
   76 LTR = LTFN

C--                Execute in reverse

   81 LTR  = LTR - 4
      LOLD = LQ(LTR)
      NW   = LQ(LTR+1) - LOLD
      IF (NW.EQ.0)                 GO TO 88
      LNEW = LOLD + LQ(LTR+2)
      CALL UCOPY2 (LQ(KQS+LOLD),LQ(KQS+LNEW),NW)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9887) LOLD, LNEW, NW
 9887 FORMAT (16X,'Copy2 From / To / N=',3I8)
#endif
   88 IF (LTR.NE.LTF)              GO TO 81
      LTF = LTFN
      IF (LTF.NE.LQTE)             GO TO 61
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
