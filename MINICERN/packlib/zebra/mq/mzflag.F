*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:32  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZFLAG (IXSTOR,LHEADP,KBITP,CHOPT)

C-    Run through d/s to set status-bit, user called

#include "zebra/zlimit.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzcwk.inc"
C--------------    End CDE                             --------------
      DIMENSION    KBITP(9),LHEADP(9)
      CHARACTER    *(*) CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZFL, 4HAG   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZFLAG /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZFLAG  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit0.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_sbit.inc"


      LHEAD = LHEADP(1)
      IF (LHEAD.EQ.0)        RETURN

#include "zebra/qtrace.inc"
#include "zebra/qstore.inc"

#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
      CALL MZCHLS (-7,LHEAD)
      IF (IQFOUL.NE.0)             GO TO 92
#endif
#if !defined(CERNLIB_QDEBUG)
      IQNS = IQ(KQS+LHEAD-2)
#endif
      LQLIML = LQSTA(KQT+21)
      LQLIMH = 0

      IQTBIT = KBITP(1)
      CALL UOPTC (CHOPT,'ZLV',IQUEST)
      IQTVAL = 1 - IQUEST(1)
      IOPTS  = 1 - IQUEST(3)
      IOPTH  = IQUEST(2)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)
     +WRITE (IQLOG,9814) LHEAD,IQTBIT,IQTVAL,IOPTS,IOPTH
 9814 FORMAT (1X/' DEVZE MZFLAG.   LHEAD,IQTBIT,IQTVAL,IOPTS,IOPTH='
     F,I6,6I4)
#endif

      LEV  = LQWKTB + 3
      LEVE = LEV + NQWKTB - 10
      LQ(LEV-2) = 0
      LQ(LEV-1) = 0
      LQ(LEV)   = LHEAD

      LCUR = LHEAD
      LX   = LHEAD - 1 + IOPTH
      LAST = LHEAD - IQNS
      IQ(KQS+LCUR) = MSBIT1 (IQ(KQS+LCUR),IQSYSX)
      GO TO 24

C--                Mark bank

   20 LAST = LCUR - IQ(KQS+LCUR-2)
      IQ(KQS+LNEW) = MSBIT0 (IQ(KQS+LNEW),IQSYSX)
      IQ(KQS+LNEW) = MSBIT  (IQTVAL,IQ(KQS+LNEW),IQTBIT)
      LQLIML = MIN (LQLIML,LNEW)
      LQLIMH = MAX (LQLIMH,LNEW)

C----              Look at next link

   24 IF (LX.LT.LAST)              GO TO 41
      LNEW = LQ(KQS+LX)
      LX   = LX - 1
      IF (LNEW.EQ.0)               GO TO 24

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LNEW)
      IF (IQFOUL.NE.0)             GO TO 94
#endif
      IF (JBIT(IQ(KQS+LNEW),IQSYSX).NE.0)   GO TO 24

C----              New bank LNEW, push down

      LQ(LEV+1) = LX
      LQ(LEV+2) = LCUR

      LEV = LEV + 3
      IF (LEV.GE.LEVE)             GO TO 91
      LQ(LEV)   = LNEW
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9831) LEV,LCUR,LX+1,LNEW
 9831 FORMAT (' DEVZE MZFLAG,  Down:   LEV,LCUR,LX+1,LNEW=',6I8)
#endif

C--                Move to end of linear structure

   32 LCUR = LNEW
      IQ(KQS+LCUR) = MSBIT1 (IQ(KQS+LCUR),IQSYSX)
      LNEW = LQ(KQS+LCUR)
      IF (LNEW.EQ.0)               GO TO 36
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9833) LCUR,LNEW
 9833 FORMAT (' DEVZE MZFLAG,  Along:  LCUR,LNEW=',6I8)
#endif
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LNEW)
      IF (IQFOUL.NE.0)             GO TO 93
      IF (JBIT(IQ(KQS+LNEW),IQSYSX).NE.0)  GO TO 36
      IF (LQ(KQS+LNEW+2).NE.LCUR)          GO TO 95
      GO TO 32

   36 CONTINUE
#endif
#if !defined(CERNLIB_QDEBUG)
      IF (JBIT(IQ(KQS+LNEW),IQSYSX).EQ.0)  GO TO 32

   36 IQNS = IQ(KQS+LCUR-2)
#endif
      LAST = LCUR - IQNS
      LX   = LCUR - 1
      GO TO 24

C----              Bank at LCUR has no further secondaries
C--                     step back in the linear structure

   41 LNEW = LCUR
      IF (LCUR.EQ.LQ(LEV))         GO TO 46
      LCUR = LQ(KQS+LCUR+2)
      LX   = LCUR - 1
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9841) LCUR,LNEW
 9841 FORMAT (' DEVZE MZFLAG,  Back:   LCUR,LNEW=',6I8)
#endif
      GO TO 20

C--                Start of linear structure reached, pop up

   46 LEV  = LEV - 3
      LX   = LQ(LEV+1)
      LCUR = LQ(LEV+2)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9846) LEV,LCUR,LX
 9846 FORMAT (' DEVZE MZFLAG,  Up:     LEV,LCUR,LX=',6I8)
#endif
      IF (LCUR.NE.0)               GO TO 20

C----              Done, mark header bank

   61 IQ(KQS+LHEAD) = MSBIT0 (IQ(KQS+LHEAD),IQSYSX)
      IF (IOPTS.EQ.0)              GO TO 999
      IQ(KQS+LHEAD) = MSBIT (IQTVAL,IQ(KQS+LHEAD),IQTBIT)
      LQLIML = MIN (LQLIML,LHEAD)
      LQLIMH = MAX (LQLIMH,LHEAD)
#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

#if defined(CERNLIB_QDEBUG)
   95 NQCASE = 2
      NQFATA = 1
      IQUEST(14) = LQ(KQS+LNEW+2)
      GO TO 93

   94 NQCASE = 1
      NQFATA = 1
      IQUEST(14) = LX+1 - LCUR
   93 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 2
      IQUEST(12) = LNEW
      IQUEST(13) = LCUR
   92 NQCASE = NQCASE + 1
#endif
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(11) = LHEAD
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
