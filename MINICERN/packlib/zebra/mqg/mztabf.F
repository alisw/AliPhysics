*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABF

C-    Finalize Link Relocation table

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZTA, 4HBF   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZTABF /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZTABF  ')
#endif

#include "zebra/q_jbit.inc"


#include "zebra/qtrace.inc"

      LMT    = LQMTA
      NCOLL  = 0
      NGARB  = 0
      NQNOOP = 0
      LFIXLO = NQLINK + 1

C----              Truncate initial inactive division entries

   21 JDIV = LQ(LMT)
      IACT = LQ(LMT+1)
      IF (IACT.EQ.4)               GO TO 26
      IF (IACT.GE.2)               GO TO 28

      IF (IACT.GE.0)  LFIXLO=LQEND(KQT+JDIV)

      LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 21
      NQNOOP = -7
      IF (NQDVMV.EQ.0)             GO TO 81
      NQNOOP = 7
      GO TO  81

C--                First active division is wipe

   26 IF (LQ(LMT+9).NE.4)          GO TO 28
      IF (LMT+8.GE.LQMTE)          GO TO 28
      LMT = LMT + 8
      GO TO 26

C--                First active division has garb. coll. or shift

   28 LFIR = LMT
      LQTA = LQRTA + LQ(LMT+5)

C--------          Loop remaining divisions

   31 IACT = LQ(LMT+1)
      NSHF = LQ(LMT+2)
      LTU  = LQRTA + LQ(LMT+5)
      IF (IACT.EQ.4)               GO TO 71
      IF (IACT.EQ.3)               GO TO 61
      IF (IACT.EQ.2)               GO TO 41
      IF (IACT.LT.0)               GO TO 79
      NCOLL = NCOLL + 1
      IF (NCOLL.NE.1)              GO TO 79
      LCOLL = LMT
      GO TO 79

C------            Inactive division with shift

   41 IF (NCOLL+NGARB.LT.2)        GO TO 49

C----              Squeeze out useless entries

   43 LCOLE = LMT - 8
      LT    = LQ(LCOLE+5)
      LTF   = LQ(LCOLL+5)
      N     = LT - LTF

C--                Shift active start of the table

      NW = LQRTA + LTF+1 - LQTA
      CALL UCOPY2 (LQ(LQTA),LQ(LQTA+N),NW)
      LQTA  = LQTA + N
      NCOLL = 0
      IF (IACT.EQ.4)               GO TO 71
      IF (IACT.EQ.3)               GO TO 61

   49 LQ(LTU+2) = NSHF
      GO TO 77

C----              Garbage collection division

   61 IF (NCOLL+NGARB.GE.2)        GO TO 43
      JDIV = LQ(LMT)
      LT   = LTU
      LTE  = LQRTA + LQ(LMT+6)
      MODE = JBIT (IQMODE(KQT+JDIV),1)
      IF (MODE.NE.0)               GO TO 65

C--                Calculate relocation constants for forward mode

      NCUM = NSHF
      GO TO 66

C--                Calculate relocation constants for reverse mode

   65 NCUM = LQ(LMT+7) + NSHF
   66 LQ(LT+2) = NCUM
      NCUM = NCUM - (LQ(LT+4)-LQ(LT+1))
      LT   = LT + 4
      IF (LT.LT.LTE)               GO TO 66
      NGARB = -64
      GO TO 77

C--                Wipe division

   71 IF (NCOLL+NGARB.GE.2)        GO TO 43
   77 NCOLL = 0

C--                End of loop, set table termination

   79 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 31
      LQ(LQTA-1) = LFIXLO

      IF (NCOLL.EQ.0)              GO TO 81
      LQTE = LQRTA + LQ(LCOLL+5)

   81 CONTINUE
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.11)            GO TO 999
      WRITE (IQLOG,9881) LQMTA,LQMTE
 9881 FORMAT (1X/' DEVZE MZTABF,  Memory Occupation Table,'
     F,' LQMTA,LQMTE= ',2I8/16X,
     F' DIV ACT     NWS    LFBK  LLBK+1     JFE     JLE    FREE')

      WRITE (IQLOG,9882) (LQ(J),J=LQMTA,LQMTE-1)
 9882 FORMAT (16X,2I4,6I8)

      WRITE (IQLOG,9883) LQRTA,LQTA,LQTE,LFIXLO
 9883 FORMAT (1H0,15X,'Link Relocation Table,'
     F,' LQRTA,LQTA,LQTE= ',3I8
     F/16X,'   LOC       L      LE    NREL  BG'
     F/30X,I8)

      IF (LQTE.LE.LQTA)            GO TO 85
      I  = LQRTA
      JA = LQTA - I
      JE = LQTE - I - 1
      WRITE (IQLOG,9884) (J,LQ(I+J),LQ(I+J+1),LQ(I+J+2),LQ(I+J+3),
     +                                           J=JA,JE,4)
 9884 FORMAT (16X,I6,3I8,I4)

   85 WRITE (IQLOG,9885) LQ(LQTE)
 9885 FORMAT (22X,I8)
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
