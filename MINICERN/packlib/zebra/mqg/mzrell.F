*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:38  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZRELL (MDESV)

C-    Relocator for links in link areas

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    MDESV(99)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZRE, 4HLL   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZRELL /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZRELL  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEVZE)
      LABS = LOCF (MDESV)
      LZEB = LABS - LQASTO
      LPRI = LABS - LQPSTO
      N    = MIN(MDESV(1),11)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9812) LABS,LZEB,LPRI,(MDESV(J),J=1,N)
 9812 FORMAT (1X/' DEVZE MZRELL.  Table of Link Areas, LABS,LZEB,LPRI='
     F,3I10/16X,'First 2 entries,  First word=',I10,
#endif
#if (defined(CERNLIB_QDEVZE))&&(!defined(CERNLIB_HEX))
     F2(/16X,2I10,2X,O11,2X,2A4))
#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_HEX))
     F2(/16X,2I10,2X,Z8,2X,2A4))
#endif

      LFIXLO = LQ(LQTA-1)
      LFIXRE = LQ(LQTA)
      LFIXHI = LQ(LQTE)
      JHIGO  = (LQTE-LQTA) / 4
      NENTR  = JHIGO - 1

      IF (NENTR.EQ.0)  THEN
          LADTB1 = LQ(LQTA+1)
          NRLTB2 = LQ(LQTA+2)
          IFLTB3 = LQ(LQTA+3)
        ENDIF

      JDESMX = MDESV(1) - 4
      JDES   = -4
      IF (MDESV(2).GE.MDESV(3))   JDES =1

C------            Next link area,  check if dead group

   17 JDES = JDES + 5
      IF (JDES.GE.JDESMX)          GO TO 999
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9818) JDES,(MDESV(J+JDES),J=1,5)
 9818 FORMAT (1X/' DEVZE MZRELL.  Do area at',I5,' with',
#endif
#if (defined(CERNLIB_QDEVZE))&&(!defined(CERNLIB_HEX))
     F2I10,2X,O11,2X,2A4)
#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_HEX))
     F2I10,2X,Z8,2X,2A4)
#endif
      LOCAR  = MDESV(JDES+1)
      LIX    = LOCAR
      LOCARE = MDESV(JDES+2)
      MODAR  = MDESV(JDES+3)
      IF (JBIT(MODAR,31).NE.0)  THEN
          IF (LQ(KQS+LOCAR).EQ.0)     GO TO 17
          LIX = LIX + 2
        ENDIF

C--                Next link area,  alive

      LIR = LOCAR + JBYT (MODAR,1,15)
      IF   (NENTR)           66, 46, 26

C--------------    2 OR MORE RELOCATION INTERVALS       -------------

C----              Next link

   24 LQ(KQS+LIX)= 0

   25 LIX = LIX + 1
      IF (LIX.EQ.LOCARE)           GO TO 17
   26 LFIRST= LQ(KQS+LIX)
   27 LINK  = LQ(KQS+LIX)
      IF (LINK.EQ.0)               GO TO 25
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9827) LINK,LIX
 9827 FORMAT (16X,'Link =',I7,' from LIX =',I7)
#endif
      IF (LINK.LT.LFIXLO)          GO TO 25
      IF (LINK.GE.LFIXHI)          GO TO 25
      IF (LINK.LT.LFIXRE)          GO TO 24

C--                Binary search in relocator table

      JLOW = 0
      JHI  = JHIGO

   29 JEX = (JHI+JLOW) / 2
      IF (JEX.EQ.JLOW)             GO TO 31
      IF (LINK.GE.LQ(LQTA+4*JEX))  GO TO 30
      JHI  = JEX
      GO TO 29

   30 JLOW = JEX
      GO TO 29

C--                Relocate

   31 JTB = LQTA + 4*JLOW
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9831) JLOW, (LQ(JTB+J-1),J=1,4)
 9831 FORMAT (50X,'Entry',I5,',',4I7)
#endif
      IF (LINK.GE.LQ(JTB+1))             GO TO 33
      LQ(KQS+LIX) = LINK + LQ(JTB+2)
      GO TO 25

C----              Link into dead area

   33 IF (LIX.GE.LIR)              GO TO 24

C--                Bridge structural link

      IF (LQ(JTB+3).LE.0)          GO TO 24

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.NE.0)             GO TO 91
#endif
      LINK = LQ(KQS+LINK)
      LQ(KQS+LIX) = LINK
      IF (LINK.NE.LFIRST)          GO TO 27
      GO TO 24

C--------------    1 RELOCATION INTERVAL ONLY           -------------

C----              Next link

   44 LQ(KQS+LIX)= 0

   45 LIX = LIX + 1
      IF (LIX.EQ.LOCARE)           GO TO 17
   46 LFIRST= LQ(KQS+LIX)
   47 LINK  = LQ(KQS+LIX)
      IF (LINK.EQ.0)               GO TO 45
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9827) LINK,LIX
#endif
      IF (LINK.LT.LFIXLO)          GO TO 45
      IF (LINK.GE.LFIXHI)          GO TO 45
      IF (LINK.LT.LFIXRE)          GO TO 44

C--                Relocate

      IF (LINK.GE.LADTB1)          GO TO 53
      LQ(KQS+LIX) = LINK + NRLTB2
      GO TO 45

C----              Link into dead area

   53 IF (LIX.GE.LIR)              GO TO 44

C--                Bridge structural link

      IF (IFLTB3.LE.0)             GO TO 44

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.NE.0)             GO TO 91
#endif
      LINK = LQ(KQS+LINK)
      LQ(KQS+LIX) = LINK
      IF (LINK.NE.LFIRST)          GO TO 47
      GO TO 44

C--------------    NO RELOCATION INTERVAL               -------------

C----              Next link

   64 LQ(KQS+LIX)= 0

   65 LIX = LIX + 1
      IF (LIX.EQ.LOCARE)           GO TO 17
   66 LINK  = LQ(KQS+LIX)
      IF (LINK.EQ.0)               GO TO 65
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9827) LINK,LIX
#endif
      IF (LINK.LT.LFIXLO)          GO TO 65
      IF (LINK.GE.LFIXHI)          GO TO 65
      GO TO 64

C------            Error conditions

#if defined(CERNLIB_QDEBUG)
   91 NQCASE = 1
      NQFATA = 5
      IQUEST(11) = LOCAR + LQSTOR
      IQUEST(12) = LIX - LOCAR + 1
      IQUEST(13) = LINK
      IQUEST(14) = MDESV(JDES+4)
      IQUEST(15) = MDESV(JDES+5)
#include "zebra/qtofatal.inc"
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
