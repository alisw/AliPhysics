*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:47  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABM

C-    Construct Memory Occupation table with one entry for
C-    each division, setting actions according to :
C-       MQDVGA    division bits for garbage collection
C-       MQDVWI                  for wipe
C-
C-    Return division bits for all active divisions in MQDVAC

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZTA, 4HBM   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZTABM /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZTABM  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_sbit1.inc"


#include "zebra/qtrace.inc"

      LQMTBR = 0
      IQTBIT = IQDROP
      IQTVAL = 0
      NQFREE = 0
      IQPART = 0
      MQDVAC = 0
      IQFLIO = 0

      IF (JQSTMV.LT.0)  THEN
          IQTNMV = 0
          JQDVM1 = 0
          JQDVM2 = 0
          NQDVMV = 0
        ENDIF

      JQGAPM = 0
      JQGAPR = 0
      LQMTE  = LQWKTB + NQWKTB - 1
      LQMTA  = LQMTE - 160
      LQMTB  = LQMTA
      LQRTE  = LQMTA - 10
      LQTC2  = LQRTE
      LQRTA  = LQWKTB
      LQTA   = LQRTA + 1
      LQTE   = LQTA
      LQTC1  = LQTA

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  WRITE (IQLOG,9828) LQRTA,LQRTE,LQMTA,LQMTE
 9828 FORMAT (1X/' DEVZE MZTABM,  LQRTA,LQRTE,LQMTA,LQMTE=',4I9)
#endif

C----              Loop all divisions

      LMT  = LQMTA
      JDIV = 1

   32 LQ(LMT)   = JDIV
      LQ(LMT+1) = 0
      LQ(LMT+2) = 0
      LQ(LMT+3) = LQSTA(KQT+JDIV)
      LQ(LMT+4) = LQEND(KQT+JDIV)
      LQ(LMT+5) = 0
      LQ(LMT+6) = 0
      LQ(LMT+7) = 0
      NW = LQ(LMT+4) - LQ(LMT+3)
      IF (NW.EQ.0)                  GO TO 37
      NQDSIZ(KQT+JDIV) = MAX (NQDSIZ(KQT+JDIV),NW)
      IF (JBIT(MQDVWI,JDIV).NE.0)   GO TO 41
      IF (JBIT(MQDVGA,JDIV).NE.0)   GO TO 44
      GO TO 48

C--                Division is empty

   37 LQ(LMT+1) = -1
      GO TO 48

C--                Division for wipe

   41 IF (JDIV.EQ.JQDVSY)          GO TO 48
      LQ(LMT+1) = 4
      GO TO 45

C--                Division with garbage collection

   44 LQ(LMT+1) = 3
   45 MQDVAC = MSBIT1 (MQDVAC,JDIV)

   48 LMT  = LMT  + 8
      JDIV = JDIV + 1
      IF (JDIV.EQ.JQDVLL+1)  JDIV=JQDVSY
      IF (JDIV.LT.21)              GO TO 32
      LQMTE  = LMT
      LQMTLU = LMT
      LQ(LQMTE) = 21
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.11)            GO TO 999
      WRITE (IQLOG,9849) LQMTE,MQDVGA,MQDVWI,MQDVAC
      WRITE (IQLOG,9850)
      WRITE (IQLOG,9851) (LQ(J),J=LQMTA,LQMTE-1)
 9850 FORMAT (16X,' DIV ACT NWS      LF    LE+1')
 9851 FORMAT (16X,3I4,2I8,3I4)
#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_HEX))
 9849 FORMAT (1X/' DEVZE MZTABM,  LQMTE,MQDVGA,MQDVWI,MQDVAC='
     F,5X,I8,3(2X,Z8))
#endif
#if (defined(CERNLIB_QDEVZE))&&(!defined(CERNLIB_HEX))
 9849 FORMAT (1X/' DEVZE MZTABM,  LQMTE,MQDVGA,MQDVWI,MQDVAC='
     F,5X,I8,3(2X,O11))
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
