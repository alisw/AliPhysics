*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:48  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABR

C-    Create link reloctation table,
C-    driven by the Memory Occupation table

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LMT,LQMTB)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZTA, 4HBR   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZTABR /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZTABR  ')
#endif

#include "zebra/q_sbit0.inc"
#include "zebra/q_sbit1.inc"


#include "zebra/qtrace.inc"


      IF (LQMTBR.NE.0)             GO TO 81
      LQTA  = LQRTA + 2
      LQTE  = LQTA
      LMT   = LQMTA

C-------           Loop for division in the Mem. Occ. table

      LQ(LQTA-1) = NQLINK + 1

   41 JDIV = LQ(LMT)
      LQ(LMT+5) = LQTE - LQRTA

   42 LQ(LQTE)   = LQ(LMT+3)
      LQ(LQTE+1) = LQ(LMT+4)
      LQ(LQTE+2) = 0
      LQ(LQTE+3) = 0

      IACT = LQ(LMT+1)
      IF (IACT.EQ.3)               GO TO 61
      IF (IACT.EQ.-1)              GO TO 78
      IF (IACT.EQ.4)               GO TO 56
      IF (LQ(LMT+6).EQ.-3)         GO TO 45
      LQTE  = LQTE  + 4
      GO TO 78

   45 LQ(LMT+6) = 0
      LQ(LMT+1) = 3
      MQDVAC = MSBIT1 (MQDVAC,JDIV)
      GO TO 42

C----              Division for wipe

   56 LQ(LMT+7) = LQ(LMT+4) - LQ(LMT+3)
      GO TO 78

C----              Division with garbage collection

   61 IF (IQPART.NE.0)             GO TO 66
      LQTC1  = LQTE
      LQTC2  = LQRTE - (LQMTE-LMT)/2
      IF (LQTC1.GE.LQTC2)          GO TO 65

      LQMTC1 = LQ(LMT+3)
      LQMTC2 = LQ(LMT+4)

      CALL MZTABC
      NQFREE = NQFREE + NQFRTC
      IF (NQLIVE.EQ.0)             GO TO 64
      IF (IQPART.NE.0)  LQMTBR=LMT
      IF (NQFRTC.EQ.0)             GO TO 67

      LQ(LMT+6) = LQTE - LQRTA
      LQ(LMT+7) = NQFRTC
      GO TO 78

C--                Division is all dead

   64 LQTE = LQTC1
      LQ(LMT+1) = 4
      GO TO 42

C--                Division has no dead banks

   65 LQMTBR = LMT
      IQPART = 7
   66 LQ(LMT+6) = -3
      LQTE  = LQTE  + 4
   67 LQ(LMT+1) = 2
      IF (LQ(LMT+2).EQ.0)  THEN
          LQ(LMT+1) = 0
          MQDVAC = MSBIT0 (MQDVAC,JDIV)
        ENDIF

   78 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 41
      JDIV   = LQ(LMT)
      LQ(LQTE) = LQSTA(KQT+JDIV)

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9879) LQMTA,LQMTE,LQRTA,LQTA,LQTE,LQRTE
 9879 FORMAT (' DEVZE MZTABR,  LQMTA,LQMTE,LQRTA,LQTA,LQTE,LQRTE= '
     F,6I8)
#endif

#include "zebra/qtrace99.inc"
      RETURN

C------            Re-start table building

   81 LMT    = LQMTBR
      LQMTBR = 0
      IQPART = 0
      JDIV   = LQ(LMT)
      MQDVAC = MSBIT1 (MQDVAC,JDIV)

#if defined(CERNLIB_QPRINT)
      WRITE (IQLOG, 9882)
 9882 FORMAT (1X/' MZTABR!! !!!!****  re-entry with LQMTBR non-zero',
     F'****!!!!'/1X)
#endif

      JWAY = LQ(LMT+6)
      IF (JWAY.EQ.-3)  THEN
          LQTE = LQ(LMT+5)
          GO TO 45
        ENDIF

      LQTE   = JWAY - 4
      LQMTC1 = LQ(LQTE)
      LQMTC2 = LQ(LQTE+1)
      LQTC1  = LQTE
      LQTC2  = LQRTE - (LQMTE-LMT)/2

      CALL MZTABC
      LQ(LMT+6) = LQTE - LQRTA
      LQ(LMT+7) = LQ(LMT+7) + NQFRTC
      GO TO 78
      END
*      ==================================================
#include "zebra/qcardl.inc"
