*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:46  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABH

C-    Help : the standard table working space is too small;
C-    shift tables to free gaps of reserve areas in dynamic stores

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZTA, 4HBH   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZTABH /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZTABH  ')
#endif


#include "zebra/qtrace.inc"

      IF (JQGAPM.NE.0)             GO TO 41

C----              On first entry : move Memory Occupation table

      CALL MZFGAP
      NW     = LQMTE+1 - LQMTA
      JQGAPM = NQGAPN
      IF (JQGAPM.LT.2)             GO TO 26
   23 LNEW  = IQGAP(2,JQGAPM)
      NSH   = LNEW - LQMTA
      CALL UCOPY (LQ(LQMTA),LQ(LNEW),NW)

      LQMTA = LNEW
      LQMTB = LQMTB + NSH
      LQMTE = LQMTE + NSH

      LQTC2 = LQTC2 + 161
      LQRTE = LQRTE + 161
      IQGAP(1,JQGAPM) = IQGAP(1,JQGAPM) - NW
      IQGAP(2,JQGAPM) = IQGAP(2,JQGAPM) + NW

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  THEN
        JDIV = IQGAP(3,JQGAPM)
        JSTO = IQGAP(4,JQGAPM)
        WRITE (IQLOG,9824) JQGAPM,JDIV,JSTO,NQOFFS(JSTO+1)
     +,           LQMTA,LQMTB,LQMTE,LQRTA,LQTC1,LQTC2,LQRTE
      ENDIF
 9824 FORMAT (1X/' DEVZE MZTABH,  Gap,Div,Store,Off-set =',3I4,I10
     F/16X,'LQMTA,LQMTB,LQMTE =     ',3I10
     F/16X,'LQRTA,LQTC1,LQTC2,LQRTE=',5I10)
#endif
#include "zebra/qtrace99.inc"
      RETURN

C--                Less than 2 gaps available

   26 IF (IQTNMV.EQ.0)  JQGAPM=NQGAP
      IF (JQGAPM.NE.0)             GO TO 23
      IF (IQTNMV.LT.0)             GO TO 31
   29 IQPART = 7
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  WRITE (IQLOG,9828)
 9828 FORMAT (1X/' DEVZE MZTABH,  no space.')
#endif
      GO TO 999

C----              Move table to forbidden gap, if IQTNMV < 0

   31 JQGAPM = NQGAP
      IF (JQGAPM.EQ.0)             GO TO 29
      IQPART = -7
      GO TO 23

   36 IF (IQTNMV.GE.0)             GO TO 29
      IF (JQGAPR.GT.NQGAPN)        GO TO 29
      JQGAPR = NQGAP
      IF (JQGAPR.EQ.0)             GO TO 29
      IF (IQGAP(1,NQGAP-1).GT.IQGAP(1,NQGAP))  JQGAPR=NQGAP-1
      NNEW = IQGAP(1,JQGAPR) - (LQRTE-LQRTA) - 10
      IF (NNEW.LT.16)              GO TO 29
      IQPART = -7
      GO TO 44

C----              Re-entry : move the Link Relocation table

   41 IF (JQGAPR.NE.0)             GO TO 36
      IF (NQGAPN.EQ.0)             GO TO 36
      NNEW = IQGAP(1,1) - NQWKTB
      IF (NNEW.LT.16)              GO TO 36
      JQGAPR = 1

   44 LNEW = IQGAP(2,JQGAPR)
      NSH  = LNEW - LQRTA
      NW   = LQTE+4 - LQRTA
      CALL UCOPY (LQ(LQRTA),LQ(LNEW),NW)

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  WRITE (IQLOG,9842) LQRTA,LQTE,LQTC2,LQRTE
 9842 FORMAT (1X/' DEVZE MZTABH,  LQRTA,LQTE,LQTC2,LQRTE=',4I10)
#endif
      LQTA   = LQTA  + NSH
      LQTE   = LQTE  + NSH
      LQTC1  = LQTC1 + NSH
      LQTC2  = LNEW + (LQTC2-LQRTA) + NNEW
      LQRTE  = LNEW + (LQRTE-LQRTA) + NNEW
      LQRTA  = LNEW

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +          WRITE (IQLOG,9844) LQRTA,LQTE,LQTC2,LQRTE,NNEW
 9844 FORMAT (16X,'LQRTA,LQTE,LQTC2,LQRTE=',4I10,' NEW=',I10)
#endif
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
