*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:04  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:16  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZEND

C-    PRINT PROCESSOR SUMMARY

#include "zebra/mzcn.inc"
#include "zebra/zbcd.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IDV(5), IDSHV(2), TIMEV(3)
      EQUIVALENCE (IDV(1),  IQUEST(1))
      EQUIVALENCE (IDSHV(1),IQUEST(1))
      EQUIVALENCE (TIMEV(1), IQUEST(6))

      DIMENSION    NUMS(10)
      EQUIVALENCE (NUMS(1),IQNUM2(2))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZEN, 4HD    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZEND  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZEND   ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      DATA  MYES   / 4H T   /

      IF (LQJZ.EQ.0)               RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"

      CALL MZCHLS (-7,LQJZ)
      IF (IQFOUL.EQ.0)             GO TO 12
      WRITE (IQLOG,9000) LQJZ
#include "zebra/qtrace99.inc"
      RETURN

   12 CALL JZROOT

      WRITE (IQLOG,9001)
      WRITE (IQLOG,9002) NUMS

      IFLAG = 0
   21 LSV = LQJZ - 3
   22 LSV = LQ(KQS+LSV)
      IF (LSV.EQ.0)  THEN
          IF (JQTIME.EQ.0)         GO TO 999
          IF (IFLAG .NE.0)         GO TO 999
          WRITE (IQLOG,9003)
          IFLAG = 7
          GO TO 21
        ENDIF

      LCR = LSV + JQNACC
      NCR = IULAST (0,IQ(KQS+LCR+1),IQ(KQS+LCR)) + 1
      LLC = LCR + IQ(KQS+LCR)

      IQUEST(1) = IQ(KQS+LLC+1)
      IQUEST(2) = IQBLAN
      IF (JBIT(IQ(KQS+LSV),17).NE.0)  IQUEST(2)=MYES
      IQUEST(3) = JBYT(IQ(KQS+LSV),1,8)
      IQUEST(4) = IQ(KQS+LSV+1)
      IQUEST(5) = IQ(KQS+LSV+2)
      IF (IFLAG .NE.0)             GO TO 24
      WRITE (IQLOG,9008) IDV,(IQ(KQS+J+LCR-1),J=1,NCR)
      GO TO 22

   24 IQUEST(1) = IQUEST(4)
      IQUEST(2) = IQUEST(5)
      IF (IQUEST(2).LE.0)          GO TO 22
      Q(KQS+LSV+8) = Q(KQS+LSV+8) + Q(KQS+LSV+5)
      Q(KQS+LSV+7) = Q(KQS+LSV+7) + Q(KQS+LSV+8)
      Q(KQS+LSV+5) = 0.
      Q(KQS+LSV+8) = 0.
      TIMEV(1) = Q(KQS+LSV+7)
      TIMEV(2) = Q(KQS+LSV+6)
      TIMEV(3) = TIMEV(1)/FLOAT(IQUEST(2))
      WRITE (IQLOG,9009) IDSHV,TIMEV
      GO TO 22

 9000 FORMAT (/' JZEND.   Bad address LQJZ = ',I6)
 9001 FORMAT (/' JZEND.   Summary of processor usage'/1X)
 9002 FORMAT ('   NAN T LV Name  Calls   Conditions'
     F/7X,1HT,20X,1HN,10(6X,A1))
 9003 FORMAT (/12X,'Name  Calls    Total  Maximum     Mean time')
 9008 FORMAT (/1X,I5,A2,I3,1X,A4,I7,I6,10I7,:/(29X,10I7))
 9009 FORMAT (/12X,A4,I7,F9.1,2F9.4)
      END
*      ==================================================
#include "zebra/qcardl.inc"
