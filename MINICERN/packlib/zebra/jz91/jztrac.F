*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:08  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZTRAC (MODEP)

C-    Dump processor trace for JZ91

#include "zebra/mzcn.inc"
#include "zebra/zbcd.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    MODEP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZTR, 4HAC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZTRAC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZTRAC  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"

C     For bits in MODE - mark critical
C                 bit 1  SV banks
C                     2  Call banks
C                     3  Banks pointed to by links
C                        in current LQDW and LQUP

      IF (LQJZ.EQ.0)               RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"
      MODE = MODEP(1)
      WRITE (IQLOG,9000) JQLEV
      NQUSED = NQUSED + 4

      LEV = MIN (JQLEV,JQMLEV)
      KHI = LQJZ - 6
      KLO = KHI - LEV

      LDW = LQ(KQS+KLO)
      IF (LQDW.EQ.0)  LDW=0

      DO 49 K=KLO,KHI
      LSV = LQ(KQS+K-JQMLEV)

C--                Check LSV valid

      CALL MZCHLS (-7,LSV)
      IF (IQFOUL.EQ.0)             GO TO 23
      ID = IQBLAN
      GO TO 24

   23 ID = IQ(KQS+LSV+1)
      IF (JBIT(MODE,1).EQ.0)       GO TO 24
      IQ(KQS+LSV) = MSBIT1 (IQ(KQS+LSV),IQCRIT)

C----              Print 1 line for this level

   24 WRITE (IQLOG,9024) LEV,ID,LSV,LDW
      NQUSED = NQUSED + 2
      IF (IQFOUL.NE.0)  WRITE (IQLOG,9025)

C----              Check down bank

      IF (LDW.EQ.0)                GO TO 36
      CALL MZCHLS (-7,LDW)
      IF (IQFOUL.EQ.0)             GO TO 31
      WRITE (IQLOG,9026)
      GO TO 36

   31 J = JBYT (MODE,2,2)
      IF (J.EQ.0)                  GO TO 36
      IQ(KQS+LDW) = MSBIT1 (IQ(KQS+LDW),IQCRIT)
      IF (J.LT.2)                  GO TO 36
      IF (LEV.LT.JQLEV-1)          GO TO 36

      DO 34  J=1,IQNL
      L = LQ(KQS+LDW-J)
      IF (L.EQ.0)                  GO TO 34
      CALL MZCHLS (-7,L)
      IF (IQFOUL.NE.0)             GO TO 34
      IQ(KQS+L) = MSBIT1 (IQ(KQS+L),IQCRIT)
   34 CONTINUE
   36 CONTINUE
      LDW = LQ(KQS+K+1)
   49 LEV = LEV - 1
#include "zebra/qtrace99.inc"
      RETURN

 9000 FORMAT (/' JZTRAC.  Processor trace at level',I3
     F/10X,'Depth Name   LQSV   LQDW')
 9024 FORMAT (/10X,I5,1X,A4,2I7)
 9025 FORMAT (' !!!',32X,'***  bank SV bad !!!')
 9026 FORMAT (' !!!',32X,'***  bank DW bad !!!')
      END
*      ==================================================
#include "zebra/qcardl.inc"
