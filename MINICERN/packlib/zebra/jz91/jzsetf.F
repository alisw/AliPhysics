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
      SUBROUTINE JZSETF (CHPA1,IPA2,IPA3)

C-    Set flag word JFL for processor IAM

#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IPA2(7),IPA3(7)
      CHARACTER    CHPA1*4
#if defined(CERNLIB_A4)
      CHARACTER    CHIAM*4
#endif
#if defined(CERNLIB_A8)
      CHARACTER    CHIAM*8
#endif
#if defined(CERNLIB_EQUHOLCH)
      EQUIVALENCE (CHIAM, IAMID)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZSE, 4HTF   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZSETF /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZSETF  ')
#endif


#include "zebra/qtraceq.inc"

#if defined(CERNLIB_QDEBUG)
#include "zebra/qstorjz.inc"
      CHIAM = CHPA1
#endif
#if (defined(CERNLIB_QDEBUG))&&(!defined(CERNLIB_EQUHOLCH))
      CALL UCTOH (CHIAM, IAMID,4,4)
#endif
#if defined(CERNLIB_QDEBUG)
      JFL   = IPA2(1)
      IVAL  = IPA3(1)

C----              Find SV bank

   21 L = LQ(KQS+LQJZ-4)
      IF (L.EQ.0)                  GO TO 24
      J = IUCOMP (IAMID,IQ(KQS+L+2),IQ(KQS+L+1))
      IF (J.EQ.0)                  GO TO 24
      LSV = LQ(KQS+L-J)
      GO TO 25

   24 LSV   = LZFIND (IXSTJZ,LQ(KQS+LQJZ-3), IAMID,1)
      IF (LSV.EQ.0)                GO TO 41
   25 LFL = LSV + JQNACC
      LFL = LFL + IQ(KQS+LFL) + 1
      LFL = LFL + IQ(KQS+LFL) + 1 + IQ(KQS+LSV+4)
      NFL = IQ(KQS+LFL)

C----              Set flag value

   31 IF (JFL.GT.NFL)              GO TO 81
      IF (JFL.LE.0)                GO TO 91
      IQUEST(3) = IQ(KQS+LFL+JFL)
      IQ(KQS+LFL+JFL) = IVAL
      GO TO 82

C----              Find flag bank if processor not yet init.

   41 LFL = LZFIND (IXSTJZ,LQ(KQS+LQJZ-2), IAMID,1)
      IF (LFL.EQ.0)                GO TO 81
      NFL = IQ(KQS+LFL-1) - 1
      LFL = LFL + 1
      GO TO 31

#endif
C----              Flag not available

   81 LFL = 0
   82 IQUEST(1) = LFL
      IQUEST(2) = NFL
#include "zebra/qtrace99.inc"
      RETURN

C----              Negative flag number

   91 NQCASE = 1
      NQFATA = 3
      IQUEST(11) = NQME(1)
      IQUEST(12) = IAMID
      IQUEST(13) = JFL
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
