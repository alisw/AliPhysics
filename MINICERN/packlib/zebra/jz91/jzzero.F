*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZZERO (NLK,NDK)

C-    Clear unused part of call bank

#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZZE, 4HRO   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZZERO /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZZERO  ')
#endif


#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"

      IF (LQDW.EQ.0)               GO TO 91
      IF (NLK.LT.0)                GO TO 92
      IF (NDK.LT.0)                GO TO 93

      DO 24  J=NLK+1,JQCBNL
   24 LQ(KQS+LQDW-J) = 0

      DO 27  J=NDK+1,JQCBND
   27 IQ(KQS+LQDW+J) = 0

#include "zebra/qtrace99.inc"
      RETURN

C----              Trouble

   93 NQCASE = 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = 3
      IQUEST(11) = LQDW
      IQUEST(12) = NLK
      IQUEST(13) = NDK
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
