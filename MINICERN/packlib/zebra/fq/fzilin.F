*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZILIN

C-    Chain banks into one linear structure, clear all links

#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------


      IQFOUL = 0
      LENTRI = 0

      K   = 0
      LMT = LQMTA
   22 IF (LQ(LMT+1).LE.0)          GO TO 29
      IQNX  = LQ(LMT+3)
      LEND  = LQ(LMT+4)

   24 CALL MZCHLN (-7,IQNX)
      IF (IQFOUL.NE.0)             GO TO 91
      IF (IQND.LT.0)               GO TO 27
      IF (K.EQ.0)  THEN
          LENTRI = IQLS
        ELSE
          LQ(KQS+K) = IQLS
        ENDIF

      L = IQLS - IQNL - 1
      DO 26  J=1, IQNL+2
   26 LQ(KQS+L+J) = 0
      LQ(KQS+IQLS+2) = K
      K = IQLS
   27 IF (IQNX.LT.LEND)            GO TO 24
      IF (IQNX.NE.LEND)            GO TO 91

   29 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 22
      IQUEST(1) = LENTRI
      RETURN

   91 IQFOUL = 7
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
