*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:59  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZCHLS (IXST,LP)

C-    Check LS is a valid bank status-word address, can be user called

#include "zebra/zmach.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXST(9), LP(9)

#include "zebra/q_jbyt.inc"


      IXSTOR = IXST(1)
      IQLS   = LP(1)

      IF (IXSTOR.EQ.-7)                 GO TO 21
#include "zebra/qstore.inc"

   21 IF (IQLS.LT.LQSTA(KQT+1))       GO TO 98
      IF (IQLS.GE.LQSTA(KQT+21))      GO TO 98

      IQNIO = JBYT (IQ(KQS+IQLS),19,4)
      IQID  = IQ(KQS+IQLS-4)
      IQNL  = IQ(KQS+IQLS-3)
      IQNS  = IQ(KQS+IQLS-2)
      IQND  = IQ(KQS+IQLS-1)
      IF (  JBYT(IQNL,IQBITW-3,4)
     +    + JBYT(IQNS,IQBITW-3,4)
     +    + JBYT(IQND,IQBITW-3,4) .NE.0)    GO TO 91
      IQNX  = IQLS + IQND + 9
      IF (IQNX.GT.LQSTA(KQT+21))      GO TO 91

      IQLN  = IQLS - IQNL - IQNIO - 1
      IF (IQLN.LT.LQSTA(KQT+1))       GO TO 91

      NST = JBYT (LQ(KQS+IQLN),1,16) - 12
      IF (NST.NE.IQNIO+IQNL)       GO TO 91
      IF (IQNS.GT.IQNL)            GO TO 91
      IQFOUL = 0
      RETURN

   91 IQFOUL = 7
      RETURN
   98 IQFOUL = -7
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
