*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:58  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZCHLN (IXST,LP)

C-    Check LN is a valid start-of-bank address, system called

#include "zebra/zmach.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXST(9), LP(9)

#include "zebra/q_jbyt.inc"


      IXSTOR = IXST(1)
      IQLN   = LP(1)

      IF (IXSTOR.EQ.-7)                 GO TO 21
#include "zebra/qstore.inc"

   21 IF (IQLN.LT.LQSTA(KQT+1))       GO TO 98
      IF (IQLN.GE.LQSTA(KQT+21))      GO TO 98

      IWD = LQ(KQS+IQLN)
      NST = JBYT (IWD,1,16) - 12
      IF (NST.LT.0)                GO TO 41

C--                TRUE BANK

      IQLS  = IQLN + NST + 1
      IF (IQLS.GE.LQSTA(KQT+21))      GO TO 91

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
      IF (IQNS.GT.IQNL)            GO TO 91
      IF (NST.NE.IQNIO+IQNL)       GO TO 91
      IQFOUL = 0
      RETURN

C--                SHORT DEAD REGION

   41 NWD  = JBYT (IWD,17,IQDROP-17)
      IQLS = IQLN - 8
      IQNX = IQLN + NWD
      IQND = -NWD
      IF (NWD.EQ.0)                GO TO 91
      IF (NWD.NE.NST+12)           GO TO 91

      NST = JBYT (IWD,IQDROP,IQBITW-IQDROP)
      IF (NST.NE.1)                GO TO 91
      IQFOUL= 0
      RETURN

   91 IQFOUL = 7
      RETURN
   98 IQFOUL = -7
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
