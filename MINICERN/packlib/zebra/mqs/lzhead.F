*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:57  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION LZHEAD (IXSTOR, LGO)

C-    Find start of the linear structure containing bank at LGO
C-    Returns : LZHEAD    L of the first bank,
C-                        = zero if trouble
C-              IQUEST(1) JBIAS : -ve for normal dependent
C-                                  1 for top structure
C-                                  2 for stand-alone
C-              IQUEST(2) adr of the supporting link if JBIAS=1
C-              IQUEST(3) error indication if trouble

#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      LH = LGO
      IQUEST(1) = 0
      IQUEST(2) = 0
      IQUEST(3) = 0
      IF (LH.EQ.0)                 GO TO 91
      LIML = LQSTA(KQT+1)
      LIMU = LQSTA(KQT+21)
      LUP  = LQ(KQS+LH+1)
      IF (LUP.EQ.0)                GO TO 41

C----              L/s normal dependent : LUP .NE. 0

      NS  = IQ(KQS+LUP-2)
      LOW = LUP - NS
   22 LNX = LQ(KQS+LH+2)
      IF (LNX.LT.LIML)             GO TO 92
      IF (LNX.GE.LIMU)             GO TO 92
      IF (LNX.LT.LOW.OR.LNX.GE.LUP) THEN
          LH = LNX
          GO TO 22
        ENDIF

      JBIAS = LNX - LUP
      IF (LQ(KQS+LUP+JBIAS).NE.LH)    GO TO 93

      IQUEST(1) = JBIAS
      LZHEAD    = LH
      GO TO 999

C----              L/s is a top structure : LUP = 0

   41 LNX = LQ(KQS+LH+2)
      IF (LNX.EQ.0)                GO TO 43
      IF (LNX.GE.LIML .AND. LNX.LT.LIMU)  THEN
          LH = LNX
          GO TO 41
        ENDIF

      IQUEST(1) = 1
      IF (LQ(KQS+LNX).NE.LH)          GO TO 94
      GO TO 44

   43 IQUEST(1) = 2
   44 IQUEST(2) = LNX
      LZHEAD    = LH
      GO TO 999

C----              Errors

   94 IQUEST(3) = 1
   93 IQUEST(3) = IQUEST(3) + 1
   92 IQUEST(3) = IQUEST(3) + 1
   91 IQUEST(3) = IQUEST(3) + 1
      LZHEAD = 0
  999 RETURN
      END
