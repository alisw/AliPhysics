*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:39  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION LZLONG (IXSTOR,LGO,NTEXT,ITEXT,JWORD)

C-    FIND BANK WITH ITEXT(1-NTEXT) IN WORDS JWORD,..,JWORD+NTEXT-1
C-    SEARCH LINEAR STRUCTURE STARTING WITH BANK AT LGO
C-    RETURN BANK-ADR AS VALUE

#include "zebra/q.inc"
C--------------    END CDE                             --------------
      DIMENSION    ITEXT(9)

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      N2 = NTEXT  - 2
      L  = LGO
      IF (L.EQ.0)                    GO TO 39
      IF (IQ(KQS+L+JWORD).EQ.ITEXT(1))  GO TO 24

   21 L = LQ(KQS+L)
      IF (L.EQ.0)                    GO TO 39
      IF (IQ(KQS+L+JWORD).NE.ITEXT(1))  GO TO 21

   24 IF  (N2)               39, 28, 26

   26 LB = L + JWORD - 1
      DO 27 J=2,NTEXT
      IF (IQ(KQS+LB+J).NE.ITEXT(J))     GO TO 21
   27 CONTINUE
      GO TO 39

   28 IF (IQ(KQS+L+JWORD+1).NE.ITEXT(2))  GO TO 21
   39 LZLONG = L
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
