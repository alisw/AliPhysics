      FUNCTION LZFVAL (IXSTOR,LGO,VAL,TOL,JWORD)

C-    FIND BANK CONTAINING IN WORD JWORD FL.P. VAL WITHIN TOL
C-    SEARCH LINEAR STRUCTURE STARTING WITH BANK AT LGO
C-    RETURN BANK-ADR AS VALUE

#include "zebra/q.inc"
C--------------    END CDE                             --------------
      EQUIVALENCE (IA,AA)

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      VLOW = VAL - TOL
      VHI  = VAL + TOL

      L = LGO
   21 IF (L.EQ.0)                  GO TO 29
      IA = IQ(KQS+L+JWORD)
      IF (IA.EQ.0)                 GO TO 22
      JEXP = JBYT (IA,8,8)
      IF (JEXP.EQ.0)               GO TO 27
      IF (JEXP.EQ.255)             GO TO 27

   22 IF (AA.LT.VLOW)              GO TO 27
      IF (AA.LE.VHI)               GO TO 29
   27 L = LQ(KQS+L)
      GO TO 21

   29 LZFVAL= L
      RETURN
      END
