*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:44  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZOAPK (LBUFCP,LBUFEP)

C-    Unpack the set of words LQ(LBUFC) to LQ(LBUFE-1) at most
C-    into LQ(LUPKA) to LQ(LUPKE-1) at most;
C-    service routine to FZOASC.

C-    Return LUPKE the adr of the first word after the unpack vector.
C-    For each word unpacked recognise the type and store 8 numbers:
C-
C-           normal number          repetition N+1 times
C-         L + 0  type            L + 0  -43 to signal repeat
C-           + 1  sub-type          + 1  N for N+2 numbers in all
C-           + 2  bits 26-30
C-           + 3  bits 21-25
C-           ...
C-           + 7  bits  1- 5

#include "zebra/zstate.inc"
#include "zebra/mqsysh.inc"
C--------------    End CDE                             --------------
      DIMENSION    LBUFCP(9), LBUFEP(9)

      EQUIVALENCE (LUPKA,IQUEST(93)), (LUPKE, IQUEST(94))


#include "zebra/q_jbyt.inc"


      LBUFC = LBUFCP(1)
      LBUFE = LBUFEP(1)
      LUPK  = LUPKA
      LUPKE = LUPKE - 10
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)  CALL VZERO (LQ(LUPK),LUPKE+8-LUPK)
#endif

   11 IWORD = LQ(LBUFC)
      LBUFC = LBUFC + 1
      JTYPS = 0
      M31   = JBYT (IWORD,31,2)

C--                Short cut if integer 0->9

      IF (M31.NE.0)                GO TO 14
      IF (IWORD.GE.10)             GO TO 14
          JTYPE = 0
          JTYPS = IWORD + 26
          LQ(LUPK+7) = IWORD
          GO TO 38
   14 CONTINUE

C--                Unpack all bytes

      LQ(LUPK+2) = JBYT (IWORD,26,5)
      LQ(LUPK+3) = JBYT (IWORD,21,5)
      LQ(LUPK+4) = JBYT (IWORD,16,5)
      LQ(LUPK+5) = JBYT (IWORD,11,5)
      LQ(LUPK+6) = JBYT (IWORD, 6,5)
      LQ(LUPK+7) = JBYT (IWORD, 1,5)

C--                Type : small +ve integer

      IF (M31.NE.0)                GO TO 21
      IF (LQ(LUPK+2).NE.0)         GO TO 31
      JTYPE = 5
      DO 16  J=3,6
      IF (LQ(LUPK+J).NE.0)         GO TO 38
   16 JTYPS = JTYPS + 1
      GO TO 38

C--                Type : small negative integer

   21 IF (M31.NE.3)                GO TO 31
      IF (LQ(LUPK+2).NE.31)        GO TO 31
      JTYPE = 6
      DO 24  J=3,6
      IF (LQ(LUPK+J).NE.31)        GO TO 38
   24 JTYPS = JTYPS + 1
      GO TO 38

C--                Type : normal words, check short mantissa

   31 JTYPE = M31 + 1
      DO 36  J=7,4,-1
      IF (LQ(LUPK+J).NE.0)         GO TO 38
   36 JTYPS = JTYPS + 1

C--                Store Main type and sub-type, check repeat

   38 LQ(LUPK)   = JTYPE
      LQ(LUPK+1) = JTYPS
      LUPK = LUPK + 8
      IF (LBUFC.EQ.LBUFE)          GO TO 49
      IF (LQ(LBUFC).EQ.IWORD)      GO TO 61
   39 IF (LUPK.LT.LUPKE)           GO TO 11

   49 LQ(LUPK) = -1
      LUPKE = LUPK
      RETURN

C----              Check set of identical words

   61 NL = LBUFE - LBUFC - 1
      NC = 0
      DO 63  J=1,NL
      IF (LQ(LBUFC+J).NE.IWORD)    GO TO 64
   63 NC = NC + 1

   64 IF     (NC.LT.2)  THEN
          IF (JTYPE.EQ.0)          GO TO 39
        ENDIF
      LQ(LUPK)   = -43
      LQ(LUPK+1) = NC
      LUPKE = LUPK + 2
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
