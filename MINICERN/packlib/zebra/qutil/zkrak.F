*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:42  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZKRAK (MFLAGP,JLP,JRP,VALP)

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zkrakc.inc"
#include "zebra/zkrakqu.inc"
C--------------    END CDE                             --------------

      COMMON /SLATE/ NDIG, JNEXT, DUMMY(2), NUMOCT(36)

      DIMENSION    MFLAGP(9), JLP(9), JRP(9), VALP(9)
      INTEGER      IVAL
      EQUIVALENCE (VAL,IVAL)


#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbyt.inc"


      CALL ZKRAKN (MFLAGP,JLP,JRP,VALP)
      IF   (MFMT+1)          21, 40, 99


C-    MODE =    3  BCD
C-              2  NUMERIC
C-              1  EMPTY FIELD
C-              0  FAULT
C-             -1  SEPARATOR IS FIRST CHAR. IN FIELD
C-             -2  $A, $Q CONTROL ITEM
C-             -3  $. COMMENT $

C-    MFMT=     3  HEX
C-              2  OCTAL
C-              1  INTEGER
C-              0  FLOATING
C-             -1  BCD IN A-FORMAT
C-             -2  BCD IN Q-FORMAT

C---  SIGNIFICANCE OF BITS IN MFLAG

C-    BIT   1 - 9  LEGAL SEPARATORS
C-             10  MULTIPLIER LEGAL
C-             11  NON-DELIMITED BCD STRING LEGAL
C-             12  $A, $Q CONTROL ITEMS LEGAL
C-             15  Q-FORMAT (ELSE A-FORMAT)
C-          16-20  'N' FOR AN/QN


C----              NON-DELIMITED BCD STRING

   21 IF (JBIT(MFLAG,11).EQ.0)     GO TO 94
      JR = JDO
   22 JR = JR + 1
      IF (JR.GT.JRRAN)             GO TO 23
      NSEP = IQCETK(JR)
      IF (NSEP.LT.40)              GO TO 22
   23 NCH = JR - JDO
      JGOTO = 1
      GO TO 41

C----              $H'----',  DELIMITED BCD STRING

   25 IT = IQCETK(JDO)
      JR = JDO
   26 JR = JR + 1
      IF (JR.GT.JRRAN)             GO TO 94
      IF (IQCETK(JR).NE.IT)        GO TO 26

      JDO = JDO + 1
      NCH = JR - JDO
      IF (NCH.EQ.0)                GO TO 94
      JGOTO = 1
      GO TO 69

   28 MODE = 3
      MFMT = -1 - JBIT(MFLAG,15)
      J    = JBYT (MFLAG,16,5)
      IF (J.EQ.0)      J=31
      IF (MFMT.EQ.-2)  J=MIN (J,4)
      CALL UTRANS (IQHOLK(JDO),VALP(1),NCH, 1,J)
      NWORDS = JNEXT
      IF (MFMT.NE.-2)        RETURN
      CALL ZHTOI (VALP(1),VALP(1),NWORDS)
   99 RETURN

C----              SEPARATOR IS THE FIRST CHARACTER

   40 IF (NSEP.EQ.43)              GO TO 51
      MODE  = -1
      JR    = JDO
      JGOTO = 2

C----              LOCALISE + VALIDATE TERMINATOR NSEP IN IQHOLK(JR)

   41 IF (JR.GT.JRRAN)             GO TO 48
      IF (NSEP.NE.45)              GO TO 44
      IF (IFLBLT.NE.0)             GO TO 45

C--                IF NSEP IS NON-TERMINATING BLANK, STEP TO TERMINATOR

   43 NSEP = IQCETK(JR)
      IF (NSEP.NE.45)              GO TO 44
      JR = JR + 1
      IF (JR.GT.JRRAN)             GO TO 48
      GO TO 43

C--                CHECK VALID NSEP

   44 NSEP = MIN  (NSEP,47)
      IF (NSEP.LT.40)              GO TO 94
      IF (JBIT(MFLAG,NSEP-39).EQ.0)  GO TO 94
      IF (NSEP.EQ.43)  JR=JR-1

   45 JNXGO = JR
   46 JNXGO = JNXGO + 1
      IF (JNXGO.GT.JRRAN)          GO TO 49
      IF (IQCETK(JNXGO).EQ.45)     GO TO 46
      GO TO 49

   48 NSEP = 0
   49 GO TO (28,99,58,89), JGOTO

C----              SEPARATOR 'DOLLAR', SPECIAL ITEMS

   51 IF (JDO+1.GE.JRRAN)          GO TO 94
      JDO = JDO + 2
      J   = IQCETK(JDO-1)
      IF (J.EQ.8)                  GO TO 25
      IF (J.EQ.47)                 GO TO 66
      IF (J.EQ.15)                 GO TO 71
      IF (J.EQ.27)                 GO TO 71

      IF (IFLBLT.EQ.0)             GO TO 54
      JTERM = JDO - 1
   53 JTERM = JTERM + 1
      IF (JTERM.GT.JRRAN)          GO TO 54
      IF (IQCETK(JTERM).EQ.45)     GO TO 53

   54 JR    = JTERM - 1
      NCH   = 2*IUFORW (IQHOLK,JDO,JR)
      JR    = JNEXT

C--                HANDLE  $A, $A7, $Q, $Q2

      IF (J.EQ.1)                  GO TO 57
      IF (J.NE.17)                 GO TO 94
      NCH = NCH + 1
   57 IF (JBIT(MFLAG,12).EQ.0)     GO TO 94
      JGOTO = 3
      NSEP = IQCETK(JR)
      GO TO 41

   58 MFLAGP(1) = MSBYT (NCH,MFLAGP(1),15,6)
      MODE = -2
      RETURN

C--                $. COMMENT $

   66 JR = JDO - 1
   67 JR = JR + 1
      IF (JR.GT.JRRAN)             GO TO 68
      IF (IQCETK(JR).NE.43)        GO TO 67

   68 MODE  = -3
      JGOTO = 2
   69 NSEP  = 0
      IF (JR.GT.JRRAN)       RETURN
      GO TO 45

C----              $O OCTAL

   71 JDO = JDO - 1
   72 JDO = JDO + 1
      IF (JDO.GT.JRRAN)            GO TO 94
      IF (IQCETK(JDO).EQ.45)       GO TO 72

      IF (IFLBLT.EQ.0)             GO TO 76
      JTERM = JDO
   74 JTERM = JTERM + 1
      IF (JTERM.GT.JRRAN)          GO TO 76
      IF (IQCETK(JTERM).NE.45)     GO TO 74

   76 NOCT = 0
      JR   = JDO - 1
   77 JR   = JR  + 1
      NSEP = IQCETK(JR)
      IF (JR.EQ.JTERM)             GO TO 81
      IF (NSEP.EQ.45)              GO TO 77
      J = NSEP - 27
      IF (J+NOCT.EQ.0)             GO TO 77
      IF (J.LT.0)                  GO TO 94
      IF (J.GE.8)                  GO TO 81
      IF (NOCT.EQ.36)              GO TO 94
      NOCT = NOCT + 1
      NUMOCT(NOCT) = J
      GO TO 77

   81 JGOTO = 4
      IVAL  = 0
      IF (NOCT.EQ.0)               GO TO 41
      NDOK  = (IQBITW-1) / 3
      N     = MIN  (NDOK,NOCT)
      JDF = 1 + NOCT - N

      DO 86  J=JDF,NOCT
   86 IVAL = 8*IVAL + NUMOCT(J)
      IF (N.EQ.NOCT)               GO TO 41
      J = 3*NDOK
      N = IQBITW - J
      IVAL = MSBYT (NUMOCT(JDF-1),IVAL,J+1,N)
      GO TO 41

   89 VALP(1) = VAL
      MODE    = 2
      NWORDS  = 1
      MFMT    = 2
      RETURN

C----              EXITS

   94 MODE = 0
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
