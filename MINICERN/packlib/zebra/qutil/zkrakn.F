*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZKRAKN (MFLAGP,JLP,JRP,VALP)

#include "zebra/zbcd.inc"
#include "zebra/zkrakc.inc"
#include "zebra/zkrakqu.inc"
C--------------    END CDE                             --------------

      COMMON /SLATE/ NDIG, JNEXT, DUMMY(38)

      DIMENSION    MFLAGP(9), JLP(9), JRP(9), VALP(9)
      INTEGER      IVAL
      EQUIVALENCE (VAL,IVAL)


#include "zebra/q_jbit.inc"


      MFLAG  = MFLAGP(1)
      JLRAN  = JLP(1)
      JRRAN  = JRP(1)
      MODE   = 0
      NWORDS = 0
      MFMT   = 0
      MULT   = JBIT(MFLAG,10) - 1
      IFLBLT = JBIT(MFLAG,6)
   12 NSEP   = 0
      JR     = JRRAN
      JTERM  = JRRAN + 1
      JNXGO  = JTERM

C-    IZBCD VAL.             1 27 36 37 38 39 40 41 42 43 44 45 46 47
C-    CHARACTER              A  0  9  +  -  *  /  (  )  $  = BL  ,  .
C-    BIT J IN MFLAG IF SEPARATOR LEGAL,   J=  1  2  3  4  5  6  7
C-    J=10   MULTIPLIER LEGAL

C--                IGNORE LEADING BLANKS

      JDO = JLRAN - 1
   13 JDO = JDO + 1
      IF (JDO.GT.JRRAN)            GO TO 93
      IFIRST = IQCETK(JDO)
      IF (IFIRST.EQ.45)            GO TO 13

      IF (IFIRST.LT.27)            GO TO 91
      IF (IFIRST.LT.39)            GO TO 14
      IF (IFIRST.NE.47)            GO TO 92

C--                FIND BLANK TERMINATOR, IF SELECTED

   14 IF (IFLBLT.EQ.0)             GO TO 17
      JTERM = JDO
   15 JTERM = JTERM + 1
      IF (JTERM.GT.JRRAN)          GO TO 16
      IF (IQCETK(JTERM).NE.45)     GO TO 15
   16 JR    = JTERM - 1

C--                STEP OVER SIGN-BIT

   17 IF (IFIRST.EQ.37)            GO TO 18
      IF (IFIRST.NE.38)            GO TO 21
   18 JDO = JDO +1
      IF (JDO.EQ.JTERM)            GO TO 94

C--                READ INITIAL INTEGER

   21 NDALL = 0
      IVAL  = 0
      JDO   = JDO - 1
   22 JDO   = JDO + 1
      IF (JDO.EQ.JTERM)            GO TO 28
      J = IQCETK(JDO) - 27
      IF (J.LT.0)                  GO TO 24
      IF (J.EQ.18)                 GO TO 22
      IF (J.GE.10)                 GO TO 24
      IVAL  = IVAL*10 + J
      NDALL = NDALL + 1
      GO TO 22

   24 JEXP = 0
      JEXF = 0
      FRAC = 0.
      NSEP = IQCETK(JDO)
      IF   (NSEP-39)         41, 81, 26
   26 IF (NSEP.EQ.47)              GO TO 31

   28 IF (IFIRST.EQ.38)  IVAL=-IVAL
      MFMT = 1
      GO TO 71

C--                READ FRACTION

   31 JEXW = 0
   32 JDO  = JDO + 1
      IF (JDO.EQ.JTERM)            GO TO 50

      NSEP = IQCETK(JDO)
      IF (NSEP-27)           41, 33, 34
   33 NDALL = NDALL + 1
      JEXW  = JEXW  + 1
      GO TO 32

   34 IF (NSEP.GE.37)              GO TO 38
   35 IF (JEXW.EQ.0)               GO TO 36
      JEXW  = JEXW  - 1
      JEXF  = JEXF  - 1
      FRAC  = FRAC * 10.
      GO TO 35

   36 NDALL = NDALL + 1
      JEXF  = JEXF  - 1
      FRAC  = FRAC*10. + REAL(NSEP-27)
      GO TO 32

   38 IF (NSEP.EQ.45)              GO TO 32
      IF (NSEP-39)           44, 94, 51

C--                READ EXPONENT

   41 IF (NSEP.NE.5)               GO TO 44
   42 JDO = JDO + 1
      IF (JDO.EQ.JTERM)            GO TO 50
      NSEP = IQCETK(JDO)
      IF (NSEP.EQ.45)              GO TO 42
      IF (NSEP-39)           44, 94, 51

   44 IF (NSEP.LT.27)              GO TO 94
      IF (NSEP.LT.37)              GO TO 46
      JDO = JDO + 1
      IF (JDO.EQ.JTERM)            GO TO 94

   46 JEXP = IUFORW (IQHOLK,JDO,JR)
      IF (JEXP.EQ.0)               GO TO 47
      IF (NSEP.EQ.38)  JEXP=-JEXP
      JEXF = JEXF + JEXP
   47 JDO  = JNEXT
      IF (JDO.EQ.JTERM)            GO TO 50
      NSEP = IQCETK(JDO)
      IF (NSEP.GE.40)              GO TO 51
      GO TO 94

C--                SET FLOATING RESULT

   50 NSEP = 0
   51 VAL  = REAL (IVAL)
      IF (VAL.EQ.0.)               GO TO 61
      IF   (JEXP)            58, 61, 52
   52 IF (JEXP.GE.9)               GO TO 58
      DO 54  J=1,JEXP
   54 VAL = VAL * 10.
      GO TO 61

   58 VAL = VAL * 10.**JEXP

C--                SET FLOATING RESULT, FRACTIONAL PART

   61 IF (FRAC.EQ.0.)              GO TO 66
      IF   (JEXF)            64, 65, 62
   62 IF (JEXF.GE.9)               GO TO 64
      DO 63  J=1,JEXF
   63 FRAC = FRAC * 10.
      GO TO 65

   64 FRAC = FRAC * 10.**JEXF

   65 VAL = VAL + FRAC
   66 IF (IFIRST.NE.38)            GO TO 71
      VAL = -VAL

C--                STORE RESULT

   71 IF (NDALL.EQ.0)              GO TO 94
      IF (NSEP.EQ.0)               GO TO 72
      NSEP = MIN  (NSEP,47)
      IF (JBIT(MFLAG,NSEP-39).EQ.0)  GO TO 94
      IF (NSEP.EQ.43)              GO TO 72
      JDO    = JDO + 1

   72 JNXGO  = JDO
   73 IF (JNXGO.GT.JRRAN)          GO TO 75
      IF (IQCETK(JNXGO).NE.45)     GO TO 75
      JNXGO  = JNXGO + 1
      GO TO 73

   75 MODE   = 2
      NWORDS = 1
      VALP(1) = VAL
      GO TO 94

C----              MULTIPLIER

   81 IF (IVAL.LE.0)               GO TO 94
      IF (MULT.NE.0)               GO TO 94
      MULT  = IVAL
      JLRAN = JDO + 1
      GO TO 12

C----              NON-NUMERIC

   91 MFMT = -1
   92 MFMT = MFMT - 1
      NSEP = IFIRST
      GO TO 94

   93 MODE = 1
   94 MULT = MAX (0,MULT-1)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
