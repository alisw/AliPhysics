*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:23  mclareni
* Kernlib
*
*
      SUBROUTINE IE3TOS (MS,MT,NWDO,JBAD)
C
C CERN PROGLIB# M220    IE3TOS          .VERSION KERNVAX  2.25  880302

C-    Convert single precision for input with copy
C-    from source in IEEE to target in native data format

      DIMENSION    MS(99), MT(99)

      DOUBLE PRECISION   THDB
      DIMENSION    THIS(2)
      EQUIVALENCE  (THDB,THIS)
      EQUIVALENCE  (ITHA,THA,THIS(1)), (ITHB,THB,THIS(2))

      PARAMETER    (JEXMIN = -127,  JEXMAX=125)
      PARAMETER    (IOVPMS = '00007F81'X)
      PARAMETER    (IOVNMS = '0000FF81'X)
      PARAMETER    (IBADMS = '00008001'X)
      PARAMETER    (IBADMD = IBADMS, IBADME = 0)
      PARAMETER    (IOVPMD = IOVPMS, IOVNMD = IOVNMS)

#include "kerngen/q_jbit.inc"

      JBAD = 0
      JMS  = 0
      JMT  = 0

  301 DO 349  JL=1,NWDO
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)         GO TO 347
      JSIGN = JBIT (ITHA,32)
      JEXP  = JBYT (ITHA,24,8)
      JMANT = JBYT (ITHA,1,23)
      IF (JEXP.EQ.0)         GO TO 331
      IF (JEXP.EQ.255)       GO TO 333
      JEXP  = JEXP - 127
      IF (JEXP.LT.JEXMIN)    GO TO 341
      IF (JEXP.GT.JEXMAX)    GO TO 332
      JPRE  = ISIGN (23,JEXP)
      JMANT = JMANT + 8 388 608
      THA   = FLOAT(JMANT) * 2.**(JPRE-23)
      THA   = THA * 2.**(JEXP-JPRE)
      IF (JSIGN.EQ.0)        GO TO 347
      THA  = -THA
      GO TO 347

  331 IF (JMANT.EQ.0)        GO TO 341
      THA  = FLOAT(JMANT) * 2.**(-63)
      THA  = THA * 2.**(-86)
      IF (ITHA .EQ.0)        GO TO 347
      IF (JSIGN.EQ.0)        GO TO 347
      THA  = -THA
      GO TO 347

  332 JMANT = 0
  333 IF (JMANT.EQ.0)    THEN
          IF (JSIGN.EQ.0)  THEN
              ITHB = IOVPMS
            ELSE
              ITHB = IOVNMS
            ENDIF
        ELSE
          ITHB = IBADMS
        ENDIF
      JBAD = JL
      ITHA = ITHB
      GO TO 347

  341 ITHA = 0
  347 MT(JMT+1) = ITHA
      JMT = JMT + 1
  349 JMS = JMS + 1

      RETURN
      END
