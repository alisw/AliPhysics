*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
*   cv IEEE -> single
*
* fzicvff.inc
*
#if defined(CERNLIB_FQIE3TSC)
C--       Default conversion from single IEEE to internal
  301 CALL IE3TOS (MS(JMS+1),MT(JMT+1),NWDO,JBAD)
      IF (JBAD.NE.0)               GO TO 331
  304 JMT = JMT + NWDO
      JMS = JMS + NWDO
      GO TO 801

  331 IFOCON(1) = 3
      IFOCON(2) = JMT + JBAD
      IFOCON(3) = MS(JMS+JBAD)
      GO TO 304
#endif
*    cv IEEE -> single
*
* fzicvff.inc
*
#if defined(CERNLIB_FQCVDFF)
C--       Default conversion from single IEEE to internal
  301 DO 349  JL=1,NWDO
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 347
      JSIGN = JBIT (ITHA,32)
      JEXP  = JBYT (ITHA,24,8)
      JMANT = JBYT (ITHA,1,23)
      IF (JEXP.EQ.0)               GO TO 331
      IF (JEXP.EQ.255)             GO TO 333
      JEXP  = JEXP - 127
*   test on exponent range
#if defined(CERNLIB_FQCVDFFT)
#include "fzicvff2.inc"
#endif
      JPRE  = ISIGN (23,JEXP)
      JMANT = JMANT + 8 388 608
      THA   = FLOAT(JMANT) * 2.**(JPRE-23)
      THA   = THA * 2.**(JEXP-JPRE)
      IF (JSIGN.EQ.0)              GO TO 347
      THA  = -THA
      GO TO 347

  331 IF (JMANT.EQ.0)              GO TO 341
      THA  = FLOAT(JMANT) * 2.**(-63)
      THA  = THA * 2.**(-86)
      IF (ITHA .EQ.0)              GO TO 347
      IF (JSIGN.EQ.0)              GO TO 347
      THA  = -THA
      GO TO 347

#if defined(CERNLIB_FQCVDFFT)
#include "fzicvff3.inc"
#endif
  333 IF (JMANT.EQ.0)    THEN
          IF (JSIGN.EQ.0)  THEN
              ITHB = IOVPMS
            ELSE
              ITHB = IOVNMS
            ENDIF
        ELSE
          ITHB = IBADMS
        ENDIF
      IFOCON(1) = 3
      IFOCON(2) = JMT
      IFOCON(3) = ITHA
      ITHA = ITHB
      GO TO 347

  341 ITHA = 0
  347 MT(JMT+1) = ITHA
      JMT = JMT + 1
  349 JMS = JMS + 1
      GO TO 801
#endif
