*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   cv single-pr. F from CRAY -> IEEE
*
* fzocxff.inc
*
#if defined(CERNLIB_QMCRY)
C--      Cray single-pr. to IEEE single
  301 DO 349  JL=1,NWDO
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 347
      JSIGN = SHIFTR (MASK(1).AND.ITHA, 32)
      JEXP  = SHIFTR (SHIFTL(ITHA,1), 49)
      JEXP  = JEXP - 40000B + 126
      JMANT = MASK(128-24) .AND. SHIFTR(ITHA,23)
      IF (JMANT.EQ.MASK(128-24))  THEN
          JMANT = 0
          JEXP  = JEXP + 1
        ELSE
          JMANT = SHIFTR(JMANT+1,1)
        ENDIF
      IF (JEXP.LE.0)               GO TO 331
      IF (JEXP.GE.255)             GO TO 333
      ITHA  = JSIGN .OR. SHIFTL(JEXP,23) .OR. JMANT
      GO TO 347

C--      make de-normalized number
  331 IF (JEXP.LT.-23)             GO TO 341
      JMANT = SHIFTR (JMANT+40000000B,1-JEXP)
      ITHA  = JSIGN .OR. JMANT
      GO TO 347

C--      overflow
  333 IFOCON(1) = 3
      IFOCON(2) = JMS
      IFOCON(3) = ITHA
      IF (JSIGN.EQ.0)    THEN
          ITHA = IOVPCS
        ELSE
          ITHA = IOVNCS
        ENDIF
      GO TO 347

  341 ITHA = 0
  347 MT(JMT+1) = ITHA
      JMT = JMT + 1
  349 JMS = JMS + 1
      GO TO 801
#endif
*   cv single-pr. F from ND -> IEEE
*
* fzocxff.inc
*
#if defined(CERNLIB_QMND3)
C--      NORD single-precision to IEEE single
  301 DO 349  JL=1,NWDO
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 347
      JSIGN = ISHFT (ITHA,-31)
      JEXP  = ISHFT (ISHFT(ITHA,1), -23)
      JEXP  = JEXP - 130
      IF (JEXP.LE.0)               GO TO 331
      IF (JEXP.GE.255)             GO TO 333
      JMANT = ISHFT (ISHFT(ITHA,10), -9)
      ITHA  = ISHFT(JSIGN,31) .OR. ISHFT(JEXP,23) .OR. JMANT
      GO TO 347

C--      make de-normalized number
  331 IF (JEXP.LE.-23)             GO TO 341
      ITHA = ITHA .OR. 20000000 B
      ITHA = ISHFT (ISHFT(ITHA,9),JEXP-9)
      ITHA = ITHA .OR. ISHFT(JSIGN,31)
      GO TO 347

C--      overflow
  333 IFOCON(1) = 3
      IFOCON(2) = JMS
      IFOCON(3) = ITHA
      IF (ITHA.NE.IBADMS)  THEN
          IF (JSIGN.EQ.0)    THEN
              ITHA = IOVPCS
            ELSE
              ITHA = IOVNCS
            ENDIF
        ELSE
          ITHA = IBADCS
        ENDIF
      GO TO 347

  341 ITHA = 0
  347 MT(JMT+1) = ITHA
      JMT = JMT + 1
  349 JMS = JMS + 1
      GO TO 801
#endif
