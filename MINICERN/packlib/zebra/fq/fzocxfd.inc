*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   cv double-pr. F from CRAY -> IEEE
*
* fzocxfd.inc
*
#if defined(CERNLIB_QMCRY)
C--      Cray single-pr. to IEEE double, ignoring 2nd word
C-                   for the moment (loss of 5 bits)
      DO 449  JL=1,NWDODB,2
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 442
      JSIGN = SHIFTR (MASK(1).AND.ITHA, 32)
      JEXP  = SHIFTR (SHIFTL(ITHA,1), 49)
      JEXP  = JEXP - 40000B + 1022
      IF (JEXP.LE.0)               GO TO 441
      IF (JEXP.GE.2047)            GO TO 433
      JMANT = MASK(128-20) .AND. SHIFTR (ITHA,27)
      ITHB  = SHIFTL (ITHA, 5) .AND. 37777777777B
      ITHA  = JSIGN .OR. SHIFTL(JEXP,20) .OR. JMANT
      GO TO 447

C--      overflow
  433 IFOCON(1) = 4
      IFOCON(2) = JMS
      IFOCON(3) = ITHA
      IF (JSIGN.EQ.0)    THEN
          ITHA = IOVPCD
        ELSE
          ITHA = IOVNCD
        ENDIF
      GO TO 442

  441 ITHA = 0
  442 ITHB = 0
  447 MT(JMT+1) = ITHA
      MT(JMT+2) = ITHB
      JMT = JMT + 2
  449 JMS = JMS + 2
#endif
*   cv double-pr. F from ND -> IEEE
*
* fzocxfd.inc
*
#if defined(CERNLIB_QMND3)
C--      NORD double-precision to IEEE double
      DO 449  JL=1,NWDODB,2
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 442
      ITHB  = MS(JMS+2)
      JSIGN = ISHFT (ITHA,-31)
      JEXP  = ISHFT (ISHFT(ITHA,1), -23)
      IF (JEXP.EQ.0)               GO TO 441
      IF (JEXP.EQ.511)             GO TO 433
      JMANT = ISHFT (ISHFT(ITHA,10), -12)
      ITHB  = ISHFT (ITHA,30) .OR. ISHFT (ITHB,-2)
      ITHA  = ISHFT(JSIGN,31) .OR. ISHFT(JEXP+766,20) .OR. JMANT
      GO TO 447

C--      overflow / NaN
  433 IFOCON(1) = 4
      IFOCON(2) = JMS
      IFOCON(3) = ITHA
      IF (ITHA.NE.IBADMS)  THEN
          IF (JSIGN.EQ.0)    THEN
              ITHA = IOVPCD
            ELSE
              ITHA = IOVNCD
            ENDIF
        ELSE
          ITHA = IBADCD
        ENDIF
      GO TO 442

  441 ITHA = 0
  442 ITHB = 0
  447 MT(JMT+1) = ITHA
      MT(JMT+2) = ITHB
      JMT = JMT + 2
  449 JMS = JMS + 2
*    --------------  sequences for input   ---------------------
#endif
*   cv double-pr. F from VAX -> IEEE
*
* fzocxfd.inc
*
#if defined(CERNLIB_QMVAX)
C--      VAX double-precision to IEEE double
      DO 449  JL=1,NWDODB,2
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 442
      ITHB  = MS(JMS+2)
      JSIGN = ITHA .AND. MSKB16
      JEXP  = JIBITS (ITHA,7,8)
      IF (JEXP.EQ.0)               GO TO 431
      IF (JEXP.EQ.255)             GO TO 433
      JMANT = JISHFT (JIBITS(ITHA, 0, 7), 13) .OR. JIBITS(ITHA,19,13)
      ITHB  = JISHFT (JIBITS(ITHA,16, 3), 29) .OR.
     +        JISHFT (JIBITS(ITHB, 0,16), 13) .OR.
     +                JIBITS(ITHB,19,13)
      ITHA  = JISHFT(JSIGN,16) .OR. JISHFT(JEXP+894,20) .OR. JMANT
      GO TO 447

C--      zero / NaN
  431 IF (JSIGN.EQ.0)              GO TO 441
      ITHB = IBADCD
      GO TO 436

C--      overflow
  433 IF (JSIGN.EQ.0)    THEN
          ITHB = IOVPCD
        ELSE
          ITHB = IOVNCD
        ENDIF
  436 IFOCON(1) = 4
      IFOCON(2) = JMS
      IFOCON(3) = ITHA
      ITHA = ITHB
      GO TO 442

  441 ITHA = 0
  442 ITHB = 0
  447 MT(JMT+1) = ITHA
      MT(JMT+2) = ITHB
      JMT = JMT + 2
  449 JMS = JMS + 2
#endif
