*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      SUBROUTINE IE3FOD (MS,MT,NDPN,JBAD)
C
C CERN PROGLIB# M220    IE3FOD          .VERSION KERNVAX  2.20  861204

C-    Convert double precision for output with copy
C-    from source in native to target in IEEE data format

      DIMENSION    MS(99), MT(99)

      DOUBLE PRECISION   THDB
      DIMENSION    THIS(2)
      EQUIVALENCE  (THDB,THIS)
      EQUIVALENCE  (ITHA,THA,THIS(1)), (ITHB,THB,THIS(2))

      PARAMETER    (MSKB16 = '00008000'X)
      PARAMETER    (IBADCS = '7F80FE00'X)
      PARAMETER    (IBADCD = '7FF01FC0'X)
      PARAMETER    (IOVPCS = '7F800000'X)
      PARAMETER    (IOVPCD = '7FF00000'X)
      PARAMETER    (IOVNCS = 'FF800000'X)
      PARAMETER    (IOVNCD = 'FFF00000'X)


      JFAI = 0
      JMS  = 0
      JMT  = 0

      DO 449  JL=1,NDPN
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)         GO TO 442
      ITHB  = MS(JMS+2)
      JSIGN = ITHA .AND. MSKB16
      JEXP  = JIBITS (ITHA,7,8)
      IF (JEXP.EQ.0)         GO TO 431
      IF (JEXP.EQ.255)       GO TO 433
      JMANT = JISHFT (JIBITS(ITHA, 0, 7), 13) .OR. JIBITS(ITHA,19,13)
      ITHB  = JISHFT (JIBITS(ITHA,16, 3), 29) .OR.
     +        JISHFT (JIBITS(ITHB, 0,16), 13) .OR.
     +                JIBITS(ITHB,19,13)
      ITHA  = JISHFT(JSIGN,16) .OR. JISHFT(JEXP+894,20) .OR. JMANT
      GO TO 447

C--      zero / NaN
  431 IF (JSIGN.EQ.0)        GO TO 441
      ITHB = IBADCD
      GO TO 436

C--      overflow
  433 IF (JSIGN.EQ.0)    THEN
          ITHB = IOVPCD
        ELSE
          ITHB = IOVNCD
        ENDIF
  436 JFAI = JL
      ITHA = ITHB
      GO TO 442

  441 ITHA = 0
  442 ITHB = 0
  447 MT(JMT+1) = ITHA
      MT(JMT+2) = ITHB
      JMT = JMT + 2
  449 JMS = JMS + 2

      JBAD = JFAI
      RETURN
      END
