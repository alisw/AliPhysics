*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:23  mclareni
* Kernlib
*
*
      SUBROUTINE IE3FOS (MS,MT,NWDO,JBAD)
C
C CERN PROGLIB# M220    IE3FOS          .VERSION KERNVAX  2.20  861204

C-    Convert single precision for output with copy
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


      JBAD = 0
      JMS  = 0
      JMT  = 0

  301 DO 349  JL=1,NWDO
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)         GO TO 347
      JSIGN = ITHA .AND. MSKB16
      JEXP  = JIBITS (ITHA,7,8)
      IF (JEXP.LE.2)         GO TO 331
      IF (JEXP.EQ.255)       GO TO 333
      JMANT = JISHFT(JIBITS(ITHA,0,7), 16) .OR. JIBITS(ITHA,16,16)
      ITHA  = JISHFT(JSIGN,16) .OR. JISHFT(JEXP-2,23) .OR. JMANT
      GO TO 347

C--      zero / NaN
  331 IF (JSIGN.EQ.0)        GO TO 341
      IF (JEXP.NE.0)         GO TO 341
      ITHB = IBADCS
      GO TO 336

C--      overflow
  333 IF (JSIGN.EQ.0)    THEN
          ITHB = IOVPCS
        ELSE
          ITHB = IOVNCS
        ENDIF
  336 JBAD = JL
      ITHA = ITHB
      GO TO 347

  341 ITHA = 0
  347 MT(JMT+1) = ITHA
      JMT = JMT + 1
  349 JMS = JMS + 1

      RETURN
      END
