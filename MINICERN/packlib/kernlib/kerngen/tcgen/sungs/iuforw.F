*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION IUFORW (IC,JL,JR)
C
C CERN PROGLIB# M410    IUFORW          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
      DIMENSION    IC(99)
      COMMON /SLATE/ ND,JEND, DUMMY(38)

      ND   = 0
      JEND = JL
      IVAL = 0

   12 IF (JEND.GT.JR)        GO TO 20
#if defined(CERNLIB_BUGLRSHFT)
      K = ishft (IC(JEND),-24)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      K = lrshft (IC(JEND),24)
#endif
      IF (K.EQ.32)           GO TO 19
      IF (K.GE.58)           GO TO 20
      K = K - 48
      IF (K.LT.0)            GO TO 20
      IVAL = 10*IVAL + K
      ND   = ND + 1
   19 JEND = JEND + 1
      GO TO 12

   20 IUFORW = IVAL
      RETURN
      END
