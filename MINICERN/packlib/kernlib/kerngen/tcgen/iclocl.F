*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICLOCL (CHI,NI,CHV,JL,JR)
C
C CERN PROGLIB# M432    ICLOCL          .VERSION KERNFOR  4.22  890913
C ORIG. 09/02/89, JZ
C
C-    Locate CHI of NI characters inside CHV(JL:JR),
C-    CHI must be given in lower case, CHV case insensitive

      CHARACTER    CHI*(*), CHV*(*), CHWK1*1
#include "kerngen/qnatch.inc"
* Ignoring t=pass

      JVV   = JL - 1
      JVVE  = JR+1 - NI
      CHWK1 = CHI(1:1)
#if defined(CERNLIB_QASCII)
      NAT1  = 0
      IVX1  = ICHAR(CHWK1) - 32
      IF (IVX1.LT.65)        GO TO 12
      IF (IVX1.GE.91)        GO TO 12
      NAT1  = 3
#endif
#if defined(CERNLIB_QEBCDIC)
      IVX1  = ICHAR(CHWK1)
      NAT1  = NATCH(IVX1+1)
      IVX1  = IVX1 + 64
#endif

C--          Scan CHV to find the first char. of CHI

   12 JVV = JVV + 1
      IF (JVV.GT.JVVE)                     GO TO 90
      IF (CHV(JVV:JVV).NE.CHWK1)  THEN
          IF (NAT1.NE.3)                   GO TO 12
          IF (ICHAR(CHV(JVV:JVV)).NE.IVX1) GO TO 12
        ENDIF

C--          Match the further characters

      JVX = JVV
      JI  = 1
   14 JI  = JI + 1
      IF (JI.GT.NI)                        GO TO 91
      JVX = JVX + 1
      IF (CHV(JVX:JVX).EQ.CHI(JI:JI))      GO TO 14
#if defined(CERNLIB_QASCII)
      IVAL = ICHAR(CHI(JI:JI))
      IF (ICHAR(CHV(JVX:JVX)).NE.IVAL-32)  GO TO 12
      IF (IVAL.LT.97)                      GO TO 12
      IF (IVAL.GE.123)                     GO TO 12
#endif
#if defined(CERNLIB_QEBCDIC)
      IVAL = ICHAR(CHI(JI:JI))
      IF (ICHAR(CHV(JVX:JVX)).NE.IVAL+64)  GO TO 12
      IF (NATCH(IVAL+1).NE.3)              GO TO 12
#endif
      GO TO 14

   90 JVV = 0
   91 ICLOCL = JVV
      RETURN
      END
