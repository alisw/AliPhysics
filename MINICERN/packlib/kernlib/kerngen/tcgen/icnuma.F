*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNUMA (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICNUMA          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Verify CHV(JL:JR) is alpha-numeric

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*), CHWK*1

#include "kerngen/qnatch.inc"
* Ignoring t=pass

      JJ = JLP(1)
      JR = JRP(1)

      NGSLAT = 0
      JPN = 0
      JPA = 0
      ND  = 0

   12 IF (JJ.GT.JR)          GO TO 99
      CHWK = CHV(JJ:JJ)
      IF (CHWK.EQ.' ')       GO TO 19
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHWK)
      IF (K.LT.48)           GO TO 91
      IF (K.LT.58)  THEN
          IF (JPN.EQ.0)  JPN = JJ
          GO TO 18
        ENDIF

      IF     (K.LT.65)       GO TO 91
      IF     (K.GE.91)  THEN
          IF (K.LT.97)       GO TO 91
          IF (K.GE.123)      GO TO 91
        ENDIF
#endif
#if defined(CERNLIB_QEBCDIC)
      J = ICHAR (CHWK)
      K = NATCH(J+1)
      IF (K.LE.1)            GO TO 91
      IF (K.EQ.2)  THEN
          IF (JPN.EQ.0)  JPN = JJ
          GO TO 18
        ENDIF
#endif
      IF (JPA.EQ.0)  JPA = JJ
   18 ND  = ND + 1
   19 JJ  = JJ + 1
      GO TO 12

   91 NGSLAT = JJ
   99 NDSLAT = ND
      NESLAT = JPN
      NFSLAT = JPA
      ICNUMA = JJ
      RETURN
      END
