*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNUMU (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICNUMU          .VERSION KERNFOR  4.40  940929
C ORIG. 19/07/93, JZ
C
C-    Verify CHV(JL:JR) is alpha-numeric or underscore

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
      IF (CHWK.EQ.'_')       GO TO 17
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
   17 IF (JPA.EQ.0)  JPA = JJ
   18 ND  = ND + 1
   19 JJ  = JJ + 1
      GO TO 12

   91 NGSLAT = JJ
   99 NDSLAT = ND
      NESLAT = JPN
      NFSLAT = JPA
      ICNUMU = JJ
      RETURN
      END
