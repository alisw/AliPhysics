*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNUM (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICNUM           .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Verify CHV(JL:JR) is numeric

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)

      JJ = JLP(1)
      JR = JRP(1)

      NGSLAT = 0
      ND  = 0

   12 IF (JJ.GT.JR)          GO TO 99
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ:JJ))
      IF (K.EQ.32)           GO TO 19
      K = K - 48
      IF (K.LT.0)            GO TO 91
      IF (K.GE.10)           GO TO 91
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ:JJ))
      IF (K.EQ.64)           GO TO 19
      K = K - 240
      IF (K.LT.0)            GO TO 91
      IF (K.GE.10)           GO TO 91
#endif
      ND  = ND + 1
   19 JJ  = JJ + 1
      GO TO 12

   91 NGSLAT = JJ
   99 NDSLAT = ND
      ICNUM  = JJ
      RETURN
      END
