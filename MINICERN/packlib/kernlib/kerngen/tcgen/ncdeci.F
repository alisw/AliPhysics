*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION NCDECI (CHV)
C
C CERN PROGLIB# M432    NCDECI          .VERSION KERNFOR  4.33  920511
C ORIG. 03/06/92, JZ
C
C-    Read decimal integer from CHV

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)

      JJ = 1
      JR = LEN(CHV)

      IVAL = 0
      NDG  = 0
      NEG  = 0
      NGSLAT = 0

   12 IF (JJ.GT.JR)          GO TO 99
      IF (CHV(JJ:JJ).EQ.' ')  THEN
          JJ = JJ + 1
          GO TO 12
        ELSEIF (CHV(JJ:JJ).EQ.'+')  THEN
          JJ = JJ + 1
        ELSEIF (CHV(JJ:JJ).EQ.'-')  THEN
          NEG = 7
          JJ  = JJ + 1
        ENDIF

   21 IF (JJ.GT.JR)          GO TO 99
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ:JJ))
      K = K - 48
      IF (K.LT.0)            GO TO 98
      IF (K.GE.10)           GO TO 98
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ:JJ))
      K = K - 240
      IF (K.LT.0)            GO TO 98
      IF (K.GE.10)           GO TO 98
#endif
      IVAL = 10*IVAL + K
      NDG  = NDG + 1
      JJ   = JJ + 1
      GO TO 21

   98 IF (CHV(JJ:JJ).NE.' ')  NGSLAT = JJ
   99 NDSLAT = NDG
      NESLAT = JJ
      IF (NEG.NE.0)  THEN
          IF (IVAL.NE.0)  IVAL = -IVAL
        ENDIF
      NCDECI = IVAL
      RETURN
      END
