*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:45  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICHEXI (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICHEXI          .VERSION KERNFOR  4.35  930308
C ORIG. 04/10/88, JZ
C
C-    Read hexadecimal integer from CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)

#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
#endif

      JJ = JLP(1)
      JR = JRP(1)

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
      IF (K.GE.10)  THEN
          K = K - 7
          IF (K.LT.10)       GO TO 98
          IF (K.GE.16)  THEN
              K = K - 32
              IF (K.LT.10)   GO TO 98
              IF (K.GE.16)   GO TO 98
            ENDIF
        ENDIF
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ:JJ))
      K = K - 119
      IF (K.LT.10)           GO TO 98
      IF (K.GE.16)  THEN
          K = K - 64
          IF (K.LT.10)       GO TO 98
          IF (K.GE.16)  THEN
              K = K - 57
              IF (K.LT.0)    GO TO 98
              IF (K.GE.10)   GO TO 98
            ENDIF
        ENDIF
#endif
#if !defined(CERNLIB_QISASTD)
      IVAL = IOR (ISHFTL(IVAL,4), K)
#endif
#if defined(CERNLIB_QISASTD)
      IVAL = IOR (ISHFT (IVAL,4), K)
#endif
      NDG  = NDG + 1
      JJ   = JJ + 1
      GO TO 21

   98 IF (CHV(JJ:JJ).NE.' ')  NGSLAT = JJ
   99 NDSLAT = NDG
      NESLAT = JJ
      IF (NEG.NE.0)  THEN
          IF (IVAL.NE.0)  IVAL = -IVAL
        ENDIF
      ICHEXI = IVAL
      RETURN
      END
