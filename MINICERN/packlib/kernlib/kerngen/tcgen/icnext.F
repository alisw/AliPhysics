*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNEXT (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICNEXT          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Find limits of next word in CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT, DUMMY(38)
      CHARACTER    CHV*(*)

C--       Ignore leading blanks

      JJ = JLP(1)
      JR = JRP(1)

   12 IF (JJ.GT.JR)          GO TO 99
      IF (CHV(JJ:JJ).EQ.' ')  THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF

C--       Find terminating blank

      JF = JJ
   14 JJ = JJ + 1
      IF (JJ.LE.JR)  THEN
          IF (CHV(JJ:JJ).NE.' ')  GO TO 14
      ENDIF

      NESLAT = JJ
      NDSLAT = JJ - JF
      ICNEXT = JF
      RETURN

   99 NDSLAT = 0
      NESLAT = JJ
      ICNEXT = JJ
      RETURN
      END
