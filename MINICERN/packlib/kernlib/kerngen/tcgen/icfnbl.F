*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:45  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICFNBL (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICFNBL          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Find first non-blank character in CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)

      JJ = JLP(1)
      JR = JRP(1)

   12 IF (JJ.GT.JR)          GO TO 19
      IF (CHV(JJ:JJ).EQ.' ')  THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF
      NGSLAT = JJ
      ICFNBL = JJ
      RETURN

   19 NGSLAT = 0
      ICFNBL = JJ
      RETURN
      END
