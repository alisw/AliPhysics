*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:45  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICFILA (CHI, CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICFILA          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Find last occurrence of CHI(1:1) in CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*), CHI*1

      JL = JLP(1)
      JR = JRP(1)
      JJ = JR

   11 IF (JJ.LT.JL)          GO TO 19
      IF (CHV(JJ:JJ).NE.CHI)  THEN
          JJ = JJ - 1
          GO TO 11
        ENDIF
      NGSLAT = JJ
      ICFILA = JJ
      RETURN

   19 NGSLAT = 0
      ICFILA = JR+1
      RETURN
      END
