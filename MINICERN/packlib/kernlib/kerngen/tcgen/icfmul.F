*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:45  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICFMUL (CHI, CHV,JLP,JRP)
C
C CERN PROGLIB# M432    ICFMUL          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Find first occurrence of any of CHI(j:j) in CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*), CHI*(*)

      JJ = JLP(1)
      JR = JRP(1)

   11 IF (JJ.GT.JR)          GO TO 19
      JTH = INDEX (CHI,CHV(JJ:JJ))
      IF (JTH.EQ.0)  THEN
          JJ = JJ + 1
          GO TO 11
        ENDIF
      NDSLAT = JTH
      NGSLAT = JJ
      ICFMUL = JJ
      RETURN

   19 NDSLAT = 0
      NGSLAT = 0
      ICFMUL = JR+1
      RETURN
      END
