*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CLEFT (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CLEFT           .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Left-justify CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT, DUMMY(38)
      CHARACTER    CHV*(*)

C--       Find the first blank

      JL = JLP(1)
      JR = JRP(1)

      JP = JR + 1
      JJ = JL
   12 IF (JJ.GE.JP)          GO TO 99
      IF (CHV(JJ:JJ).NE.' ')   THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF

C--       Copy shifted

      JP = JJ
      JJ = JJ + 1
   21 IF (JJ.GT.JR)          GO TO 91
      IF (CHV(JJ:JJ).NE.' ')   THEN
          CHV(JP:JP) = CHV(JJ:JJ)
          JP = JP+1
        ENDIF
      JJ = JJ+1
      GO TO 21

   91 CHV(JP:JR) = ' '
   99 NESLAT = JP
      NDSLAT = JP - JL
      RETURN
      END
