*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CRIGHT (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CRIGHT          .VERSION KERNFOR  4.21  890323
C ORIG. 04/10/88, JZ
C
C-    Right-justify CHV(JL:JR)

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT, DUMMY(38)
      CHARACTER    CHV*(*)

C--       Find first blank

      JL = JLP(1)
      JR = JRP(1)

      JP = JL - 1
      JJ = JR
   12 IF (JJ.LE.JP)          GO TO 99
      IF (CHV(JJ:JJ).NE.' ')   THEN
          JJ = JJ - 1
          GO TO 12
        ENDIF

C--       Copy shifted

      JP = JJ
      JJ = JJ - 1
   21 IF (JJ.LT.JL)          GO TO 91
      IF (CHV(JJ:JJ).NE.' ')   THEN
          CHV(JP:JP) = CHV(JJ:JJ)
          JP = JP - 1
        ENDIF
      JJ = JJ - 1
      GO TO 21

   91 CHV(JL:JP) = ' '
   99 NESLAT = JP
      NDSLAT = JR - JP
      RETURN
      END
