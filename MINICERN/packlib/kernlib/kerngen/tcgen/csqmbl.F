*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSQMBL (CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CSQMBL          .VERSION KERNFOR  4.22  890913
C ORIG. 28/06/89, JZ
C
C-    Squeeze multiple blanks in CHV(JL:JR), shifting left

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT, DUMMY(38)
      CHARACTER    CHV*(*)

C----              Find the first blank

      JL = JLP(1)
      JR = JRP(1)

      JP = JR + 1
      JJ = JL
   12 IF (JJ.GE.JP)                GO TO 99
      IF (CHV(JJ:JJ).NE.' ')   THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF

C--       is it multiple ?

      JJ = JJ + 1
      IF (JJ.GE.JP)                GO TO 99
      IF (CHV(JJ:JJ).NE.' ')       GO TO 12

      JP = JJ
      JJ = JJ + 1
   16 IF (JJ.GT.JR)                GO TO 99
      IF (CHV(JJ:JJ).EQ.' ')   THEN
          JJ = JJ + 1
          GO TO 16
        ENDIF

C----              Copy shifted

   24 CHV(JP:JP) = CHV(JJ:JJ)
      JP = JP + 1
      JJ = JJ + 1
      IF (JJ.GT.JR)                GO TO 91
      IF (CHV(JJ:JJ).NE.' ')       GO TO 24

      IF (JJ.EQ.JR)                GO TO 91
      IF (CHV(JJ+1:JJ+1).NE.' ')   GO TO 24
      CHV(JP:JP) = ' '
      JP = JP + 1
      JJ = JJ + 1

   26 IF (JJ.GT.JR)                GO TO 91
      IF (CHV(JJ:JJ).NE.' ')       GO TO 24
      JJ = JJ + 1
      GO TO 26

   91 CHV(JP:JR) = ' '
   99 IF (CHV(JP-1:JP-1).EQ.' ')  JP = JP - 1
      NESLAT = JP
      NDSLAT = JP - JL
      RETURN
      END
