*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSQMCH (CHAPAR,CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CSQMCH          .VERSION KERNFOR  4.24  900424
C ORIG. 13/03/90, Jamie Shiers
C
C-    Squeeze multiple occurrences of the character CHA in CHV(JL:JR),
C-    shifting left

      DIMENSION    JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT, DUMMY(38)
      CHARACTER    CHV*(*)
      CHARACTER*1  CHAPAR, CHA

C----              Find the first occurrence of CHA

      CHA = CHAPAR
      JL  = JLP(1)
      JR  = JRP(1)

      JP = JR + 1
      JJ = JL
   12 IF (JJ.GE.JP)                GO TO 99
      IF (CHV(JJ:JJ).NE.CHA)   THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF

C--       is it multiple ?

      JJ = JJ + 1
      IF (JJ.GE.JP)                GO TO 99
      IF (CHV(JJ:JJ).NE.CHA)       GO TO 12

      JP = JJ
      JJ = JJ + 1
   16 IF (JJ.GT.JR)                GO TO 91
      IF (CHV(JJ:JJ).EQ.CHA)   THEN
          JJ = JJ + 1
          GO TO 16
        ENDIF

C----              Copy shifted

   24 CHV(JP:JP) = CHV(JJ:JJ)
      JP = JP + 1
      JJ = JJ + 1
      IF (JJ.GT.JR)                GO TO 91
      IF (CHV(JJ-1:JJ-1).NE.CHA)   GO TO 24

   26 IF (CHV(JJ:JJ).NE.CHA)       GO TO 24
      JJ = JJ + 1
      IF (JJ.LE.JR)                GO TO 26

   91 CHV(JP:JR) = ' '
   99 NESLAT = JP
      NDSLAT = JP - JL
      RETURN
      END
