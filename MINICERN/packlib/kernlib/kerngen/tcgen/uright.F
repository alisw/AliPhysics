*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE URIGHT (IC,JL,JR)
C
C CERN PROGLIB# M410    URIGHT          .VERSION KERNFOR  4.08  840613
C ORIG. 27/05/67
C
      COMMON /SLATE/NC,JP, DUMMY(38)
      DIMENSION IC(80)
      DATA  IBLANK / 4H     /
C
C
      JP= JR
      JT= JR
      NC= 0
C
   12 IF (JT.LT.JL)  RETURN
      IF (IC(JT).EQ.IBLANK)  GO TO 19
      IF (JP.EQ.JT)  GO TO 17
      IC (JP)= IC(JT)
      IC(JT)= IBLANK
   17 NC= NC+1
      JP= JP-1
   19 JT= JT-1
      GO TO 12
      END
