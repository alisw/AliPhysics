*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUFILA (IT,ICH,JL,JR)
C
C CERN PROGLIB# V304    IUFILA          .VERSION KERNFOR  4.18  880425
C ORIG. 25/04/88 JZ
C
      DIMENSION ICH(9)
C
      DO 12   JF=JR,JL,-1
      IF (ICH(JF).EQ.IT)     GO TO 24
   12 CONTINUE
      JF = JR+1
C
   24 IUFILA = JF
      RETURN
      END
