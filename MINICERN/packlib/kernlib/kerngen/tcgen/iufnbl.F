*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUFNBL (IV,JA,JE)
C
C CERN PROGLIB# M410    IUFNBL          .VERSION KERNFOR  3.06  780915
C ORIG. 15/09/78 JZ
C
      DIMENSION    IV(80)
      DATA  IBLANK /4H    /
C
C
      J  = JA
   12 IF (J.GT.JE)           GO TO 19
      IF (IV(J).NE.IBLANK)   GO TO 19
      J  = J + 1
      GO TO 12
C
   19 IUFNBL = J
      RETURN
      END
