*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUNEXT (IVEC,N)
C
C CERN PROGLIB# M410    IUNEXT          .VERSION KERNFOR  4.08  840613
C ORIG. 27/05/67
C
      DIMENSION IVEC(80)
      DATA  IBLANK /4H     /
C
C
      J= N
   10 IF (IVEC(J).NE.IBLANK)  GO TO 12
      J= J + 1
      GO TO 10
C
   12 IUNEXT= J
      RETURN
      END
