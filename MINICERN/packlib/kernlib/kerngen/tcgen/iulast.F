*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IULAST (ITEXT,IC,N)
C
C CERN PROGLIB# V305    IULAST          .VERSION KERNFOR  4.19  880728
C ORIG. 01/08/69
C
      DIMENSION IC(80)
C
      J =    N
      IF (J.EQ.0) GO TO 20
C
      DO 17  J=N,1,-1
      IF (IC(J).NE.ITEXT)  GO TO 20
   17 CONTINUE
      J = 0
C
   20 IULAST= J
      RETURN
      END
