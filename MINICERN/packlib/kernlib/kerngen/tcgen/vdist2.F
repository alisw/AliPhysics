*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VDIST2 (X,Y,N)
C
C CERN PROGLIB# F121    VDIST2          .VERSION KERNFOR  3.01  761204
C ORIG. 26/09/76
C
      DIMENSION X(*), Y(*)
C
      XX = 0.
      IF (N.LE.0)            GO TO 24
C
      DO 9 J=1,N
    9 XX = XX + (X(J)-Y(J)) **2
   24 VDIST2 = XX
      RETURN
      END
