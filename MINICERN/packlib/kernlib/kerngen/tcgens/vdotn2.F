*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VDOTN2 (X,Y,N)
C
C CERN PROGLIB# F121    VDOTN2          .VERSION KERNFOR  1.0   720620
C ORIG. 20/06/72
C
      DIMENSION X(*),Y(*)
C
      IF (N.LE.0)  GO TO 100
      XX = 0.
      XY = 0.
      YY = 0.
C
      DO 9 J=1,N
      XX = XX + X(J)*X(J)
      XY = XY + X(J)*Y(J)
    9 YY = YY + Y(J)*Y(J)
C
      VDOTN2= XY**2/(XX*YY)
      RETURN
  100 VDOTN2= 0.
      RETURN
      END
