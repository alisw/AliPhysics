*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VDOTN (X,Y,N)
C
C CERN PROGLIB# F121    VDOTN           .VERSION KERNFOR  2.09  751101
C ORIG. 01/07/71
C
      DIMENSION X(*), Y(*)
C
      IF (N.LE.0) GO TO 100
      XX= 0.
      XY= 0.
      YY= 0.
C
      DO 9 J=1,N
      XX = XX + X(J)*X(J)
      XY = XY + X(J)*Y(J)
    9 YY = YY + Y(J)*Y(J)
C
      VDOTN= XY / SQRT(XX*YY)
      IF (ABS (VDOTN).LT.1.) RETURN
      VDOTN= SIGN (1.,VDOTN)
      RETURN
  100 VDOTN = 0.
      RETURN
      END
