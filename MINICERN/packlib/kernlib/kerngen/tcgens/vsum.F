*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VSUM (A,N)
C
C CERN PROGLIB# F121    VSUM            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      XX= 0.
      IF (N.LE.0)  GO TO 100
      DO 9 I= 1,N
    9 XX= XX + A(I)
C
  100 VSUM= XX
      RETURN
      END
