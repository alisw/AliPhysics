*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VMAX (A,N)
C
C CERN PROGLIB# F121    VMAX            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      FMAX = A(1)
      DO 9 I= 1,N
      FMAX = MAX (A(I),FMAX)
    9 CONTINUE
      VMAX = FMAX
      RETURN
  100 VMAX = 0.
      RETURN
      END
