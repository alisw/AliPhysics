*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VMIN (A,N)
C
C CERN PROGLIB# F121    VMIN            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      FMIN= A(1)
      DO 9 I= 1,N
      FMIN = MIN (A(I),FMIN)
    9 CONTINUE
      VMIN= FMIN
      RETURN
  100 VMIN = 0.
      RETURN
      END
