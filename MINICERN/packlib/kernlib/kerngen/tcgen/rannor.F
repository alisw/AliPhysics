*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE RANNOR (A,B)
C
C CERN PROGLIB# V100    RANNOR          .VERSION KERNFOR  4.18  880425
C ORIG. 18/10/77
C
      Y = RNDM(-1.)
      IF (Y.EQ.0.)  Y = RNDM(-1.)
      Z = RNDM(-1.)

      X = 6.283185*Z
      R = SQRT (-2.0*LOG(Y))
      A = R*SIN (X)
      B = R*COS (X)
      RETURN
      END
