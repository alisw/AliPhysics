*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE NORMCO(A,B,UA,UB,SIGMAA,SIGMAB,RO)
      Y=RNDM(-1.)
      Z=RNDM(-1.)
      X=6.2831853*Z
      AA=SQRT(-2.*LOG(Y))
      A=UA+SIGMAA*AA*(SQRT(1.-RO**2)*COS(X)+RO*SIN(X))
      B=UB+SIGMAB*AA*SIN(X)
      RETURN
      END
