*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:45  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE U100DV (J,K,L)
      DIMENSION  J(40), K(40), L(40)
      DO 1 LL=1,40
    1 L(LL)=J(LL)-K(LL)
      RETURN
      END
