*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION RGQUAD(F,A,B,N)

      DIMENSION X(1),W(1)
      EXTERNAL F

      CALL D107R1(1,F,A,B,N,X,W)
      RGQUAD=X(1)
      RETURN
      END
