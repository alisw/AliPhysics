*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      COMPLEX FUNCTION RGAM(Z,A,B)

      COMPLEX Z,A,B,CLOGAM

      RGAM=EXP(CLOGAM(Z+A)-CLOGAM(Z+B))
      RETURN
      END
