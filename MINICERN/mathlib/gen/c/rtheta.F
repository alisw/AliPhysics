*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:12  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RTHETA( K, X, Q )
      DOUBLE PRECISION D,DTHETA
C
C  RTHETA is set equal to the correctly rounded value of DTHETA.
C  (Version 9.04.1992)
C
C--Start.
      D = DTHETA( K, DBLE(X), DBLE(Q) )
      RTHETA = D + ( D - SNGL(D) )
      RETURN
      END
#endif
