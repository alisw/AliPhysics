*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (!defined(CERNLIB_CDC))&&(!defined(CERNLIB_VAX)||defined(CERNLIB_FORTRAN))
      SUBROUTINE NRAN(VECTOR,N)
      DIMENSION VECTOR(N)
      DO 100 I=1,N
      VECTOR(I) = RNDM(I)
  100 CONTINUE
      RETURN
      ENTRY NRANIN (V)
      CALL RDMIN(V)
      RETURN
      ENTRY NRANUT (V)
      CALL RDMOUT(V)
      RETURN
      END
#endif
