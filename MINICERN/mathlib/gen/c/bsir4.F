*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION BSIR4(X,NU)
      IMPLICIT DOUBLE PRECISION (D)

      SROUND(D)=D+(D-SNGL(D))
      BSIR4=SROUND(DBSIR4(DBLE(X),NU))
      RETURN

      ENTRY BSKR4(X,NU)
      BSKR4=SROUND(DBSKR4(DBLE(X),NU))
      RETURN

      ENTRY EBSIR4(X,NU)
      EBSIR4=SROUND(DEBIR4(DBLE(X),NU))
      RETURN

      ENTRY EBSKR4(X,NU)
      EBSKR4=SROUND(DEBKR4(DBLE(X),NU))
      RETURN
      END
#endif
