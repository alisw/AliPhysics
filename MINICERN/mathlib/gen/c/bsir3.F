*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION BSIR3(X,NU)
C     IMPLICIT DOUBLE PRECISION (D)
#include "gen/def64.inc"
     +    D,DBSIR3,DBLE,DBSKR3,DEBIR3,DEBKR3

      SROUND(D)=D+(D-SNGL(D))
      BSIR3=SROUND(DBSIR3(DBLE(X),NU))
      RETURN

      ENTRY BSKR3(X,NU)
      BSKR3=SROUND(DBSKR3(DBLE(X),NU))
      RETURN

      ENTRY EBSIR3(X,NU)
      EBSIR3=SROUND(DEBIR3(DBLE(X),NU))
      RETURN

      ENTRY EBSKR3(X,NU)
      EBSKR3=SROUND(DEBKR3(DBLE(X),NU))
      RETURN
      END
#endif
