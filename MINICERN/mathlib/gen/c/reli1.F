*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:10  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RELI1(X,AKP)
      IMPLICIT DOUBLE PRECISION (D)

      SROUND(D)=D+(D-SNGL(D))
      RELI1=SROUND(DELI1(DBLE(X),DBLE(AKP)))
      RETURN

      ENTRY RELI2(X,AKP,A,B,MODE)
      RELI2=SROUND(DELI2(DBLE(X),DBLE(AKP),DBLE(A),DBLE(B),MODE))
      RETURN

      ENTRY RELI3(X,AKP,P)
      RELI3=SROUND(DELI3(DBLE(X),DBLE(AKP),DBLE(P)))
      RETURN
      END
#endif
