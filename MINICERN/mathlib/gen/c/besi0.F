*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION BESI0(X)
      IMPLICIT DOUBLE PRECISION (D)

      SROUND(D)=D+(D-SNGL(D))
      BESI0=SROUND(DBESI0(DBLE(X)))
      RETURN

      ENTRY BESI1(X)
      BESI1=SROUND(DBESI1(DBLE(X)))
      RETURN

      ENTRY EBESI0(X)
      EBESI0=SROUND(DEBSI0(DBLE(X)))
      RETURN

      ENTRY EBESI1(X)
      EBESI1=SROUND(DEBSI1(DBLE(X)))
      RETURN

      ENTRY BESK0(X)
      BESK0=SROUND(DBESK0(DBLE(X)))
      RETURN

      ENTRY BESK1(X)
      BESK1=SROUND(DBESK1(DBLE(X)))
      RETURN

      ENTRY EBESK0(X)
      EBESK0=SROUND(DEBSK0(DBLE(X)))
      RETURN

      ENTRY EBESK1(X)
      EBESK1=SROUND(DEBSK1(DBLE(X)))
      RETURN
      END
#endif
