#include "isajet/pilot.h"
      FUNCTION EVOLMS(J,FUDGE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     Set evolution mass scale for parton J
C-
C-   Returned value  : maximum mass
C-
C-   Inputs  : 
C-     J    = index to PJSET array
C-     FUDGE= fudge factor
C-
C-   Created  16-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL    EVOLMS,FUDGE
      INTEGER J
#include "isajet/limevl.inc"
#include "isajet/jetset.inc"
#include "isajet/jetpar.inc"
C----------------------------------------------------------------------
C
      IF ( USELIM ) THEN
        EVOLMS=SQRT(PJSET(1,J)**2+PJSET(2,J)**2)*CONCUT
      ELSE
        EVOLMS=FUDGE*SQRT(QSQ)
      ENDIF
  999 RETURN
      END
