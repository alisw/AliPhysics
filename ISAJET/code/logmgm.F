#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGMGM(IERR)
C
C       Set and check limits for multijet mass
C
C       Note we use the convention that not setting an upper limit
C       gives a fixed variable, even though that currently is not
C       implemented in N-jet phase space.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
#include "isajet/mglims.inc"
C
      REAL UNDEF
      INTEGER IERR
      DATA UNDEF/-.9E9/
C
      LOGMGM=.TRUE.
C
      IF(EHMGMN.LT.UNDEF.OR.EHMGMX.LT.UNDEF) THEN
        LOGMGM=.FALSE.
      ENDIF
C
      RETURN
      END
