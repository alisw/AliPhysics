#include "isajet/pilot.h"
      SUBROUTINE RANFST(SEED)
C
C          Set seed for RANF() from real or double precision SEED
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#if defined(CERNLIB_SINGLE)
      REAL SEED
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION SEED
#endif
#if defined(CERNLIB_RANFCALL)
      CALL RANSET(SEED)
#endif
#if defined(CERNLIB_CRAY)
      INTEGER ISEED
      ISEED=SEED
      CALL RANSET(ISEED)
#endif
      RETURN
      END
