#include "isajet/pilot.h"
      SUBROUTINE RANFGT(SEED)
C
C          Get seed for RANF() in real or double precision SEED.
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
      CALL RANGET(SEED)
#endif
#if defined(CERNLIB_CRAY)
      INTEGER ISEED,RANGET,IDUMMY
      ISEED=RANGET(IDUMMY)
      SEED=ISEED
#endif
      RETURN
      END
