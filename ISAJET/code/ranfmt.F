#include "isajet/pilot.h"
      SUBROUTINE RANFMT
C
C          Get RANF seed and translate it to a character variable
C          to ensure exactly the same seed with a formatted read.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/seed.inc"
#if defined(CERNLIB_SINGLE)
*      REAL SEED
#endif
#if defined(CERNLIB_DOUBLE)
*      DOUBLE PRECISION SEED
#endif
*      CALL RANFGT(SEED)
*      WRITE(XSEED,'(E24.15)') SEED
*      READ(XSEED,'(E24.15)') SEED
*      CALL RANFST(SEED)
      RETURN
      END
