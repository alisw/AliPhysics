#include "isajet/pilot.h"
      COMPLEX*16 FUNCTION SSB0(QSQ,M1,M2)
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/ssinf.inc"
      COMPLEX*16 SSF0
      REAL QSQ,M1,M2
      SSB0=XLAM*(1.D0,0.D0)-SSF0(QSQ,M1,M2)
      RETURN
      END
