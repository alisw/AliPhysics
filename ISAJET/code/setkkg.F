#include "isajet/pilot.h"
      SUBROUTINE SETKKG
C
C          Set the standard KKG parameters in /KKGRAVI/.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/kkgrav.inc"
#include "isajet/const.inc"
C
      REAL DIM2,GMMA,GAMMA
      EXTERNAL GAMMA
C          Calculate D-surface:
      DIM2 = (NEXTRAD*1.0)/2.
      GMMA = GAMMA(DIM2)
      SURFD = (2.*PI**DIM2) / GMMA
      KKGSD = SURFD / (MASSD**(NEXTRAD+2))
      RETURN
      END
