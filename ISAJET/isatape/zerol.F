#include "isajet/pilot.h"
      SUBROUTINE ZEROL(Z,N)
C          SET N VALUES OF Z IN LEVEL2 TO ZERO
#include "isajet/itapes.inc"
      DIMENSION Z(N)
#if defined(CERNLIB_LEVEL2)
      LEVEL2,Z
#endif
      DO 1 I=1,N
    1 Z(I)=0
      RETURN
      END
