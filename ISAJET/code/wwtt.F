#include "isajet/pilot.h"
      FUNCTION WWTT(T,U,T1,U1,T3,U3)
C          DECAY DISTRIBUTION FOR W+ W- PAIRS FROM SCHOONSCHIP(1980).
C          SQUARE OF T GRAPH.
#include "isajet/itapes.inc"
#include "isajet/wwpar.inc"
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION WWTT
      DOUBLE PRECISION T,U,T1,U1,T3,U3
#endif
      WWTT=
     1(+CQ**2*(-512.*T1**2*T3**2)
     1 +CQ**2*T*(-512.*T1*U1*T3-512.*T1*T3*S13+512.*T1*T3**2)
     1 +CQ**2*T*WM2*(-512.*U1*T3-512.*T3*S13+512.*T3**2)
     1 +CQ**2*T**2*(512.*U1*T3+512.*T3*S13)
     1 +CQ**2*WM2*(-1024.*T1*T3**2)
     1 +CQ**2*WM2**2*(-512.*T3**2))/T**2
      WWTT=2.*WWTT
      RETURN
      END
