#include "isajet/pilot.h"
      REAL FUNCTION KKGF3(S,T,M2)
      REAL S,T,M2
      REAL XG,YG
      XG = T/S
      YG = M2/S
      KKGF3 = 
     &   ( 1. + 2.*XG + 3.*XG**2 + 2.*XG**3 + XG**4 - 2.*YG*(1.+XG**3)
     &     + 3.*YG**2*(1.+XG**2) - 2.*YG**3*(1.+XG) + YG**4 ) /
     &   ( XG*(YG-1.-XG) )
      RETURN
      END
