#include "isajet/pilot.h"
      REAL FUNCTION KKGF2(S,T,M2)
      REAL S,T,M2
      REAL XG,YG
      XG = T/S
      YG = M2/S
      KKGF2 =
     &   ( -4.*XG*(1.+XG**2) + YG*(1.+XG)*(1.+8.*XG+XG**2) - 
     &      3.*YG**2*(1.+4.*XG+XG**2) + 4.*YG**3*(1.+XG) - 2.*YG**4 ) /
     &   ( XG*(YG-1.-XG) )
      RETURN
      END
