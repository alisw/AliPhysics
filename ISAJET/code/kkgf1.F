#include "isajet/pilot.h"
      REAL FUNCTION KKGF1(S,T,M2)
      REAL S,T,M2
      REAL XG,YG
      XG = T/S
      YG = M2/S
      KKGF1 =
     &   ( -4.*XG*(1.+XG)*(1.+2.*XG+2.*XG**2) 
     &     + YG*(1.+6.*XG+18.*XG**2+16.*XG**3) 
     &     - 6.*YG**2*XG*(1.+2.*XG) + YG**3*(1.+4.*XG) ) / 
     &   ( XG*(YG-1.-XG) )
      RETURN
      END
