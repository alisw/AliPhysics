#include "pdf/pilot.h"
C-----------------------------------------------------------------------
#if defined(CERNLIB_SINGLE)
      function grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
#endif
#if defined(CERNLIB_DOUBLE)
      double precision function grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,
     +                                ge,gep)
#endif
C
#include "pdf/impdp.inc"
C
      grsf2=(s*x**a*(ga+gb*sqrt(x)+gc*x**b)+
     +      s**alp*exp(-ge+sqrt(gep*s**bet*log(1.d0/x))))*
     +      (1.d0-x)**gd
      return
      end

