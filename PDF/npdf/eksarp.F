#include "pdf/pilot.h"
      function eksarp(x,v2,x1,x116,ptm,indpi)
#include "pdf/impdp.inc"
      dimension ptm(10)

      z=v2
      xx=x-x1
      if (x.le.x116) then
         eksarp=ptm(1)+ptm(2)*z+ptm(3)*z**2+ptm(4)*xx+ptm(5)*xx**2
         else
            z1= log(x116/x1)
            xx1=x116-x1
            ff0=ptm(1)+ptm(2)*z1+ptm(3)*z1**2+ptm(4)*xx1+ptm(5)*xx1**2
            xx16=x-x116
            z16= log(x/x116)
            qexp=19.d0+(indpi-1)*(indpi-2)*8.d0/2.d0
            eksarp=ff0+ptm(6)*xx16+ptm(7)*xx16**2+ptm(8)*xx16**3
     #             +ptm(9)*xx16**qexp+ptm(10)*z16
      endif
11    return
      end


