*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:31  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:58  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION AFCPLU(X,Q2)
      CMS=1.5**2
      BETS=1-4.*CMS*X/(1.-X)/Q2
      IF(BETS.LE..0) THEN
         AFCPLU=.0
         RETURN
      ENDIF
      BETA=SQRT(BETS)
      CPLU=(8.*X*(1.-X)-1.-4.*CMS*X*(1.-X)/Q2)*BETA
      CAU=X**2+(1.-X)**2+4.*CMS*X*(1.-3.*X)/Q2-8.*CMS**2*X**2/Q2**2
      CPLU=CPLU+CAU* LOG((1.+BETA)/(1.-BETA))
      AFCPLU=3.*(4./9.)*CPLU*X/(3.1415*137.)
  1   RETURN
      END
