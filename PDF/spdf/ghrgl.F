*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:29  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION GHRGL(Y,SQ)
C
C  STRUCTURE FUNCTION GLUON FROM A PROTON
      ASQ=ABS(SQ)
      SI=1.
      IF (ASQ.GT.0.) SI=SQ/ABS(SQ)
      A=0.9243+2.51*SI*ASQ**0.5
      B=8.558-9.227*SI*ASQ**0.3-0.655*SI*ASQ**1.5
      C=53.57-68.78*SI*ASQ**0.3+19.3*SQ
      D=6.0+1.454*SQ
      E=11.29*SQ**2
      F=41.24+50.71*SQ
      GHRGL =A*(1.0+B*Y+C*Y**2)*(1.0-Y)**D+E* EXP(-F*Y)
      RETURN
      END
