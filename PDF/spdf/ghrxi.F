*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:30  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION GHRXI(Y,SQ)
C
C  STRUCTURE FUNCTION ANTIQUARK FROM A PROTON
      ASQ=ABS(SQ)
      SI=1.
      IF (ASQ.GT.0.) SI=SQ/ABS(SQ)
      A=0.25+0.088*SI*ASQ**1.3
      B=0.8128*SQ-2.003*SI*ASQ**1.8+0.0831*SQ**2
      C=3.97*SQ
      D=7.0+1.666*SQ
      E=0.2487*SI*ASQ**2.5
      F=27.80+59.68*SQ
      GHRXI=A*(1.0+B*Y+C*Y**2)*(1.0-Y)**D+E* EXP(-F*Y)
      RETURN
      END
