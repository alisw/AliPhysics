*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:31  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION GHRXS(Y,SQ)
C
C  STRUCTURE FUNCTION  S-QUARK FROM A PROTON
      ASQ=ABS(SQ)
      SI=1.
      IF ( ASQ.GT.0.) SI=SQ/ABS(SQ)
      A=0.0625+0.1132*SI*ASQ**1.3
      B=12.64*SQ-51.70*SI*ASQ**1.8+38.02*SQ**2
      C=4.448*SQ
      D=7.0+1.562*SQ
      E=0.3081*SI*ASQ**2.5
      F=47.24+67.91*SQ
      GHRXS=A*(1.0+B*Y+C*Y**2)*(1.0-Y)**D+E* EXP(-F*Y)
      RETURN
      END
