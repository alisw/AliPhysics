*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:33  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
*
C      DOUBLE PRECISION FUNCTION GR94FWS (X,S,AL,BE,AK,AG,B,D,E,ES)
       FUNCTION GR94FWS (X,S,AL,BE,AK,AG,B,D,E,ES)
C
C+SEQ, IMPDP.
       IMPLICIT REAL (A - Z)
C
       DX =  SQRT (X)
       LX =  LOG (1./X)
       GR94FWS = S**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1       *  EXP (-E +  SQRT (ES * S**BE * LX))
       RETURN
       END
