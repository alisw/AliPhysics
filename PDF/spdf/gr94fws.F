*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:28:09  mclareni
* Import version 8.04
*
* Revision 1.2  2000/04/15 14:10:46  plothow
* Version 8.03: Nuclear PDFs introduced
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
       if (X.LE.0.) then
         GR94FWS = 0.
         return
       endif
       DX =  SQRT (X)
       LX =  LOG (1./X)
       GR94FWS = S**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1       *  EXP (-E +  SQRT (ES * S**BE * LX))
       RETURN
       END
