*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:32  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
*
C      DOUBLE PRECISION FUNCTION GR94FW (X,S,AL,BE,AK,BK,A,B,C,D,E,ES)
       FUNCTION GR94FW (X,S,AL,BE,AK,BK,A,B,C,D,E,ES)
C
C+SEQ, IMPDP.
       IMPLICIT REAL (A - Z)
C
       LX =  LOG (1./X)
       GR94FW = (X**AK * (A + X * (B + X*C)) * LX**BK + S**AL
     1      *  EXP (-E +  SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
