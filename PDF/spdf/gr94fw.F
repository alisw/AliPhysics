*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:28:08  mclareni
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
C      DOUBLE PRECISION FUNCTION GR94FW (X,S,AL,BE,AK,BK,A,B,C,D,E,ES)
       FUNCTION GR94FW (X,S,AL,BE,AK,BK,A,B,C,D,E,ES)
C
C+SEQ, IMPDP.
       IMPLICIT REAL (A - Z)
C
       if (X.LE.0.) then
          GR94FW=0.
          return
       endif
       LX =  LOG (1./X)
       GR94FW = (X**AK * (A + X * (B + X*C)) * LX**BK + S**AL
     1      *  EXP (-E +  SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
