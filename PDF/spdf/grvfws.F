*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:39  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C                                                                       GRV02120
       FUNCTION GRVFWS (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       LX = ALOG (1./X)
       IF (S .LE. ST) THEN
          GRVFWS = 0.0
       ELSE
          GRVFWS = (S-ST)**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1          * EXP (-E + SQRT (ES * S**BE * LX))
       END IF
       RETURN
       END
