*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:41  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       FUNCTION GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT REAL (A - Z)
       IF (S .LE. SF) THEN
          GRVGFS = 0.0
       ELSE
          SX = SQRT (X)
          LX = ALOG (1./X)
          DS = S - SF
          GRVGFS = (DS * X**AK * (AG + BG * SX + C * X**BK) + DS**AL
     1         * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       END IF
       RETURN
       END
