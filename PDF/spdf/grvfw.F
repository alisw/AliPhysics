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
C
       FUNCTION GRVFW (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT REAL (A - Z)
       LX = ALOG (1./X)
       GRVFW = (X**AK * (AG + X * (BG + X*C)) * LX**BK + S**AL
     1      * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
