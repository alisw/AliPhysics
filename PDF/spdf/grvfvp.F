*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:38  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       FUNCTION GRVFVP (X, N, AK, AG, D)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       GRVFVP = N * X**AK * (1.+ AG*DX) * (1.- X)**D
       RETURN
       END
