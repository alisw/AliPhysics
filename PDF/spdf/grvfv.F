*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:37  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
C...FUNCTIONAL FORMS FOR HO AND LO PARAMETRIZATIONS :
       FUNCTION GRVFV (X, N, AK, AG, B, D)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       GRVFV = N * X**AK * (1.+ AG*DX + B*X) * (1.- X)**D
       RETURN
       END
