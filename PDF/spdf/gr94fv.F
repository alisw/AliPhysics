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
*...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
*
C      DOUBLE PRECISION FUNCTION GR94FV (X, N, AK, BK, A, B, C, D)
       FUNCTION GR94FV (X, N, AK, BK, A, B, C, D)
C
C+SEQ, IMPDP.
       IMPLICIT REAL (A - Z)
C
       DX =  SQRT (X)
       GR94FV = N * X**AK * (1.+ A*X**BK + X * (B + C*DX)) * (1.- X)**D
       RETURN
       END
