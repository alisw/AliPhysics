*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION NIC311(X)
C     X=1/2 OF 4N+1 SET NIC311=1
C     X=1/2 OF 4N+2 SET NIC311=2
C     X=1/2 OF 4N+3 SET NIC311=3
C     X=1/2 OF 4N+4 SET NIC311=4
C     OTHERWISE SET NIC311=5
      Y=X+X
      N=INT(Y)
      IF(Y-N)9,10,9
    9 NIC311=5
      RETURN
   10 NIC311=MOD(N,4)
      IF(NIC311.LE.0) NIC311=4+NIC311
      RETURN
      END
