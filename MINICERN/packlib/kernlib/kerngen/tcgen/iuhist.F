*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUHIST (X,PAR,SPILL)
C
C CERN PROGLIB# Y203    IUHIST          .VERSION KERNFOR  3.15  820113
C ORIG. 21/08/72 P.RASTL
C
      DIMENSION PAR(3)
      EQUIVALENCE (NX,A)
      LOGICAL SPILL
C
C
      A = PAR(1)
      SPILL=.FALSE.
#if defined(CERNLIB_B36M)
      AN = (X-PAR(3))*PAR(2) + 1.000001
#endif
#if !defined(CERNLIB_B36M)
      AN = (X-PAR(3))*PAR(2) + 1.00001
#endif
      IF (AN.LT.1.) GO TO 20
      IF (AN.GE.FLOAT(NX+1))  GO TO 18
      N=AN
      IF (N.LT.0)            GO TO 18
      IUHIST= N
      RETURN
C
C--                OVERFLOW
   18 SPILL= .TRUE.
      IUHIST=NX+1
      RETURN
C
C--                UNDERFLOW
   20 SPILL=.TRUE.
      IUHIST=0
      RETURN
      END
