*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:55  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRALT (A,U,B,M,N)
C
C CERN PROGLIB# F112    TRALT           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),U(*),B(*)
C
      IB = M*N
      INDU = (N*N+N)/2
C
   20 IU = INDU
C
      DO 40 J=1,N
      IA = IB
      SUM = 0.
C
      DO 30 K=J,N
      SUM = SUM + A(IA)*U(IU)
      IA = IA - 1
   30 IU = IU - 1
      B(IB) = SUM
   40 IB = IB - 1
      IF (IB.GT.0)            GO TO 20
C
      RETURN
      END
