*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRLA (U,A,B,M,N)
C
C CERN PROGLIB# F112    TRLA            .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),U(*),B(*)
C
      IB = M*N
      IPIV = (M*M+M)/2
C
   20 IA = IB
      IU = IPIV
C
      SUM = 0.
   30 SUM = SUM + A(IA)*U(IU)
      IU = IU - 1
      IA = IA - N
      IF (IA.GT.0)            GO TO 30
C
      B(IB) = SUM
      IB = IB - 1
      IF (IA.GT.1-N)          GO TO 20
C
      IPIV = IU
      IF (IU.GT.0)            GO TO 20
C
      RETURN
      END
