*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRAL (A,U,B,M,N)
C
C CERN PROGLIB# F112    TRAL            .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),U(*),B(*)
C
      IB = 1
C
      DO 40 I=1,M
      INDU = 0
C
      DO 40 J=1,N
      INDU = INDU + J
      IA = IB
      IU = INDU
      SUM = 0.
C
      DO 30 K=J,N
      SUM = SUM + A(IA)*U(IU)
      IA = IA + 1
   30 IU = IU + K
      B(IB) = SUM
   40 IB = IB + 1
C
      RETURN
      END
