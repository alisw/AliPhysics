*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:55  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRATS (A,S,B,M,N)
C
C CERN PROGLIB# F112    TRATS           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),S(*),B(*)
C
      IB = 0
      INDS = 0
      I = 0
C
    5 INDS = INDS + I
      IB = I + 1
C
      DO 40 J=1,M
      IA = J
      IS = INDS
      SUM = 0.
      K = 0
C
   10 IF (K.GT.I)             GO TO 20
      IS = IS + 1
        GO  TO                30
   20 IS = IS + K
   30 SUM = SUM + A(IA)*S(IS)
      IA = IA + M
      K = K + 1
      IF (K.LT.N)             GO TO 10
C
      B(IB) = SUM
   40 IB = IB + N
      I = I + 1
      IF (I.LT.N)             GO TO 5
C
      RETURN
      END
