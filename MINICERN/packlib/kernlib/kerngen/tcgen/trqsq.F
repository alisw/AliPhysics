*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRQSQ (Q,S,R,M)
C
C CERN PROGLIB# F112    TRQSQ           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION S(*),Q(*),R(*)
C
      IMAX = (M*M+M)/2
      CALL VZERO (R,IMAX)
      INDS = 0
      I = 0
C
    5 INDS = INDS + I
      IR = 0
      INDQ = 0
      J = 0
C
   10 INDQ = INDQ + J
      IS = INDS
      IQ = INDQ
      SUM = 0.
      K = 0
C
   15 IF (K.GT.I)             GO TO 20
      IS = IS + 1
        GO  TO                30
   20 IS = IS + K
   30 IF (K.GT.J)             GO TO 40
      IQ = IQ + 1
        GO  TO                50
   40 IQ = IQ + K
   50 SUM = SUM + S(IS)*Q(IQ)
      K = K + 1
      IF (K.LT.M)             GO TO 15
      IQQ = INDS
      L = 0
C
   60 IR = IR + 1
      IF (L.GT.I)             GO TO 70
      IQQ = IQQ + 1
        GO  TO                80
   70 IQQ = IQQ + L
   80 R(IR) = R(IR) + Q(IQQ)*SUM
      L = L + 1
      IF (L.LE.J)             GO TO 60
      J = J + 1
      IF (J.LT.M)             GO TO 10
      I = I + 1
      IF (I.LT.M)             GO TO 5
C
      RETURN
      END
