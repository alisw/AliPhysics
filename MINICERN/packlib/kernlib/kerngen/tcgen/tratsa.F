*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRATSA (A,S,R,M,N)
C
C CERN PROGLIB# F112    TRATSA          .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION S(*),A(*),R(*)
C
      IMAX = (M*M+M)/2
      CALL VZERO (R,IMAX)
      MN = M*N
      IND = 0
      I = 0
C
    5 IND = IND + I
      IR = 0
C
      DO 40 J=1,M
      IS = IND
      IA = J
      SUM = 0.
      K = 0
C
   15 IF (K.GT.I)             GO TO 20
      IS = IS + 1
        GO  TO                30
   20 IS = IS + K
   30 SUM = SUM + S(IS)*A(IA)
      IA = IA + M
      K = K + 1
      IF (K.LT.N)             GO TO 15
      IAA = I*M
C
      DO 40 K=1,J
      IAA = IAA + 1
      IR = IR + 1
   40 R(IR) = R(IR) + SUM*A(IAA)
      I = I + 1
      IF (I.LT.N)             GO TO 5
C
      RETURN
      END
