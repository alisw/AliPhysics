*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRAAT  (A,S,M,N)
C
C CERN PROGLIB# F112    TRAAT           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),S(*)
C
      IA = 0
      IS = 0
C
      DO 30 I=1,M
      IPIV = IA
      IPIVN = IPIV + N
      IAT = 0
C
      DO 30 J=1,I
      IA = IPIV
C
      SUM = 0.
   20 IA = IA + 1
      IAT = IAT + 1
      SUM = SUM + A(IA)*A(IAT)
      IF (IA.LT.IPIVN)        GO TO 20
C
      IS = IS + 1
   30 S(IS) = SUM
C
      RETURN
      END
