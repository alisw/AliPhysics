*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:55  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRATA  (A,R,M,N)
C
C CERN PROGLIB# F112    TRATA           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),R(*)
C
      MN = M*N
      IR = 0
C
      DO 30 I=1,M
C
      DO 30 J=1,I
      IA = I
      IAT = J
C
      SUM = 0.
   20 SUM = SUM + A(IA)*A(IAT)
      IA = IA + M
      IAT = IAT + M
      IF (IA.LE.MN)           GO TO 20
C
      IR = IR + 1
   30 R(IR) = SUM
C
      RETURN
      END
