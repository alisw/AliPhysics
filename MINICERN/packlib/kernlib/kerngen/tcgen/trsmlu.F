*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRSMLU (U,S,N)
C
C CERN PROGLIB# F112    TRSMLU          .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION U(*),S(*)
C
      IND = (N*N+N)/2
C
      DO 30 I=1,N
      LVER = IND
C
      DO 30 K=I,N
      LHOR = IND
      SUM = 0.
C
      DO 20 L=K,N
      SUM = SUM + U(LVER)*U(LHOR)
      LVER = LVER - 1
   20 LHOR = LHOR - 1
      S(IND) = SUM
   30 IND = IND - 1
C
      RETURN
      END
