*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRSMUL (G,GI,N)
C
C CERN PROGLIB# F112    TRSMUL          .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION G(*),GI(*)
C
      IND = 1
      LPIV = 0
      DO 30 I=1,N
      LPIV = LPIV + I
      DO 30 J=1,I
      LVER = LPIV
      LHOR = IND
      SUM = 0.
      DO 20 K=I,N
      SUM = SUM + G(LVER)*G(LHOR)
      LHOR = LHOR + K
   20 LVER = LVER + K
      GI(IND) = SUM
   30 IND = IND + 1
C
      RETURN
      END
