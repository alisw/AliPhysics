*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:14  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE SBYTPK (I,X,JTH,NBITS)
C
C CERN PROGLIB# M422    SBYTPK          .VERSION KERNFOR  4.06  830526
C ORIG. 10/02/71
C
      DIMENSION X(99), NBITS(3)
C
C
      NBIT  = NBITS(1)
      IF (NBIT.GT.0)  GO TO 7
      NBIT  = 1
#if defined(CERNLIB_B32)
      INWORD = 32
#endif
#if defined(CERNLIB_B36)
      INWORD = 36
#endif
#if defined(CERNLIB_B48)
      INWORD = 48
#endif
#if defined(CERNLIB_B60)
      INWORD = 60
#endif
#if defined(CERNLIB_B64)
      INWORD = 64
#endif
      GO TO 11
C
    7 INWORD= NBITS(2)
   11 JWORD= (JTH-1)/INWORD + 1
      JPOS = JTH - (JWORD-1)*INWORD
      LOC  = (JPOS-1)*NBIT + 1
      CALL SBYT (I,X(JWORD),LOC,NBIT)
      RETURN
      END
