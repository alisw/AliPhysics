#include "isajet/pilot.h"
      FUNCTION DBLPCM(A,B,C)
C          Calculate com momentum for A-->B+C with double precision.
C          Needed to fix bug on 32-bit machines at high energy.
C          Ver. 7.27: Rewrite order and then take abs value to be sure.
#include "isajet/itapes.inc"
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION DA,DB,DC,DVAL
#endif
C          Convert to double precision
      DA=A
      DB=B
      DC=C
      DVAL=(DA-(DB+DC))*(DA+(DB+DC))*(DA-(DB-DC))*(DA+(DB-DC))
C          Convert back to single precision
      VAL=DVAL
      DBLPCM=SQRT(ABS(VAL))/(2.*A)
      RETURN
      END
