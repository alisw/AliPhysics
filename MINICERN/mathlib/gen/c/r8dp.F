*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R8(Z,ACC)
      COMPLEX*16 C309R8,Z
      DOUBLE PRECISION ACC
      DOUBLE PRECISION X,Y,AX,AY,A

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      X=DREAL(Z)
      Y=DIMAG(Z)
      AX=ABS(X)
      AY=ABS(Y)
      A=5*ACC*(AX+AY)
      IF(AX .LT. A) X=0
      IF(AY .LT. A) Y=0
      C309R8=DCMPLX(X,Y)
      RETURN
      END
#endif
