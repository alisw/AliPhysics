
*
* $Id$
*
* $Log$
* Revision 1.2  1997/04/08 14:39:00  mclareni
* A fourth term is needed for part of the range of X, Fred James
*
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C     This corresponds to PROBKL,IF=DOUBLE and PROBKL64,IF=-DOUBLE
      FUNCTION PROBKL(X)

      DIMENSION FJ(4),R(4)

      PARAMETER (PI = 3.14159 265D0)
      PARAMETER (W  = 2.50662 827D0)
      PARAMETER (C1 = -PI**2/8, C2 = 9*C1, C3 = 25*C1)

      DATA FJ /-2,-8,-18,-32/

      U=ABS(X)
      IF(U .LT. 0.2) THEN
       P=1.
      ELSEIF(U .LT. 0.755) THEN
       V=1/U**2
       P=1-W*(EXP(C1*V)+EXP(C2*V)+EXP(C3*V))/U
      ELSEIF(U .LT. 6.8116) THEN
       R(2)=0.
       R(3)=0.
       R(4)=0.
       V=U**2
       DO 1 J = 1,MAX(1,NINT(3/U))
    1  R(J)=EXP(FJ(J)*V)
       P=2*(R(1)-R(2)+R(3)-R(4))
CCC         PRINT '(35x,4e10.2)', (R(JJ),JJ=1,4)
      ELSE
       P=0
      ENDIF
      PROBKL=P
      RETURN
      END






