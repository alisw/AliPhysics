*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION PROD2(A,B,MD,M,N,I,VAL)
      DIMENSION A(MD,1),B(1)
      DOUBLE PRECISION SUM
C
      SUM=VAL
      DO 5  J=M,N
 5    SUM=SUM+DBLE(A(I,J))*DBLE(B(J))
      PROD2=SUM
      RETURN
      END
