*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION PROD1(A,B,M,N)
      DIMENSION A(1),B(1)
      DOUBLE PRECISION SUM
C
      SUM=0.
      DO 5  I=M,N
 5    SUM=SUM+DBLE(A(I))*DBLE(B(I))
      PROD1=SUM
      RETURN
      END
