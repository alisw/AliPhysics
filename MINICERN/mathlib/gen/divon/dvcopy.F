*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DVCOPY(N,A,B)
      INTEGER N
      DOUBLE PRECISION A(N), B(N)
      INTEGER I
      DO 10 I=1,N
      B(I)=A(I)
 10   CONTINUE
      RETURN
      END
