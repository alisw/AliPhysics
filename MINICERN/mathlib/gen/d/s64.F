*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE S64(A,B)
      REAL A(32),B(32),C(32)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      DO 10 I=3,15,2
      A1 = A(I-1)-A(I+1)
      A2 = A(33-I)-A(35-I)
      C1R = A1+A2
      C2R = A1-A2
      C1I = A(I)-A(34-I)
      C2I = A(I)+A(34-I)
      C1 = D(8*I-9)
      C2 = D(8*I-8)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(34-I) = C1R-C3R
   10 C(35-I) = C3I-C1I
      C(1) = 2.*(A(32)-A(2))
      C(2) = -2.*(A(32)+A(2))
      C(17) = 2.*(A(16)-A(18))
      C(18) = 2.*A(17)
      CALL I16(C,C(2),2,C,C(2),2)
      B(17) = 0.5*C(17)
      DO 20 I=2,16
      E1 = 0.5*(C(I)-C(34-I))
      E2 = E(8*I-8)*(C(I)+C(34-I))
      B(I) = E2-E1
   20 B(34-I) = E2+E1
      RETURN
      END
