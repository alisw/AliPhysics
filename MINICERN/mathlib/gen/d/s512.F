*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE S512(A,B)
      REAL A(256),B(256),C(256)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      DO 10 I=3,127,2
      A1 = A(I-1)-A(I+1)
      A2 = A(257-I)-A(259-I)
      C1R = A1+A2
      C2R = A1-A2
      C1I = A(I)-A(258-I)
      C2I = A(I)+A(258-I)
      C1 = D(I-2)
      C2 = D(I-1)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(258-I) = C1R-C3R
   10 C(259-I) = C3I-C1I
      C(1) = 2.*(A(256)-A(2))
      C(2) = -2.*(A(256)+A(2))
      C(129) = 2.*(A(128)-A(130))
      C(130) = 2.*A(129)
      CALL I128(C,C(2),2,C,C(2),2)
      B(129) = 0.5*C(129)
      DO 20 I=2,128
      E1 = 0.5*(C(I)-C(258-I))
      E2 = E(I-1)*(C(I)+C(258-I))
      B(I) = E2-E1
   20 B(258-I) = E2+E1
      RETURN
      END
