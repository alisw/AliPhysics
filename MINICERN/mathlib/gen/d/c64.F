*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE C64(A,B)
      REAL A(33),B(33),C(32)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      S2 = A(16)+A(18)
      DO 10 I=3,15,2
      S2 = S2+A(I-1)+A(35-I)
      A1 = A(I+1)-A(I-1)
      A2 = A(35-I)-A(33-I)
      C1R = A(I)+A(34-I)
      C2R = A(I)-A(34-I)
      C1I = A1-A2
      C2I = A1+A2
      C1 = -D(8*I-9)
      C2 = D(8*I-8)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(34-I) = C1R-C3R
   10 C(35-I) = C3I-C1I
      C(1) = A(1)+A(33)
      C(2) = A(1)-A(33)
      C(17) = 2.*A(17)
      C(18) = 2.*(A(16)-A(18))
      CALL D16(C,C(2),2,C,C(2),2)
      S2 = 2.*S2
      B(1) = C(1)+S2
      B(33) = C(1)-S2
      B(17) = C(17)
      DO 20 I=2,16
      E1 = 0.5*(C(I)+C(34-I))
      E2 = E(8*I-8)*(C(I)-C(34-I))
      B(I) = E1+E2
   20 B(34-I) = E1-E2
      RETURN
      END
