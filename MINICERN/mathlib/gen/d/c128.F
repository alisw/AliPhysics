*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE C128(A,B)
      REAL A(65),B(65),C(64)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      S2 = A(32)+A(34)
      DO 10 I=3,31,2
      S2 = S2+A(I-1)+A(67-I)
      A1 = A(I+1)-A(I-1)
      A2 = A(67-I)-A(65-I)
      C1R = A(I)+A(66-I)
      C2R = A(I)-A(66-I)
      C1I = A1-A2
      C2I = A1+A2
      C1 = -D(4*I-5)
      C2 = D(4*I-4)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(66-I) = C1R-C3R
   10 C(67-I) = C3I-C1I
      C(1) = A(1)+A(65)
      C(2) = A(1)-A(65)
      C(33) = 2.*A(33)
      C(34) = 2.*(A(32)-A(34))
      CALL D32(C,C(2),2,C,C(2),2)
      S2 = 2.*S2
      B(1) = C(1)+S2
      B(65) = C(1)-S2
      B(33) = C(33)
      DO 20 I=2,32
      E1 = 0.5*(C(I)+C(66-I))
      E2 = E(4*I-4)*(C(I)-C(66-I))
      B(I) = E1+E2
   20 B(66-I) = E1-E2
      RETURN
      END
