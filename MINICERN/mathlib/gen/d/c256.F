*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE C256(A,B)
      REAL A(129),B(129),C(128)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      S2 = A(64)+A(66)
      DO 10 I=3,63,2
      S2 = S2+A(I-1)+A(131-I)
      A1 = A(I+1)-A(I-1)
      A2 = A(131-I)-A(129-I)
      C1R = A(I)+A(130-I)
      C2R = A(I)-A(130-I)
      C1I = A1-A2
      C2I = A1+A2
      C1 = -D(2*I-3)
      C2 = D(2*I-2)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(130-I) = C1R-C3R
   10 C(131-I) = C3I-C1I
      C(1) = A(1)+A(129)
      C(2) = A(1)-A(129)
      C(65) = 2.*A(65)
      C(66) = 2.*(A(64)-A(66))
      CALL D64(C,C(2),2,C,C(2),2)
      S2 = 2.*S2
      B(1) = C(1)+S2
      B(129) = C(1)-S2
      B(65) = C(65)
      DO 20 I=2,64
      E1 = 0.5*(C(I)+C(130-I))
      E2 = E(2*I-2)*(C(I)-C(130-I))
      B(I) = E1+E2
   20 B(130-I) = E1-E2
      RETURN
      END
