*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE S256(A,B)
      REAL A(128),B(128),C(128)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      DO 10 I=3,63,2
      A1 = A(I-1)-A(I+1)
      A2 = A(129-I)-A(131-I)
      C1R = A1+A2
      C2R = A1-A2
      C1I = A(I)-A(130-I)
      C2I = A(I)+A(130-I)
      C1 = D(2*I-3)
      C2 = D(2*I-2)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(130-I) = C1R-C3R
   10 C(131-I) = C3I-C1I
      C(1) = 2.*(A(128)-A(2))
      C(2) = -2.*(A(128)+A(2))
      C(65) = 2.*(A(64)-A(66))
      C(66) = 2.*A(65)
      CALL I64(C,C(2),2,C,C(2),2)
      B(65) = 0.5*C(65)
      DO 20 I=2,64
      E1 = 0.5*(C(I)-C(130-I))
      E2 = E(2*I-2)*(C(I)+C(130-I))
      B(I) = E2-E1
   20 B(130-I) = E2+E1
      RETURN
      END
