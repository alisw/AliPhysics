*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE S128(A,B)
      REAL A(64),B(64),C(64)
      COMMON /TAB/D(126),E(127)
      CALL D701BD
      DO 10 I=3,31,2
      A1 = A(I-1)-A(I+1)
      A2 = A(65-I)-A(67-I)
      C1R = A1+A2
      C2R = A1-A2
      C1I = A(I)-A(66-I)
      C2I = A(I)+A(66-I)
      C1 = D(4*I-5)
      C2 = D(4*I-4)
      C3R = C1*C2R-C2*C2I
      C3I = C1*C2I+C2*C2R
      C(I) = C1R+C3R
      C(I+1) = C1I+C3I
      C(66-I) = C1R-C3R
   10 C(67-I) = C3I-C1I
      C(1) = 2.*(A(64)-A(2))
      C(2) = -2.*(A(64)+A(2))
      C(33) = 2.*(A(32)-A(34))
      C(34) = 2.*A(33)
      CALL I32(C,C(2),2,C,C(2),2)
      B(33) = 0.5*C(33)
      DO 20 I=2,32
      E1 = 0.5*(C(I)-C(66-I))
      E2 = E(4*I-4)*(C(I)+C(66-I))
      B(I) = E2-E1
   20 B(66-I) = E2+E1
      RETURN
      END
