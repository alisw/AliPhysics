*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:38  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
      INTEGER I,J,K,L,M,N,NM
      REAL A(NM,N),E(N),Z(NM,M)
      REAL H,S
      IF (N .EQ. 1) GO TO 200
      DO 140 I = 2, N
         L = I - 1
         H = E(I) * A(I,L)
         IF (H .EQ. 0.0) GO TO 140
         DO 130 J = 1, M
            S = 0.0
            DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
            S = S / H
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
