*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:36  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
      INTEGER I,J,K,L,M,N,NM
      REAL AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      REAL H,S,SI
      DO 50 K = 1, N
         DO 50 J = 1, M
            ZI(K,J) = - ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
      IF (N .EQ. 1) GO TO 200
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0) GO TO 140
         DO 130 J = 1, M
            S = 0.0
            SI = 0.0
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
            S = S / H
            SI = SI / H
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
