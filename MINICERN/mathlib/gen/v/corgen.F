*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE CORGEN(C,X,NP)
C
C         CORSET sets up the generation by calculating C from V.
C         CORGEN generates a set of NP random numbers
C                Gaussian-distributed with covariance matrix V
C                (V = C*C') and mean values zero.
C
      DIMENSION C(NP,NP), X(NP)
      PARAMETER (NMAX=100)
      DIMENSION Z(NMAX)
C
      IF (NP .GT. NMAX)  GO TO 120
C
      CALL RNORML(Z,NP)
C
      DO 100 I= 1, NP
         X(I) = 0.
         DO 90 J= 1, I
         X(I) = X(I) + C(I,J)*Z(J)
   90    CONTINUE
  100 CONTINUE
C
      RETURN
C                Error return
  120 CONTINUE
      WRITE (6,121) NP,NMAX
  121 FORMAT (' ERROR IN CORGEN. VECTOR LENGTH NP=',I5,
     1   ', BUT MAXIMUM ALLOWED IS',I5)
      RETURN
      END
