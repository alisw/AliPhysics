*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRINV (T,S,N)
C
C CERN PROGLIB# F112    TRINV           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM, R
#endif
      DIMENSION T(*),S(*)
C
      MX = (N*N+N)/2
      IPIV = MX
C
      I = N
   10 R = 0.
      IF (T(IPIV).GT.0.)  R=1./T(IPIV)
      S(IPIV) = R
      NDSTEP = N
      IND = MX - N + I
C
   20 IF (IND.EQ.IPIV)        GO TO 50
      SUM = 0.
      IF (R.EQ.0.)            GO TO 40
      LHOR = IPIV
      LVER = IND
      J = I
C
   30 LHOR = LHOR + J
      LVER = LVER + 1
      SUM = SUM + T(LHOR)*S(LVER)
      J = J + 1
      IF (LHOR.LT.IND)        GO TO 30
C
   40 S(IND) = - SUM*R
      NDSTEP = NDSTEP - 1
      IND = IND - NDSTEP
        GO  TO                20
C
   50 IPIV = IPIV - I
      I = I - 1
      IF (I.GT.0)             GO TO 10
C
      RETURN
      END
