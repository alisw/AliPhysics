*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION RANGAM(P)
C     THIS FUNCTION GENERATES A GAMMA DISTRIBUTED RANDOM NUMBER
C     WITH PARAMETER P .GT.0.
C     THE JOHNK S ALGORITHM IS USED
      DIMENSION STOR(20)
      RANGAM=0.
      IF(P.GT.15.) GO TO 40
      M=INT(P)
      F=P-M
      IF(M.EQ.0) GO TO 20
      X=1.
      CALL NRAN(STOR,M)
      DO 1 I=1,M
    1 X = X * STOR(I)
      RANGAM=-LOG(X)
   20 IF ( F .LT. 1.0E-5)  RETURN
C
      X1=-LOG(RNDM(3))
      IF (F .LT. 0.9999)  GO TO 25
   22 RANGAM = RANGAM + X1
      RETURN
C         NON-INTEGER VALUE OF P
   25 CONTINUE
C     ....W1=R1**(1/F)
      WLOG = LOG (RNDM (1) ) / F
      IF (WLOG .LT. -100.)  RETURN
      W1 = EXP(WLOG)
C     ....W2=R2**(1/(1-F))
      WLOG = LOG(RNDM(2)) / (1.-F)
      IF (WLOG .LT. -100.)  GO TO 22
      W2 = EXP(WLOG)
      W=W1+W2
      IF(W.GT.1.) GO TO 25
      X2=W1/W
      RANGAM=RANGAM+X1*X2
      RETURN
C
C         WILSON - HILFERTY APPROXIMATION
#if (defined(CERNLIB_IBM)||defined(CERNLIB_CDC))&&(!defined(CERNLIB_FORTRAN))
   40 CALL NORRAN(A)
#endif
#if (!defined(CERNLIB_IBM)||defined(CERNLIB_FORTRAN))&&(!defined(CERNLIB_CDC)||defined(CERNLIB_FORTRAN))
   40 CALL RANNOR(A,B)
#endif
      Q=1.-1./(9.*P)+A/(3.*SQRT(P))
      RANGAM=P*Q*Q*Q
      IF (RANGAM.LE.0.) GO TO 40
      RETURN
      END
