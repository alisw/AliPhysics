*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:21  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MCONV (N)
C-----MCONV INVERTS THE POSITIVE DEFINITE PACKED SYMMETRIC MATRIX Z
C-----BY THE SQUARE-ROOT METHOD
#include "d510pl.inc"
#include "d510si.inc"
#include "d510uo.inc"
C-----MAXIMUM REAL NUMBER AND 10.*MAXIMUM RELATIVE PRECISION ON CDC6000
#if defined(CERNLIB_CDC)
      DATA AM,RP/1.E300,1.E-14/
#endif
#if defined(CERNLIB_IBM)
      DATA AM,RP/1.E75,1.E-14/
#endif
#if (!defined(CERNLIB_CDC))&&(!defined(CERNLIB_IBM))
      DATA AM,RP / 1.0E37, 1.0E-14/
#endif
      IF (N.LT.1) RETURN
      APS=SQRT(AM/N)
      AP=1./(APS*APS)
      IR=0
      DO 11 I=1,N
 1    IR=IR+1
      IF (PL(IR)) 1,1,2
 2    NI=I*(I-1)/2
      II=NI+I
      K=N+1
      IF (Z(II).LE.RP*ABS(R(IR)).OR.Z(II).LE.AP) GO TO 19
      Z(II)=1./SQRT(Z(II))
      NL=II-1
 3    IF (NL-NI) 5,5,4
 4    Z(NL)=Z(NL)*Z(II)
      IF (ABS(Z(NL)).GE.APS) GO TO 16
      NL=NL-1
      GO TO 3
 5    IF (I-N) 6,12,12
 6    K=K-1
      NK=K*(K-1)/2
      NL=NK
      KK=NK+I
      D=Z(KK)*Z(II)
      C=D*Z(II)
      L=K
 7    LL=NK+L
      LI=NL+I
      Z(LL)=Z(LL)-Z(LI)*C
      L=L-1
      NL=NL-L
      IF (L-I) 9,9,7
 8    LL=NK+L
      LI=NI+L
      Z(LL)=Z(LL)-Z(LI)*D
 9    L=L-1
      IF (L) 10,10,8
 10   Z(KK)=-C
      IF (K-I-1) 11,11,6
 11   CONTINUE
 12   DO 14 I=1,N
      DO 14 K=I,N
      NL=K*(K-1)/2
      KI=NL+I
      D=0.
      DO 13 L=K,N
      LI=NL+I
      LK=NL+K
      D=D+Z(LI)*Z(LK)
      NL=NL+L
 13   CONTINUE
      KI=K*(K-1)/2+I
      Z(KI)=D
 14   CONTINUE
 15   RETURN
 16   K=I+NL-II
      IR=0
      DO 18 I=1,K
 17   IR=IR+1
      IF (PL(IR)) 17,17,18
 18   CONTINUE
 19   PL(IR)=-2.
      R(IR)=0.
      INDFLG(1)=IR
      GO TO 15
      END
