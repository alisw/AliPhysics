*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CWERF(Z)
#include "gen/defc64.inc"
     +     CWERF
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION WWERF(Z)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +     WWERF
#endif
#include "gen/defc64.inc"
C    +     Z,ZH,R(37),S,T,V,W,GCMPLX,GCONJG
     +     Z,ZH,R(37),S,T,V,W

      PARAMETER (Z1 = 1, HF = Z1/2, Z10 = 10)
      PARAMETER (C1 = 74/Z10, C2 = 83/Z10, C3 = Z10/32, C4 = 16/Z10)
      PARAMETER (C = 1.12837 91670 95512 57D0, P = (2*C4)**33)

#include "gen/gcmpfun.inc"

C     GREAL(W)=DREAL(W)
C     GIMAG(W)=DIMAG(W)
C     GCONJG(W)=DCONJG(W)
C     GCMPLX(X,Y)=DCMPLX(X,Y)

      X=Z
      Y=GIMAG(Z)
      XA=ABS(X)
      YA=ABS(Y)
      IF(YA .LT. C1 .AND. XA .LT. C2) THEN
       ZH=GCMPLX(YA+C4,XA)
       R(37)=0
       DO 1 N = 36,1,-1
       T=ZH+N*GCONJG(R(N+1))
    1  R(N)=HF*T/(GREAL(T)**2+GIMAG(T)**2)
       XL=P
       S=0
       DO 2 N = 33,1,-1
       XL=C3*XL
    2  S=R(N)*(S+XL)
       V=C*S
      ELSE
       ZH=GCMPLX(YA,XA)
       R(1)=0
       DO 3 N = 9,1,-1
       T=ZH+N*GCONJG(R(1))
    3  R(1)=HF*T/(GREAL(T)**2+GIMAG(T)**2)
       V=C*R(1)
      END IF
      IF(YA .EQ. 0) V=GCMPLX(EXP(-XA**2),GIMAG(V))
      IF(Y .LT. 0) THEN
       V=2*EXP(-GCMPLX(XA,YA)**2)-V
       IF(X .GT. 0) V=GCONJG(V)
      ELSE
       IF(X .LT. 0) V=GCONJG(V)
      END IF
#if !defined(CERNLIB_DOUBLE)
      CWERF=V
#endif
#if defined(CERNLIB_DOUBLE)
      WWERF=V
#endif
      RETURN
      END
