*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CDIGAM(Z)
#include "gen/defc64.inc"
     +  CDIGAM
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION WDIGAM(Z)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +   WDIGAM
#endif
#include "gen/defc64.inc"
     + Z,U,V,H,R,P
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CDIGAM')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CDIGAM/WDIGAM')
#endif
      DIMENSION C(6)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)

#include "gen/gcmpfun.inc"
CSEQ,GCMPLX.

      DATA C(1) / 8.33333 33333 33333 33D-2/
      DATA C(2) /-8.33333 33333 33333 33D-3/
      DATA C(3) / 3.96825 39682 53968 25D-3/
      DATA C(4) /-4.16666 66666 66666 67D-3/
      DATA C(5) / 7.57575 75757 57575 76D-3/
      DATA C(6) /-2.10927 96092 79609 28D-2/

      U=Z
      X=U
      A=ABS(X)
      IF(GIMAG(U) .EQ. 0 .AND. -A .EQ. INT(X)) THEN
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C307.1',ERRTXT)
      ELSE
       IF(X .LT. 0) U=-U
       V=U
       H=0
       IF(A .LT. 15) THEN
        H=1/V
        DO 1 I = 1,14-INT(A)
        V=V+1
    1   H=H+1/V
        V=V+1
       END IF
       R=1/V**2
       P=R*C(1)
       DO 2 I = 6,1,-1
    2  P=R*(C(I)+P)
       H=LOG(V)-HF/V-P-H
       IF(X .LT. 0) THEN
        V=PI*U
        X=V
        A=SIN(X)
        X=COS(X)
        Y=TANH(GIMAG(V))
        H=H+1/U+PI*GCMPLX(X,-A*Y)/GCMPLX(A,X*Y)
       END IF
      ENDIF
#if defined(CERNLIB_DOUBLE)
      WDIGAM=H
#endif
#if !defined(CERNLIB_DOUBLE)
      CDIGAM=H
#endif
      RETURN
  101 FORMAT(1X,'ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
