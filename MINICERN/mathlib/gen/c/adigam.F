*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#if defined(CERLIB_OBSOLETE)
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION ADIGAM(RX)
#include "gen/def64.inc"
     +      DDIGAM,DX,X,C(6),ZERO,ONE,HALF,C0,PI,A,V,H,R,P
      CHARACTER NAME*6

      DATA ZERO /0D0/, ONE /1D0/, HALF /0.5D0/, C0 /15D0/

      DATA PI /3.14159 26535 89793 24D0/
      DATA C(1) / 8.33333 33333 33333 33D-2/
      DATA C(2) /-8.33333 33333 33333 33D-3/
      DATA C(3) / 3.96825 39682 53968 25D-3/
      DATA C(4) /-4.16666 66666 66666 67D-3/
      DATA C(5) / 7.57575 75757 57575 76D-3/
      DATA C(6) /-2.10927 96092 79609 28D-2/

      X=RX
      NAME='ADIGAM'
      ADIGAM = ZERO
    9 A=ABS(X)
      IF(-A .EQ. INT(X)) THEN
       WRITE(6,100) NAME,X
       RETURN
      END IF
      V=A
      H=ZERO
      IF(A .LT. C0) THEN
       H=ONE/V
       DO 1 I = 1,14-INT(A)
       V=V+ONE
    1  H=H+ONE/V
       V=V+ONE
      END IF
      R=ONE/V**2
      P=R*C(1)
      DO 2 I = 6,1,-1
    2 P=R*(C(I)+P)
      H=LOG(V)-HALF/V-P-H
      IF(X .LT. ZERO) H=H+ONE/A+PI/TAN(PI*A)
      IF(NAME .EQ. 'DDIGAM') THEN
      GO TO 10
      ELSE
       ADIGAM=SNGL(H+(H-DBLE(SNGL(H))))
      END IF
      RETURN

      ENTRY DDIGAM(DX)
      X=DX
      NAME='DDIGAM'
      DDIGAM = ZERO
      GO TO 9
  10  DDIGAM = H


  100 FORMAT(7X,'***** CERN C317 ',A6,' ... ARGUMENT EQUALS ',
     1          'NON-POSITIVE INTEGER =',E10.1)
      END
#endif
#if defined(CERNLIB_SINGLE)
      FUNCTION ADIGAM(X)
      DIMENSION B(6)
      DATA PI /3.14159 26535 89793/
      DATA B /+8.33333 33333 33333E-2, -8.33333 33333 33333E-3,
     1        +3.96825 39682 53968E-3, -4.16666 66666 66667E-3,
     2        +7.57575 75757 57576E-3, -2.10927 96092 79609E-2/
      A=ABS(X)
      IF(-A .EQ. AINT(X)) GO TO 4
      V=A
      H=0.
      IF(A .GE. 15.0) GO TO 3
      N=14-INT(A)
      H=1.0/V
      IF(N .EQ. 0) GO TO 2
      DO 1 I = 1,N
      V=V+1.0
    1 H=H+1.0/V
    2 V=V+1.0
    3 R=1.0/V**2
      ADIGAM=LOG(V)-0.5/V-R*(B(1)+R*(B(2)+R*(B(3)+R*(B(4)+R*(B(5)+R*
     1 (B(6)+R*B(1)))))))-H
      IF(X .GE. 0.0) RETURN
      H=PI*A
      ADIGAM=ADIGAM+1.0/A+PI*COS(H)/SIN(H)
      RETURN
    4 WRITE(6,100)X
      ADIGAM=0.
      RETURN
  100 FORMAT(' ADIGAM ... ARGUMENT IS NON-POSITIVE INTEGER =',F20.2)
      END
#endif
#endif
