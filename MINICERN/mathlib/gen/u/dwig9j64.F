*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DWIG9J(A,B,C,P,Q,R,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      PARAMETER (R1 = 1, HF = R1/2)
 
      IA=NINT(2*A)
      IB=NINT(2*B)
      IC=NINT(2*C)
      IP=NINT(2*P)
      IQ=NINT(2*Q)
      IR=NINT(2*R)
      IX=NINT(2*X)
      IY=NINT(2*Y)
      IZ=NINT(2*Z)
 
      H=0
      IF(IA .LT. 0 .OR. IB .LT. 0 .OR. IC .LT. 0 .OR.
     1   IP .LT. 0 .OR. IQ .LT. 0 .OR. IR .LT. 0 .OR.
     2   IX .LT. 0 .OR. IY .LT. 0 .OR. IZ .LT. 0) GO TO 99
      J0=MAX(ABS(IA-IZ),ABS(IY-IP),ABS(IB-IR))
      J1=MIN(IA+IZ,IY+IP,IB+IR)
      S=0
      V=(-1)**J0
      DO 1 J = J0,J1
      AJ=HF*J
      H=DWIG6J(A,P,X,Y,Z,AJ)
      IF(H .NE. 0) H=H*DWIG6J(B,Q,Y,P,AJ,R)
      IF(H .NE. 0) H=H*DWIG6J(C,R,Z,AJ,A,B)
      S=S+V*(AJ+HF)*H
    1 V=-V
      H=2*S
   99 DWIG9J=H
      RETURN
      END
#endif
