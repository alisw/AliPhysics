*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION RTRGSM(X,A,N,B,M,IOP)
 
      DIMENSION A(0:*),B(*)
 
      PARAMETER (R1 = 1, HF = R1/2)
 
      C=2*COS(X)
      H=0
      IF(IOP .NE. 3) THEN
       A2=0
       A1=0
       A0=A(N)
       DO 1 K = N-1,0,-1
       A2=A1
       A1=A0
       A0=A(K)+C*A1-A2
    1  CONTINUE
       H=HF*(A(0)+A0-A2)
      ENDIF
      IF(IOP .NE. 2) THEN
       A2=0
       A1=0
       A0=B(M)
       DO 2 K = M-1,1,-1
       A2=A1
       A1=A0
       A0=B(K)+C*A1-A2
    2  CONTINUE
       H=H+A0*SIN(X)
      ENDIF
      RTRGSM=H
      RETURN
      END
