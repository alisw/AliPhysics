*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE PRIPAR (A1,A2,B)
      COMMON /SLATE/XQX(9), V(4), G1(9), G2(9), G(9)
      DIMENSION A1(3),A2(3),B(3)
      CALL UCOPY (A1,G,3)
      V(3)=0.
      I= 0
    1 CONTINUE
      GINV= 1./G(1)
      G2(7)=SIN(G(2))
      G2(8)=COS(G(2))
      G2(2)=COS (G(3))
      G2(5)= SIN(G(3))
      G2(3)=-G2(8)*G2(5)
      G2(6)=G2(8)*G2(2)
      G2(2)=-G2(7)*G2(2)
      G2(5)=-G2(7)*G2(5)
      G2(1)=-G2(6)*GINV
      G2(4)= G2(3)*GINV
      G2(7)=-G2(7)*GINV
      G2(9)=0.
      V(3)=V(3)-G2(7)
      CALL VSCALE (G2,GINV,G2,9)
      IF (I.EQ.2) GO TO 5
      CALL UCOPY (G2,G1,9)
      CALL UCOPY (A2,G,3)
      I= 2
      GO TO 1
    5 V(1)=G1(6)+G2(6)
      V(2)=-G1(3)-G2(3)
      B(1)= V(1)**2+V(2)**2
      B(2)=ATAN(V(3)/SQRT(B(1)))
      V(4)=SQRT(V(3)**2+ B(1))
      B(3)= ATAN2(V(2),V(1))
      IF(B(3).LT.0.)  B(3)=B(3)+6.28318
      B(1)=1./V(4)
      RETURN
      END
