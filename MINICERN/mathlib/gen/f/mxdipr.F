*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:32  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MXDIPR(A,B,C,IAD,JBD,IJD,IA,KA,JB,LB)
C
C     THIS ROUTINE COMPUTES THE DIRECT, OR TENSOR, OR KRONECKER, PRODUCT
C     OF TWO MATRICES A AND B, AND STORES THE RESULT IN C.
C
      DIMENSION A(IAD,1),B(JBD,1),C(IJD,1)
C
      IAJB=IA*JB
      KALB=KA*LB
      IF(IAJB*KALB .EQ. 0) RETURN
      JB1=JB-1
      LB1=LB-1
      DO 1 I = 1,IAJB
      IIA=(I+JB1)/JB
      IIB=MOD(I-1,JB)+1
      DO 1 K = 1,KALB
      KKA=(K+LB1)/LB
      KKB=MOD(K-1,LB)+1
    1 C(I,K)=A(IIA,KKA)*B(IIB,KKB)
      RETURN
      END
