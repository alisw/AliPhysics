*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RRTEQ4(A,B,C,D,Z,DC,MT)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DRTEQ4(A,B,C,D,Z,DC,MT)
#endif
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +  I,Z,Z0,W1,W2,W3
#include "gen/def128.inc"
     +      ZQ1,AA,PP,QQ,RR,Q1,Q2,Q3,Q4,Q8
      DIMENSION Z(*),Z0(4),U(3),V(3)

      PARAMETER(I = (0,1), ZD1 = 1, ZQ1 = 1)
      PARAMETER(R4 = ZD1/4, R12 = ZD1/12)
      PARAMETER(Q2 = ZQ1/2, Q4 = ZQ1/4, Q8 = ZQ1/8)
      PARAMETER(Q1 = 3*ZQ1/8, Q3 = 3*ZQ1/16)

      IF(B .EQ. 0 .AND. C .EQ. 0) THEN
       IF(D .EQ. 0) THEN
        MT=1
        Z(1)=-A
        Z(2)=0
        Z(3)=0
        Z(4)=0
        DC=0
        RETURN
       ELSEIF(A .EQ. 0) THEN
        IF(D .GT. 0) THEN
         MT=2
         Z(1)=SQRT(I*SQRT(D))
         Z(2)=-Z(1)
         Z(4)=SQRT(-Z(1)**2)
         Z(3)=-Z(4)
        ELSE
         MT=3
         Z(1)=SQRT(SQRT(-D))
         Z(2)=-Z(1)
         Z(3)=SQRT(-Z(1)**2)
         Z(4)=-Z(3)
        ENDIF
        DC=(-R12*D)**3
        RETURN
       ENDIF
      ENDIF
      AA=A**2
      PP=B-Q1*AA
      QQ=C-Q2*A*(B-Q4*AA)
      RR=D-Q4*(A*C-Q4*AA*(B-Q3*AA))
      RC=Q2*PP
      SC=Q4*(Q4*PP**2-RR)
      TC=-(Q8*QQ)**2
#if defined(CERNLIB_DOUBLE)
      CALL DRTEQ3(RC,SC,TC,U,DC)
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL RRTEQ3(RC,SC,TC,U,DC)
#endif
      Q=QQ
      H=R4*A
      IF(DC .EQ. 0) U(3)=U(2)
      IF(DC .LE. 0) THEN
       MT=2
       V(1)=ABS(U(1))
       V(2)=ABS(U(2))
       V(3)=ABS(U(3))
       V1=MAX(V(1),V(2),V(3))
       IF(V1 .EQ. V(1)) THEN
        K1=1
        V2=MAX(V(2),V(3))
       ELSEIF(V1 .EQ. V(2)) THEN
        K1=2
        V2=MAX(V(1),V(3))
       ELSE
        K1=3
        V2=MAX(V(1),V(2))
       ENDIF
       IF(V2 .EQ. V(1)) THEN
        K2=1
       ELSEIF(V2 .EQ. V(2)) THEN
        K2=2
       ELSE
        K2=3
       ENDIF
       W1=SQRT(U(K1)+I*0)
       W2=SQRT(U(K2)+I*0)
      ELSE
       MT=3
       W1=SQRT(U(2)+I*U(3))
       W2=SQRT(U(2)-I*U(3))
      ENDIF
      W3=0
      IF(W1*W2 .NE. 0) W3=-Q/(8*W1*W2)
      Z0(1)=W1+W2+W3-H
      Z0(2)=-W1-W2+W3-H
      Z0(3)=-W1+W2-W3-H
      Z0(4)=W1-W2-W3-H
      IF(MT .EQ. 2) THEN
       IF(U(K1) .GE. 0 .AND. U(K2) .GE. 0) THEN
        MT=1
        DO 1 J = 1,4
        RZ0=Z0(J)
    1   Z(J)=RZ0
       ELSEIF(U(K1) .GE. 0 .AND. U(K2) .LT. 0) THEN
        Z(1)=Z0(1)
        Z(2)=Z0(4)
        Z(3)=Z0(3)
        Z(4)=Z0(2)
       ELSEIF(U(K1) .LT. 0 .AND. U(K2) .GE. 0) THEN
        Z(1)=Z0(1)
        Z(2)=Z0(3)
        Z(3)=Z0(4)
        Z(4)=Z0(2)
       ELSEIF(U(K1) .LT. 0 .AND. U(K2) .LT. 0) THEN
        Z(1)=Z0(1)
        Z(2)=Z0(2)
        Z(3)=Z0(4)
        Z(4)=Z0(3)
       ENDIF
      ELSEIF(MT .EQ. 3) THEN
       DO 2 J = 1,2
       RZ0=Z0(J)
    2  Z(J)=RZ0
       Z(3)=Z0(4)
       Z(4)=Z0(3)
      ENDIF
      RETURN
      END
