*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RRTEQ3(R,S,T,X,D)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DRTEQ3(R,S,T,X,D)
#endif
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +  I,Z(0:2)
#include "gen/def128.inc"
     + ZQ1,QQ,PP,Q1,Q2,Q3
      DIMENSION X(*),Y(0:2)

      PARAMETER(EPS = 1D-6, DELTA = 1D-15)
      PARAMETER(I = (0,1), ZD1 = 1, ZQ1 = 1)
      PARAMETER(R1 = 2*ZD1/27, R2 = ZD1/2, R3 = ZD1/3)
      PARAMETER(W3 = 1.73205 08075 68877 294D0, R4 = W3/2)
      PARAMETER(Q1 = 2*ZQ1/27, Q2 = ZQ1/2, Q3 = ZQ1/3)

      IF(S .EQ. 0 .AND. T .EQ. 0) THEN
       X(1)=-R
       X(2)=0
       X(3)=0
       D=0
       RETURN
      ENDIF
      P=S-R3*R**2
      Q=(R1*R**2-R3*S)*R+T
      D=(R2*Q)**2+(R3*P)**3
      IF(ABS(D) .LE. EPS) THEN
       PP=S-Q3*R**2
       QQ=(Q1*R**2-Q3*S)*R+T
       D=(Q2*QQ)**2+(Q3*PP)**3
       P=PP
       Q=QQ
      ENDIF
      H=R3*R
      H1=R2*Q
      IF(D .GT. DELTA) THEN
       H2=SQRT(D)
       U0=-H1+H2
       V0=-H1-H2
       U=SIGN(ABS(U0)**R3,U0)
       V=SIGN(ABS(V0)**R3,V0)
       X(1)=U+V-H
       X(2)=-R2*(U+V)-H
       X(3)=R4*ABS(U-V)
       IF(ABS(U0) .LE. EPS .OR. ABS(V0) .LE. EPS) THEN
        Y(0)=X(1)
        DO 1 K = 0,1
    1   Y(K+1)=Y(K)-(((Y(K)+R)*Y(K)+S)*Y(K)+T)/((3*Y(K)+2*R)*Y(K)+S)
        X(1)=Y(2)
        Z(0)=X(2)+I*X(3)
        DO 2 K = 0,1
    2   Z(K+1)=Z(K)-(((Z(K)+R)*Z(K)+S)*Z(K)+T)/((3*Z(K)+2*R)*Z(K)+S)
        X(2)=Z(2)
        X(3)=-I*Z(2)
       ENDIF
      ELSEIF(ABS(D) .LE. DELTA) THEN
       D=0
       U=SIGN(ABS(H1)**R3,-H1)
       X(1)=U+U-H
       X(2)=-U-H
       X(3)=X(2)
       IF(ABS(H1) .LE. EPS) THEN
        Y(0)=X(1)
        DO 5 K = 0,1
        H1=(3*Y(K)+2*R)*Y(K)+S
        IF(ABS(H1) .GT. DELTA) THEN
         Y(K+1)=Y(K)-(((Y(K)+R)*Y(K)+S)*Y(K)+T)/H1
        ELSE
         X(1)=-R3*R
         X(2)=X(1)
         X(3)=X(1)
         RETURN
        ENDIF
    5   CONTINUE
        X(1)=Y(2)
        X(2)=-R2*(R+X(1))
        X(3)=X(2)
       ENDIF
      ELSE
       H3=SQRT(ABS(R3*P)**3)
       H2=R3*ACOS(-H1/H3)
       H1=H3**R3
       U=H1*COS(H2)
       V=W3*H1*SIN(H2)
       X(1)=U+U-H
       X(2)=-U-V-H
       X(3)=-U+V-H
       IF(H3 .LE. EPS .OR. X(1) .LE. EPS .OR. X(2) .LE .EPS .OR.
     1                     X(3) .LE. EPS) THEN
        DO 3 J = 1,3
        Y(0)=X(J)
        DO 4 K = 0,1
    4   Y(K+1)=Y(K)-(((Y(K)+R)*Y(K)+S)*Y(K)+T)/((3*Y(K)+2*R)*Y(K)+S)
    3   X(J)=Y(2)
       ENDIF
      ENDIF
      RETURN
      END
