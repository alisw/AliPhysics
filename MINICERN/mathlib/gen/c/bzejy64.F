*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RBZEJY(A,N,MODE,REL,X)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DBZEJY(A,N,MODE,REL,X)
#include "gen/imp64.inc"
#endif
C     Computes the first n positive (in the case Jo'(x) the first n
C     non-negative) zeros of the Bessel functions
C               Ja(x), Ya(x), Ja'(x), Ya'(x),
C     where a >= 0 and ' = d/dx.
C
C     Based on Algol procedures published in
C
C     N.M. TEMME, An algorithm with Algol 60 program for the compu-
C     tation of the zeros of ordinary Bessel functions and those of
C     their derivatives, J. Comput. Phys. 32 (1979) 270-279, and
C
C     N.M. TEMME, On the numerical evaluation of the ordinary Bessel
C     function of the second kind, J. Comput. Phys. 21 (1976) 343-350.

      DIMENSION X(*),FCO(4),ICO(4),H(13)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RBZEJY')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RBZEJY/DBZEJY')
#endif

      PARAMETER (Z1 = 1, Z2 = 2, Z5 = 5, Z6 = 6, Z7 = 7)
      PARAMETER (C1 = Z5/48, C2 = -Z5/36, C3 = -Z7/48)
      PARAMETER (C4 = 5*Z7/288, C5 = Z1/64, C6 = Z6/5, C7 = 3*Z6/7)
      PARAMETER (E1 = Z1/3, E2 = Z2/3, HALF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PIH = PI/2, PI1 = 3*PI/8)

      DATA FCO /-2.33811,-1.17371,-1.01879,-2.29444/, ICO /1,3,3,1/

#if !defined(CERNLIB_DOUBLE)
      ENTRY BZEJY(A,N,MODE,REL,X)
#endif

      IF(N .LE. 0) RETURN
      IF(A .LT. 0) THEN
       WRITE(ERRTXT,101) A
       CALL MTLPRT(NAME,'C345.1',ERRTXT)
       RETURN
      ENDIF
      NA=NINT(A)
      A0=A-NA
      H(1)=A**2
      H(2)=4*H(1)
      H(3)=H(2)**2
      H(4)=H(3)+H(3)
      H(5)=12*A+6
      H(6)=H(5)-4
      H(7)=3*A-8
      H(8)=HALF*(A-HALF*ICO(MODE))
      H(9)=COS(A0*PI)/(PI*REL)
      H(10)=PIH*(A+HALF)
      H(11)=HALF+A0
      H(12)=H(11)*(HALF-A0)
      H(13)=1
      IF(A .GE. 3) H(13)=A**(-E2)
      P1=0
      Q1=0
      IF(MODE .LT. 3) THEN
       IF(N .GE. H(7)) THEN
        P=7*H(2)-31
        PP=H(2)-1
        IF(1+P .NE. P) THEN
         P1=PP*((1012*H(2)-14888)*H(2)+71476)/(15*P)
         Q1=((664*H(2)-7856)*H(2)+30232)/(5*P)
        ENDIF
       ENDIF
       P2=C1
       Q2=C2
       R2=C6
      ELSE
       IF(N .GE. H(7)) THEN
        P=(7*H(2)+82)*H(2)-9
        PP=H(2)+3
        IF(1+P .NE. P) THEN
         P1=((((1012*H(2)+32816)*H(2)-55496)*H(2)-104400)*H(2)+
     1         253044)/(15*P)
         Q1=(((664*H(2)+16600)*H(2)-24312)*H(2)+28296)/(5*P)
        ENDIF
       ENDIF
       P2=C3
       Q2=C4
       R2=C7
      ENDIF
      X(1)=0
      N0=1
      IF(A .EQ. 0 .AND. MODE .EQ. 3) N0=2

      DO 1 L = N0,N
      IF(L .GE. H(7)) THEN
       B=(L+H(8))*PI
       C=C5/B**2
       X0=B-(PP-P1*C)/(8*B*(1-Q1*C))
      ELSE
       IF(L .EQ. 1) THEN
        X0=FCO(MODE)
       ELSE
        X0=PI1*(4*L-ICO(MODE))
        V=1/X0**2
        X0=-X0**E2*(1+V*(P2+Q2*V))
       ENDIF
       U=H(13)*X0
       Y=E2*SQRT(ABS(U)**3)

       IF(Y .LT. 1) THEN
        G0=(3*Y)**E1
        G1=G0**2
        G0=G0*(1-G1*(210+G1*(G1+G1-27))/1575)
       ELSE
        G0=1/(Y+PIH)
        G1=G0**2
        G0=PIH-G0*(1+G1*(2310+G1*(3003+G1*(4818+G1*
     1          (8591+G1*16328))))/3465)
       ENDIF
       G2=Y+G0
       G1=G2**2
       G=(G0-ATAN(G2))/G1
       W=1/COS(G0-(1+G1)*(1+G/G2)*G)
       G=1/(1-W**2)
       C=SQRT(U*G)
       X0=W*(A-P2*C*(1/U-C*(G+G-R2))/(A*U))
      ENDIF

      DO 2 J = 1,5
      XX=X0**2
      AX=H(1)-XX

      B=X0+X0
      E=(H(9)*X0)**2
      P=1
      Q=-X0
      S=1+XX
      R=S
      DO 3 K = 2,500
      IF(R*(K**2) .GE. E) GO TO 9
      RK=Z1/(K+1)
      D=((K-1)+H(12)/K)/S
      P=RK*((K+K)-P*D)
      Q=RK*(Q*D-B)
      S=P**2+Q**2
      R=R*S
    3 CONTINUE
    9 G=1/S
      P=G*P
      F=P
      Q=-G*Q
      G=Q
      DO 4 M = K+1,2,-1
      R=M*(2-P)-2
      S=B+M*Q
      D=((M-2)+H(12)/(M-1))/(R**2+S**2)
      P=D*R
      Q=D*S
      E=F+1
      F=P*E-G*Q
      G=Q*E+P*G
    4 CONTINUE
      F=F+1
      D=1/(F**2+G**2)
      PA=F*D
      QA=-G*D
      D=H(11)-P
      Q=Q+X0
      RX=1/X0
      PA1=RX*(PA*Q-QA*D)
      QA1=RX*(QA*Q+PA*D)
      XR=RX+RX
      B=(1+A0)*XR
      DO 5 M = 1,NA
      P0=PA-QA1*B
      Q0=QA+PA1*B
      PA=PA1
      PA1=P0
      QA=QA1
      QA1=Q0
      B=B+XR
    5 CONTINUE

      SX=SIN(X0-H(10))
      CX=COS(X0-H(10))
      GO TO (11,12,13,14), MODE
   11 R0=(PA*CX-QA*SX)/(PA1*SX+QA1*CX)
      V=R0/(H(5)*X0)
      W=R0*(1+V*(1-4*AX))/(1+V*(4*(XX-H(2))-H(6)))
      GO TO 10
   12 R0=(PA*SX+QA*CX)/(QA1*SX-PA1*CX)
      V=R0/(H(5)*X0)
      W=R0*(1+V*(1-4*AX))/(1+V*(4*(XX-H(2))-H(6)))
      GO TO 10
   13 R0=A/X0-(PA1*SX+QA1*CX)/(PA*CX-QA*SX)
      GO TO 15
   14 R0=A/X0-(QA1*SX-PA1*CX)/(PA*SX+QA*CX)
   15 V=E2*R0*X0*AX/(H(1)+XX)
      W=C5/AX**3
      W=-XX*R0*(1+V*(1-W*(H(3)-XX*(40*H(2)+48*XX))))/
     1     (AX*(1+V*(1+W*(H(4)+XX*(64*H(2)+96*XX)))))
   10 X0=X0+W
      IF(ABS(W) .LE. REL*X0) GO TO 6
    2 CONTINUE
    6 X(L)=X0
    1 CONTINUE
      RETURN
  101 FORMAT('NEGATIVE INDEX A = ',1P,D15.8)
      END
