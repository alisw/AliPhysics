*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:11  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DELI3(X,AKP,P)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELI3/DELI3')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RELI3(X,AKP,P)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELI3')
#endif
C
C     Translation of Algol procedure el3(x,kc,p) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions III.,  Numer. Math. 13 (1969) 305-315
C
      LOGICAL LBO,LBK

      PARAMETER (ID = 16, IB = 4)
      PARAMETER (PI  = 3.14159 26535 89793 24D0)
      PARAMETER (AL2 = 0.69314 71805 59945 31D0)
      PARAMETER (ALB = IB*AL2, RLB = 1/ALB)
      PARAMETER (Z1 = 1, Z10 = 10, HF = Z1/2, C1 = Z1/10)
      PARAMETER (ND = ID-2, CA = Z10**(-ID/2), CB = Z10**(-(ID+2)))
      PARAMETER (ZD = HF/(ND+1))

      CHARACTER*80 ERRTXT

      DIMENSION RA(2:ND),RB(2:ND),RR(2:ND)

      DATA (RB(K),RA(K),K=2,ND)
     2/2.50000000000000000D-01, 7.50000000000000000D-01,
     3 1.66666666666666667D-01, 8.33333333333333333D-01,
     4 1.25000000000000000D-01, 8.75000000000000000D-01,
     5 1.00000000000000000D-01, 9.00000000000000000D-01,
     6 8.33333333333333333D-02, 9.16666666666666667D-01,
     7 7.14285714285714286D-02, 9.28571428571428571D-01,
     8 6.25000000000000000D-02, 9.37500000000000000D-01,
     9 5.55555555555555556D-02, 9.44444444444444445D-01,
     A 5.00000000000000000D-02, 9.50000000000000000D-01,
     B 4.54545454545454545D-02, 9.54545454545454545D-01,
     C 4.16666666666666667D-02, 9.58333333333333333D-01,
     D 3.84615384615384615D-02, 9.61538461538461538D-01,
     E 3.57142857142857143D-02, 9.64285714285714286D-01/

      IF(X .EQ. 0) THEN
       H=0
       GO TO 9
      ENDIF
      HH=X**2
      F=P*HH
      S=AKP
      IF(S .EQ. 0) S=CA/(1+ABS(X))
      T=S**2
      PM=HF*T
      E=HH*T
      Z=ABS(F)
      R=ABS(P)
      H1=HH+1
      IF(E .LT. C1 .AND. Z .LT. C1 .AND. T .LT. 1 .AND. R .LT. 1) THEN
       S=P+PM
       DO 1 K = 2,ND
       RR(K)=S
       PM=PM*T*RA(K)
    1  S=S*P+PM
       U=S*ZD
       S=U
       LBO=.FALSE.
       DO 2 K = ND,2,-1
       U=U+(RR(K)-U)*RB(K)
       LBO=.NOT.LBO
       V=U
       IF(LBO) V=-U
    2  S=S*HH+V
       IF(LBO) S=-S
       U=HF*(U+1)
#if defined(CERNLIB_DOUBLE)
       H=(U-S*H1)*SQRT(H1)*X+U*DASINH(X)
#endif
#if !defined(CERNLIB_DOUBLE)
       H=(U-S*H1)*SQRT(H1)*X+U*ASINH(X)
#endif
      ELSE
       W=1+F
       IF(W .EQ. 0) THEN
        H=0
        WRITE(ERRTXT,103) X,P
        CALL MTLPRT(NAME,'C346.3',ERRTXT)
        GO TO 9
       ENDIF
       PP=P
       IF(PP .EQ. 0) PP=CB/HH
       S=ABS(S)
       Y=ABS(X)
       G=PP-1
       IF(G .EQ. 0) G=CB
       F=PP-T
       IF(F .EQ. 0) F=CB*T
       AM=1-T
       AP=1+E
       R=PP*H1
       FA=G/(F*PP)
       LBO=FA .GT. 0
       FA=ABS(FA)
       PZ=ABS(G*F)
       DE=SQRT(PZ)
       Q=SQRT(ABS(PP))
       PM=PP-MIN(PM,HF)
       IF(PM .GE. 0) THEN
        U=SQRT(R*AP)
        V=Y*DE
        IF(G .LT. 0) V=-V
        D=1/Q
        C=1
       ELSE
        U=SQRT(H1*AP*PZ)
        YE=Y*Q
        V=AM*YE
        Q=-DE/G
        D=-AM/DE
        C=0
        PZ=AP-R
       ENDIF
       IF(LBO) THEN
        R=V/U
        Z=1
        K=1
        IF(PM .LT. 0) THEN
         H1=Y*SQRT(H1/(AP*FA))
         H1=1/H1-H1
         Z=H1-2*R
         R=2+R*H1
         IF(R .EQ. 0) R=CB
         IF(Z .EQ. 0) Z=H1*CB
         R=R/Z
         Z=R
         W=PZ
        ENDIF
        U=U/W
        V=V/W
       ELSE
        T=U+ABS(V)
        LBK=.TRUE.
        IF(PP .LT. 0) THEN
         DE=V/PZ
         YE=2*U*YE
         U=T/PZ
         V=-(F+G*E)/T
         T=PZ*ABS(W)
         Z=(HH*R*F-G*AP+YE)/T
         YE=YE/T
        ELSE
         DE=V/W
         YE=0
         U=(E+PP)/T
         V=T/W
         Z=1
        ENDIF
        IF(S .GT. 1) THEN
         H1=U
         U=V
         V=H1
        ENDIF
       ENDIF
       Y=1/Y
       E=S
       T=1
       N=1
       L=0
       M=0
    3  Y=Y-E/Y
       IF(Y .EQ. 0) Y =SQRT(E)*CB
       F=C
       C=D/Q+C
       G=E/Q
       D=2*(F*G+D)
       Q=G+Q
       G=T
       T=S+T
       N=2*N
       M=2*M
       IF(LBO) THEN
        IF(Z .LT. 0) M=K+M
        K=0
        IF(R .GT. 0) K=1
        IF(R .LT. 0) K=-1
        H1=E/(U**2+V**2)
        U=U*(1+H1)
        V=V*(1-H1)
       ELSE
        R=U/V
        H1=Z*R
        Z=H1*Z
        HH=E/V
        IF(LBK) THEN
         DE=DE/U
         YE=YE*(H1+1/H1)+DE*(1+R)
         DE=DE*(U-HH)
         LBK=ABS(YE) .LT. 1
        ELSE
         A=LOG(Z)
         K=RLB*A
         IF(Z .GE. 1) K=K+1
         Z=EXP(A-K*ALB)
         M=M+K
        ENDIF
       ENDIF
       IF(ABS(G-S) .GT. CA*G) THEN
        IF(LBO) THEN
         G=HF*(1/R-R)
         HH=U+V*G
         H1=G*U-V
         IF(HH .EQ. 0) HH=U*CB
         IF(H1 .EQ. 0) H1=V*CB
         Z=R*H1
         R=HH/H1
        ELSE
         U=U+E/U
         V=V+HH
        ENDIF
        S=2*SQRT(E)
        E=S*T
        L=2*L
        IF(Y .LT. 0) L=L+1
        GO TO 3
       ENDIF
       IF(Y .LT. 0) L=L+1
       E=(ATAN(T/Y)+PI*L)*(C*T+D)/(T*(T+Q))
       IF(LBO) THEN
        H1=V/(T+U)
        Z=1-R*H1
        H1=R+H1
        IF(Z .EQ. 0) Z=CB
        IF(Z .LT. 0) THEN
         IF(H1 .GT. 0) M=M+1
         IF(H1 .LT. 0) M=M-1
        ENDIF
        S=ATAN(H1/Z)+M*PI
       ELSE
        IF(LBK) THEN
#if defined(CERNLIB_DOUBLE)
         S=HF*DASINH(YE)
#endif
#if !defined(CERNLIB_DOUBLE)
         S=HF*ASINH(YE)
#endif
        ELSE
         S=HF*(LOG(Z)+M*ALB)
        ENDIF
       ENDIF
       H=(E+SQRT(FA)*S)/N
       IF(X .LT. 0) H=-H
      ENDIF
#if defined(CERNLIB_DOUBLE)
    9 DELI3=H
#endif
#if !defined(CERNLIB_DOUBLE)
    9 RELI3=H
#endif
      RETURN
  103 FORMAT('FUNCTION SINGULAR FOR ',1P,
     1 'X = ',D15.8,3X,'P = ',D15.8,4X,'(P*X**2 = -1)')
      END
