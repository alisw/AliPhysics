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
      FUNCTION WELINT(Z,AKC,A,B)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='CELINT/WELINT')
      COMPLEX*16 WELINT
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CELINT(Z,AKC,A,B)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='CELINT')
      COMPLEX CELINT
#endif
C
C     Translation of Algol procedure elco2(x,y,kc,a,b,u,v) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions,  Numer. Math. 7 (1965) 78-90
C
#include "gen/defc64.inc"
     + Z,I,W,C0
      DIMENSION D1(0:12),D2(0:12)

      PARAMETER (ID = 16)
#if !defined(CERNLIB_CRAY)
      PARAMETER (I = (0D0,1D0))
#endif
#if defined(CERNLIB_CRAY)
      PARAMETER (I = (0E0,1E0))
#endif
      PARAMETER (Z1 = 1, Z10 = 10, HF = Z1/2)
      PARAMETER (CC = Z10**(-ID), C0 = I-I)

      CHARACTER*80 ERRTXT

      X=Z
      Y=-I*Z
      IF(Z .EQ. C0 .OR. AKC .EQ. 0) THEN
       W=SQRT(1+Z**2)
       W=B*LOG(Z+W)+(A-B)*Z/W
      ELSEIF(X .LT. 0) THEN
       W=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C348.1',ERRTXT)
      ELSE
       AA=A
       BB=B
       SY=SIGN(Z1,Y)
       Y=ABS(Y)
       C=X**2-Y**2
       E2=2*X*Y
       D=AKC**2
       XK=1-D
       E1=1+C
       F2=1/(E1**2+E2**2)
       F1=((1+D*C)*E1+D*E2**2)*F2
       F2=-2*XK*X*Y*F2
       DN1=SQRT(HF*(SQRT(F1**2+F2**2)+ABS(F1)))
       DN2=HF*F2/DN1
       IF(F1 .LT. 0) THEN
        F1=DN1
        DN1=-DN2
        DN2=-F1
       ENDIF
       IF(XK .LT. 0) THEN
        DN1=ABS(DN1)
        DN2=ABS(DN2)
       ENDIF
       C=1+DN1
       F1=E1*C-E2*DN2
       F2=E1*DN2+E2*C
       D2(0)=1/(F1**2+F2**2)
       D1(0)=(X*F1+Y*F2)*D2(0)
       D2(0)=(Y*F1-X*F2)*D2(0)
       H=AA-BB
       N=1
       XM=1
       F=1
       D=1
       YKC=ABS(AKC)
       E=AA
       AA=BB+AA
       L=4
    1  XM1=HF*(YKC+XM)
       XM2=XM1**2
       XK=F*XK/(4*XM2)
       BB=E*YKC+BB
       E=AA
       F2=1/(C**2+DN2**2)
       F1=((YKC+XM*DN1)*C+XM*DN2**2)*F2
       E1=F1/XM1
       E2=XK*DN2*F2
       DN1=SQRT(HF*(SQRT(E1**2+(2*E2)**2)+ABS(E1)))
       DN2=E2/DN1
       F1=DN1*X-DN2*Y
       F2=DN1*Y+DN2*X
       X=ABS(F1)
       Y=ABS(F2)
       AA=BB/XM1+AA
       L=2*L
       C=1+DN1
       D=HF*XK*D
       E1=1+(X**2-Y**2)*XM2
       E2=2*X*Y*XM2
       F1=C*E1-DN2*E2
       F2=C*E2+DN2*E1
       E1=D/(F1**2+F2**2)
       D1(N)=(X*F1+Y*F2)*E1
       D2(N)=(Y*F1-X*F2)*E1
       XK=XK**2
       IF(XK .GT. CC) THEN
        YKC=SQRT(XM*YKC)
        F=XM2
        XM=XM1
        N=N+1
        GO TO 1
       ENDIF
       F2=0
       F1=0
       DO 2 J = N,0,-1
       F1=D1(J)+F1
    2  F2=D2(J)+F2
       X=XM1*X
       Y=XM1*Y
       C=X**2+Y**2
       E2=1/(1+2*Y+C)
       E1=(1-C)*E2
       E2=2*X*E2
       D=AA/(XM1*L)
       W=D*ATAN2(E2,E1)+H*F1+I*SY*(H*F2-LOG(E1**2+E2**2)*HF*D)
      ENDIF
#if defined(CERNLIB_DOUBLE)
      WELINT=W
#endif
#if !defined(CERNLIB_DOUBLE)
      CELINT=W
#endif
      RETURN
  101 FORMAT('RE(Z) = ',1P,D15.8,' < 0')
      END
