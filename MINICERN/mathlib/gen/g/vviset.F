*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE VVISET(RKA,BE2,MODE,XL,XU)

      COMMON /G116C1/ H(7),T0,T1,T,OMEGA,A(155),B(155),X0
      EXTERNAL G116F1,G116F2
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'VVISET')
      DIMENSION XP(8),XQ(6)

C     H1 = 5 ln 10 - ln(2/PI**2), h2 = 3 ln 10, h0 = Ln 0.0005...

      PARAMETER (PI = 3.14159 265, EU = 0.57721 566)
      PARAMETER (PI2 = 2*PI, RPI = 1/PI)
      PARAMETER (H1 = 5*2.30258 509-1.59631 259, H2 = 3*2.30258 509)
      PARAMETER (H0 =-7.6, PIH = PI/2, EPS = 1E-5)

      DATA  XP /9.29,  2.47, 0.89, 0.36, 0.15, 0.07, 0.03, 0.02/
      DATA  XQ /0.012, 0.03, 0.08, 0.26, 0.87, 3.83/

      IF(RKA .LT. 0.01 .OR. RKA .GT. 10) THEN
       WRITE(ERRTXT,101) RKA
       CALL MTLPRT(NAME,'G116.1',ERRTXT)
      ELSEIF(BE2 .LT. 0 .OR. BE2 .GT. 1) THEN
       WRITE(ERRTXT,102) BE2
       CALL MTLPRT(NAME,'G116.2',ERRTXT)
      ELSE
       H(5)=1-BE2*(1-EU)-H0/RKA
       H(6)=BE2
       H(7)=1-BE2
       H4=H0/RKA-(1+BE2*EU)
       H5=LOG(RKA)
       H6=1/RKA
       T0=(H4-H(5)*H5-(H(5)+BE2)*(LOG(H(5))+REXPIN(H(5)))+EXP(-H(5)))/
     1     H(5)
       DO 1 LP = 1,8
       IF(RKA .GE. XP(LP)) GO TO 11
    1  CONTINUE
       LP=9
   11  DO 2 LQ = 1,6
       IF(RKA .LE. XQ(LQ)) GO TO 22
    2  CONTINUE
       LQ=7
   22  CALL RZERO(-LP-0.5,LQ-7.5,U,XX,EPS,1000,G116F2)
       Q=1/U
       T1=H4*Q-H5-(1+BE2*Q)*(LOG(ABS(U))+REXPIN(U))+EXP(-U)*Q

       T=T1-T0
       OMEGA=PI2/T
       H(1)=RKA*(2+BE2*EU)+H1
       IF(RKA .GE. 0.07) H(1)=H(1)+H2
       H(2)=BE2*RKA
       H(3)=H6*OMEGA
       H(4)=PIH*OMEGA
       CALL RZERO(5.,155.,X0,XX,EPS,1000,G116F1)
       N=X0+1

       D=RPI*EXP(RKA*(1+BE2*(EU-H5)))
       A(N)=0
       IF(MODE .EQ. 0) A(N)=RPI*OMEGA
       Q=-1
       Q2=2
       DO 3 K = 1,N-1
       L=N-K
       X=OMEGA*K
       X1=H6*X
       C1=LOG(X)-COSINT(X1)
       C2=SININT(X1)
       C3=SIN(X1)
       C4=COS(X1)
       XF1=RKA*(BE2*C1-C4)-X*C2
       XF2=X*C1+RKA*(C3+BE2*C2)+T0*X
       IF(MODE .EQ. 0) THEN
        D1=Q*D*OMEGA*EXP(XF1)
        A(L)=D1*COS(XF2)
        B(L)=-D1*SIN(XF2)
       ELSE
        D1=Q*D*EXP(XF1)/K
        A(L)=D1*SIN(XF2)
        B(L)=D1*COS(XF2)
        A(N)=A(N)+Q2*A(L)
       ENDIF
       Q=-Q
    3  Q2=-Q2
      ENDIF
      XL=T0
      XU=T1
      RETURN
  101 FORMAT('KAPPA =   ',E10.3,' - OUT OF RANGE')
  102 FORMAT('BETA**2 = ',E10.3,' - OUT OF RANGE')
      END
