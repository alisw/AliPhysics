*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:04  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RFCONC(X,TAU,M)
#include "gen/defc64.inc"
     +           CGAMMA,CLOGAM
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DFCONC(X,TAU,M)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +           WGAMMA,WLOGAM
#endif
C
C     Computes the (real) Conical Function of the first kind
C                 P M (-1/2 + I*TAU) (X)
C     for M = 0 and M = 1. P m nu (x) is the Legendre function of
C     the first kind.
C     Based on K.S. Koelbig, A program for computing the conical
C     functions of the first kind P ... for m = 0 and m = 1,
C     Computer Phys. Comm. 23 (1981) 51-61.
C
#include "gen/defc64.inc"
     +     CGM,CLG,CRG,I,A,B,C,TI,R,RR,U(0:3),V(0:3),W(19)
      LOGICAL LM0,LM1,LTA
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RFCONC')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RFCONC/DFCONC')
#endif
      DIMENSION T(7),H(9),S(5),P(11),D(-1:6)


      PARAMETER (PI  = 3.14159 26535 89793 24D+0)
      PARAMETER (RPI = 1.77245 38509 05516 03D+0)
      PARAMETER (I = (0,1))
      PARAMETER (Z1 = 1, HF = Z1/2, TH = 1+HF, C1 = Z1/10, C2 = Z1/5)
      PARAMETER (RPH = 2/PI, RPW = 2/RPI, TW = 20, NMAX = 200)

      DATA EPS /1D-14/

#if !defined(CERNLIB_DOUBLE)
      ENTRY FCONC(X,TAU,M)
#endif

#if defined(CERNLIB_DOUBLE)
      FEK(ARG)=DELIKC(ARG)
      FEE(ARG)=DELIEC(ARG)
      FJ0(ARG)=DBESJ0(ARG)
      FJ1(ARG)=DBESJ1(ARG)
      FI0(ARG)=DBESI0(ARG)
      FI1(ARG)=DBESI1(ARG)
      CGM(CRG)=WGAMMA(CRG)
      CLG(CRG)=WLOGAM(CRG)
#endif
#if !defined(CERNLIB_DOUBLE)
      FEK(ARG)=RELIKC(ARG)
      FEE(ARG)=RELIEC(ARG)
      FJ0(ARG)= BESJ0(ARG)
      FJ1(ARG)= BESJ1(ARG)
      FI0(ARG)= BESI0(ARG)
      FI1(ARG)= BESI1(ARG)
      CGM(CRG)=CGAMMA(CRG)
      CLG(CRG)=CLOGAM(CRG)
#endif
      LM0=M .EQ. 0
      LM1=M .EQ. 1
      IF(X .LE. -1 .OR. TAU .LT. 0 .OR. .NOT.(LM0 .OR. LM1)) THEN
       FC=0
       WRITE(ERRTXT,101) X,TAU,M
       CALL MTLPRT(NAME,'C331.1',ERRTXT)
       GO TO 99
      END IF
      FC=1-M
      IF(X .EQ. 1) GO TO 99

      IF(TAU .EQ. 0) THEN
       IF(X .LT. 1) THEN
        Y=SQRT(HF*(1-X))
        IF(LM0) THEN
         FC=RPH*FEK(Y)
        ELSE
         FC=RPH*(FEE(Y)-HF*(1+X)*FEK(Y))/SQRT(1-X**2)
        END IF
       ELSE
        Y=SQRT((X-1)/(X+1))
        IF(LM0) THEN
         FC=RPH*FEK(Y)/SQRT(HF*(X+1))
        ELSE
         FC=RPH*SQRT(HF/(X-1))*(FEE(Y)-FEK(Y))
        END IF
       END IF
       GO TO 99
      END IF

      TI=I*TAU
      FM=M
      IF(-1 .LT. X .AND. X .LE. 0) GO TO 11
      IF(0  .LT. X .AND. X .LE. C1 .AND. TAU .LE. 17) GO TO 11
      IF(C1 .LT. X .AND. X .LE. C2 .AND. TAU .LE.  5) GO TO 11
      IF(C1 .LT. X .AND. X .LE. C2 .AND. TAU .LE. 17) GO TO 12
      IF(C2 .LT. X .AND. X .LE. TH .AND. TAU .LE. 20) GO TO 12
      IF(TH .LT. X .AND. TAU .LE. MAX(TW,X)) GO TO 13

      IF(X .LT. 1) THEN
       Y=SQRT(1-X**2)
       T(1)=ACOS(X)
       H(1)=TAU*T(1)
       B0=FI0(H(1))
       B1=FI1(H(1))
       Z=-1
      ELSE
       Y=SQRT(X**2-1)
       T(1)=LOG(X+Y)
       H(1)=TAU*T(1)
       B0=FJ0(H(1))
       B1=FJ1(H(1))
       Z=1
      END IF
      H(1)=T(1)*X/Y
      P(1)=1/TAU
      S(1)=1/T(1)
      DO 5 J = 2,7
      T(J)=T(J-1)*T(1)
    5 H(J)=H(J-1)*H(1)
      DO 6 J = 2,11
    6 P(J)=P(J-1)*P(1)
      DO 7 J = 2,5
    7 S(J)=S(J-1)*S(1)
C
      IF(LM0) THEN
       D(-1)=0
       D(0)=1
       D(1)=(H(1)-1)/(8*T(1))
       D(2)=(9*H(2)+6*H(1)-15-Z*8*T(2))/(128*T(2))
       D(3)=5*(15*H(3)+27*H(2)+21*H(1)-63-Z*T(2)*(16*H(1)+24))/
     1  (1024*T(3))
       D(4)=7*(525*H(4)+1500*H(3)+2430*H(2)+1980*H(1)-6435
     1 +192*T(4)-Z*T(2)*(720*H(2)+1600*H(1)+2160))/(32768*T(4))
       D(5)=21*(2835*H(5)+11025*H(4)+24750*H(3)+38610*H(2)
     1  +32175*H(1)-109395+T(4)*(1984*H(1)+4032)
     2  -Z*T(2)*(4800*H(3)+15120*H(2)+26400*H(1)+34320))/
     3  (262144*T(5))
       D(6)=11*(218295*H(6)+1071630*H(5)+3009825*H(4)+6142500*
     1  H(3)+9398025*H(2)+7936110*H(1)-27776385+T(4)*(254016*H(2)
     2  +749952*H(1)+1100736)-Z*T(2)*(441000*H(4)+1814400*H(3)
     3  +4127760*H(2)+6552000*H(1)+8353800+31232*T(4)))/
     4  (4194304*T(6))
      ELSE
       D(-1)=-1
       D(0)=3*(1-H(1))/(8*T(1))
       D(1)=(-15*H(2)+6*H(1)+9+Z*8*T(2))/(128*T(2))
       D(2)=3*(-35*H(3)-15*H(2)+15*H(1)+35+Z*T(2)*(32*H(1)+8))/
     1  (1024*T(3))
       D(3)=(-4725*H(4)-6300*H(3)-3150*H(2)+3780*H(1)+10395
     1   -1216*T(4)+Z*T(2)*(6000*H(2)+5760*H(1)+1680))/(32768*T(4))
       D(4)=7*(-10395*H(5)-23625*H(4)-28350*H(3)-14850*H(2)
     1  +19305*H(1)+57915-T(4)*(6336*H(1)+6080)+Z*T(2)*
     2  (16800*H(3)+30000*H(2)+25920*H(1)+7920))/(262144*T(5))
       D(5)=(-2837835*H(6)-9168390*H(5)-16372125*H(4)-18918900*H(3)
     1  -10135125*H(2)+13783770*H(1)+43648605-T(4)*(3044160*H(2)
     2  +5588352*H(1)+4213440)+Z*T(2)*(5556600*H(4)+14817600*H(3)
     3  +20790000*H(2)+17297280*H(1)+5405400+323072*T(4)))/
     4  (4194304*T(6))
       D(6)=0
      END IF
      S0=D(0)+(-4*D(3)*S(1)+D(4))*P(4)+
     1 (-192*D(5)*S(3)+144*D(6)*S(2))*P(8)
     2 +Z*(-D(2)*P(2)+(-24*D(4)*S(2)+12*D(5)*S(1)-D(6))*P(6)
     3 +(-1920*D(6)*S(4))*P(10))
      S1=D(1)*P(1)+(8*(D(3)*S(2)-D(4)*S(1))+D(5))*P(5)
     1 +(384*D(5)*S(4)-768*D(6)*S(3))*P(9)
     2 +Z*(D(-1)*TAU+(2*D(2)*S(1)-D(3))*P(3)+(48*D(4)*S(3)
     3 -72*D(5)*S(2)+18*D(6)*S(1))*P(7)+(3840*D(6)*S(5))*P(11))
      FC=SQRT(T(1)/Y)*(B0*S0+B1*S1)
      GO TO 99

   11 LTA=TAU .LE. 10
      X1=X**2
      A=HF*((HF-FM)-TI)
      B=HF*((HF-FM)+TI)
      C=HF
      ASSIGN 1 TO JP
      GO TO 20
    1 R1=R
      R1=R1/ABS(CGM(A+HF))**2
      A=HF*((TH-FM)-TI)
      B=HF*((TH-FM)+TI)
      C=TH
      ASSIGN 2 TO JP
      GO TO 20
    2 R2=R
      FC=RPI*(R1-2*X*R2/ABS(CGM(A-HF))**2)
      IF(LM1) FC=2*FC/SQRT(1-X1)
      GO TO 99

   12 LTA=X .GT. 1 .OR. X .LE. 1 .AND. TAU .LE. 5
      X1=HF*(1-X)
      A=(HF+FM)-TI
      B=(HF+FM)+TI
      C=FM+1
      ASSIGN 3 TO JP
      GO TO 20
    3 FC=R
      IF(LM1) FC=SIGN(HF,1-X)*(TAU**2+HF**2)*SQRT(ABS(X**2-1))*FC
      GO TO 99

   13 LTA=.TRUE.
      X1=1/X**2
      A=HF*((HF-FM)-TI)
      B=HF*((TH-FM)-TI)
      C=1-TI
      ASSIGN 4 TO JP
      GO TO 20
    4 R1=EXP((TI-HF)*LOG(X+X)+CLG(1+TI)-CLG((TH-FM)+TI))*
     1        R*((HF-FM)+TI)/TI
      FC=RPW*R1
      IF(LM1) FC=FC/SQRT(1-X1)
      GO TO 99

   20 IF(LTA) THEN
       Y=-X1
       Y2=Y**2
       Y3=Y*Y2
       W(1)=A+1
       W(2)=A+2
       W(3)=B+1
       W(4)=B+2
       W(5)=C+1
       W(6)=C*W(5)
       W(7)=A+B
       W(8)=A*B
       W(9)=(W(8)/C)*Y
       W(10)=W(1)*W(3)
       W(11)=W(2)*W(4)
       W(12)=1+(W(11)/(W(5)+W(5)))*Y
       W(13)=W(7)-6
       W(14)=W(7)+6
       W(15)=2-W(8)
       W(16)=W(15)-W(7)-W(7)

       V(0)=1
       V(1)=1+(W(10)/(C+C))*Y
       V(2)=W(12)+(W(10)*W(11)/(12*W(6)))*Y2
       U(0)=1
       U(1)=V(1)-W(9)
       U(2)=V(2)-W(9)*W(12)+(W(8)*W(10)/(W(6)+W(6)))*Y2

       R=1
       DO 21 N = 3,NMAX
       FN=N
       RR=R
       H(1)=FN-1
       H(2)=FN-2
       H(3)=FN-3
       H(4)=FN+FN
       H(5)=H(4)-3
       H(6)=H(5)+H(5)
       H(7)=4*(H(4)-1)*H(5)
       H(8)=8*H(5)**2*(H(4)-5)
       H(9)=3*FN**2
       W(1)=A+H(1)
       W(2)=A+H(2)
       W(3)=B+H(1)
       W(4)=B+H(2)
       W(5)=C+H(1)
       W(6)=C+H(2)
       W(7)=C+H(3)
       W(8)=H(2)-A
       W(9)=H(2)-B
       W(10)=H(1)-C
       W(11)=W(1)*W(3)
       W(12)=W(5)*W(6)

       W(17)=1+((H(9)+W(13)*FN+W(16))/(H(6)*W(5)))*Y
       W(18)=-((W(11)*W(10)/H(6)+(H(9)-W(14)*FN+W(15))*W(11)*Y/H(7))/
     1          W(12))*Y
       W(19)=(W(2)*W(11)*W(4)*W(8)*W(9)/(H(8)*W(7)*W(12)))*Y3
       V(3)=W(17)*V(2)+W(18)*V(1)+W(19)*V(0)
       U(3)=W(17)*U(2)+W(18)*U(1)+W(19)*U(0)
       R=U(3)/V(3)
       IF(ABS(R-RR) .LT. EPS) GO TO JP, (1,2,3,4)
       DO 22 J = 1,3
       V(J-1)=V(J)
   22  U(J-1)=U(J)
   21  CONTINUE
      ELSE
       W(1)=X1*A*B/C
       R=1+W(1)
       DO 23 N = 1,NMAX
       FN=N
       RR=R
       W(1)=W(1)*X1*(A+FN)*(B+FN)/((C+FN)*(FN+1))
       R=R+W(1)
       IF(ABS(R-RR) .LT. EPS) GO TO JP, (1,2,3,4)
   23  CONTINUE
      END IF
      FC=0
      WRITE(ERRTXT,102) X
      CALL MTLPRT(NAME,'C331.2',ERRTXT)
#if defined(CERNLIB_DOUBLE)
   99 DFCONC=FC
#endif
#if !defined(CERNLIB_DOUBLE)
   99 RFCONC=FC
#endif
      RETURN

  101 FORMAT('ILLEGAL ARGUMENT(S)  X = ',D15.8,' TAU = ',D15.8,
     1       ' M = ',I3)
  102 FORMAT('CONVERGENCE PROBLEM FOR HYPERGEOMETRIC FUNCTION, X = ',
     1        D15.8)
      END
