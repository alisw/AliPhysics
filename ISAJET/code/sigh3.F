#include "isajet/pilot.h"
      SUBROUTINE SIGH3
C
C          Calculate angular distributions for W decays from Higgs,
C          d(sigma)/d(qmw**2)d(yw)d(omega)d(omega1)d(omega2)
C
C          Ver 7.14: Only modification needed for MSSM is to check
C          GOMSSM flag instead of INITYP
C
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/wsig.inc"
#include "isajet/wwsig.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/wwpar.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
C
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      DIMENSION IDADDR(4),IW(2),LAM(3),LISTJ(29)
     $,T12(3,3),T34(3,3),FTERM(4),FR(3,3),FI(3,3)
     $,CPHI12(3),SPHI12(3),CPHI34(3),SPHI34(3)
      DIMENSION PFCM(5,4),PWCM(5,2)
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION TERM,FTERM,ZCM
#endif
      DATA LAM/0,1,-1/
      DATA LISTJ/
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     $10,80,-80,90/
C
C          FUNCTIONS
      DOTP(I,J)=PPAIR(4,I)*PPAIR(4,J)-PPAIR(1,I)*PPAIR(1,J)
     $-PPAIR(2,I)*PPAIR(2,J)-PPAIR(3,I)*PPAIR(3,J)
C
C          ENTRY
      IF(NPAIR.NE.4) RETURN
C
C          RECONSTRUCT W-->FF DECAY ANGLES
C
C          INITIALIZE PFCM AND PWCM
      DO 10 I=1,4
      DO 10 K=1,5
10    PFCM(K,I)=PPAIR(K,I)
      DO 11 I=1,2
      DO 11 K=1,5
11    PWCM(K,I)=PJETS(K,I)
C
C          Z BOOST TO WW CENTER OF MASS
      CHWW=QWJET(4)/QWJET(5)
      SHWW=QWJET(3)/QWJET(5)
      DO 20 I=1,4
      TMP=CHWW*PFCM(4,I)-SHWW*PFCM(3,I)
      PFCM(3,I)=-SHWW*PFCM(4,I)+CHWW*PFCM(3,I)
20    PFCM(4,I)=TMP
      DO 21 I=1,2
      TMP=CHWW*PWCM(4,I)-SHWW*PWCM(3,I)
      PWCM(3,I)=-SHWW*PWCM(4,I)+CHWW*PWCM(3,I)
21    PWCM(4,I)=TMP
C
C          ROTATE W1 TO +Z AXIS
      PTW1=SQRT(PWCM(1,1)**2+PWCM(2,1)**2)
      CPHIW1=PWCM(1,1)/PTW1
      SPHIW1=PWCM(2,1)/PTW1
      PW1=SQRT(PTW1**2+PWCM(3,1)**2)
      CTHW1=PWCM(3,1)/PW1
      STHW1=PTW1/PW1
C          Z ROTATION
      DO 30 I=1,4
      TMP=CPHIW1*PFCM(1,I)+SPHIW1*PFCM(2,I)
      PFCM(2,I)=-SPHIW1*PFCM(1,I)+CPHIW1*PFCM(2,I)
30    PFCM(1,I)=TMP
C          Y ROTATION
      DO 31 I=1,4
      TMP=CTHW1*PFCM(1,I)-STHW1*PFCM(3,I)
      PFCM(3,I)=STHW1*PFCM(1,I)+CTHW1*PFCM(3,I)
31    PFCM(1,I)=TMP
C
C          BOOST TO W REST FRAMES
      CHW1=PWCM(4,1)/PWCM(5,1)
      SHW1=PW1/PWCM(5,1)
      DO 40 I=1,4
      IF(I.LE.2) THEN
        SHWI=SHW1
      ELSE
        SHWI=-SHW1
      ENDIF
      TMP=CHW1*PFCM(4,I)-SHWI*PFCM(3,I)
      PFCM(3,I)=-SHWI*PFCM(4,I)+CHW1*PFCM(3,I)
40    PFCM(4,I)=TMP
C
C          COMPUTE ANGLES
      TH12=ACOS(PFCM(3,1)/SQRT(PFCM(1,1)**2+PFCM(2,1)**2+PFCM(3,1)**2))
      PHI12=ATAN2(PFCM(2,1),PFCM(1,1))
      TH34=ACOS(PFCM(3,3)/SQRT(PFCM(1,3)**2+PFCM(2,3)**2+PFCM(3,3)**2))
      PHI34=ATAN2(PFCM(2,3),PFCM(1,3))
C
C          COMPUTE DECAY ANGULAR DISTRIBUTIONS.
C
      DO 100 I=1,4
      IDADDR(I)=IABS(IDPAIR(I))
100   IF(IDADDR(I).GE.11) IDADDR(I)=IDADDR(I)-4
      IF(GOMSSM) THEN
        IW(1)=JETTYP(1)-76
        IW(2)=JETTYP(2)-76
      ELSE
        IW(1)=JETTYP(1)-25
        IW(2)=JETTYP(2)-25
      ENDIF
C
      AMV=PJETS(5,1)
      GAMV=WGAM(IW(1))
      QMH=QMW
C          COUPLINGS
      A12=AQ(IDADDR(1),IW(1))
      B12=BQ(IDADDR(1),IW(1))
      A34=AQ(IDADDR(3),IW(2))
      B34=BQ(IDADDR(3),IW(2))
C          DECAY DISTRIBUTIONS
      TVV12=8.*PI*ALFA*(A12**2+B12**2)
      TVA12=16.*PI*ALFA*A12*B12
      COS12=COS(TH12)
      SIN12=SIN(TH12)
      T12(1,1)=TVV12*SIN12**2
      T12(1,2)=TVV12*SIN12*COS12/SQRT2+TVA12*SIN12/SQRT2
      T12(1,3)=-TVV12*SIN12*COS12/SQRT2+TVA12*SIN12/SQRT2
      T12(2,1)=T12(1,2)
      T12(2,2)=TVV12*(.5+.5*COS12**2)+TVA12*COS12
      T12(2,3)=TVV12*.5*SIN12**2
      T12(3,1)=T12(1,3)
      T12(3,2)=T12(2,3)
      T12(3,3)=TVV12*(.5+.5*COS12**2)-TVA12*COS12
C
      TVV34=8.*PI*ALFA*(A34**2+B34**2)
      TVA34=16.*PI*ALFA*A34*B34
      COS34=COS(TH34)
      SIN34=SIN(TH34)
      T34(1,1)=TVV34*SIN34**2
      T34(1,2)=TVV34*SIN34*COS34/SQRT2+TVA34*SIN34/SQRT2
      T34(1,3)=-TVV34*SIN34*COS34/SQRT2+TVA34*SIN34/SQRT2
      T34(2,1)=T34(1,2)
      T34(2,2)=TVV34*(.5+.5*COS34**2)+TVA34*COS34
      T34(2,3)=TVV34*.5*SIN34**2
      T34(3,1)=T34(1,3)
      T34(3,2)=T34(2,3)
      T34(3,3)=TVV34*(.5+.5*COS34**2)-TVA34*COS34
C
      CPHI12(1)=1.
      CPHI12(2)=COS(PHI12)
      CPHI12(3)=COS(2.*PHI12)
      SPHI12(1)=0.
      SPHI12(2)=SIN(PHI12)
      SPHI12(3)=SIN(2.*PHI12)
      CPHI34(1)=1.
      CPHI34(2)=COS(PHI34)
      CPHI34(3)=COS(2.*PHI34)
      SPHI34(1)=0.
      SPHI34(2)=SIN(PHI34)
      SPHI34(3)=SIN(2.*PHI34)
C
      TCPHI=CPHI12(2)*CPHI34(2)-SPHI12(2)*SPHI34(2)
      TSPHI=SPHI12(2)*CPHI34(2)+CPHI12(2)*SPHI34(2)
      TC2PHI=CPHI12(3)*CPHI34(3)-SPHI12(3)*SPHI34(3)
      TS2PHI=SPHI12(3)*CPHI34(3)+CPHI12(3)*SPHI34(3)
C
C          PURE HIGGS --> W W. CALCULATE ANGULAR DISTRIBUTION FOR
C          HIGGS DECAY AND MULTIPLY BY CROSS SECTION.
C
      IF(INITYP(1).LE.25.OR.GOMSSM) THEN
        F0=.5*QMH**2/AMV**2-1.
        F1=1.
        TOTAL=(8.*PI/3.)**2*TVV12*TVV34*(F0**2+2.*F1**2)
        DIFF=F0**2*T12(1,1)*T34(1,1)
     $  +F0*F1*(2.*T12(1,2)*T34(1,2)+2.*T12(1,3)*T34(1,3))*TCPHI
     $  +F1**2*(T12(2,2)*T34(1,2)+T12(3,3)*T34(3,3)
     $    +2.*T12(2,3)*T34(2,3)*TC2PHI)
        WWSIG=SIGLLQ*DIFF/TOTAL
        RETURN
      ENDIF
C
C          W W FUSION. CALCULATE ANGULAR DISTRIUBTION FOR DECAY
C          INCLUDING ALL GRAPHS.
C
C          KINEMATICS
      IFL1=LISTJ(JETTYP(1))
      IFL2=LISTJ(JETTYP(2))
      IFIN1=LISTJ(INITYP(1))
      IFIN2=LISTJ(INITYP(2))
      WMF=AMASS(IFL1)
      WMI=AMASS(IFIN1)
      PINPF=SQRT((S-4.*WMI**2)*(S-4.*WMF**2))
      ZCM=(.5*S+T-WMI**2-WMF**2)/(.5*PINPF)
C          PRODUCTION AMPLITUDES. REMEMBER MISSING SIN(THETA)/SQRT(2)
      DO 110 L=1,4
      FTERM(L)=0.
      DO 120 J=1,4
      TERM=0.
      DO 130 I=1,4
130   TERM=TERM+ANWWWW(I,J,L)*ZCM**(I-1)
      TERM=TERM/(ADWWWW(1,J)+ADWWWW(2,J)*ZCM)
120   FTERM(L)=FTERM(L)+TERM
110   CONTINUE
      FTERM(4)=FTERM(4)*SQRT(ABS(1.-ZCM**2))/SQRT2
C          HELICITY AMPLITUDES. NOTATION IS 0,+,-
      FR(1,1)=FTERM(1)
      FI(1,1)=AIWWWW(1)
      FR(1,2)=FTERM(4)
      FI(1,2)=AIWWWW(4)
      FR(2,2)=FTERM(3)
      FI(2,2)=AIWWWW(3)
      FR(2,3)=FTERM(2)
      FI(2,3)=AIWWWW(2)
C
      FR(1,3)=FR(1,2)
      FI(1,3)=FI(1,2)
      FR(3,1)=FR(1,3)
      FI(3,1)=FI(1,3)
      FR(2,1)=FR(1,2)
      FI(2,1)=FI(1,2)
C
      FR(3,3)=FR(2,2)
      FI(3,3)=FI(2,2)
      FR(3,2)=FR(2,3)
      FI(3,2)=FI(2,3)
C
C          DIFFERENTIAL DISTRIBUTION FROM DENSITY MATRIX
      DIFF=0.
      DO 140 I1=1,3
      L1=LAM(I1)
      DO 140 I2=1,3
      L2=LAM(I2)
      DO 140 I3=1,3
      L3=LAM(I3)
      DO 140 I4=1,3
      L4=LAM(I4)
      L12=L1-L2
      I12=IABS(L12)+1
      IF(I12.EQ.0) I12=3
      L34=L3-L4
      I34=IABS(L34)+1
      IF(I34.EQ.0) I34=3
      C1234=CPHI12(I12)*CPHI34(I34)
     $-SPHI12(I12)*ISIGN(1,L12)*SPHI34(I34)*ISIGN(1,L34)
      S1234=SPHI12(I12)*ISIGN(1,L12)*CPHI34(I34)
     $+CPHI12(I12)*SPHI34(I34)*ISIGN(1,L34)
      DIFF=DIFF+(FR(I1,I2)*FR(I3,I4)+FI(I1,I2)*FI(I3,I4))
     $*T12(I3,I1)*T34(I4,I2)*C1234
     $+(FR(I1,I2)*FI(I3,I4)-FI(I1,I2)*FR(I3,I4))
     $*T12(I3,I1)*T34(I4,I2)*S1234
140   CONTINUE
C          INTEGRATED DISTRIBUTION
      TOTAL=0.
      DO 150 I1=1,3
      DO 150 I2=1,3
      TOTAL=TOTAL+FR(I1,I2)**2+FI(I1,I2)**2
150   CONTINUE
      FAC=(16.*PI/3.*4.*PI*ALFA)**2
      FAC=FAC*(A12**2+B12**2)*(A34**2+B34**2)
      TOTAL=TOTAL*FAC
      WWSIG=DIFF/TOTAL*SIGLLQ
      RETURN
      END
