#include "isajet/pilot.h"
      SUBROUTINE SIGTC3
C
C          Calculate angular distributions for W decays from technirho:
C          d(sigma)/d(qmw**2)d(yw)d(omega)d(omega1)d(omega2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
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
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/wwpar.inc"
#include "isajet/tcpar.inc"
C
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      INTEGER I,K,IDADDR(4),IW(2)
      REAL T12(3,3),T34(3,3),FR(3,3),FI(3,3),CPHI12(3),SPHI12(3),
     $CPHI34(3),SPHI34(3),PFCM(5,4),PWCM(5,2),CHWW,SHWW,TMP,PTW1,
     $CPHIW1,SPHIW1,PW1,CTHW1,STHW1,CHW1,SHW1,SHWI,TH12,PHI12,TH34,
     $PHI34,AMV,GAMV,QMH,A12,B12,A34,B34,TVV12,TVA12,COS12,SIN12,
     $TVV34,TVA34,COS34,SIN34,TCPHI,TSPHI,TC2PHI,TS2PHI,F0,F1,TOTAL,
     $DIFF,T,U,S
C
      IF(NPAIR.NE.4) RETURN
C
C          Reconstruct W-->FF decay angles
C
C          Initialize PFCM and PWCM
      DO 10 I=1,4
      DO 10 K=1,5
        PFCM(K,I)=PPAIR(K,I)
10    CONTINUE
      DO 11 I=1,2
      DO 11 K=1,5
        PWCM(K,I)=PJETS(K,I)
11    CONTINUE
C
C          Z boost to WW center of mass
      CHWW=QWJET(4)/QWJET(5)
      SHWW=QWJET(3)/QWJET(5)
      DO 20 I=1,4
        TMP=CHWW*PFCM(4,I)-SHWW*PFCM(3,I)
        PFCM(3,I)=-SHWW*PFCM(4,I)+CHWW*PFCM(3,I)
        PFCM(4,I)=TMP
20    CONTINUE
      DO 21 I=1,2
        TMP=CHWW*PWCM(4,I)-SHWW*PWCM(3,I)
        PWCM(3,I)=-SHWW*PWCM(4,I)+CHWW*PWCM(3,I)
        PWCM(4,I)=TMP
21    CONTINUE
C
C          Rotate W1 to +z axis
      PTW1=SQRT(PWCM(1,1)**2+PWCM(2,1)**2)
      CPHIW1=PWCM(1,1)/PTW1
      SPHIW1=PWCM(2,1)/PTW1
      PW1=SQRT(PTW1**2+PWCM(3,1)**2)
      CTHW1=PWCM(3,1)/PW1
      STHW1=PTW1/PW1
C          Z rotation
      DO 30 I=1,4
        TMP=CPHIW1*PFCM(1,I)+SPHIW1*PFCM(2,I)
        PFCM(2,I)=-SPHIW1*PFCM(1,I)+CPHIW1*PFCM(2,I)
        PFCM(1,I)=TMP
30    CONTINUE
C          Y rotation
      DO 31 I=1,4
        TMP=CTHW1*PFCM(1,I)-STHW1*PFCM(3,I)
        PFCM(3,I)=STHW1*PFCM(1,I)+CTHW1*PFCM(3,I)
        PFCM(1,I)=TMP
31    CONTINUE
C
C          Boost to W rest frames
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
        PFCM(4,I)=TMP
40    CONTINUE
C
C          Compute angles
      TH12=ACOS(PFCM(3,1)/SQRT(PFCM(1,1)**2+PFCM(2,1)**2+PFCM(3,1)**2))
      PHI12=ATAN2(PFCM(2,1),PFCM(1,1))
      TH34=ACOS(PFCM(3,3)/SQRT(PFCM(1,3)**2+PFCM(2,3)**2+PFCM(3,3)**2))
      PHI34=ATAN2(PFCM(2,3),PFCM(1,3))
C
C          Compute decay angular distributions.
C
      DO 100 I=1,4
        IDADDR(I)=IABS(IDPAIR(I))
        IF(IDADDR(I).GE.11) IDADDR(I)=IDADDR(I)-4
100   CONTINUE
      IW(1)=JETTYP(1)-25
      IW(2)=JETTYP(2)-25
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
C          Pure technirho --> WW. Calculate angular distribution for
C          decay and multiply by cross section.
C
      F0=.5*QMH**2/AMV**2-1.
      F1=1.
      TOTAL=(8.*PI/3.)**2*TVV12*TVV34*(F0**2+2.*F1**2)
      DIFF=F0**2*T12(1,1)*T34(1,1)
     $+F0*F1*(2.*T12(1,2)*T34(1,2)+2.*T12(1,3)*T34(1,3))*TCPHI
     $+F1**2*(T12(2,2)*T34(1,2)+T12(3,3)*T34(3,3)
     $  +2.*T12(2,3)*T34(2,3)*TC2PHI)
      WWSIG=SIGLLQ*DIFF/TOTAL
      RETURN
      END
