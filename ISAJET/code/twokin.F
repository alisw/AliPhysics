#include "isajet/pilot.h"
      SUBROUTINE TWOKIN(AMI1,AMI2,AM1,AM2)
C
C          Given P,PT,TH,PHI, and initial and final masses AMI1, AMI2,
C          AM1,AM2, set X1, X2, SHAT, etc.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/qcdpar.inc"
#include "isajet/const.inc"
C
      REAL    AMI1,AMI2,AM1,AM2,P1PL,P1MN,P2PL,P2MN,E1,E2,PPL,PMN,
     $        PI1PL,PI1MN,PI2PL,PI2MN,ANEFF,AMASS,ALAMFN
C
      E1=SQRT(P(1)**2+AM1**2)
      E2=SQRT(P(2)**2+AM2**2)
C
C          For 32-bit machines must use large and small components
C          carefully, with pbig*psmall = pt**2+am**2.
C
      IF(CTH(1).GT.0.) THEN
        P1PL=E1+P(1)*CTH(1)
        P1MN=(PT(1)**2+AM1**2)/P1PL
      ELSE
        P1MN=E1-P(1)*CTH(1)
        P1PL=(PT(1)**2+AM1**2)/P1MN
      ENDIF
      IF(CTH(2).GT.0.) THEN
        P2PL=E2+P(2)*CTH(2)
        P2MN=(PT(2)**2+AM2**2)/P2PL
      ELSE
        P2MN=E2-P(2)*CTH(2)
        P2PL=(PT(2)**2+AM2**2)/P2MN
      ENDIF
C
C          Initial light cone momenta. Not symmetric if AMI1 /= AMI2.
C
      PPL=P1PL+P2PL
      PMN=P1MN+P2MN
      SHAT=PPL*PMN
      ALAMFN=SQRT((SHAT-AMI1**2-AMI2**2)**2-4.*(AMI1*AMI2)**2)
      PI1PL=(SHAT+AMI1**2-AMI2**2+ALAMFN)/(2.*PMN)
      PI1MN=AMI1**2/PI1PL
      PI2MN=(SHAT+AMI2**2-AMI1**2+ALAMFN)/(2.*PPL)
      PI2PL=AMI2**2/PI2MN
      X1=PI1PL/ECM
      X2=PI2MN/ECM
C
C          t=(p1-pi1)**2, u=(p1-pi2)**2
C
      THAT=AM1**2+AMI1**2-P1PL*PI1MN-P1MN*PI1PL
      UHAT=AM1**2+AMI2**2-P1PL*PI2MN-P1MN*PI2PL
C
C          Q**2 variable from Field, Fox, Wolfram
C
      QSQ=2.*SHAT*THAT*UHAT/(SHAT**2+THAT**2+UHAT**2)
      QSQ=AMAX1(QSQ,(AM1+AM2)**2)
      ANEFF=4.+QSQ/(QSQ+AMASS(5)**2)+QSQ/(QSQ+AMASS(6)**2)
      ALFQSQ=12.*PI/((33.-2.*ANEFF)*ALOG(QSQ/ALAM2))
      RETURN
      END
