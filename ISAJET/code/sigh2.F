#include "isajet/pilot.h"
      SUBROUTINE SIGH2
C
C          COMPUTE THE WEINBERG-SALAM HIGGS CROSS SECTION
C          D(SIGMA)/D(QMW**2)D(YW)D(OMEGA)
C          FOR THE SPECIFIED JET TYPES. TRIVIAL EXCEPT FOR W W FUSION,
C          WHICH HAS INTERFERENCE WITH W W SCATTERING.
C
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/wsig.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/hcon.inc"
C
      DIMENSION X(2),LISTJ(29),WTHELI(4)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT),(X(1),X1)
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION C,TERM,SUM,DENOM,ZCM
#endif
C
C          WTHELI ARE WEIGHTS OF HELICITY AMPLITUDES IN SIGMA.
      DATA WTHELI/1.,2.,2.,4./
      DATA LISTJ/
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     $10,80,-80,90/
C
C          QUARK OR GLUON FUSION TO HIGGS
C
      IF(INITYP(1).LE.25) THEN
        SIGLLQ=SIGEVT/(4.*PI)
        RETURN
      ENDIF
C
C          W+W FUSION AND W+W->W+W IN EFFECTIVE W APPROXIMATION.
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
C          RESET COEFFICIENTS FOR SELECTED PROCESS
      IABSI=IABS(IFIN1)
      IABSF=IABS(IFL1)
      IF(IABSI.EQ.80) THEN
        IF(IABSF.EQ.80) THEN
          CALL XWWWW
        ELSE
          CALL XWWZZ
        ENDIF
      ELSE
        IF(IABSF.EQ.80) THEN
          CALL XZZWW
        ELSE
          CALL XZZZZ
        ENDIF
      ENDIF
C          SUM CROSS SECTION TERMS. I,J RUN OVER AMPLITUDE TERMS.
C          L RUNS OVER HELICITY STATES. N RUNS OVER POWERS.
C          REMEMBER THAT L=4 IS MISSING SIN(THETA)/SQRT(2)
      SUM=0.
      DO 111 I=1,4
      DO 111 J=I,4
      DENOM=1./((ADWWWW(1,I)+ADWWWW(2,I)*ZCM)
     $*(ADWWWW(1,J)+ADWWWW(2,J)*ZCM))
        DO 112 L=1,4
        TERM=0.
          DO 113 N=0,6
          C=0.
          N1=MAX(N-3,0)
          N2=MIN(3,N)
          DO 114 K=N1,N2
114       C=C+ANWWWW(K+1,I,L)*ANWWWW(N-K+1,J,L)
          C=C*WTHELI(L)
          IF(J.NE.I) C=2.*C
          TERM=TERM+C*ZCM**N
113       CONTINUE
        IF(L.EQ.4) TERM=TERM*(1.-ZCM**2)/2.
        TERM=TERM*DENOM
        SUM=SUM+TERM
112     CONTINUE
111   CONTINUE
C          ADD IMAGINARY PART SQUARED.
      SUM=SUM+WTHELI(1)*AIWWWW(1)**2+WTHELI(2)*AIWWWW(2)**2
     $+WTHELI(3)*AIWWWW(3)**2+WTHELI(4)*AIWWWW(4)**2
C          CROSS SECTION. NOTE D(OMEGA)=2.*PI*D(Z)
      SIG0=SUM/(64.*PI**2*S*SCM)*UNITS
      SIG0=SIG0*TBRWW(JETTYP(1)-25,1)*TBRWW(JETTYP(2)-25,2)
C          SYMMETRY FACTOR
      IF(IABSF.EQ.90) SIG0=.5*SIG0
      SIGLLQ=SIG0*QSAVE(INITYP(1),1)*QSAVE(INITYP(2),2)
      RETURN
      END
