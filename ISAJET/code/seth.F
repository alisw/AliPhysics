#include "isajet/pilot.h"
      SUBROUTINE SETH
C
C          Set the standard Weinberg-Salam Higgs parameters in /HCON/.
C          HMASS  = Higgs mass
C          HGAM   = Higgs width
C          HGAMS  = Higgs partial width
C          ZSTARS = minimum allowed mass for Z*
C
C          IQ = 1  2  3  4  5  6  7  8  9  10 11 12 13
C               GL UP UB DN DB ST SB CH CB BT BB TP TB
C          IQ = 14  15   16 17 18   19   20  21  22  23   24   25
C               NUE ANUE E- E+ NUMU ANUM MU- MU+ NUT ANUT TAU- TAU+
C          IQ = 26 27 28 29
C               GM W+ W- Z0
C
C          Ver 6.25: Added H -> GM GM.
C          Ver 6.26: Added H -> Z0 Z* from Keung and Marciano, Phys. 
C                    Rev. D30, 248 (1984).
C          Ver 6.30: Fixed sign of FFR in H -> GM GM for TAU<1. Added
C                    H -> W W* to total width but not to partial widths
C                    to get right branching ratios.
C          Ver 7.38: Add H_SM decay modes to SSSAVE for use in WHIGGS
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/wcon.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/nodcay.inc"
#include "isajet/const.inc"
#include "isajet/hcon.inc"
C
      REAL GAMFCN,X,AMASS,AMQ,GAMQ,AML,WM,GAMWW,TAU,FFR,FFI,FR,FI,
     $ROOT,ROOTLN,TM,SUMBR,TERM,ETAR,ETAI,RQ,RQLOG,PHIR,PHII
      REAL EPS,FEPS,AM12
      INTEGER IQ,IQ1,IQ2,I,IW
      INTEGER LISTJ(25),LISTW(4)
      DATA LISTJ/
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
      DATA LISTW/10,80,-80,90/
C
      GAMFCN(X)=SQRT(1.-4*X**2)*(1.-4.*X**2+12.*X**4)
C
C          Calculate Higgs mass and width
C
      HMASS=AMASS(81)
      HGAM=0.
      DO 100 IQ=1,29
100   HGAMS(IQ)=0.
C
C          Quarks and leptons
      DO 110 IQ=1,6
        AMQ=AMASS(IQ)
        IF(AMQ.GT.0..AND.AMQ.LT..5*HMASS) THEN
          GAMQ=3.*GF*AMQ**2*HMASS/(4.*PI*SQRT2)
     $    *(SQRT(1.-4.*AMQ**2/HMASS**2))**3
          HGAM=HGAM+GAMQ
          HGAMS(2*IQ)=.5*GAMQ
          HGAMS(2*IQ+1)=.5*GAMQ
          CALL SSSAVE(81,GAMQ,IQ,-IQ,0,0,0)
        ENDIF
        AML=AMASS(IQ+10)
        IF(AML.GT.0..AND.AML.LT..5*HMASS) THEN
          GAMQ=GF*AML**2*HMASS/(4.*PI*SQRT2)
     $    *(SQRT(1.-4.*AML**2/HMASS**2))**3
          HGAM=HGAM+GAMQ
          HGAMS(2*IQ+12)=.5*GAMQ
          HGAMS(2*IQ+13)=.5*GAMQ
          CALL SSSAVE(81,GAMQ,IQ+10,-(IQ+10),0,0,0)
        ENDIF
110   CONTINUE
C
C          W+ W- and Z0 Z0, including W W* and Z Z*.
      WM=WMASS(2)
      IF(HMASS.GT.2.*WM) THEN
        GAMWW=GF*HMASS**3*GAMFCN(WM/HMASS)/(8.*PI*SQRT2)
        HGAM=HGAM+GAMWW
        HGAMS(27)=.5*GAMWW
        HGAMS(28)=.5*GAMWW
        CALL SSSAVE(81,GAMWW,80,-80,0,0,0)
      ELSEIF(HMASS.GT.WM) THEN
        EPS=WM/HMASS
        FEPS=3.*(1.-8.*EPS**2+20.*EPS**4)/SQRT(4.*EPS**2-1.)
     $  *ACOS((3.*EPS**2-1.)/(2.*EPS**3))
     $  -(1.-EPS**2)*(47./2.*EPS**2-13./2.+1./EPS**2)
     $  -3.*(1.-6.*EPS**2+4.*EPS**4)*ALOG(EPS)
        GAMWW=3.*ALFA**2*HMASS/(32.*PI*SIN2W**2)*FEPS
        HGAM=HGAM+GAMWW
        HGAMS(27)=.5*GAMWW
        HGAMS(28)=.5*GAMWW
        CALL SSSAVE(81,GAMWW/18.,80,12,-11,0,0)
        CALL SSSAVE(81,GAMWW/18.,-80,-12,11,0,0)
        CALL SSSAVE(81,GAMWW/18.,80,14,-13,0,0)
        CALL SSSAVE(81,GAMWW/18.,-80,-14,13,0,0)
        CALL SSSAVE(81,GAMWW/18.,80,16,-15,0,0)
        CALL SSSAVE(81,GAMWW/18.,-80,-16,15,0,0)
        CALL SSSAVE(81,GAMWW/6.,80,-1,2,0,0)
        CALL SSSAVE(81,GAMWW/6.,-80,1,-2,0,0)
        CALL SSSAVE(81,GAMWW/6.,80,-4,3,0,0)
        CALL SSSAVE(81,GAMWW/6.,-80,4,-3,0,0)
      ENDIF
      WM=WMASS(4)
      IF(HMASS.GT.2.*WM) THEN
        GAMWW=GF*HMASS**3*GAMFCN(WM/HMASS)/(16.*PI*SQRT2)
        HGAM=HGAM+GAMWW
        HGAMS(29)=GAMWW
        CALL SSSAVE(81,GAMWW,90,90,0,0,0)
      ELSEIF(HMASS.GT.WM) THEN
        EPS=WM/HMASS
        FEPS=3.*(1.-8.*EPS**2+20.*EPS**4)/SQRT(4.*EPS**2-1.)
     $  *ACOS((3.*EPS**2-1.)/(2.*EPS**3))
     $  -(1.-EPS**2)*(47./2.*EPS**2-13./2.+1./EPS**2)
     $  -3.*(1.-6.*EPS**2+4.*EPS**4)*ALOG(EPS)
        GAMWW=ALFA**2*HMASS/(128.*PI*SIN2W**2*(1.-SIN2W)**2)
     $  *(7.-40./3.*SIN2W+160./9.*SIN2W**2)*FEPS
        HGAM=HGAM+GAMWW
        HGAMS(29)=GAMWW
        CALL SSSAVE(81,.11922*GAMWW,90,-1,1,0,0)
        CALL SSSAVE(81,.15375*GAMWW,90,-2,2,0,0)
        CALL SSSAVE(81,.15375*GAMWW,90,-3,3,0,0)
        CALL SSSAVE(81,.11922*GAMWW,90,-4,4,0,0)
        CALL SSSAVE(81,.15375*GAMWW,90,-5,5,0,0)
        CALL SSSAVE(81,.06668*GAMWW,90,-11,11,0,0)
        CALL SSSAVE(81,.03343*GAMWW,90,-12,12,0,0)
        CALL SSSAVE(81,.06668*GAMWW,90,-13,13,0,0)
        CALL SSSAVE(81,.03343*GAMWW,90,-14,14,0,0)
        CALL SSSAVE(81,.06668*GAMWW,90,-15,15,0,0)
        CALL SSSAVE(81,.03343*GAMWW,90,-16,16,0,0)
      ENDIF
C          W* and Z* mass limits
      DO 120 I=1,2
        ZSTARS(1,I)=0.
        DO 130 IW=2,4
          ZSTARS(IW,I)=AMASS(LISTW(IW))
          DO 140 IQ1=2,25
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 140
            IF(GOWW(IQ1,1).AND.GOWW(IQ2,2)) THEN
              AM12=AMASS(LISTJ(IQ1))+AMASS(LISTJ(IQ2))
              ZSTARS(IW,I)=MIN(ZSTARS(IW,I),AM12)
            ENDIF
140       CONTINUE
130     CONTINUE
120   CONTINUE
C
C          GM GM -- W loop term
      WM=WMASS(2)
      TAU=4.*WM**2/HMASS**2
      IF(TAU.GE.1.0) THEN
        FFR=(ASIN(1./SQRT(TAU)))**2
        FFI=0.
      ELSE
        ROOT=SQRT(1.-TAU)
        ROOTLN=ALOG((1.+ROOT)/(1.-ROOT))
        FFR=-0.25*(ROOTLN**2-PI**2)
        FFI=0.5*PI*ROOTLN
      ENDIF
      FR=2.+3.*TAU+3.*TAU*(2.-TAU)*FFR
      FI=3.*TAU*(2.-TAU)*FFI
C          Top loop term
      TM=AMASS(6)
      TAU=4.*TM**2/HMASS**2
      IF(TAU.GE.1.0) THEN
        FFR=(ASIN(1./SQRT(TAU)))**2
        FFI=0.
      ELSE
        ROOT=SQRT(1.-TAU)
        ROOTLN=ALOG((1.+ROOT)/(1.-ROOT))
        FFR=-0.25*(ROOTLN**2-PI**2)
        FFI=0.5*PI*ROOTLN
      ENDIF
      FR=FR-8./3.*TAU*(1.+(1.-TAU)*FFR)
      FI=FI-8./3.*TAU*(1.-TAU)*FFI
C          Total GM GM
      HGAMS(26)=ALFA**3/(256.*PI**2*SIN2W)*HMASS**3/WM**2*(FR**2+FI**2)
      HGAM=HGAM+HGAMS(26)
      CALL SSSAVE(81,HGAMS(26),10,10,0,0,0)
C
C          Calculate Higgs-gluon-gluon coupling
C
      ETAR=0.
      ETAI=0.
      DO 300 IQ=1,8
        AMQ=AMASS(IQ)
        IF(AMQ.LE.0.) GO TO 300
        RQ=(2.*AMQ/HMASS)**2
        IF(RQ.GE.1.) THEN
          ETAR=ETAR+.5*RQ*(1.+(1.-RQ)*ASIN(1./SQRT(RQ))**2)
        ELSE
          RQLOG=ALOG((1.+SQRT(1.-RQ))/(1.-SQRT(1.-RQ)))
          PHIR=.25*(RQLOG**2-PI**2)
          ETAR=ETAR+.5*RQ*(1.+(RQ-1.)*PHIR)
          PHII=.5*PI*RQLOG
          ETAI=ETAI+.5*RQ*(1.+(RQ-1.)*PHII)
        ENDIF
300   CONTINUE
      ETAHGG=ETAR**2+ETAI**2
C
      RETURN
      END
