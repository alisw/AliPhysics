#include "isajet/pilot.h"
      SUBROUTINE QCDINT(J0)
C
C          AUXILIARY ROUTINE FOR QCDINI.  GENERATE A NEW MASS FOR
C          SPACELIKE PARTON J0.
C
#include "isajet/itapes.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/jwork2.inc"
#include "isajet/qcdpar.inc"
#include "isajet/primar.inc"
C
      DIMENSION GAMS(13),FX0S(13)
      DATA CA/3./,CF/1.333333333/
C
C          FUNCTIONS -- USE DZMAX FOR PRECISION
      GQQ(Z,DZ)=CF*(-2.*ALOG(DZ)+Z*(-1.-.5*Z))
      GQG(Z)=CF*(+2.*ALOG(Z)+Z*(-2.+.5*Z))
      GGQ(Z)=(Z**3-(1.-Z)**3)/6.
      GGG(Z,DZ)=2.*CA*(ALOG(Z/DZ)+Z*(-2.+Z*(.5-Z/3.)))
      GBQQ(RZ,DZ)=CF*(2.*ALOG((1.+RZ)**2/DZ)+RZ*(-2.-2./3.*RZ**2))
      GBQG(RZ)=CF*(-4./RZ+RZ*(-4.+2./3.*RZ**2))
C
      GLFORC(JET-10)=.FALSE.
      IDABS=IABS(JTYPE(J0))
      IF(JTYPE(J0).EQ.9) THEN
        ITYP=1
      ELSEIF(JTYPE(J0).GT.0) THEN
        ITYP=2*IDABS
      ELSE
        ITYP=2.*IDABS+1
      ENDIF
      IBEAM=JET-10
      AM0=ABS(PJSET(5,J0))
1     T0=AM0**2
      X0=ZMIN
      ANF=3
      DO 110 I=4,6
      AMQ2=AMASS(I)**2
110   ANF=ANF+T0/(AMQ2+T0)
      B0=11.-2.*ANF/3.
C
C          SET UP ANOMALOUS DIMENSIONS. ALSO USE THESE TO DETERMINE TYPE
C          OF INCOMING PARTON (TO BE USED IN QCDINZ).
C
C          GLUON
      IF(IDABS.EQ.9) THEN
        AMQ=0.
        GAMG=GGG(ZMAX,DZMAX)-GGG(ZMIN,1.-ZMIN)
        GAMS(1)=GAMG
        FX0=STRUC(X0,T0,1,IDIN(IBEAM))
        FX0S(1)=FX0
        GAMFAC=(GBQG(SQRT(ZMAX))-GBQG(SQRT(ZMIN)))/FX0
        GAMQ=0.
        DO 210 IQ=2,13
        FX0S(IQ)=STRUC(X0,T0,IQ,IDIN(IBEAM))
        GAMS(IQ)=GAMFAC*FX0S(IQ)
210     GAMQ=GAMQ+GAMS(IQ)
        GAM=GAMG+GAMQ
        AM1=CUTJET
C
        TRY=RANF()
        SUM=0.
        DO 220 IQ=1,13
        SUM=SUM+GAMS(IQ)/GAM
        IF(SUM.LT.TRY) GO TO 220
        JIN(J0)=IQ
        FXTEST(J0)=FX0S(IQ)
        GO TO 300
220     CONTINUE
C
C          LIGHT QUARK
      ELSEIF(IDABS.LE.3) THEN
        AMQ=AMASS(IDABS)
        GAMQ=GBQQ(SQRT(ZMAX),DZMAX)-GBQQ(SQRT(ZMIN),1.-ZMIN)
        FX0=STRUC(X0,T0,ITYP,IDIN(IBEAM))
        FXG=STRUC(X0,T0,1,IDIN(IBEAM))
        GAMFAC=FXG/FX0
        GAMG=GAMFAC*(GGQ(ZMAX)-GGQ(ZMIN))
        GAM=GAMQ+GAMG
        AM1=AMQ+CUTJET
C
        IF(GAMQ/GAM.GT.RANF()) THEN
          JIN(J0)=ITYP
          FXTEST(J0)=FX0
        ELSE
          JIN(J0)=1
          FXTEST(J0)=FXG
        ENDIF
C
C          HEAVY QUARK -- SPECIAL TREATMENT NEEDED TO ALWAYS FORCE
C          GL-->QK+QB BEFORE END OF EVOLUTION.
C          USE SMALLER MASS FOR FORCED DECAYS TO PREVENT INFINITE LOOP.
      ELSE
        AMQ=AMASS(IDABS)
        THRESH=4.*AMQ**2*X0/(1.-X0)
        THRESH=(SQRT(THRESH)+CUTJET)**2
        IF(STRUC(X0,T0,ITYP,IDIN(IBEAM)).LE.0..OR.
     $  T0.LE.THRESH) THEN
          PJSET(5,J0)=-AM0*SQRT(RANF())-ALAM
          GLFORC(JET-10)=.TRUE.
          JDCAY(J0)=-2
          JIN(J0)=1
          FXTEST(J0)=1.
          RETURN
        ENDIF
        T1=SQRT(T0*THRESH)
230     AM1=SQRT(T1)
        FX0=STRUC(X0,T1,ITYP,IDIN(IBEAM))
        IF(FX0.LE.0.) THEN
          T1=SQRT(T1*T0)
          GO TO 230
        ENDIF
        FXG=STRUC(X0,T1,1,IDIN(IBEAM))
        GAMFAC=FXG/FX0
        GAMQ=GQQ(ZMAX,DZMAX)-GQQ(ZMIN,1.-ZMIN)
        GAMG=GAMFAC*(GGQ(ZMAX)-GGQ(ZMIN))
        GAM=GAMQ+GAMG
C
        IF(GAMQ/GAM.GT.RANF()) THEN
          JIN(J0)=ITYP
          FXTEST(J0)=FX0
        ELSE
          JIN(J0)=1
          FXTEST(J0)=FXG
        ENDIF
      ENDIF
C
C          LEADING-LOG MASS GENERATION.
C
300   GB=2.*GAM/B0
      IF(AM1.GT.ALAM.AND.AM0.GT.ALAM) THEN
        PROBL=GB*ALOG(ALOG(AM1/ALAM)/ALOG(AM0/ALAM))
      ELSE
        PROBL=0.
      ENDIF
      IF(PROBL.GT.0.) THEN
        PROB=1.
      ELSEIF(PROBL.GT.-50.) THEN
        PROB=EXP(PROBL)
      ELSE
        PROB=0.
      ENDIF
      IF(PROB.GT.RANF()) THEN
        IF(IDABS.LE.3.OR.IDABS.EQ.9) THEN
          PJSET(5,J0)=AMQ
          JDCAY(J0)=JPACK*J0+J0
          RETURN
        ELSEIF(AM0.LT.AM1+CUTJET) THEN
          PJSET(5,J0)=-SQRT(T0)
          GLFORC(JET-10)=.TRUE.
          JDCAY(J0)=-2
          JIN(J0)=1
          FXTEST(J0)=1
          RETURN
        ELSE
          AM0=AM1
          GO TO 1
        ENDIF
      ELSE
        POW=(1.-(1.-PROB)*RANF())**(1./GB)
        AMNEW=ALAM*(AM0/ALAM)**POW
        IF(AMNEW.GE.AM1) THEN
          PJSET(5,J0)=-AMNEW
          JDCAY(J0)=-2
          RETURN
        ELSEIF(IDABS.LE.3.OR.IDABS.EQ.9) THEN
          PJSET(5,J0)=AMQ
          JDCAY(J0)=JPACK*J0+J0
          RETURN
        ELSEIF(AM0.LT.AM1+CUTJET) THEN
          PJSET(5,J0)=-AM0*SQRT(RANF())-ALAM
          GLFORC(JET-10)=.TRUE.
          JDCAY(J0)=-2
          JIN(J0)=1
          FXTEST(J0)=1
          RETURN
        ELSE
          AM0=AM1
          GO TO 1
        ENDIF
      ENDIF
      END
