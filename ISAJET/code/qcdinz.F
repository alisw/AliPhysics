#include "isajet/pilot.h"
      SUBROUTINE QCDINZ(J0)
C
C          AUXILIARY ROUTINE FOR QCDINI. GENERATE A Z AND TWO DAUGHTER
C          PARTONS FOR SPACELIKE PARTON J0.
C
#include "isajet/itapes.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/jwork2.inc"
#include "isajet/qcdpar.inc"
#include "isajet/primar.inc"
C
      DATA CA/3./,CF/1.333333333/
C          FUNCTIONS.
      PQQ(Z)=CF*(1.+Z**2)/(1.-Z)
      PQG(Z)=CF*(1.+(1.-Z)**2)/Z
      PGQ(Z)=.5*(Z**2+(1.-Z)**2)
      PGG(Z)=2.*CA*(1.-Z*(1.-Z))**2/(Z*(1.-Z))
C
C          INITIALIZE
      IDABS=IABS(JTYPE(J0))
      AM0=ABS(PJSET(5,J0))
      T0=AM0**2
      JIN0=JIN(J0)
      X0=(PJSET(4,J0)+SGN*PJSET(3,J0))/ECM
      ZGOOD=.FALSE.
      IF(ZMIN.GE.ZMAX) RETURN
C
C          SELECT BRANCHING AND GENERATE Z ACCORDING TO ALTARELLI-PARISI
C          FUNCTIONS.  THEN CHECK WITH STRUCTURE FUNCTIONS
C
C          GLUON
C
      IF(IDABS.EQ.9) THEN
C
C          GL->GL+GL
        IF(JIN0.EQ.1) THEN
110       ZGEN=DZMAX/ZMAX*(ZMAX*(1.-ZMIN)/(ZMIN*DZMAX))**RANF()
          Z=1./(1.+ZGEN)
          DZ=ZGEN/(1.+ZGEN)
          GZ=2.*CA/(Z*DZ)
          PGGZ=2.*CA*(1.-Z*(1.-Z))**2/(Z*DZ)
          IF(PGGZ.LT.GZ*RANF()) GO TO 110
          JTYPE(NJSET+1)=9
          JTYPE(NJSET+2)=9
          ZZC(J0)=Z
C
          X1=X0/Z
          FX1=STRUC(X1,T0,1,IDIN(JET-10))
          FX0=FXTEST(J0)
          IF(FX1/FX0.GT.RANF()) ZGOOD=.TRUE.
C
C          QK->GL+QK
        ELSE
120       RZMAX=SQRT(ZMAX)
          RZMIN=SQRT(ZMIN)
          ZGEN=1./RZMAX-RANF()*(1./RZMAX-1./RZMIN)
          Z=1./ZGEN**2
          RZ=SQRT(Z)
          GZ=2.*CF/RZ**3
          IF(PQG(Z)/RZ.LT.GZ*RANF()) GO TO 120
          IFL=JIN0/2
          IF(JIN0.NE.2*IFL) IFL=-IFL
          JTYPE(NJSET+1)=IFL
          JTYPE(NJSET+2)=IFL
          ZZC(J0)=Z
C
          X1=X0/Z
          FX1=STRUC(X1,T0,JIN0,IDIN(JET-10))
          FX0=FXTEST(J0)
          IF(RZ*FX1/FX0.GT.RANF()) ZGOOD=.TRUE.
        ENDIF
C
C          QUARK
C
      ELSE
C
C          GL->QK+QB
        IF(JIN0.EQ.1) THEN
130       Z=ZMIN+(ZMAX-ZMIN)*RANF()
          IF(PGQ(Z).LT..5*RANF()) GO TO 130
          JTYPE(NJSET+1)=9
          JTYPE(NJSET+2)=-JTYPE(J0)
          ZZC(J0)=Z
C
          X1=X0/Z
          FX1=STRUC(X1,T0,1,IDIN(JET-10))
          FX0=FXTEST(J0)
          IF(FX1/FX0.GT.RANF().OR.GLFORC(JET-10)) ZGOOD=.TRUE.
C
C          QK->QK+GL
        ELSE
140       DZ=DZMAX*((1.-ZMIN)/DZMAX)**RANF()
          Z=1.-DZ
          GZ=2.*CF/DZ
          RZ=1.
          IF(IDABS.LE.3) RZ=SQRT(Z)
          PQQZ=CF*(1.+Z**2)/DZ
          IF(PQQZ/RZ.LT.GZ*RANF()) GO TO 140
          JTYPE(NJSET+1)=JTYPE(J0)
          JTYPE(NJSET+2)=9
          ZZC(J0)=Z
C
          X1=X0/Z
          FX1=STRUC(X1,T0,JIN0,IDIN(JET-10))
          FX0=FXTEST(J0)
          IF(RZ*FX1/FX0.GT.RANF()) ZGOOD=.TRUE.
        ENDIF
      ENDIF
      JMATCH(NJSET+1)=0
      JMATCH(NJSET+2)=0
      RETURN
      END
