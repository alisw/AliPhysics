#include "isajet/pilot.h"
      SUBROUTINE QCDZ(J)
C
C          Auxiliary routine for QCDJET.  Generate Z for parton J and 
C          store in ZZC(J). Add possible new partons to /JETSET/.
C
C          Include GM, W+, W-, and Z0 radiation.
C
C          Ver 7.20: Anomalous dimensions were coded incorrectly!
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/qcdpar.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/q1q2.inc"
C
      REAL PQQ,PGQ,PQG,PGG,Z,PGSGS,PQSQS,ALFAS,QQ,AM0,ZC,AMASS
      REAL GAMQQ,GAMGG,PROBG,PROBQ,RND,RANF,ZGEN,GZ
      REAL GAMZC,GAMSUM,AM1W,AM2W,T1W,T2W,ZCW,T0,GAMZCW,TERM,SUM
      REAL SUMBR,BRMODE,TRY,HELPL,HELMN,HEL,PZ
      INTEGER NF,J,JTABS,IQ,IFL,IW,JT0,JT1,IFL1,IFL2
      INTEGER IWTYPE,JET,JW,IQ1,IQ2,JPAR,IFLPAR,NJ1,NJ2,IDABS1,IDABS2
      REAL GAMSAV(5),ZCSAV(5),BRANCH(25)
      INTEGER JSAV(5),LISTW(5),LISTJ(25)
      DATA LISTW/9,10,80,-80,90/
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
C
C          Altarelli-Parisi functions.
      PQQ(Z)=4*(1+Z**2)/(3*(1-Z))
      PGQ(Z)=.5*(Z**2+(1-Z)**2)
      PGG(Z)=6*(1-Z*(1-Z))**2/(Z*(1-Z))
      PGSGS(Z)=3.*(1.+Z**2)/(1.-Z)
      PQSQS(Z)=8./3.*Z/(1.-Z)
      ALFAS(QQ)=12.*PI/((33.-2.*NF)*ALOG(QQ/ALAM2))
C
C          Initialize.
C
      AM0=PJSET(5,J)
      ZC=ZZC(J)
      JTABS=IABS(JTYPE(J))
      NF=3
      DO 110 IQ=4,6
        IF(AM0.LT.2*AMASS(IQ)) GO TO 120
        NF=NF+1
110   CONTINUE
120   CONTINUE
      NJ1=NJSET+1
      NJ2=NJSET+2        
C
C          Initial gluon
C
      IF (JTABS.EQ.9) THEN
        GAMQQ=(1-2*ZC)*(1-ZC*(1-ZC))/3.
        GAMGG=12*ALOG((1-ZC)/ZC)-9*(1-2*ZC)-6*GAMQQ
        PROBG=GAMGG/(GAMGG+2*NF*GAMQQ)
        PROBQ=GAMQQ/(GAMGG+2*NF*GAMQQ)
        RND=RANF()
C          GL--->GL+GL
        IF(PROBG.GT.RND) THEN
130       ZGEN=(ZC/(1-ZC))**(1-2*RANF())
          Z=ZGEN/(1.+ZGEN)
          GZ=6./(Z*(1.-Z))
          IF(PGG(Z).LT.GZ*RANF()) GO TO 130
          JTYPE(NJ1)=9
          JTYPE(NJ2)=9
          ZZC(J)=Z
C          GL--->QK+QB
        ELSE 
140       Z=RANF()
          IF(PGQ(Z).LT.0.5*RANF()) GO TO 140
          IFL=(RND-PROBG)/PROBQ+1.
          IF(IFL.GT.NF) IFL=NF-IFL
          JTYPE(NJ1)=IFL
          JTYPE(NJ2)=-IFL
          ZZC(J)=Z
        ENDIF
C
C          Initial quark - may radiate GL, GM, W+-, Z0
C
      ELSEIF(JTABS.LT.9) THEN
C          Gluon
        GAMZC=2.*ALOG((1-ZC)/ZC)-1.5*(1.-2.*ZC)
        GAMSAV(1)=4./3.*ALFAS(AM0**2)*GAMZC
        ZCSAV(1)=ZC
        JSAV(1)=JTYPE(J)
C          Photon
        GAMSAV(2)=ALFA*AQ(JTABS,1)**2*GAMZC
        ZCSAV(2)=ZC
        GAMSUM=GAMSAV(1)+GAMSAV(2)
        JSAV(2)=JTYPE(J)
C          W+- and Z0 
        IF(AM0.GT.WMASS(4)) THEN
          DO 200 IW=2,4
            GAMSAV(IW+1)=0.
            ZCSAV(IW+1)=.5
            JSAV(IW+1)=0
            JT0=2*IABS(JTYPE(J))
            IF(JTYPE(J).LT.0) JT0=JT0+1
            JT1=MATCH(JT0,IW)
            IF(JT1.EQ.0) GO TO 200
            JT1=MATCH(JT1,4)
            IFL1=JT1/2
            AM1W=AMASS(IFL1)
            AM2W=AMASS(LISTW(IW+1))
            IF(AM1W+AM2W.GE.AM0) GO TO 200
            T0=AM0**2
            T1W=AM1W**2
            T2W=AM2W**2
            ZCW=(T0-T1W+T2W-SQRT((T0-T1W-T2W)**2-4*T1W*T2W))/(2*T0)
            GAMZCW=2.*ALOG((1-ZCW)/ZCW)-2.*(1.-2.*ZCW)
            TERM=(AQ(JTABS,IW)**2+BQ(JTABS,IW)**2)*ALFA*GAMZCW
            GAMSAV(IW+1)=TERM
            ZCSAV(IW+1)=ZCW
            JSAV(IW+1)=IFL1*ISIGN(1,JTYPE(J))
            GAMSUM=GAMSUM+TERM
200       CONTINUE
        ELSE
          DO 210 IW=2,4
            GAMSAV(IW+1)=0.
            ZCSAV(IW+1)=.5
            JSAV(IW+1)=0
210       CONTINUE
        ENDIF
C          Select decay mode
        RND=RANF()
        SUM=0.
        DO 220 IW=1,5
          IWTYPE=IW
          SUM=SUM+GAMSAV(IW)/GAMSUM
          IF(RND.LT.SUM) GO TO 230
220     CONTINUE
C          Generate Z
230     CONTINUE
        Z=1-(ZC/(1-ZC))**RANF()*(1-ZC)
        GZ=8./(3.*(1-Z))
        IF(PQQ(Z).LT.GZ*RANF()) GO TO 230
        IF(Z.LT.ZCSAV(IWTYPE).OR.Z.GT.1.-ZCSAV(IWTYPE)) GO TO 230
        JTYPE(NJ1)=JSAV(IWTYPE)
        JTYPE(NJ2)=LISTW(IWTYPE)
        ZZC(J)=Z
C
C          Initial diquark
C
      ELSEIF(MOD(JTABS,100).EQ.0) THEN
300     CONTINUE
        Z=1-(ZC/(1-ZC))**RANF()*(1-ZC)
        GZ=8./(3.*(1-Z))
        IF(PQQ(Z).LT.GZ*RANF()) GO TO 300
        JTYPE(NJ1)=JTYPE(J)
        JTYPE(NJ2)=9
        ZZC(J)=Z
C
C          Initial gluino
C
       ELSEIF (JTABS.EQ.29) THEN
400      Z=1.-(ZC/(1.-ZC))**RANF()*(1.-ZC)
         GZ=6./(1.-Z)
         IF(PGSGS(Z) .LT. GZ*RANF()) GOTO 400
         JTYPE(NJ1)=JTYPE(J)
         JTYPE(NJ2)=9
         ZZC(J)=Z
C
C          Initial squark
C
      ELSEIF(JTABS.GT.20.AND.JTABS.LT.29) THEN
500     CONTINUE
        Z=1-(ZC/(1-ZC))**RANF()*(1-ZC)
        GZ=8./(3.*(1-Z))
        IF(PQSQS(Z).LT.GZ*RANF()) GO TO 500
        JTYPE(NJ1)=JTYPE(J)
        JTYPE(NJ2)=9
        ZZC(J)=Z
C
C          Initial W+, W-, or Z0
C
      ELSEIF(JTABS.EQ.80.OR.JTABS.EQ.90) THEN
C          Select decay mode
        IF(JTYPE(J).EQ.+80) JW=2
        IF(JTYPE(J).EQ.-80) JW=3
        IF(JTYPE(J).EQ.+90) JW=4
        TRY=RANF()
        DO 610 IQ=2,25
        IF(TRY.LT.CUMWBR(IQ,JW-1)) THEN
          IQ1=IQ
          IQ2=MATCH(IQ,JW)
          GO TO 620
        ENDIF
610     CONTINUE
620     JTYPE(NJ1)=LISTJ(IQ1)
        JTYPE(NJ2)=LISTJ(IQ2)
C          Select W helicity
        JPAR=MOD(JORIG(J),JPACK)
        IFLPAR=IABS(JTYPE(JPAR))
        HELPL=(AQ(IFLPAR,JW)-BQ(IFLPAR,JW))**2
        HELMN=(AQ(IFLPAR,JW)+BQ(IFLPAR,JW))**2
        IF(RANF()*(HELPL+HELMN).LT.HELMN) THEN
          HEL=-ISIGN(1,JTYPE(NJ1))
        ELSE
          HEL=+ISIGN(1,JTYPE(NJ1))
        ENDIF
630     Z=RANF()
        PZ=(1.+HEL*(2.*Z-1.))**2
        IF(PZ.LT.4.*RANF()) GO TO 630
        ZZC(J)=Z
      ENDIF
C
C          Set masses and flags.
C
      JET=IABS(JORIG(J))/JPACK
      JORIG(NJ1)=JPACK*JET+J
      JORIG(NJ2)=JPACK*JET+J
      IDABS1=IABS(JTYPE(NJ1))
      IDABS2=IABS(JTYPE(NJ2))
      JMATCH(NJ1)=NJ2
      JMATCH(NJ2)=NJ1
C          JDCAY=-1 implies further decay
      IF(IDABS1.LE.9.OR.(IDABS1.GT.20.AND.IDABS1.LT.30.).OR.
     $MOD(IDABS1,100).EQ.0) THEN
        PJSET(5,NJ1)=Z*AM0
        JDCAY(NJ1)=-1
      ELSEIF(IDABS1.GE.80.OR.IDABS1.LE.90) THEN
        PJSET(5,NJ1)=AMASS(IDABS1)
        JDCAY(NJ1)=-1
      ELSE
        PJSET(5,NJ1)=AMASS(IDABS1)
        JDCAY(NJ1)=0
      ENDIF
      IF(IDABS2.LE.9.OR.(IDABS2.GT.20.AND.IDABS2.LT.30.).OR.
     $MOD(IDABS2,100).EQ.0) THEN
        PJSET(5,NJ2)=(1.-Z)*AM0
        JDCAY(NJ2)=-1
      ELSEIF(IDABS2.EQ.80.OR.IDABS2.EQ.90) THEN
        PJSET(5,NJ2)=AMASS(IDABS2)
        JDCAY(NJ2)=-1
      ELSE
        PJSET(5,NJ2)=AMASS(IDABS2)
        JDCAY(NJ2)=0
      ENDIF
      RETURN
      END
