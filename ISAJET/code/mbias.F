#include "isajet/pilot.h"
      SUBROUTINE MBIAS
C
C          Generate minbias event or beam jets for high-pt event using
C          parameters set in MBSET:
C
C          (1) Select number NPOM of cut pomerons -- cf cut Reggeon
C              field theory of Abramovskii, Kanchelli, and Gribov.
C          (2) Generate xf for leading baryons including 1/(1-xf)
C              diffractive term and guessed NPOM dependence,
C              F(XF)=(1-XF)**(A+B/NPOM)
C          (3) Select xf for each half of each Pomeron. then fragment
C              each half Pomeron into mesons and baryons independently
C              in the Pomeron-Pomeron center of mass. This avoids
C              making xf=0 a singular point.
C
C          Note that multiple cut Pomerons give approximate KNO scaling.
C          The only short-range correlations are from resonances.
C
C          Ver. 7.09: Add traps on free loops and IMPLICIT NONE.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/mbgen.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/partcl.inc"
#include "isajet/mbpar.inc"
C
      DIMENSION IFL(3),IFLEXC(2),PXEXC(2),PYEXC(2),SIGN(2)
      DIMENSION PSUM(5)
      DIMENSION LDIFFR(2)
      LOGICAL LDIFFR
      REAL RANF,AMASS
      REAL RND,XX,XSUM,P0,PPOM,PXEXC,PSUM,SIGN,PYEXC,GAM,BETA,X,
     $AM,PPLUS,EPSDIF,PEND0,TRY,PX,PY,PX2,PY2,QMINUS,PZ,QPLUS,PT1,
     $PHI1,XBGEN,PT2,PHI2,PX1,PY1,AMT2
      INTEGER ID1,ID2,IFL1,IFL2,IMOD1,IMOD2,ITWIST,IPOM,LOOP,NFIRST,
     $ID3,IFLEXC,IFL,I,NP2,IP,NP1,IFAIL,NPTLV1,IDHAD,IB,NEW,JSPIN,
     $INDEX,NBEGIN,IPASS,MXPASS,N,IDIFF,IPASSB,IFLNEW,ISWAP
      DATA SIGN/1.,-1./,PEND0/.14/
      DATA PSUM/5*0./
      DATA MXPASS/200/
C
C          Start
      NBEGIN=NPTCL+1
      IPASS=1
      IPASSB=1
C
C          Select number of cut Pomerons.
C
1     CONTINUE
      TRY=RANF()
      DO 10 N=MNPOM,MXPOM
        NPOM=N
        IF(POMGEN(N).GT.TRY) GO TO 20
10    CONTINUE
20    CONTINUE
C
C          Decide if diffractive event
      IF(RANF().LT.PDIFFR) THEN
        IDIFF=INT(1.99999*RANF())+1
        LDIFFR(IDIFF)=.TRUE.
        LDIFFR(3-IDIFF)=.FALSE.
      ELSE
        LDIFFR(1)=.FALSE.
        LDIFFR(2)=.FALSE.
      ENDIF
C
C          Generate leading baryons.
C
      DO 100 IB=1,2
        PPLUS=2.*PBEAM(IB)
C
C          Special treatment for diffractive beam.
        IF(LDIFFR(IB)) THEN
          IDHAD=IDIN(IB)
          AM=AMASS(IDHAD)
          CALL FLAVOR(IDIN(IB),IFL(1),IFL(2),IFL(3),JSPIN,INDEX)
          NEW=INT(3.*RANF())+1
          IFLEXC(1)=+IFL(NEW)
          IFLEXC(2)=-IFL(NEW)
          EPSDIF=2./SCM
          DXBARY(IB)=EPSDIF**RANF()
          XBARY(IB)=1.-DXBARY(IB)
          GO TO 115
        ENDIF
C
C          If not diffractive, construct IDENT of leading baryon
        CALL FLAVOR(IDIN(IB),IFL(1),IFL(2),IFL(3),JSPIN,INDEX)
        NEW=INT(3.*RANF())+1
        IFLNEW=ISIGN(INT(RANF()/PUD0)+1,IDIN(IB))
        IFLEXC(1)=IFL(NEW)
        IFLEXC(2)=-IFLNEW
        IFL(NEW)=IFLNEW
        IF(IABS(IFL(1)).GT.IABS(IFL(2))) THEN
          ISWAP=IFL(1)
          IFL(1)=IFL(2)
          IFL(2)=ISWAP
        ENDIF
        IF(IABS(IFL(2)).GT.IABS(IFL(3))) THEN
          ISWAP=IFL(2)
          IFL(2)=IFL(3)
          IFL(3)=ISWAP
        ENDIF
        IF(IABS(IFL(1)).GT.IABS(IFL(2))) THEN
          ISWAP=IFL(1)
          IFL(1)=IFL(2)
          IFL(2)=ISWAP
        ENDIF
        JSPIN=1
        IF(IFL(1).EQ.IFL(2).AND.IFL(2).EQ.IFL(3)) THEN
          JSPIN=1
        ELSE
          JSPIN=INT(RANF()+PJSPN)
        ENDIF
        IF(JSPIN.EQ.0.AND.IFL(1).NE.IFL(2).AND.IFL(2).NE.IFL(3)) THEN
          IF(RANF().GT.PISPN) THEN
            ISWAP=IFL(1)
            IFL(1)=IFL(2)
            IFL(2)=ISWAP
          ENDIF
        ENDIF
        IDHAD=1000*IABS(IFL(1))+100*IABS(IFL(2))+10*IABS(IFL(3))+JSPIN
        IDHAD=ISIGN(IDHAD,IDIN(IB))
        AM=AMASS(IDHAD)
C
C          Select xf for nondiffractive baryon, flat for NPOM=1 and
C          like mesons for NPOM=infinity.
110     XBGEN=XGEN0(2)*(1.-1./NPOM)
        DXBARY(IB)=RANF()**(1./(XBGEN+1.))
        XBARY(IB)=1.-DXBARY(IB)
C
C          Select transverse momentum of baryon
115     CALL GETPT(PT1,SIGQT0)
        PHI1=2.*PI*RANF()
        PX1=PT1*COS(PHI1)
        PY1=PT1*SIN(PHI1)
        PXEXC(1)=PX1
        PYEXC(1)=PY1
        CALL GETPT(PT2,SIGQT0)
        PHI2=2.*PI*RANF()
        PX2=PT2*COS(PHI2)
        PY2=PT2*SIN(PHI2)
        PXEXC(2)=PX2
        PYEXC(2)=PY2
        PX=-PX1-PX2
        PY=-PY1-PY2
        AMT2=PX**2+PY**2+AM**2
C
        QPLUS=XBARY(IB)*PPLUS
        QPLUS=AMAX1(QPLUS,1.E-6)
        QMINUS=AMT2/QPLUS
        PZ=.5*(QPLUS-QMINUS)
        P0=.5*(QPLUS+QMINUS)
C
C          Add baryon to /PARTCL/ if PZ>0.
        IF(NPTCL.GE.MXPTCL) GO TO 9999
        IF(PZ.GE.0.) THEN
          NPTCL=NPTCL+1
          PPTCL(1,NPTCL)=PX
          PPTCL(2,NPTCL)=PY
          PPTCL(3,NPTCL)=PZ*SIGN(IB)
          PPTCL(4,NPTCL)=P0
          PPTCL(5,NPTCL)=AM
          IORIG(NPTCL)=0
          IDCAY(NPTCL)=0
          IDENT(NPTCL)=IDHAD
        ELSE
          IPASSB=IPASSB+1
          IF(IPASSB.LT.MXPASS) GO TO 110
C          Just give up if it fails MXPASS times
          WRITE(ITLIS,998)
998       FORMAT(//5X,'ERROR IN MBIAS ... COULD NOT MAKE BARYON')
          XBARY(IB)=0.
          DXBARY(IB)=1.
        ENDIF
C
C          Having accepted baryon, set up XPOM array for cut Pomerons,
C          rescaling to 1.-XBARY(IB).
        XSUM=0.
        DO 120 N=1,NPOM
          XX=RANF()
          XPOM(N,IB)=XX
          XSUM=XSUM+XX
120     CONTINUE
        XSUM=1./XSUM
        DO 130 N=1,NPOM
          XPOM(N,IB)=XSUM*XPOM(N,IB)*DXBARY(IB)
130     CONTINUE
100   CONTINUE
C
C          Fragment each Pomeron into mesons and baryon pairs in the
C          Pomeron-Pomeron center of mass.
C
      DO 1000 IB=1,2
        DO 2000 IPOM=1,NPOM
          PPOM=SQRT(PBEAM(1)*XPOM(IPOM,1)*PBEAM(2)*XPOM(IPOM,2))
          PPLUS=2.*PPOM
          NFIRST=NPTCL+1
          LOOP=0
C
200       CONTINUE
          ITWIST=INT(1.99999*RANF())+1
          LOOP=LOOP+1
C
C          Select new quark or diquark. Old diquark implies new quark.
C          Old quark implies new diquark with probability PBARY0.
          IFL1=IFLEXC(ITWIST)
          IF(MOD(IFL1,100).EQ.0) THEN
            IFL2=ISIGN(INT(RANF()/PUD0)+1,+IFL1)
          ELSEIF(RANF().GT.PBARY0) THEN
            IFL2=ISIGN(INT(RANF()/PUD0)+1,-IFL1)
          ELSE
            ID1=INT(RANF()/PUD0)+1
            ID2=INT(RANF()/PUD0)+1
            IF(IABS(ID1).GT.IABS(ID2)) THEN
              ISWAP=ID1
              ID1=ID2
              ID2=ISWAP
            ENDIF
            IFL2=ISIGN(1000*ID1+100*ID2,IFL1)
          ENDIF
          IFLEXC(ITWIST)=-IFL2
C          Construct meson from quark+antiquark. Else, construct baryon
C          IDENT from quark+diquark.
          IMOD1=MOD(IFL1,100)
          IMOD2=MOD(IFL2,100)
          IF(IMOD1.NE.0.AND.IMOD2.NE.0) THEN
            JSPIN=INT(RANF()+PJSPN)
            ID1=IFL1
            ID2=IFL2
            IF(ID1+ID2.EQ.0) THEN
              RND=RANF()
              ID1=IABS(ID1)
              ID1=INT(PMIX01(ID1,JSPIN+1)+RND)
     $        +INT(PMIX02(ID1,JSPIN+1)+RND)+1
              ID2=-ID1
            ELSEIF(IABS(ID1).GT.IABS(ID2)) THEN
              ISWAP=ID1
              ID1=ID2
              ID2=ISWAP
            ENDIF
            IDHAD=ISIGN(100*IABS(ID1)+10*IABS(ID2)+JSPIN,ID1)
          ELSE
            IF(IMOD1.EQ.0) THEN
              ID3=MOD(IFL1/100,10)
              ID2=IFL1/1000
              ID1=IFL2
            ELSE
              ID3=MOD(IFL2/100,10)
              ID2=IFL2/1000
              ID1=IFL1
            ENDIF
            IF(IABS(ID1).GT.IABS(ID2)) THEN
              ISWAP=ID1
              ID1=ID2
              ID2=ISWAP
            ENDIF
            IF(IABS(ID2).GT.IABS(ID3)) THEN
              ISWAP=ID2
              ID2=ID3
              ID3=ISWAP
            ENDIF
            IF(IABS(ID1).GT.IABS(ID2)) THEN
              ISWAP=ID1
              ID1=ID2
              ID2=ISWAP
            ENDIF
            IF(ID1.EQ.ID2.AND.ID2.EQ.ID3) THEN
              JSPIN=1
            ELSE
              JSPIN=INT(RANF()+PJSPN)
            ENDIF
            IF(JSPIN.EQ.0.AND.ID1.NE.ID2.AND.ID2.NE.ID3) THEN
              IF(RANF().LT.PISPN) THEN
                ISWAP=ID1
                ID1=ID2
                ID2=ISWAP
              ENDIF
            ENDIF
            IDHAD=1000*IABS(ID1)+100*IABS(ID2)+10*IABS(ID3)+JSPIN
            IDHAD=ISIGN(IDHAD,IFL1)
          ENDIF
C
          AM=AMASS(IDHAD)
          PX1=PXEXC(ITWIST)
          PY1=PYEXC(ITWIST)
          CALL GETPT(PT2,SIGQT0)
          PHI2=2.*PI*RANF()
          PX2=PT2*COS(PHI2)
          PY2=PT2*SIN(PHI2)
          PXEXC(ITWIST)=PX2
          PYEXC(ITWIST)=PY2
          PX=PX1-PX2
          PY=PY1-PY2
          AMT2=PX**2+PY**2+AM**2
C
C          Select x -- same distribution for all particles.
          X=RANF()
          IF(RANF().LT.XGEN0(1)) X=1.-X**(1./(XGEN0(2)+1.))
          QPLUS=X*PPLUS
          QPLUS=AMAX1(QPLUS,1.E-6)
          QMINUS=AMT2/QPLUS
          P0=.5*(QPLUS+QMINUS)
          PZ=.5*(QPLUS-QMINUS)
C
C          Add particle to /PARTCL/ if PZ>0.
          IF(NPTCL.GE.MXPTCL) GO TO 9999
          IF(PZ.GE.0.) THEN
            NPTCL=NPTCL+1
            PPTCL(1,NPTCL)=PX
            PPTCL(2,NPTCL)=PY
            PPTCL(3,NPTCL)=PZ*SIGN(IB)
            PPTCL(4,NPTCL)=P0
            PPTCL(5,NPTCL)=AM
            IORIG(NPTCL)=0
            IDCAY(NPTCL)=0
            IDENT(NPTCL)=IDHAD
          ENDIF
C
C          Continue if sufficient pplus
          PPLUS=(1.-X)*PPLUS
          IF(PPLUS.GT.PEND0.AND.LOOP.LT.MXPTCL) GO TO 200
C
C          Boost hadrons to lab frame.
          IF(NPTCL.LT.NFIRST) GO TO 2000
          BETA=(XPOM(IPOM,1)*PBEAM(1)-XPOM(IPOM,2)*PBEAM(2))/(2.*PPOM)
          GAM=(XPOM(IPOM,1)*PBEAM(1)+XPOM(IPOM,2)*PBEAM(2))/(2.*PPOM)
          DO 400 IP=NFIRST,NPTCL
            P0=GAM*PPTCL(4,IP)+BETA*PPTCL(3,IP)
            PZ=BETA*PPTCL(4,IP)+GAM*PPTCL(3,IP)
            PPTCL(3,IP)=PZ
            PPTCL(4,IP)=P0
400       CONTINUE
C
2000    CONTINUE
1000  CONTINUE
C
C          Rescale hadron momenta for correct four-momentum.
C
      NPTLV1=NPTCL
      IF(KEYS(4)) THEN
        PSUM(4)=ECM
        PSUM(5)=ECM
        CALL RESCAL(NBEGIN,NPTLV1,PSUM,IFAIL)
      ELSE
        CALL RESCAL(NBEGIN,NPTLV1,PBEAMS,IFAIL)
      ENDIF
      IF(IFAIL.NE.0.AND.IPASS.LT.MXPASS) THEN
        IPASS=IPASS+1
        NPTCL=NBEGIN-1
        GO TO 1
      ENDIF
C
C          Decay hadrons
C
      NP1=NBEGIN
500   NP2=NPTCL
      DO 510 I=NP1,NP2
      CALL DECAY(I)
510   CONTINUE
      NP1=NP2+1
      IF(NP1.LE.NPTCL) GO TO 500
      RETURN
C
9999  CALL PRTEVT(0)
      WRITE(ITLIS,999) NPTCL
999   FORMAT(//5X,'ERROR IN MBIAS...NPTCL >',I5)
      RETURN
      END
