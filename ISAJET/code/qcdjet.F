#include "isajet/pilot.h"
      SUBROUTINE QCDJET(NJMIN)
C
C          Carry out final state QCD jet evolution using the algorithm
C          of Fox and Wolfram.  Evolve each parton in T with fixed ZC
C          and iterate as follows--
C
C          (0) Evolve initial partons.
C          (1) Pick I and find matching J>I.
C          (2) Solve kinematics.
C          (3) For K=I,J, generate Z(K) and evolve T(K1), T(K2).  If no
C              good, evolve T(K). Otherwise, add K1 and K2 to /JETSET/.
C          (4) If I or J no good, then (2).
C          (5) Then (1).
C
C          Use Z=(E+P)/(E0+P0) and a large TCUT.
C          JMATCH(J1)=J2 if J1 and J2 match.
C          JMATCH(J)=JPACK*J1+J2 if J1,...,J2 match. Used for multiple
C          initial partons.
C          JMATCH(J)=0 for initial jet partons
C
C          Include W+- and Z0 radiation.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/partcl.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/const.inc"
C
      INTEGER J,NJMIN,JPRNT,JI1,JI2,NJI,JI,NJ1,NJ2,L,K,NPTLV1,IFAIL,J0
      REAL AM0,AM1,AM2,RANF,AMSUM,PCM2,POLD2,RATIO,PSUM,P12CM,E0,P0,Z1,
     $E1MAX,P1MAX,ZMAX,E1MIN,P1MIN,ZMIN,ZEP,E1,P1,CTHCMZ,Z2,E2MAX,P2MAX,
     $E2MIN,P2MIN,P2,E2,CTHCM,STHCM,PHICM,CPHICM,SPHICM,PT0,CTH0,STH0,
     $CPHI0,SPHI0,SGN,BP,ZLIM,ZLIM1
      DIMENSION PSUM(5)
      DATA PSUM/5*0./
C
C          (0) Evolve initial parton masses.
C
      DO 100 J=NJMIN,NJSET
        J1=J
        J2=JMATCH(J)
        IF(J2.GT.JPACK) GO TO 150
        IF(J2.LE.J1) GO TO 100
C          Two partons
        IF(JDCAY(J1).EQ.-1) CALL QCDT(J1)
        IF(JDCAY(J2).EQ.-1) CALL QCDT(J2)
        JPRNT=MOD(JORIG(J),JPACK)
        IF(JPRNT.EQ.0) THEN
          AM0=PJSET(4,J1)+PJSET(4,J2)
        ELSE
          AM0=PJSET(5,JPRNT)
        ENDIF
110     AM1=PJSET(5,J1)
        AM2=PJSET(5,J2)
        IF(AM0.LE.AM1+AM2) THEN
          J3=J1
          IF(RANF().GT..5) J3=J2
          IF(JDCAY(J3).EQ.-1) CALL QCDT(J3)
          GO TO 110
        ENDIF
        GO TO 100
C          More than two partons
150     JI1=JMATCH(J)/JPACK
        IF(J.NE.JI1) GO TO 100
        JI2=JMATCH(J)-JPACK*JI1
        NJI=JI2-JI1+1
        AM0=0.
        AMSUM=0.
        DO 160 JI=JI1,JI2
          IF(JDCAY(JI).EQ.-1) CALL QCDT(JI)
          AM0=AM0+PJSET(4,JI)
          AMSUM=AMSUM+PJSET(5,JI)
160     CONTINUE
170     IF(AM0.LT.AMSUM) THEN
          J3=NJI*RANF()+JI1
          AMSUM=AMSUM-PJSET(5,J3)
          IF(JDCAY(J3).EQ.-1) CALL QCDT(J3)
          AMSUM=AMSUM+PJSET(5,J3)
          GO TO 170
        ENDIF
100   CONTINUE
C
C          (1) Loop over active partons
C
      NJ1=NJMIN
1     NJ2=NJSET
      DO 200 J=NJ1,NJ2
        J1=J
        J2=JMATCH(J1)
        NJI=2
        IF(J2.LE.J1) GO TO 200
C
C          (2) Solve kinematics.
C
C          Initial partons--keep directions fixed.
210     IF(MOD(JORIG(J),JPACK).NE.0) GO TO 230
        IF(JMATCH(J).GT.JPACK) GO TO 400
        AM0=PJSET(4,J1)+PJSET(4,J2)
        AM1=PJSET(5,J1)
        AM2=PJSET(5,J2)
        PJSET(4,J1)=(AM0**2+AM1**2-AM2**2)/(2*AM0)
        PJSET(4,J2)=(AM0**2+AM2**2-AM1**2)/(2*AM0)
        PCM2=((AM0**2-AM1**2-AM2**2)**2-(2*AM1*AM2)**2)/(4*AM0**2)
        DO 220 L=1,2
          POLD2=PJSET(1,JJ(L))**2+PJSET(2,JJ(L))**2+PJSET(3,JJ(L))**2
          RATIO=SQRT(PCM2/POLD2)
          DO 225 K=1,3
225       PJSET(K,JJ(L))=RATIO*PJSET(K,JJ(L))
220     CONTINUE
        GO TO 300
C
C          NJI.LE.5 initial partons
400     CONTINUE
        JI1=JMATCH(J)/JPACK
        IF(J.NE.JI1) GO TO 200
        JI2=JMATCH(J)-JPACK*JI1
        NJI=JI2-JI1+1
        AM0=0.
        DO 410 JI=JI1,JI2
          AM0=AM0+PJSET(4,JI)
          JJ(JI-JI1+1)=JI
          PJSET(4,JI)=SQRT(PJSET(1,JI)**2+PJSET(2,JI)**2+PJSET(3,JI)**2
     1    +PJSET(5,JI)**2)
          DO 420 K=1,5
420       PPTCL(K,NPTCL+JI-JI1+1)=PJSET(K,JI)
410     CONTINUE
        PSUM(4)=AM0
        PSUM(5)=PSUM(4)
        NPTLV1=NPTCL
        CALL RESCAL(NPTLV1+1,NPTLV1+NJI,PSUM,IFAIL)
        DO 430 JI=JI1,JI2
        DO 430 K=1,5
          PJSET(K,JI)=PPTCL(K,NPTCL+JI-JI1+1)
430     CONTINUE
        GO TO 300
C
C          Solve kinematics for general partons.
C
230     J0=MOD(JORIG(J),JPACK)
        AM0=PJSET(5,J0)
        AM1=PJSET(5,J1)
        AM2=PJSET(5,J2)
        E1CM=(AM0**2+AM1**2-AM2**2)/(2*AM0)
        E2CM=(AM0**2+AM2**2-AM1**2)/(2*AM0)
        P12CM=SQRT((AM0**2-AM1**2-AM2**2)**2-(2*AM1*AM2)**2)/(2*AM0)
        NJI=2
C          Determine E1, P1, and COS(THCM) from Z(J0).
C          Occasionally COS(TH)>1.  If so then reset Z.
        E0=PJSET(4,J0)
        P0=SQRT(PJSET(1,J0)**2+PJSET(2,J0)**2+PJSET(3,J0)**2)
        Z1=ZZC(J0)
        IF(Z1.GT.0.5) THEN
          E1MAX=(E0*E1CM+P0*P12CM)/AM0
          P1MAX=(P0*E1CM+E0*P12CM)/AM0
          ZMAX=(E1MAX+P1MAX)/(E0+P0)
          E1MIN=(E0*E1CM-P0*P12CM)/AM0
          P1MIN=(P0*E1CM-E0*P12CM)/AM0
          P1MIN=ABS(P1MIN)
          ZMIN=(E1MIN+P1MIN)/(E0+P0)
          IF(Z1.LT.ZMIN.OR.Z1.GT.ZMAX) Z1=ZMIN+Z1*(ZMAX-ZMIN)
          ZZC(J0)=Z1
          ZEP=Z1*(E0+P0)
          P1=(ZEP**2-AM1**2)/(2.*ZEP)
          E1=(ZEP**2+AM1**2)/(2.*ZEP)
          CTHCM=(E1*AM0-E0*E1CM)/(P0*P12CM)
        ELSE
          Z2=1.-Z1
          E2MAX=(E0*E2CM+P0*P12CM)/AM0
          P2MAX=(P0*E2CM+E0*P12CM)/AM0
          ZMAX=(E2MAX+P2MAX)/(E0+P0)
          E2MIN=(E0*E2CM-P0*P12CM)/AM0
          P2MIN=(P0*E2CM-E0*P12CM)/AM0
          P2MIN=ABS(P2MIN)
          ZMIN=(E2MIN+P2MIN)/(E0+P0)
          IF(Z2.LT.ZMIN.OR.Z2.GT.ZMAX) Z2=ZMIN+Z2*(ZMAX-ZMIN)
          ZZC(J0)=Z2
          ZEP=Z2*(E0+P0)
          P2=(ZEP**2-AM2**2)/(2.*ZEP)
          E2=(ZEP**2+AM2**2)/(2.*ZEP)
          CTHCM=-(E2*AM0-E0*E2CM)/(P0*P12CM)
        ENDIF
C          Avoid disaster
        IF(ABS(CTHCM).GT.1.) CTHCM=SIGN(RANF(),CTHCM)
        STHCM=SQRT(1.-CTHCM**2)
        PHICM=2*PI*RANF()
        CPHICM=COS(PHICM)
        SPHICM=SIN(PHICM)
C
C          Construct cm momenta.
        PT0=SQRT(PJSET(1,J0)**2+PJSET(2,J0)**2)
        CTH0=PJSET(3,J0)/P0
        STH0=PT0/P0
        CPHI0=PJSET(1,J0)/PT0
        SPHI0=PJSET(2,J0)/PT0
        P1CM(1)=P12CM*(CPHI0*(CTH0*CPHICM*STHCM+STH0*CTHCM)
     1  -SPHI0*SPHICM*STHCM)
        P1CM(2)=P12CM*(SPHI0*(CTH0*CPHICM*STHCM+STH0*CTHCM)
     1  +CPHI0*SPHICM*STHCM)
        P1CM(3)=P12CM*(-STH0*CPHICM*STHCM+CTH0*CTHCM)
C          Boost with P0 to get lab momenta
        DO 240 L=1,2
          SGN=3-2*L
          BP=0
          DO 241 K=1,3
241       BP=BP+PJSET(K,J0)*SGN*P1CM(K)
          BP=BP/AM0
          PJSET(4,JJ(L))=PJSET(4,J0)*EE(L)/PJSET(5,J0)+BP
          DO 242 K=1,3
242       PJSET(K,JJ(L))=SGN*P1CM(K)+PJSET(K,J0)*EE(L)/PJSET(5,J0)
     1    +PJSET(K,J0)*BP/(PJSET(4,J0)+PJSET(5,J0))
240     CONTINUE
C
C          (3) Pick Z and decay partons. Check.
C
300     CONTINUE
        TNEW=.FALSE.
        DO 310 L=1,NJI
          IF(JDCAY(JJ(L)).GE.0) GO TO 310
          IF(NJSET+2.GT.MXJSET) GO TO 9999
          CALL QCDZ(JJ(L))
          CALL QCDT(NJSET+1)
          CALL QCDT(NJSET+2)
C
C          Check whether masses allowed.
          AM0=PJSET(5,JJ(L))
          AM1=PJSET(5,NJSET+1)
          AM2=PJSET(5,NJSET+2)
          IF(AM1+AM2.GE.AM0) GO TO 320
C
C          Check whether Z allowed.
          E1CM=(AM0**2+AM1**2-AM2**2)/(2*AM0)
          E2CM=(AM0**2+AM2**2-AM1**2)/(2.*AM0)
          P12CM=SQRT((AM0**2-AM1**2-AM2**2)**2-(2*AM1*AM2)**2)/(2*AM0)
          E0=PJSET(4,JJ(L))
          P0=SQRT(PJSET(1,JJ(L))**2+PJSET(2,JJ(L))**2+PJSET(3,JJ(L))**2)
          IF(ZZC(JJ(L)).GT.0.5) THEN
            ZEP=ZZC(JJ(L))*(E0+P0)
            P1=(ZEP**2-AM1**2)/(2.*ZEP)
            E1=(ZEP**2+AM1**2)/(2.*ZEP)
            CTHCM=(E1*AM0-E0*E1CM)/(P0*P12CM)
            IF((ABS(CTHCM).GE.1..OR.P1.LE.0.).AND.IABS(JTYPE(JJ(L)))
     $      .LT.80) GO TO 320
          ELSE
            ZEP=(1.-ZZC(JJ(L)))*(E0+P0)
            P2=(ZEP**2-AM2**2)/(2.*ZEP)
            E2=(ZEP**2+AM2**2)/(2.*ZEP)
            CTHCM=-(E2*AM0-E0*E2CM)/(P0*P12CM)
            IF((ABS(CTHCM).GE.1..OR.P2.LE.0.).AND.IABS(JTYPE(JJ(L)))
     $      .LT.80) GO TO 320
          ENDIF
C
C          Require Z and 1-Z within kinematic limits.
C
          ZLIM=(AM0/(E0+P0))**2
          ZLIM1=CUTJET/(E0+P0)
          ZLIM=AMAX1(ZLIM,ZLIM1)
          IF((ZZC(JJ(L)).GT.ZLIM.AND.ZZC(JJ(L)).LT.(1.-ZLIM)).OR.
     $    IABS(JTYPE(JJ(L))).GE.80) THEN
C          Add new partons to /JETSET/.
            JDCAY(JJ(L))=JPACK*(NJSET+1)+(NJSET+2)
            NJSET=NJSET+2
            GO TO 310
          ENDIF
C          Discard partons and evolve JJ(L) again.
320       TNEW=.TRUE.
          CALL QCDT(JJ(L))
310     CONTINUE
C
C          (4) Resolve kinematics if any parton mass is changed.
C
        IF(TNEW) GO TO 210
200   CONTINUE
C
C          (5) Iterate entire proceedure.
C
      NJ1=NJ2+1
      IF(NJ1.LE.NJSET) GO TO 1
      RETURN
C
C          Error message
C
9999  CALL PRTEVT(0)
      WRITE(ITLIS,10) NJSET
10    FORMAT(//' ERROR IN QCDJET...NJSET > ',I4)
      RETURN
      END
