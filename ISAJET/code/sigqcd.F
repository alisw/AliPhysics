#include "isajet/pilot.h"
      SUBROUTINE SIGQCD
C
C          Compute D(SIGMA)/D(PT**2)D(Y1)D(Y2)
C          Include quark masses for ch, bt, and tp and 4th generation.
C          Note ch is now treated as heavy.
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE card.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
C
C          Cross sections from Feynman, Field and Fox, P.R. D18, 3320
C          Massive cross sections from B. Combridge, N.P. B151, 429.
C          Extra factor of 1/2 needed for  non-identical jets since all
C          all jets are treated as identical.
C
C          Ver 6.35: Fix kinematics for gl + tp -> gl + tp, etc.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/const.inc"
C
      REAL    X(2),QSAVE(13,2),EBT(2)
      EQUIVALENCE (X(1),X1),(S,SHAT),(T,THAT),(U,UHAT)
      REAL    FFF1,FFF2,FFF3,FFF4,FFF5,FFF6,FFF7,S,T,U,FGQ,AM2,FQQ,
     $        QFCN,STRUC,FJAC,SIG,AMASS,SIG1,AMQ,FJACBT,SIG2,QQ,XQMIN,
     $        E1,E2
      INTEGER IQ,IH,I,J,IFL,JTYP1,JTYP2,IQ1,IQ2
C
C          Elementary cross sections from Feynman, Field, and Fox.
C
      FFF1(S,T,U)=4./9.*(S**2+U**2)/T**2
      FFF2(S,T,U)=4./9.*((S**2+U**2)/T**2+(S**2+T**2)/U**2)
     1-8./27.*S**2/(U*T)
      FFF3(S,T,U)=4./9.*((S**2+U**2)/T**2+(T**2+U**2)/S**2)
     1-8./27.*U**2/(S*T)
      FFF4(S,T,U)=32./27.*(U**2+T**2)/(U*T)-8./3.*(U**2+T**2)/S**2
      FFF5(S,T,U)=1./6.*(U**2+T**2)/(U*T)-3./8.*(U**2+T**2)/S**2
      FFF6(S,T,U)=-4./9.*(U**2+S**2)/(U*S)+(U**2+S**2)/T**2
      FFF7(S,T,U)=9./2.*(3.-U*T/S**2-U*S/T**2-S*T/U**2)
C          Heavy quark cross sections from Combridge
      FGQ(S,T,U)=2.*(S-AM2)*(AM2-U)/T**2
     1+4./9.*((S-AM2)*(AM2-U)+2.*AM2*(S+AM2))/(S-AM2)**2
     2+4./9.*((S-AM2)*(AM2-U)+2.*AM2*(AM2+U))/(AM2-U)**2
     3+1./9.*AM2*(4.*AM2-T)/((S-AM2)*(AM2-U))
     4+((S-AM2)*(AM2-U)+AM2*(S-U))/(T*(S-AM2))
     5-((S-AM2)*(AM2-U)-AM2*(S-U))/(T*(AM2-U))
      FQQ(S,T,U)=4./9.*((AM2-U)**2+(S-AM2)**2+2.*AM2*T)/T**2
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C
C          Use massless kinematics for ch and lighter quarks.
C
      CALL TWOKIN(0.,0.,0.,0.)
      FJAC=SHAT/SCM*UNITS
      FJAC=FJAC*PI*ALFQSQ**2/SHAT**2
C
C          Initialize cross sections.
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
        SIGS(I)=0.
100   CONTINUE
      IF(X1.GE.1.0.OR.X2.GE.1.0) RETURN
C          Compute structure functions
      DO 110 IH=1,2
      DO 110 IQ=1,7
        QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
110   CONTINUE
C
C          Compute cross sections summed over quark types allowed by
C          JETTYPE card.
C
C          Gluon-gluon
      IF(.NOT.(GOQ(1,1).AND.GOQ(1,2))) GO TO 210
      SIG=.5*FJAC*QSAVE(1,1)*QSAVE(1,2)*FFF7(S,T,U)
      CALL SIGFIL(SIG,1,1,1,1)
C
      DO 201 I=1,3
        SIG=.5*FJAC*QSAVE(2*I,1)*QSAVE(2*I+1,2)*FFF4(S,T,U)
        CALL SIGFIL(SIG,2*I,2*I+1,1,1)
        SIG=.5*FJAC*QSAVE(2*I+1,1)*QSAVE(2*I,2)*FFF4(S,U,T)
        CALL SIGFIL(SIG,2*I+1,2*I,1,1)
201   CONTINUE
C
C          Quark-gluon
210   CONTINUE
      DO 211 I=2,7
        IF(.NOT.(GOQ(I,1).AND.GOQ(1,2))) GO TO 212
        SIG=.5*FJAC*QSAVE(I,1)*QSAVE(1,2)*FFF6(S,T,U)
        CALL SIGFIL(SIG,I,1,I,1)
        SIG=.5*FJAC*QSAVE(1,1)*QSAVE(I,2)*FFF6(S,U,T)
        CALL SIGFIL(SIG,1,I,I,1)
212     CONTINUE
        IF(.NOT.(GOQ(1,1).AND.GOQ(I,2))) GO TO 211
        SIG=.5*FJAC*QSAVE(1,1)*QSAVE(I,2)*FFF6(S,T,U)
        CALL SIGFIL(SIG,1,I,1,I)
        SIG=.5*FJAC*QSAVE(I,1)*QSAVE(1,2)*FFF6(S,U,T)
        CALL SIGFIL(SIG,I,1,1,I)
211   CONTINUE
C
C          Identical quark-quark
      DO 220 I=2,7
        IF(.NOT.(GOQ(I,1).AND.GOQ(I,2))) GO TO 220
        SIG=.5*FJAC*QSAVE(I,1)*QSAVE(I,2)*FFF2(S,T,U)
        CALL SIGFIL(SIG,I,I,I,I)
220   CONTINUE
C
C          Identical quark-antiquark
      DO 230 I=1,3
        IF(SHAT.LT.4.*AMASS(I)**2) GO TO 230
        IF(.NOT.(GOQ(2*I,1).AND.GOQ(2*I+1,2))) GO TO 235
        SIG=.5*FJAC*QSAVE(1,1)*QSAVE(1,2)*FFF5(S,T,U)
        CALL SIGFIL(SIG,1,1,2*I,2*I+1)
        DO 231 J=1,3
          IF(J.EQ.I) GO TO 231
          SIG=.5*FJAC*QSAVE(2*J,1)*QSAVE(2*J+1,2)*FFF1(T,S,U)
          CALL SIGFIL(SIG,2*J,2*J+1,2*I,2*I+1)
          SIG=.5*FJAC*QSAVE(2*J+1,1)*QSAVE(2*J,2)*FFF1(T,S,U)
          CALL SIGFIL(SIG,2*J+1,2*J,2*I,2*I+1)
231     CONTINUE
        SIG=.5*FJAC*QSAVE(2*I,1)*QSAVE(2*I+1,2)*FFF3(S,T,U)
        CALL SIGFIL(SIG,2*I,2*I+1,2*I,2*I+1)
        SIG=.5*FJAC*QSAVE(2*I+1,1)*QSAVE(2*I,2)*FFF3(S,U,T)
        CALL SIGFIL(SIG,2*I+1,2*I,2*I,2*I+1)
C
235     CONTINUE
        IF(.NOT.(GOQ(2*I+1,1).AND.GOQ(2*I,2))) GO TO 230
        SIG=.5*FJAC*QSAVE(1,1)*QSAVE(1,2)*FFF5(S,T,U)
        CALL SIGFIL(SIG,1,1,2*I+1,2*I)
        DO 236 J=1,3
          IF(J.EQ.I) GO TO 236
          SIG=.5*FJAC*QSAVE(2*J,1)*QSAVE(2*J+1,2)*FFF1(T,S,U)
          CALL SIGFIL(SIG,2*J,2*J+1,2*I+1,2*I)
          SIG=.5*FJAC*QSAVE(2*J+1,1)*QSAVE(2*J,2)*FFF1(T,S,U)
          CALL SIGFIL(SIG,2*J+1,2*J,2*I+1,2*I)
236     CONTINUE
        SIG1=.5*FJAC*QSAVE(2*I,1)*QSAVE(2*I+1,2)*FFF3(S,U,T)
        CALL SIGFIL(SIG1,2*I,2*I+1,2*I+1,2*I)
        SIG=.5*FJAC*QSAVE(2*I+1,1)*QSAVE(2*I,2)*FFF3(S,T,U)
        CALL SIGFIL(SIG,2*I+1,2*I,2*I+1,2*I)
230   CONTINUE
C
C          General massless quark-quark
      DO 240 I=2,7
        DO 241 J=2,7
          IF(.NOT.(GOQ(I,1).AND.GOQ(J,2))) GO TO 241
          IF((I/2).EQ.(J/2)) GO TO 241
          SIG=.5*FJAC*QSAVE(I,1)*QSAVE(J,2)*FFF1(S,T,U)
          CALL SIGFIL(SIG,I,J,I,J)
          SIG=.5*FJAC*QSAVE(J,1)*QSAVE(I,2)*FFF1(S,U,T)
          CALL SIGFIL(SIG,I,J,J,I)
241     CONTINUE
240   CONTINUE
C
C          CH+CB, BT+BB, and TP+TB cross sections.
C          Y=-log(tan(theta/2)), so Jacobean contains P1*P2/E1*E2.
C          Also fourth generation.
C
      DO 250 IQ=1,5
        IFL=IQ+3
        JTYP1=2*IFL
        JTYP2=JTYP1+1
        IF(.NOT.((GOQ(JTYP1,1).AND.GOQ(JTYP2,2)).OR.
     1  (GOQ(JTYP2,1).AND.GOQ(JTYP1,2)))) GO TO 250
        AMQ=AMASS(IFL)
        IF(AMQ.LT.0.) GO TO 250
        AM2=AMQ**2
        CALL TWOKIN(0.,0.,AMQ,AMQ)
        IF(X(1).GE.1..OR.X(2).GE.1.) GO TO 250
        EBT(1)=SQRT(P(1)**2+AM2)
        EBT(2)=SQRT(P(2)**2+AM2)
        FJACBT=.5*S/SCM*UNITS*P(1)*P(2)/(EBT(1)*EBT(2))
        SIG1=12.*(AM2-T)*(AM2-U)/S**2
     1  +8./3.*((AM2-T)*(AM2-U)-2.*AM2*(AM2+T))/(AM2-T)**2
     2  +8./3.*((AM2-T)*(AM2-U)-2.*AM2*(AM2+U))/(AM2-U)**2
     3  -2./3.*AM2*(S-4.*AM2)/((AM2-T)*(AM2-U))
     4  -6.*((AM2-T)*(AM2-U)+AM2*(U-T))/(S*(AM2-T))
     5  -6.*((AM2-T)*(AM2-U)+AM2*(T-U))/(S*(AM2-U))
        SIG1=SIG1*PI**2*ALFQSQ**2/(16.*PI*S**2)
        SIG=FJACBT*SIG1*STRUC(X(1),QSQ,1,IDIN(1))/X(1)
     1  *STRUC(X(2),QSQ,1,IDIN(2))/X(2)
        IF(GOQ(JTYP1,1).AND.GOQ(JTYP2,2)) 
     $  CALL SIGFIL(SIG,1,1,JTYP1,JTYP2)
        IF(GOQ(JTYP2,1).AND.GOQ(JTYP1,2)) 
     $  CALL SIGFIL(SIG,1,1,JTYP2,JTYP1)
C
        SIG2=((AM2-T)**2+(AM2-U)**2+2.*S*AM2)/S**2
        SIG2=FJACBT*SIG2*64.*PI**2*ALFQSQ**2/(9.*16.*PI*S**2)
        DO 255 I=1,3
          QQ=STRUC(X(1),QSQ,2*I,IDIN(1))*STRUC(X(2),QSQ,2*I+1,IDIN(2))
          SIG=SIG2*QQ/(X(1)*X(2))
          IF(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))
     $    CALL SIGFIL(SIG,2*I,2*I+1,JTYP1,JTYP2)
          IF(GOQ(JTYP2,1).AND.GOQ(JTYP1,2))
     $    CALL SIGFIL(SIG,2*I,2*I+1,JTYP2,JTYP1)
          QQ=STRUC(X(1),QSQ,2*I+1,IDIN(1))*STRUC(X(2),QSQ,2*I,IDIN(2))
          SIG=SIG2*QQ/(X(1)*X(2))
          IF(GOQ(JTYP1,1).AND.GOQ(JTYP2,2))
     $    CALL SIGFIL(SIG,2*I+1,2*I,JTYP1,JTYP2)
          IF(GOQ(JTYP2,1).AND.GOQ(JTYP1,2))
     $    CALL SIGFIL(SIG,2*I+1,2*I,JTYP2,JTYP1)
255     CONTINUE
250   CONTINUE
C
C          Gluon + heavy quark
      DO 300 IQ=8,13.
        IF(.NOT.(GOQ(1,1).AND.GOQ(IQ,2))) GO TO 310
        AMQ=AMASS(IQ/2)
        AM2=AMQ**2
        XQMIN=AMQ/ECM
        E1=P(1)
        E2=SQRT(P(2)**2+AM2)
        FJAC=.5*S/SCM*UNITS*PI*ALFQSQ**2/S**2
        CALL TWOKIN(0.,AMQ,0.,AMQ)
        IF(X(1).LT.1..AND.X(2).LT.1..AND.X(2).GT.XQMIN) THEN
          SIG=FJAC*P(1)*P(2)/(E1*E2)*FGQ(S,T,U)*QFCN(1,1)*QFCN(IQ,2)
          CALL SIGFIL(SIG,1,IQ,1,IQ)
        ENDIF
        CALL TWOKIN(AMQ,0.,0.,AMQ)        
        IF(X(1).LT.1..AND.X(2).LT.1..AND.X(1).GT.XQMIN) THEN
          SIG=FJAC*P(1)*P(2)/(E1*E2)*FGQ(S,U,T)*QFCN(IQ,1)*QFCN(1,2)
          CALL SIGFIL(SIG,IQ,1,1,IQ)
        ENDIF
C
310     IF(.NOT.(GOQ(IQ,1).AND.GOQ(1,2))) GO TO 300
        AMQ=AMASS(IQ/2)
        AM2=AMQ**2
        XQMIN=AMQ/ECM
        E1=SQRT(P(1)**2+AM2)
        E2=P(2)
        FJAC=.5*S/SCM*UNITS*PI*ALFQSQ**2/S**2
        CALL TWOKIN(0.,AMQ,AMQ,0.)
        IF(X(1).LT.1..AND.X(2).LT.1..AND.X(2).GT.XQMIN) THEN
          SIG=FJAC*P(1)*P(2)/(E1*E2)*FGQ(S,U,T)*QFCN(1,1)*QFCN(IQ,2)
          CALL SIGFIL(SIG,1,IQ,IQ,1)
        ENDIF
        CALL TWOKIN(AMQ,0.,AMQ,0.)
        IF(X(1).LT.1..AND.X(2).LT.1..AND.X(1).GT.XQMIN) THEN
          SIG=FJAC*P(1)*P(2)/(E1*E2)*FGQ(S,T,U)*QFCN(IQ,1)*QFCN(1,2)
          CALL SIGFIL(SIG,IQ,1,IQ,1)
        ENDIF
300   CONTINUE
C
C          Light quark + heavy quark
      DO 320 IQ1=2,7
        DO 330 IQ2=8,13
          IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 340
          AMQ=AMASS(IQ2/2)
          AM2=AMQ**2
          XQMIN=AMQ/ECM
          E1=P(1)
          E2=SQRT(P(2)**2+AM2)
          FJAC=.5*S/SCM*UNITS*PI*ALFQSQ**2/S**2
          CALL TWOKIN(0.,AMQ,0.,AMQ)
          IF(X(1).LT.1..AND.X(2).LT.1..AND.X(2).GT.XQMIN) THEN
            SIG=FJAC*P(1)*P(2)/(E1*E2)*FQQ(S,T,U)*QFCN(IQ1,1)
     $      *QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IQ1,IQ2)
          ENDIF
          CALL TWOKIN(AMQ,0.,0.,AMQ)
          IF(X(1).LT.1..AND.X(2).LT.1..AND.X(1).GT.XQMIN) THEN
            SIG=FJAC*P(1)*P(2)/(E1*E2)*FQQ(S,U,T)*QFCN(IQ1,2)
     $      *QFCN(IQ2,1)
            CALL SIGFIL(SIG,IQ2,IQ1,IQ1,IQ2)
          ENDIF
C
340       IF(.NOT.(GOQ(IQ1,2).AND.GOQ(IQ2,1))) GO TO 330
          AMQ=AMASS(IQ2/2)
          AM2=AMQ**2
          XQMIN=AMQ/ECM
          E1=SQRT(P(1)**2+AM2)
          E2=P(2)
          FJAC=.5*S/SCM*UNITS*PI*ALFQSQ**2/S**2
          CALL TWOKIN(0.,AMQ,AMQ,0.)
          IF(X(1).LT.1..AND.X(2).LT.1..AND.X(2).GT.XQMIN) THEN
            SIG=FJAC*P(1)*P(2)/(E1*E2)*FQQ(S,U,T)*QFCN(IQ1,1)
     $      *QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IQ2,IQ1)
          ENDIF
          CALL TWOKIN(AMQ,0.,AMQ,0.)      
          IF(X(1).LT.1..AND.X(2).LT.1..AND.X(1).GT.XQMIN) THEN
            SIG=FJAC*P(1)*P(2)/(E1*E2)*FQQ(S,T,U)*QFCN(IQ1,2)
     $      *QFCN(IQ2,1)
            CALL SIGFIL(SIG,IQ2,IQ1,IQ2,IQ1)
          ENDIF
330     CONTINUE
320   CONTINUE
C
      RETURN
      END
