#include "isajet/pilot.h"
      SUBROUTINE SIGWW2
C
C          Calculate WPAIR decay distribution
C          D(SIGMA)/D(PT**2)D(Y1)D(Y2)D(OMEGA1)D(OMEGA2)
C          for modes selected in WPAIR.
C
C          Also fix the initial parton types to those selected.
C
C          Cross sections from SCHOONSCHIP (1980) neglecting W width
C          and quark masses. Hence use zero-mass vectors PZERO from
C          WPAIR to define kinematics.
C          QK(P1) + QB(P2) --> W1(P3) + W2(P4)
C                   W1(P3) --> QK(Q1) + QB(Q2)
C                   W2(P4) --> QK(Q3) + QB(Q4)
C          S=(P3+P4)**2,  T=(P3-P1)**2,  U=(P3-P2)**2
C          S1=(Q1+P4)**2, T1=(Q1-P1)**2, U1=(Q1-P2)**2
C          S3=(Q3+P3)**2, T3=(Q3-P2)**2, U3=(Q3-P1)**2
C          S13=(Q1+Q3)**2
C          Note that the W+- final couplings have been set equal to 1.
C          in the SCHOONSCHIP formulas and must be restored.
C
C          Need double precision for 32-bit machines.
C
C          Ver. 5.35 - correct symmetrization for DN DB -> W+ W-.
C          Ver. 6.22 - use W + GM decay distributions from
C                      Cortes, Hagiwara, and Herzog, NP B278, 26 (1986)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/const.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/wwsig.inc"
#include "isajet/wwpar.inc"
C
      DIMENSION P1(5),P2(5),QSGN(6),PP1(4),PP2(4)
      EQUIVALENCE (S,SWW),(T,TWW),(U,UWW)
      EQUIVALENCE (P1(1),P1WW(1)),(P2(1),P2WW(1))
C          Double precision kinematics for 32-bit.
#if defined(CERNLIB_SINGLE)
      REAL S,T,U,T1,U1,T3,U3,P1,P2
     1,TX,UX,TT,UU,TT1,UU1,TT3,UU3,PP1,PP2
      REAL TERM,WWSS,WWST,WWTT,ZZALL,WZSS,WZST,WZSU,WZTU
     1,WGSS,WGST,WGSU,WGTU
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION S,T,U,T1,U1,T3,U3,P1,P2
     1,TX,UX,TT,UU,TT1,UU1,TT3,UU3,PP1,PP2
      DOUBLE PRECISION TERM,WWSS,WWST,WWTT,ZZALL,WZSS,WZST,WZSU,WZTU
     1,WGSS,WGST,WGSU,WGTU
#endif
      REAL P3IS3,P3IS4,FJAC,AMW1,AMW2,GAM1,GAM2,SGN,QSGN,AMASS3
      REAL P1DQ2,P2DQ1
      REAL A1,B1,A2,B2,ES,SMS,SMSZG,EQ3(12)
      REAL Q(5),QB(5),KK(5),E(5),EB(5)
      INTEGER K,JQ1,JQ3,JW1,JW2,IW1,IW2,IQ1,IQ2,IQ,ISWAPQ,JW,JZ,ISGN
      INTEGER IFLI,IFLJ,JG,IL,IW
      LOGICAL LQK1
C
      DATA QSGN/1.,-1.,-1.,1.,-1.,1./
      DATA EQ3/2.,-1.,-1.,2.,-1.,2.,0.,-3.,0.,-3.,0.,-3./
C
C          Entry
C
      ES=4*PI*ALFA
      WWSIG=0.
      IF(IDJETS(1).EQ.10.OR.IDJETS(2).EQ.10) GO TO 2
C          Normal case
      IF((IDJETS(1).EQ.80.AND.IDJETS(2).EQ.-80).OR.
     $(IDJETS(1).EQ.90.AND.IDJETS(2).EQ.90).OR.
     $(IABS(IDJETS(1)).EQ.80.AND.IDJETS(2).EQ.90)) THEN
        DO 10 K=1,4
          P3(K)=P3WW(K)
          Q1(K)=PZERO(K,1)
          Q3(K)=PZERO(K,3)
10      CONTINUE
        P3IS3=1.
        P3IS4=0.
        JQ1=1
        JQ3=3
        JW1=1
        JW2=2
        TX=T
        UX=U
C          Crossed case
      ELSE
        DO 20 K=1,4
          P3(K)=P4WW(K)
          Q1(K)=PZERO(K,3)
          Q3(K)=PZERO(K,1)
20      CONTINUE
        P3IS3=0.
        P3IS4=1.
        JQ1=3
        JQ3=1
        JW1=2
        JW2=1
        TX=U
        UX=T
      ENDIF
C          Variables
      T1=-2.*(Q1(4)*P1(4)-Q1(1)*P1(1)-Q1(2)*P1(2)-Q1(3)*P1(3))
      U1=-2.*(Q1(4)*P2(4)-Q1(1)*P2(1)-Q1(2)*P2(2)-Q1(3)*P2(3))
      T3=-2.*(Q3(4)*P2(4)-Q3(1)*P2(1)-Q3(2)*P2(2)-Q3(3)*P2(3))
      U3=-2.*(Q3(4)*P1(4)-Q3(1)*P1(1)-Q3(2)*P1(2)-Q3(3)*P1(3))
      S13=2.*(Q1(4)*Q3(4)-Q1(1)*Q3(1)-Q1(2)*Q3(2)-Q1(3)*Q3(3))
C          Jacobean for 4-body cross section in terms of squared
C          matrix exement in narrow resonance approximation--
C          1/((P**2-M**2)**2+M**2*GAM**2)=1/(2*M*GAM)*DELTA(P**2-M**2)
      FJAC=S/SCM*UNITS
      FJAC=FJAC*ALFA**4/(256.*PI*3.*S**2)
      AMW1=PJETS(5,1)
      AMW2=PJETS(5,2)
      GAM1=WGAM(JETTYP(1))
      GAM2=WGAM(JETTYP(2))
      FJAC=FJAC/(AMW1*GAM1*AMW2*GAM2)
      FJAC=FJAC*P(1)*P(2)/SQRT((P(1)**2+AMW1**2)*(P(2)**2+AMW2**2))
C          Color factor
      IF(IABS(IDPAIR(1)).LT.10) FJAC=3.*FJAC
      IF(IABS(IDPAIR(3)).LT.10) FJAC=3.*FJAC
C
C          W+ W- pair decays
C          Standard order is UP + UB --> W+ + W-
C
      IF(.NOT.((JETTYP(1).EQ.2.AND.JETTYP(2).EQ.3).OR.(JETTYP(1).EQ.3
     1.AND.JETTYP(2).EQ.2))) GO TO 200
      FJAC=.5*FJAC*AQ(2,2)**4
C
C          Select W+ W- OR W- W+, swapping T and U for latter.
      IW1=JETTYP(1)
      IW2=JETTYP(2)
C
C          Select quarks
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
      CQ=AQDP(IQ,2)**2
      CV=AQDP(IQ,1)/S+EZDP*AQDP(IQ,4)/(S-ZM2)
      CA=EZDP*BQDP(IQ,4)/(S-ZM2)
      SGN=QSGN(IQ)
      ISWAPQ=1
      IF(SGN.LT.0.) ISWAPQ=-1
      IF(ISWAPQ.GT.0) THEN
        TT=TX
        UU=UX
        TT1=T1
        UU1=U1
        TT3=T3
        UU3=U3
        DO 122 K=1,4
          PP1(K)=P1(K)
          PP2(K)=P2(K)
          P3(K)=P3IS3*P3WW(K)+P3IS4*P4WW(K)
          Q1(K)=PZERO(K,JQ1)
          Q3(K)=PZERO(K,JQ3)
122     CONTINUE
      ELSE
        TT=UX
        UU=TX
        TT1=U3
        UU1=T3
        TT3=U1
        UU3=T1
        DO 123 K=1,4
          PP1(K)=P1(K)
          PP2(K)=P2(K)
          P3(K)=P3IS4*P3WW(K)+P3IS3*P4WW(K)
          Q1(K)=PZERO(K,JQ3)
          Q3(K)=PZERO(K,JQ1)
123     CONTINUE
      ENDIF
C
      IF(IQ1.EQ.2*IQ) THEN
        TERM=WWTT(TT,UU,TT1,UU1,TT3,UU3)
        TERM=TERM-SGN*WWST(TT,UU,TT1,UU1,TT3,UU3,PP1,PP2)
        TERM=TERM+WWSS(TT,UU,TT1,UU1,TT3,UU3)
        WWSIG=TERM*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)*FJAC
      ELSE
        TERM=WWTT(UU,TT,UU1,TT1,UU3,TT3)
        TERM=TERM-SGN*WWST(UU,TT,UU1,TT1,UU3,TT3,PP2,PP1)
        TERM=TERM+WWSS(UU,TT,UU1,TT1,UU3,TT3)
        WWSIG=TERM*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)*FJAC
      ENDIF
C
      RETURN
C
C          Z0 Z0 pair decays
C          Standard order is UP + UB --> Z0 + Z0
C
200   IF(.NOT.(JETTYP(1).EQ.4.AND.JETTYP(2).EQ.4)) GO TO 300
      FJAC=.5*FJAC
C
C          Select quarks
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
      CV=AQDP(IQ,4)**2+BQDP(IQ,4)**2
      CA=2.*AQDP(IQ,4)*BQDP(IQ,4)
      CV1=AQDP(JQWW(1),4)**2+BQDP(JQWW(1),4)**2
      CA1=2.*AQDP(JQWW(1),4)*BQDP(JQWW(1),4)
      CV3=AQDP(JQWW(2),4)**2+BQDP(JQWW(2),4)**2
      CA3=2.*AQDP(JQWW(2),4)*BQDP(JQWW(2),4)
C
      TERM=ZZALL(TX,UX,T1,U1,T3,U3,P1,P2)
      IF(INITYP(1).EQ.2*IQ) THEN
        WWSIG=TERM*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)*FJAC
      ELSE
        WWSIG=TERM*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)*FJAC
      ENDIF
C
      RETURN
C
C          W+- Z0 pair decays
C          Standard order is DN + UB --> W- + Z0
C
300   JW=JW1
      JZ=JW2
      ISGN=-ISIGN(1,IDJETS(JW))
      SGN=ISGN
      CV3=AQDP(JQWW(JZ),4)**2+BQDP(JQWW(JZ),4)**2
      CA3=2.*AQDP(JQWW(JZ),4)*BQDP(JQWW(JZ),4)
      FJAC=.5*FJAC*AQ(1,2)**2
C
C          Select quarks. Formulas are for DN UB --> W- Z0.
C          Use symmetry for other cases.
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
C          Find whether IQ1 should be fermion or antifermion.
      IF(IQ1.EQ.2*(IQ1/2)) THEN
        ISWAPQ=+1
        IFLI=IQ1/2
        IFLJ=IQ2/2
      ELSE
        ISWAPQ=-1
        IFLI=IQ2/2
        IFLJ=IQ1/2
      ENDIF
C
      CS=AQDP(IQ,JETTYP(JW))*EZDP/(S-WM2)
      CT=AQDP(IQ,JETTYP(JW))*(AQDP(IFLJ,4)+BQDP(IFLJ,4))
      CU=AQDP(IQ,JETTYP(JW))*(AQDP(IFLI,4)+BQDP(IFLI,4))
C
C          SWAP T AND U AS NEEDED
      IF(ISWAPQ*ISGN.GT.0) THEN
        TT=TX
        UU=UX
        TT1=T1
        UU1=U1
        TT3=T3
        UU3=U3
        DO 321 K=1,4
          PP1(K)=P1(K)
          PP2(K)=P2(K)
321     CONTINUE
      ELSE
        TT=UX
        UU=TX
        TT1=U1
        UU1=T1
        TT3=U3
        UU3=T3
        DO 323 K=1,4
          PP1(K)=P2(K)
          PP2(K)=P1(K)
323     CONTINUE
      ENDIF
C
      TERM=WZSS(TT,UU,TT1,UU1,TT3,UU3,PP1,PP2)
      TERM=TERM-SGN*WZST(TT,UU,TT1,UU1,TT3,UU3,PP1,PP2)
      TERM=TERM-SGN*WZSU(TT,UU,TT1,UU1,TT3,UU3,PP1,PP2)
      TERM=TERM+WZTU(TT,UU,TT1,UU1,TT3,UU3,PP1,PP2)
      WWSIG=TERM*QSAVE(IQ1,1)*QSAVE(IQ2,2)*FJAC
C
      RETURN
C
C     Do Z+gamma or W+gamma 3-body subprocesses
C
2     CONTINUE
C
C          Z+gamma
C          Standard order is UP + UB --> Z0 + gamma
C
      IF(.NOT.(JETTYP(1).EQ.4.AND.JETTYP(2).EQ.1)) GO TO 505
      FJAC=S/SCM*P(1)/SQRT(P(1)**2+WMASS(4)**2)*UNITS
C
C          Select quarks
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
      A1=-AQ(IQ,4)
      B1=BQ(IQ,4)
      A2=-AQ(JQWW(1),4)
      B2=BQ(JQWW(1),4)
      DO K=1,5
        Q(K)=SNGL(P1WW(K))
        QB(K)=SNGL(P2WW(K))
        KK(K)=SNGL(P4WW(K))
        E(K)=SNGL(PZERO(K,1))
        EB(K)=SNGL(PZERO(K,2))
      END DO
C
      IF(INITYP(1).EQ.2*IQ) THEN
        SMS=SMSZG(Q,QB,KK,E,EB,A1,B1,A2,B2)
        TERM=ES**3*(EQ3(IQ)/3.)**2*SMS/192./PI**4/WMASS(4)/WGAM(4)/S**2
        WWSIG=TERM*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)*FJAC/2.
      ELSE
        SMS=SMSZG(QB,Q,KK,E,EB,A1,B1,A2,B2)
        TERM=ES**3*(EQ3(IQ)/3.)**2*SMS/192./PI**4/WMASS(4)/WGAM(4)/S**2
        WWSIG=TERM*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)*FJAC/2.
      ENDIF
505   IF(.NOT.(JETTYP(1).EQ.1.AND.JETTYP(2).EQ.4)) GO TO 509
      FJAC=S/SCM*P(2)/SQRT(P(2)**2+WMASS(4)**2)*UNITS
C
C          Select quarks
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
      A1=-AQ(IQ,4)
      B1=BQ(IQ,4)
      A2=-AQ(JQWW(2),4)
      B2=BQ(JQWW(2),4)
      DO K=1,5
        Q(K)=SNGL(P1WW(K))
        QB(K)=SNGL(P2WW(K))
        KK(K)=SNGL(P3WW(K))
        E(K)=SNGL(PZERO(K,1))
        EB(K)=SNGL(PZERO(K,2))
      END DO
C
      IF(INITYP(1).EQ.2*IQ) THEN
        SMS=SMSZG(Q,QB,KK,E,EB,A1,B1,A2,B2)
        TERM=ES**3*(EQ3(IQ)/3.)**2*SMS/192./PI**4/WMASS(4)/WGAM(4)/S**2
        WWSIG=TERM*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)*FJAC/2.
      ELSE
        SMS=SMSZG(QB,Q,KK,E,EB,A1,B1,A2,B2)
        TERM=ES**3*(EQ3(IQ)/3.)**2*SMS/192./PI**4/WMASS(4)/WGAM(4)/S**2
        WWSIG=TERM*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)*FJAC/2.
      ENDIF

C          W+- GM pair decays
C          Standard order is DN + UB --> W- + GM
C
C          Swap if W is jet 2
509   IF (ABS(IDJETS(1)).EQ.80.OR.ABS(IDJETS(2)).EQ.80) THEN
      IF(IDJETS(2).EQ.10) THEN
        DO 510 K=1,4
          P3(K)=P3WW(K)
          Q1(K)=PZERO(K,1)
510     CONTINUE
        AMASS3=PJETS(5,1)
        JW=1
        JG=2
        TX=T
        UX=U
      ELSE
        DO 520 K=1,4
          P3(K)=P4WW(K)
          Q1(K)=PZERO(K,1)
520     CONTINUE
        AMASS3=PJETS(5,2)
        JW=2
        JG=1
        TX=U
        UX=T
      ENDIF
      IF(IDJETS(JW).EQ.80) THEN
        IW=2
      ELSE
        IW=3
      ENDIF
C
      T1=-2.*(Q1(4)*P1(4)-Q1(1)*P1(1)-Q1(2)*P1(2)-Q1(3)*P1(3))
      U1=-2.*(Q1(4)*P2(4)-Q1(1)*P2(1)-Q1(2)*P2(2)-Q1(3)*P2(3))
C          Jacobean
      FJAC=S/SCM*UNITS
      FJAC=FJAC*P(JW)/SQRT(P(JW)**2+WM2)
C          Sum over quarks. Formulas are for DN UB --> W- GM.
C          Use symmetry for other cases.
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IQ=IQ1/2
      IF(2*IQ.EQ.IQ1) THEN
        LQK1=.TRUE.
      ELSE
        LQK1=.FALSE.
      ENDIF
C          Swap t and u as necessary
      IF((LQK1.AND.IW.EQ.3).OR.(.NOT.LQK1.AND.IW.EQ.2)) THEN
        TT=TX
        UU=UX
        TT1=T1
        UU1=U1
      ELSE
        TT=UX
        UU=TX
        TT1=U1
        UU1=T1
      ENDIF
C          Lepton or quark pointer
      IL=IABS(IDPAIR(1))
      IF(IL.GT.6) IL=IL-4
C
C          Matrix element - properly crossed variables.
C          Remember PZERO(K,1) is always the fermion.
      IF(LQK1) THEN
        P1DQ2=P1(4)*PZERO(4,2)-P1(1)*PZERO(1,2)-P1(2)*PZERO(2,2)
     $  -P1(3)*PZERO(3,2)
        P2DQ1=P2(4)*PZERO(4,1)-P2(1)*PZERO(1,1)-P2(2)*PZERO(2,1)
     $  -P2(3)*PZERO(3,1)
      ELSE
        P1DQ2=P2(4)*PZERO(4,2)-P2(1)*PZERO(1,2)-P2(2)*PZERO(2,2)
     $  -P2(3)*PZERO(3,2)
        P2DQ1=P1(4)*PZERO(4,1)-P1(1)*PZERO(1,1)-P1(2)*PZERO(2,1)
     $  -P1(3)*PZERO(3,1)
      ENDIF
      TERM=ALFA**2/(8.*SIN2W*S**2)*TBRWW(IW,JW)*RBRWW(IL,IW,JW)
     $*(-1./3.+UU/(TT+UU))**2/(TT*UU)*(4.*P2DQ1**2+4.*P1DQ2**2)
      WWSIG=TERM*QSAVE(IQ1,1)*QSAVE(IQ2,2)*FJAC
      END IF
C
      RETURN
      END
