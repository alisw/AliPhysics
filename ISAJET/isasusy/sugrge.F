#include "isajet/pilot.h"
      SUBROUTINE SUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2
     $,NSTEP,IMODEL)
C
C          Make one complete iteration of the renormalization group
C          equations from MZ to MGUT and back, setting the boundary
C          conditions on each end.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sugpas.inc"
#include "isajet/sugnu.inc"
#include "isajet/sugxin.inc"
#include "isajet/sugmg.inc"
C
      EXTERNAL SURG26
      DOUBLE PRECISION DDILOG,XLM
      REAL M0,MHF,A0,TANB,SGNMU,MT,G(29),G0(29),W2(87)
      INTEGER IG(29),NSTEP,IMODEL
      REAL PI,TZ,A1I,A2I,A3I,GGUT,AGUTI,SIG1,SIG2,
     $MH1S,MH2S,MUS,T,MZ,TGUT,DT,AGUT,Q,ASMT,MTMT,SINB,
     $BETA,QOLD,XLAMGM,XMESGM,XN5GM,XC,G3GUT,THRF,THRG,DY,
     $BLHAT,BBHAT,BTHAT
      INTEGER I,II
      DATA MZ/91.187/
C
C          Re-initialize weak scale parameters
C
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      PI=4.*ATAN(1.)
      BETA=ATAN(XTANB)
      SINB=SIN(BETA)
      ASMZ=0.118
C      ASMT=G3MT**2/4./PI
C      MTMT=MT/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/MT))*(ASMT/PI)**2)
C      FTMT=MTMT/SINB/VEV
      G(1)=SQRT(4*PI*A1MZ)
      G(2)=SQRT(4*PI*A2MZ)
      G(3)=SQRT(4*PI*ASMZ)
      G(4)=FTAMZ
      G(5)=FBMZ
      G(6)=G(6)
      G(25)=MU
      G(26)=B
      G(27)=0.
      G(28)=0.
      G(29)=0.
C          Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF
C
C          Run back up to mgut with approximate susy spectra
C
      IF (IMODEL.EQ.1) THEN
        IF (XSUGIN(7).EQ.0.) THEN 
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=LOG(MZ/MGUT)
      TGUT=0.
      DT=(TGUT-TZ)/FLOAT(NSTEP)
      DO 250 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/FLOAT(NSTEP)
        QOLD=Q
        Q=MGUT*EXP(T)
        IF (QOLD.LE.MT.AND.Q.GT.MT) G(6)=FTMT
        IF (QOLD.LE.XNRIN(2).AND.Q.GT.XNRIN(2)) THEN
          G(27)=FNMZ
          G(28)=G0(28)
          G(29)=G0(29)
        END IF
        CALL RKSTP(29,DT,T,G,SURG26,W2)
        A1I=4*PI/G(1)**2
        A2I=4*PI/G(2)**2
        A3I=4*PI/G(3)**2
        IF (G(5).GT.10..OR.G(6).GT.10..OR.G(27).GT.10.) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 30
250   CONTINUE
      IF (IMODEL.EQ.1.AND.XSUGIN(7).EQ.0.) THEN
        WRITE(LOUT,*) 'SUGRGE ERROR: NO UNIFICATION FOUND'
        NOGOOD=1
        GO TO 100
      END IF
30    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=Q
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=(G(1)**2/4./PI+G(2)**2/4./PI)/2.
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=G(4)
      FBGUT=G(5)
      FTGUT=G(6)
      IF (XNRIN(2).LT.1.E19.AND.XNRIN(1).EQ.0.) THEN
C     IMPOSE FN-FT UNIFICATION
        FNGUT=G(6)
      ELSE
        FNGUT=G(27)
      END IF
      G3GUT=G(3)
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C
C          Set GUT boundary condition
C
      DO 260 I=1,3
        IF (IMODEL.EQ.1) THEN
          G(I)=G(I)
          G(I+6)=MHF
          G(I+9)=A0
        ELSE IF (IMODEL.EQ.2) THEN
          G(I)=G(I)
          G(I+6)=XGMIN(11+I)*XGMIN(8)*THRG*(G(I)/4./PI)**2*XLAMGM
          G(I+9)=0.
        END IF
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=FNGUT
        G(28)=XNRIN(4)**2
        G(29)=XNRIN(3)
      ELSE
        G(27)=0.
        G(28)=0.
        G(29)=0.
      END IF
260   CONTINUE
C     OVERWRITE ALFA_3 UNIFICATION TO GET ALFA_3(MZ) RIGHT
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=GGUT
      IF (IMODEL.EQ.1) THEN
        DO 270 I=13,24
          G(I)=M0**2
270     CONTINUE
C          Set possible non-universal GUT scale boundary conditions
      DO 280 I=1,6
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=XNUSUG(I)
        END IF
280   CONTINUE
      DO 281 I=7,18
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=XNUSUG(I)**2
        END IF
281   CONTINUE
      ELSE IF (IMODEL.EQ.2) THEN
       XC=2*THRF*XLAMGM**2
       DY=SQRT(3./5.)*G(1)*XGMIN(11)
       G(13)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     , XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(9)-DY
       G(14)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     , XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(10)+DY
       G(15)=XC*(.6*XGMIN(12)*(G(1)/4./PI)**4)+2*DY
       G(16)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     , XGMIN(12)*(G(1)/4./PI)**4)-DY
       G(17)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*XGMIN(12)*
     , (G(1)/4./PI)**4/9.)+2*DY/3.
       G(18)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*4*XGMIN(12)*
     , (G(1)/4./PI)**4/9.)-4*DY/3.
       G(19)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.75*XGMIN(13)*
     , (G(2)/4./PI)**4+.6*XGMIN(12)*(G(1)/4./PI)**4/36.)+DY/3.
       G(20)=G(15)
       G(21)=G(16)
       G(22)=G(17)
       G(23)=G(18)
       G(24)=G(19)
      ELSE IF (IMODEL.EQ.7) THEN
       G(1)=G(1)
       G(2)=G(2)
       G(3)=G(3)
       BLHAT=G(4)*(-9*G(1)**2/5.-3*G(2)**2+3*G(5)**2+4*G(4)**2)
       BBHAT=G(5)*(-7*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
       BTHAT=G(6)*(-13*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             6*G(6)**2+G(5)**2)
       G(7)=-33*MHF*G(1)**2/5./16./PI**2
       G(8)=-MHF*G(2)**2/16./PI**2
       G(9)=3*MHF*G(3)**2/16./PI**2
       G(10)=BLHAT*MHF/G(4)/16./PI**2
       G(11)=BBHAT*MHF/G(5)/16./PI**2
       G(12)=BTHAT*MHF/G(6)/16./PI**2
       G(13)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(5)*BBHAT+G(4)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(14)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(6)*BTHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(15)=(-198*G(1)**4/25.)*MHF**2/(16*PI**2)**2
       G(16)=(-99*G(1)**4/50.-3*G(2)**4/2.)*MHF**2/(16*PI**2)**2
       G(17)=(-22*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2
       G(18)=(-88*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2
       G(19)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4)*
     ,        MHF**2/(16*PI**2)**2
       G(20)=(-198*G(1)**4/25.+2*G(4)*BLHAT)*MHF**2/(16*PI**2)**2
       G(21)=(-99*G(1)**4/50.-3*G(2)**4/2.+G(4)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2
       G(22)=(-22*G(1)**4/25.+8*G(3)**4+2*G(5)*BBHAT)*
     , MHF**2/(16*PI**2)**2
       G(23)=(-88*G(1)**4/25.+8*G(3)**4+2*G(6)*BTHAT)*
     , MHF**2/(16*PI**2)**2
       G(24)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*PI**2)**2
       DO 284 I=13,24
284      G(I)=G(I)+M0**2
      END IF
      DO 285 I=1,29
        IG(I)=0
285   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0..OR.G(16).LT.0.) THEN
        ITACHY=2
      ELSE
        ITACHY=0
      END IF
C
C          Run back down to weak scale
C
      TZ=LOG(MZ/MGUT)
      TGUT=0.
      DT=(TZ-TGUT)/FLOAT(NSTEP)
      DO 290 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/FLOAT(NSTEP)
        QOLD=Q
        Q=MGUT*EXP(T)
        CALL RKSTP(29,DT,T,G,SURG26,W2)
        CALL SUGFRZ(Q,G,G0,IG)
        IF (QOLD.GE.AMNRMJ.AND.Q.LT.AMNRMJ.AND.XNRIN(1).EQ.0.) THEN
          FNMZ=G(27)
        END IF
        IF (Q.LT.AMNRMJ) THEN
          G(27)=0.
          G(28)=0.
          G(29)=0.
        END IF
        IF (NOGOOD.NE.0) GO TO 100
        IF (Q.LT.MZ) GO TO 40
290   CONTINUE
40    CONTINUE
C
C          Electroweak breaking constraints; tree level
C
      MUS=(G0(13)-G0(14)*TANB**2)/(TANB**2-1.)-MZ**2/2.
      IF (MUS.LT.0.) THEN
        NOGOOD=2
        GO TO 100
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2B/MU/2.
      CALL SUGMAS(G0,0,IMODEL)
      IF (NOGOOD.NE.0) GO TO 100
C
C           Electroweak breaking constraints; loop level
C
      CALL SUGEFF(G0,SIG1,SIG2)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANB**2)/(TANB**2-1.)-MZ**2/2.
      IF (MUS.LT.0.) THEN
        NOGOOD=2
        GO TO 100
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2B/MU/2.
      CALL SUGMAS(G0,1,IMODEL)
C
100   RETURN
      END
