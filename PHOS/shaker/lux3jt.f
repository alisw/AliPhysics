*CMZ :          17/07/98  15.44.35  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUX3JT(NJET,CUT,KFL,ECM,X1,X2)

C...Purpose: to select the kinematical variables of three-jet events.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.
      DIMENSION ZHUP(5,12)

C...Coefficients of Zhu second order parametrization.
      DATA ((ZHUP(IC1,IC2),IC2=1,12),IC1=1,5)/
     &    18.29,    89.56,    4.541,   -52.09,   -109.8,    24.90,
     &    11.63,    3.683,    17.50, 0.002440,   -1.362,  -0.3537,
     &    11.42,    6.299,   -22.55,   -8.915,    59.25,   -5.855,
     &   -32.85,   -1.054,   -16.90, 0.006489,  -0.8156,  0.01095,
     &    7.847,   -3.964,   -35.83,    1.178,    29.39,   0.2806,
     &    47.82,   -12.36,   -56.72,  0.04054,  -0.4365,   0.6062,
     &    5.441,   -56.89,   -50.27,    15.13,    114.3,   -18.19,
     &    97.05,   -1.890,   -139.9,  0.08153,  -0.4984,   0.9439,
     &   -17.65,    51.44,   -58.32,    70.95,   -255.7,   -78.99,
     &    476.9,    29.65,   -239.3,   0.4745,   -1.174,    6.081/

C...Dilogarithm of x for x<0.5 (x>0.5 obtained by analytic trick).
      DILOG(X)=X+X**2/4.+X**3/9.+X**4/16.+X**5/25.+X**6/36.+X**7/49.

C...Event type. Mass effect factors and other common constants.
      MSTJ(120)=2
      MSTJ(121)=0
      PMQ=ULMASS(KFL)
      QME=(2.*PMQ/ECM)**2
      IF(MSTJ(109).NE.1) THEN
        CUTL=LOG(CUT)
        CUTD=LOG(1./CUT-2.)
        IF(MSTJ(109).EQ.0) THEN
          CF=4./3.
          CN=3.
          TR=2.
          WTMX=MIN(20.,37.-6.*CUTD)
          IF(MSTJ(110).EQ.2) WTMX=2.*(7.5+80.*CUT)
        ELSE
          CF=1.
          CN=0.
          TR=12.
          WTMX=0.
        ENDIF

C...Alpha_strong and effects of optimized Q^2 scale. Maximum weight.
        ALS2PI=PARU(118)/PARU(2)
        WTOPT=0.
        IF(MSTJ(111).EQ.1) WTOPT=(33.-2.*MSTU(112))/6.*LOG(PARJ(169))*
     &  ALS2PI
        WTMAX=MAX(0.,1.+WTOPT+ALS2PI*WTMX)

C...Choose three-jet events in allowed region.
  100   NJET=3
  110   Y13L=CUTL+CUTD*RLU(0)
        Y23L=CUTL+CUTD*RLU(0)
        Y13=EXP(Y13L)
        Y23=EXP(Y23L)
        Y12=1.-Y13-Y23
        IF(Y12.LE.CUT) GOTO 110
        IF(Y13**2+Y23**2+2.*Y12.LE.2.*RLU(0)) GOTO 110

C...Second order corrections.
        IF(MSTJ(101).EQ.2.AND.MSTJ(110).LE.1) THEN
          Y12L=LOG(Y12)
          Y13M=LOG(1.-Y13)
          Y23M=LOG(1.-Y23)
          Y12M=LOG(1.-Y12)
          IF(Y13.LE.0.5) Y13I=DILOG(Y13)
          IF(Y13.GE.0.5) Y13I=1.644934-Y13L*Y13M-DILOG(1.-Y13)
          IF(Y23.LE.0.5) Y23I=DILOG(Y23)
          IF(Y23.GE.0.5) Y23I=1.644934-Y23L*Y23M-DILOG(1.-Y23)
          IF(Y12.LE.0.5) Y12I=DILOG(Y12)
          IF(Y12.GE.0.5) Y12I=1.644934-Y12L*Y12M-DILOG(1.-Y12)
          WT1=(Y13**2+Y23**2+2.*Y12)/(Y13*Y23)
          WT2=CF*(-2.*(CUTL-Y12L)**2-3.*CUTL-1.+3.289868+
     &    2.*(2.*CUTL-Y12L)*CUT/Y12)+
     &    CN*((CUTL-Y12L)**2-(CUTL-Y13L)**2-(CUTL-Y23L)**2-11.*CUTL/6.+
     &    67./18.+1.644934-(2.*CUTL-Y12L)*CUT/Y12+(2.*CUTL-Y13L)*
     &    CUT/Y13+(2.*CUTL-Y23L)*CUT/Y23)+
     &    TR*(2.*CUTL/3.-10./9.)+
     &    CF*(Y12/(Y12+Y13)+Y12/(Y12+Y23)+(Y12+Y23)/Y13+(Y12+Y13)/Y23+
     &    Y13L*(4.*Y12**2+2.*Y12*Y13+4.*Y12*Y23+Y13*Y23)/(Y12+Y23)**2+
     &    Y23L*(4.*Y12**2+2.*Y12*Y23+4.*Y12*Y13+Y13*Y23)/(Y12+Y13)**2)/
     &    WT1+
     &    CN*(Y13L*Y13/(Y12+Y23)+Y23L*Y23/(Y12+Y13))/WT1+
     &    (CN-2.*CF)*((Y12**2+(Y12+Y13)**2)*(Y12L*Y23L-Y12L*Y12M-Y23L*
     &    Y23M+1.644934-Y12I-Y23I)/(Y13*Y23)+(Y12**2+(Y12+Y23)**2)*
     &    (Y12L*Y13L-Y12L*Y12M-Y13L*Y13M+1.644934-Y12I-Y13I)/
     &    (Y13*Y23)+(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))-
     &    2.*Y12L*Y12**2/(Y13+Y23)**2-4.*Y12L*Y12/(Y13+Y23))/WT1-
     &    CN*(Y13L*Y23L-Y13L*Y13M-Y23L*Y23M+1.644934-Y13I-Y23I)
          IF(1.+WTOPT+ALS2PI*WT2.LE.0.) MSTJ(121)=1
          IF(1.+WTOPT+ALS2PI*WT2.LE.WTMAX*RLU(0)) GOTO 110
          PARJ(156)=(WTOPT+ALS2PI*WT2)/(1.+WTOPT+ALS2PI*WT2)

        ELSEIF(MSTJ(101).EQ.2.AND.MSTJ(110).EQ.2) THEN
C...Second order corrections; Zhu parametrization of ERT.
          ZX=(Y23-Y13)**2
          ZY=1.-Y12
          IZA=0
          DO 120 IY=1,5
  120     IF(ABS(CUT-0.01*IY).LT.0.0001) IZA=IY
          IF(IZA.NE.0) THEN
            IZ=IZA
            WT2=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY
          ELSE
            IZ=100.*CUT
            WTL=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY
            IZ=IZ+1
            WTU=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY
            WT2=WTL+(WTU-WTL)*(100.*CUT+1.-IZ)
          ENDIF
          IF(1.+WTOPT+2.*ALS2PI*WT2.LE.0.) MSTJ(121)=1
          IF(1.+WTOPT+2.*ALS2PI*WT2.LE.WTMAX*RLU(0)) GOTO 110
          PARJ(156)=(WTOPT+2.*ALS2PI*WT2)/(1.+WTOPT+2.*ALS2PI*WT2)
        ENDIF

C...Impose mass cuts (gives two jets). For fixed jet number new try.
        X1=1.-Y23
        X2=1.-Y13
        X3=1.-Y12
        IF(4.*Y23*Y13*Y12/X3**2.LE.QME) NJET=2
        IF(MOD(MSTJ(103),4).GE.2.AND.IABS(MSTJ(101)).LE.1.AND.QME*X3+
     &  0.5*QME**2+(0.5*QME+0.25*QME**2)*((1.-X2)/(1.-X1)+
     &  (1.-X1)/(1.-X2)).GT.(X1**2+X2**2)*RLU(0)) NJET=2
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 100

C...Scalar gluon model (first order only, no mass effects).
      ELSE
  130   NJET=3
  140   Y12=SQRT(4.*CUT**2+RLU(0)*((1.-CUT)**2-4.*CUT**2))
        IF(LOG((Y12-CUT)/CUT).LE.RLU(0)*LOG((1.-2.*CUT)/CUT)) GOTO 140
        YD=SIGN(2.*CUT*((Y12-CUT)/CUT)**RLU(0)-Y12,RLU(0)-0.5)
        X1=1.-0.5*(Y12+YD)
        X2=1.-0.5*(Y12-YD)
        IF(4.*(1.-X1)*(1.-X2)*Y12/(1.-Y12)**2.LE.QME) NJET=2
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 130
      ENDIF

      RETURN
      END
