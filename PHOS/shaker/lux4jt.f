*CMZ :          17/07/98  15.44.35  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUX4JT(NJET,CUT,KFL,ECM,KFLN,X1,X2,X4,X12,X14)

C...Purpose: to select the kinematical variables of four-jet events.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.
      DIMENSION WTA(4),WTB(4),WTC(4),WTD(4),WTE(4)

C...Common constants. Colour factors for QCD and Abelian gluon theory.
      PMQ=ULMASS(KFL)
      QME=(2.*PMQ/ECM)**2
      CT=LOG(1./CUT-5.)
      IF(MSTJ(109).EQ.0) THEN
        CF=4./3.
        CN=3.
        TR=2.5
      ELSE
        CF=1.
        CN=0.
        TR=15.
      ENDIF

C...Choice of process (qqbargg or qqbarqqbar).
  100 NJET=4
      IT=1
      IF(PARJ(155).GT.RLU(0)) IT=2
      IF(MSTJ(101).LE.-3) IT=-MSTJ(101)-2
      IF(IT.EQ.1) WTMX=0.7/CUT**2
      IF(IT.EQ.1.AND.MSTJ(109).EQ.2) WTMX=0.6/CUT**2
      IF(IT.EQ.2) WTMX=0.1125*CF*TR/CUT**2
      ID=1

C...Sample the five kinematical variables (for qqgg preweighted in y34).
  110 Y134=3.*CUT+(1.-6.*CUT)*RLU(0)
      Y234=3.*CUT+(1.-6.*CUT)*RLU(0)
      IF(IT.EQ.1) Y34=(1.-5.*CUT)*EXP(-CT*RLU(0))
      IF(IT.EQ.2) Y34=CUT+(1.-6.*CUT)*RLU(0)
      IF(Y34.LE.Y134+Y234-1..OR.Y34.GE.Y134*Y234) GOTO 110
      VT=RLU(0)
      CP=COS(PARU(1)*RLU(0))
      Y14=(Y134-Y34)*VT
      Y13=Y134-Y14-Y34
      VB=Y34*(1.-Y134-Y234+Y34)/((Y134-Y34)*(Y234-Y34))
      Y24=0.5*(Y234-Y34)*(1.-4.*SQRT(MAX(0.,VT*(1.-VT)*VB*(1.-VB)))*
     &CP-(1.-2.*VT)*(1.-2.*VB))
      Y23=Y234-Y34-Y24
      Y12=1.-Y134-Y23-Y24
      IF(MIN(Y12,Y13,Y14,Y23,Y24).LE.CUT) GOTO 110
      Y123=Y12+Y13+Y23
      Y124=Y12+Y14+Y24

C...Calculate matrix elements for qqgg or qqqq process.
      IC=0
      WTTOT=0.
  120 IC=IC+1
      IF(IT.EQ.1) THEN
        WTA(IC)=(Y12*Y34**2-Y13*Y24*Y34+Y14*Y23*Y34+3.*Y12*Y23*Y34+
     &  3.*Y12*Y14*Y34+4.*Y12**2*Y34-Y13*Y23*Y24+2.*Y12*Y23*Y24-
     &  Y13*Y14*Y24-2.*Y12*Y13*Y24+2.*Y12**2*Y24+Y14*Y23**2+2.*Y12*
     &  Y23**2+Y14**2*Y23+4.*Y12*Y14*Y23+4.*Y12**2*Y23+2.*Y12*Y14**2+
     &  2.*Y12*Y13*Y14+4.*Y12**2*Y14+2.*Y12**2*Y13+2.*Y12**3)/(2.*Y13*
     &  Y134*Y234*Y24)+(Y24*Y34+Y12*Y34+Y13*Y24-Y14*Y23+Y12*Y13)/(Y13*
     &  Y134**2)+2.*Y23*(1.-Y13)/(Y13*Y134*Y24)+Y34/(2.*Y13*Y24)
        WTB(IC)=(Y12*Y24*Y34+Y12*Y14*Y34-Y13*Y24**2+Y13*Y14*Y24+2.*Y12*
     &  Y14*Y24)/(Y13*Y134*Y23*Y14)+Y12*(1.+Y34)*Y124/(Y134*Y234*Y14*
     &  Y24)-(2.*Y13*Y24+Y14**2+Y13*Y23+2.*Y12*Y13)/(Y13*Y134*Y14)+
     &  Y12*Y123*Y124/(2.*Y13*Y14*Y23*Y24)
        WTC(IC)=-(5.*Y12*Y34**2+2.*Y12*Y24*Y34+2.*Y12*Y23*Y34+2.*Y12*
     &  Y14*Y34+2.*Y12*Y13*Y34+4.*Y12**2*Y34-Y13*Y24**2+Y14*Y23*Y24+
     &  Y13*Y23*Y24+Y13*Y14*Y24-Y12*Y14*Y24-Y13**2*Y24-3.*Y12*Y13*Y24-
     &  Y14*Y23**2-Y14**2*Y23+Y13*Y14*Y23-3.*Y12*Y14*Y23-Y12*Y13*Y23)/
     &  (4.*Y134*Y234*Y34**2)+(3.*Y12*Y34**2-3.*Y13*Y24*Y34+3.*Y12*Y24*
     &  Y34+3.*Y14*Y23*Y34-Y13*Y24**2-Y12*Y23*Y34+6.*Y12*Y14*Y34+2.*Y12*
     &  Y13*Y34-2.*Y12**2*Y34+Y14*Y23*Y24-3.*Y13*Y23*Y24-2.*Y13*Y14*
     &  Y24+4.*Y12*Y14*Y24+2.*Y12*Y13*Y24+3.*Y14*Y23**2+2.*Y14**2*Y23+
     &  2.*Y14**2*Y12+2.*Y12**2*Y14+6.*Y12*Y14*Y23-2.*Y12*Y13**2-
     &  2.*Y12**2*Y13)/(4.*Y13*Y134*Y234*Y34)
        WTC(IC)=WTC(IC)+(2.*Y12*Y34**2-2.*Y13*Y24*Y34+Y12*Y24*Y34+
     &  4.*Y13*Y23*Y34+4.*Y12*Y14*Y34+2.*Y12*Y13*Y34+2.*Y12**2*Y34-
     &  Y13*Y24**2+3.*Y14*Y23*Y24+4.*Y13*Y23*Y24-2.*Y13*Y14*Y24+
     &  4.*Y12*Y14*Y24+2.*Y12*Y13*Y24+2.*Y14*Y23**2+4.*Y13*Y23**2+
     &  2.*Y13*Y14*Y23+2.*Y12*Y14*Y23+4.*Y12*Y13*Y23+2.*Y12*Y14**2+4.*
     &  Y12**2*Y13+4.*Y12*Y13*Y14+2.*Y12**2*Y14)/(4.*Y13*Y134*Y24*Y34)-
     &  (Y12*Y34**2-2.*Y14*Y24*Y34-2.*Y13*Y24*Y34-Y14*Y23*Y34+Y13*Y23*
     &  Y34+Y12*Y14*Y34+2.*Y12*Y13*Y34-2.*Y14**2*Y24-4.*Y13*Y14*Y24-
     &  4.*Y13**2*Y24-Y14**2*Y23-Y13**2*Y23+Y12*Y13*Y14-Y12*Y13**2)/
     &  (2.*Y13*Y34*Y134**2)+(Y12*Y34**2-4.*Y14*Y24*Y34-2.*Y13*Y24*Y34-
     &  2.*Y14*Y23*Y34-4.*Y13*Y23*Y34-4.*Y12*Y14*Y34-4.*Y12*Y13*Y34-
     &  2.*Y13*Y14*Y24+2.*Y13**2*Y24+2.*Y14**2*Y23-2.*Y13*Y14*Y23-
     &  Y12*Y14**2-6.*Y12*Y13*Y14-Y12*Y13**2)/(4.*Y34**2*Y134**2)
        WTTOT=WTTOT+Y34*CF*(CF*WTA(IC)+(CF-0.5*CN)*WTB(IC)+CN*WTC(IC))/
     &  8.
      ELSE
        WTD(IC)=(Y13*Y23*Y34+Y12*Y23*Y34-Y12**2*Y34+Y13*Y23*Y24+2.*Y12*
     &  Y23*Y24-Y14*Y23**2+Y12*Y13*Y24+Y12*Y14*Y23+Y12*Y13*Y14)/(Y13**2*
     &  Y123**2)-(Y12*Y34**2-Y13*Y24*Y34+Y12*Y24*Y34-Y14*Y23*Y34-Y12*
     &  Y23*Y34-Y13*Y24**2+Y14*Y23*Y24-Y13*Y23*Y24-Y13**2*Y24+Y14*
     &  Y23**2)/(Y13**2*Y123*Y134)+(Y13*Y14*Y12+Y34*Y14*Y12-Y34**2*Y12+
     &  Y13*Y14*Y24+2.*Y34*Y14*Y24-Y23*Y14**2+Y34*Y13*Y24+Y34*Y23*Y14+
     &  Y34*Y13*Y23)/(Y13**2*Y134**2)-(Y34*Y12**2-Y13*Y24*Y12+Y34*Y24*
     &  Y12-Y23*Y14*Y12-Y34*Y14*Y12-Y13*Y24**2+Y23*Y14*Y24-Y13*Y14*Y24-
     &  Y13**2*Y24+Y23*Y14**2)/(Y13**2*Y134*Y123)
        WTE(IC)=(Y12*Y34*(Y23-Y24+Y14+Y13)+Y13*Y24**2-Y14*Y23*Y24+Y13*
     &  Y23*Y24+Y13*Y14*Y24+Y13**2*Y24-Y14*Y23*(Y14+Y23+Y13))/(Y13*Y23*
     &  Y123*Y134)-Y12*(Y12*Y34-Y23*Y24-Y13*Y24-Y14*Y23-Y14*Y13)/(Y13*
     &  Y23*Y123**2)-(Y14+Y13)*(Y24+Y23)*Y34/(Y13*Y23*Y134*Y234)+
     &  (Y12*Y34*(Y14-Y24+Y23+Y13)+Y13*Y24**2-Y23*Y14*Y24+Y13*Y14*Y24+
     &  Y13*Y23*Y24+Y13**2*Y24-Y23*Y14*(Y14+Y23+Y13))/(Y13*Y14*Y134*
     &  Y123)-Y34*(Y34*Y12-Y14*Y24-Y13*Y24-Y23*Y14-Y23*Y13)/(Y13*Y14*
     &  Y134**2)-(Y23+Y13)*(Y24+Y14)*Y12/(Y13*Y14*Y123*Y124)
        WTTOT=WTTOT+CF*(TR*WTD(IC)+(CF-0.5*CN)*WTE(IC))/16.
      ENDIF

C...Permutations of momenta in matrix element. Weighting.
  130 IF(IC.EQ.1.OR.IC.EQ.3.OR.ID.EQ.2.OR.ID.EQ.3) THEN
        YSAV=Y13
        Y13=Y14
        Y14=YSAV
        YSAV=Y23
        Y23=Y24
        Y24=YSAV
        YSAV=Y123
        Y123=Y124
        Y124=YSAV
      ENDIF
      IF(IC.EQ.2.OR.IC.EQ.4.OR.ID.EQ.3.OR.ID.EQ.4) THEN
        YSAV=Y13
        Y13=Y23
        Y23=YSAV
        YSAV=Y14
        Y14=Y24
        Y24=YSAV
        YSAV=Y134
        Y134=Y234
        Y234=YSAV
      ENDIF
      IF(IC.LE.3) GOTO 120
      IF(ID.EQ.1.AND.WTTOT.LT.RLU(0)*WTMX) GOTO 110
      IC=5

C...qqgg events: string configuration and event type.
      IF(IT.EQ.1) THEN
        IF(MSTJ(109).EQ.0.AND.ID.EQ.1) THEN
          PARJ(156)=Y34*(2.*(WTA(1)+WTA(2)+WTA(3)+WTA(4))+4.*(WTC(1)+
     &    WTC(2)+WTC(3)+WTC(4)))/(9.*WTTOT)
          IF(WTA(2)+WTA(4)+2.*(WTC(2)+WTC(4)).GT.RLU(0)*(WTA(1)+WTA(2)+
     &    WTA(3)+WTA(4)+2.*(WTC(1)+WTC(2)+WTC(3)+WTC(4)))) ID=2
          IF(ID.EQ.2) GOTO 130
        ELSEIF(MSTJ(109).EQ.2.AND.ID.EQ.1) THEN
          PARJ(156)=Y34*(WTA(1)+WTA(2)+WTA(3)+WTA(4))/(8.*WTTOT)
          IF(WTA(2)+WTA(4).GT.RLU(0)*(WTA(1)+WTA(2)+WTA(3)+WTA(4))) ID=2
          IF(ID.EQ.2) GOTO 130
        ENDIF
        MSTJ(120)=3
        IF(MSTJ(109).EQ.0.AND.0.5*Y34*(WTC(1)+WTC(2)+WTC(3)+WTC(4)).GT.
     &  RLU(0)*WTTOT) MSTJ(120)=4
        KFLN=21

C...Mass cuts. Kinematical variables out.
        IF(Y12.LE.CUT+QME) NJET=2
        IF(NJET.EQ.2) GOTO 150
        Q12=0.5*(1.-SQRT(1.-QME/Y12))
        X1=1.-(1.-Q12)*Y234-Q12*Y134
        X4=1.-(1.-Q12)*Y134-Q12*Y234
        X2=1.-Y124
        X12=(1.-Q12)*Y13+Q12*Y23
        X14=Y12-0.5*QME
        IF(Y134*Y234/((1.-X1)*(1.-X4)).LE.RLU(0)) NJET=2

C...qqbarqqbar events: string configuration, choose new flavour.
      ELSE
        IF(ID.EQ.1) THEN
          WTR=RLU(0)*(WTD(1)+WTD(2)+WTD(3)+WTD(4))
          IF(WTR.LT.WTD(2)+WTD(3)+WTD(4)) ID=2
          IF(WTR.LT.WTD(3)+WTD(4)) ID=3
          IF(WTR.LT.WTD(4)) ID=4
          IF(ID.GE.2) GOTO 130
        ENDIF
        MSTJ(120)=5
        PARJ(156)=CF*TR*(WTD(1)+WTD(2)+WTD(3)+WTD(4))/(16.*WTTOT)
  140   KFLN=1+INT(5.*RLU(0))
        IF(KFLN.NE.KFL.AND.0.2*PARJ(156).LE.RLU(0)) GOTO 140
        IF(KFLN.EQ.KFL.AND.1.-0.8*PARJ(156).LE.RLU(0)) GOTO 140
        IF(KFLN.GT.MSTJ(104)) NJET=2
        PMQN=ULMASS(KFLN)
        QMEN=(2.*PMQN/ECM)**2

C...Mass cuts. Kinematical variables out.
        IF(Y24.LE.CUT+QME.OR.Y13.LE.1.1*QMEN) NJET=2
        IF(NJET.EQ.2) GOTO 150
        Q24=0.5*(1.-SQRT(1.-QME/Y24))
        Q13=0.5*(1.-SQRT(1.-QMEN/Y13))
        X1=1.-(1.-Q24)*Y123-Q24*Y134
        X4=1.-(1.-Q24)*Y134-Q24*Y123
        X2=1.-(1.-Q13)*Y234-Q13*Y124
        X12=(1.-Q24)*((1.-Q13)*Y14+Q13*Y34)+Q24*((1.-Q13)*Y12+Q13*Y23)
        X14=Y24-0.5*QME
        X34=(1.-Q24)*((1.-Q13)*Y23+Q13*Y12)+Q24*((1.-Q13)*Y34+Q13*Y14)
        IF(PMQ**2+PMQN**2+MIN(X12,X34)*ECM**2.LE.
     &  (PARJ(127)+PMQ+PMQN)**2) NJET=2
        IF(Y123*Y134/((1.-X1)*(1.-X4)).LE.RLU(0)) NJET=2
      ENDIF
  150 IF(MSTJ(101).LE.-2.AND.NJET.EQ.2) GOTO 100

      RETURN
      END
