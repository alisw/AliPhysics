 
C*********************************************************************
 
      SUBROUTINE PYSTFL(KF,X,Q2,XPQ)
 
C...Give proton structure function at small x and/or Q^2 according to
C...correct limiting behaviour.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /PYPARS/,/PYINT1/
      DIMENSION XPQ(-25:25),XPA(-25:25),XPB(-25:25),WTSB(-3:3)
      DATA RMR/0.92/,RMP/0.38/,WTSB/0.5,1.,1.,5.,1.,1.,0.5/
 
C...Send everything but protons/neutrons/VMD pions directly to PYSTFU.
      MINT(92)=0
      KFA=IABS(KF)
      IACC=0
      IF((KFA.EQ.2212.OR.KFA.EQ.2112).AND.MSTP(57).GE.2) IACC=1
      IF(KFA.EQ.211.AND.MSTP(57).GE.3) IACC=1
      IF(KFA.EQ.22.AND.MINT(109).EQ.2.AND.MSTP(57).GE.3) IACC=1
      IF(IACC.EQ.0) THEN
        CALL PYSTFU(KF,X,Q2,XPQ)
        RETURN
      ENDIF
 
C...Reset. Check x.
      DO 100 KFL=-25,25
      XPQ(KFL)=0.
  100 CONTINUE
      IF(X.LE.0..OR.X.GE.1.) THEN
        WRITE(MSTU(11),5000) X
        RETURN
      ENDIF
 
C...Define valence content.
      KFC=KF
      NV1=2
      NV2=1
      IF(KF.EQ.2212) THEN
        KFV1=2
        KFV2=1
      ELSEIF(KF.EQ.-2212) THEN
        KFV1=-2
        KFV2=-1
      ELSEIF(KF.EQ.2112) THEN
        KFV1=1
        KFV2=2
      ELSEIF(KF.EQ.-2112) THEN
        KFV1=-1
        KFV2=-2
      ELSEIF(KF.EQ.211) THEN
        NV1=1
        KFV1=2
        KFV2=-1
      ELSEIF(KF.EQ.-211) THEN
        NV1=1
        KFV1=-2
        KFV2=1
      ELSEIF(MINT(105).LE.223) THEN
        KFV1=1
        WTV1=0.2
        KFV2=2
        WTV2=0.8
      ELSEIF(MINT(105).EQ.333) THEN
        KFV1=3
        WTV1=1.0
        KFV2=1
        WTV2=0.0
      ELSEIF(MINT(105).EQ.443) THEN
        KFV1=4
        WTV1=1.0
        KFV2=1
        WTV2=0.0
      ENDIF
 
C...Do naive evaluation and find min Q^2, boundary Q^2 and x_0.
      CALL PYSTFU(KFC,X,Q2,XPA)
      Q2MN=MAX(3.,VINT(231))
      Q2B=2.+0.052**2*EXP(3.56*SQRT(MAX(0.,-LOG(3.*X))))
      XMN=EXP(-(LOG((Q2MN-2.)/0.052**2)/3.56)**2)/3.
 
C...Large Q2 and large x: naive call is enough.
      IF(Q2.GT.Q2MN.AND.Q2.GT.Q2B) THEN
        DO 110 KFL=-25,25
        XPQ(KFL)=XPA(KFL)
  110   CONTINUE
        MINT(92)=1
 
C...Small Q2 and large x: dampen boundary value.
      ELSEIF(X.GT.XMN) THEN
 
C...Evaluate at boundary and define dampening factors.
        CALL PYSTFU(KFC,X,Q2MN,XPA)
        FV=(Q2*(Q2MN+RMR)/(Q2MN*(Q2+RMR)))**(0.55*(1.-X)/(1.-XMN))
        FS=(Q2*(Q2MN+RMP)/(Q2MN*(Q2+RMP)))**1.08
 
C...Separate valence and sea parts of structure function.
        IF(KFA.NE.22) THEN
          XFV1=XPA(KFV1)-XPA(-KFV1)
          XPA(KFV1)=XPA(-KFV1)
          XFV2=XPA(KFV2)-XPA(-KFV2)
          XPA(KFV2)=XPA(-KFV2)
        ELSE
          XPA(KFV1)=XPA(KFV1)-WTV1*VINT(232)
          XPA(-KFV1)=XPA(-KFV1)-WTV1*VINT(232)
          XPA(KFV2)=XPA(KFV2)-WTV2*VINT(232)
          XPA(-KFV2)=XPA(-KFV2)-WTV2*VINT(232)
        ENDIF
 
C...Dampen valence and sea separately. Put back together.
        DO 120 KFL=-25,25
        XPQ(KFL)=FS*XPA(KFL)
  120   CONTINUE
        IF(KFA.NE.22) THEN
          XPQ(KFV1)=XPQ(KFV1)+FV*XFV1
          XPQ(KFV2)=XPQ(KFV2)+FV*XFV2
        ELSE
          XPQ(KFV1)=XPQ(KFV1)+FV*WTV1*VINT(232)
          XPQ(-KFV1)=XPQ(-KFV1)+FV*WTV1*VINT(232)
          XPQ(KFV2)=XPQ(KFV2)+FV*WTV2*VINT(232)
          XPQ(-KFV2)=XPQ(-KFV2)+FV*WTV2*VINT(232)          
        ENDIF
        MINT(92)=2
 
C...Large Q2 and small x: interpolate behaviour.
      ELSEIF(Q2.GT.Q2MN) THEN
 
C...Evaluate at extremes and define coefficients for interpolation.
        CALL PYSTFU(KFC,XMN,Q2MN,XPA)
        VI232A=VINT(232)
        CALL PYSTFU(KFC,X,Q2B,XPB)
        VI232B=VINT(232) 
        FLA=LOG(Q2B/Q2)/LOG(Q2B/Q2MN)
        FVA=(X/XMN)**0.45*FLA
        FSA=(X/XMN)**(-0.08)*FLA
        FB=1.-FLA
 
C...Separate valence and sea parts of structure function.
        IF(KFA.NE.22) THEN
          XFVA1=XPA(KFV1)-XPA(-KFV1)
          XPA(KFV1)=XPA(-KFV1)
          XFVA2=XPA(KFV2)-XPA(-KFV2)
          XPA(KFV2)=XPA(-KFV2)
          XFVB1=XPB(KFV1)-XPB(-KFV1)
          XPB(KFV1)=XPB(-KFV1)
          XFVB2=XPB(KFV2)-XPB(-KFV2)
          XPB(KFV2)=XPB(-KFV2)
        ELSE
          XPA(KFV1)=XPA(KFV1)-WTV1*VI232A
          XPA(-KFV1)=XPA(-KFV1)-WTV1*VI232A
          XPA(KFV2)=XPA(KFV2)-WTV2*VI232A
          XPA(-KFV2)=XPA(-KFV2)-WTV2*VI232A
          XPB(KFV1)=XPB(KFV1)-WTV1*VI232B
          XPB(-KFV1)=XPB(-KFV1)-WTV1*VI232B
          XPB(KFV2)=XPB(KFV2)-WTV2*VI232B
          XPB(-KFV2)=XPB(-KFV2)-WTV2*VI232B
        ENDIF
 
C...Interpolate for valence and sea. Put back together.
        DO 130 KFL=-25,25
        XPQ(KFL)=FSA*XPA(KFL)+FB*XPB(KFL)
  130   CONTINUE
        IF(KFA.NE.22) THEN
          XPQ(KFV1)=XPQ(KFV1)+(FVA*XFVA1+FB*XFVB1)
          XPQ(KFV2)=XPQ(KFV2)+(FVA*XFVA2+FB*XFVB2)
        ELSE
          XPQ(KFV1)=XPQ(KFV1)+WTV1*(FVA*VI232A+FB*VI232B)
          XPQ(-KFV1)=XPQ(-KFV1)+WTV1*(FVA*VI232A+FB*VI232B)
          XPQ(KFV2)=XPQ(KFV2)+WTV2*(FVA*VI232A+FB*VI232B)
          XPQ(-KFV2)=XPQ(-KFV2)+WTV2*(FVA*VI232A+FB*VI232B)          
        ENDIF
        MINT(92)=3
 
C...Small Q2 and small x: dampen boundary value and add term.
      ELSE
 
C...Evaluate at boundary and define dampening factors.
        CALL PYSTFU(KFC,XMN,Q2MN,XPA)
        FB=(XMN-X)*(Q2MN-Q2)/(XMN*Q2MN)
        FA=1.-FB
        FVC=(X/XMN)**0.45*(Q2/(Q2+RMR))**0.55
        FVA=FVC*FA*((Q2MN+RMR)/Q2MN)**0.55
        FVB=FVC*FB*1.10*XMN**0.45*0.11
        FSC=(X/XMN)**(-0.08)*(Q2/(Q2+RMP))**1.08
        FSA=FSC*FA*((Q2MN+RMP)/Q2MN)**1.08
        FSB=FSC*FB*0.21*XMN**(-0.08)*0.21
 
C...Separate valence and sea parts of structure function.
        IF(KFA.NE.22) THEN
          XFV1=XPA(KFV1)-XPA(-KFV1)
          XPA(KFV1)=XPA(-KFV1)
          XFV2=XPA(KFV2)-XPA(-KFV2)
          XPA(KFV2)=XPA(-KFV2)
        ELSE
          XPA(KFV1)=XPA(KFV1)-WTV1*VINT(232)
          XPA(-KFV1)=XPA(-KFV1)-WTV1*VINT(232)
          XPA(KFV2)=XPA(KFV2)-WTV2*VINT(232)
          XPA(-KFV2)=XPA(-KFV2)-WTV2*VINT(232)
        ENDIF
 
C...Dampen valence and sea separately. Add constant terms.
C...Put back together.
        DO 140 KFL=-25,25
        XPQ(KFL)=FSA*XPA(KFL)
  140   CONTINUE
        IF(KFA.NE.22) THEN
          DO 150 KFL=-3,3
          XPQ(KFL)=XPQ(KFL)+FSB*WTSB(KFL)
  150     CONTINUE
          XPQ(KFV1)=XPQ(KFV1)+(FVA*XFV1+FVB*NV1)
          XPQ(KFV2)=XPQ(KFV2)+(FVA*XFV2+FVB*NV2)
        ELSE
          DO 160 KFL=-3,3
          XPQ(KFL)=XPQ(KFL)+VINT(281)*FSB*WTSB(KFL)
  160     CONTINUE
          XPQ(KFV1)=XPQ(KFV1)+WTV1*(FVA*VINT(232)+FVB*VINT(281))
          XPQ(-KFV1)=XPQ(-KFV1)+WTV1*(FVA*VINT(232)+FVB*VINT(281))
          XPQ(KFV2)=XPQ(KFV2)+WTV2*(FVA*VINT(232)+FVB*VINT(281))
          XPQ(-KFV2)=XPQ(-KFV2)+WTV2*(FVA*VINT(232)+FVB*VINT(281))          
        ENDIF
        XPQ(21)=XPQ(0)
        MINT(92)=4
      ENDIF
 
C...Format for error printout.
 5000 FORMAT(' Error: x value outside physical range; x =',1P,E12.3)
 
      RETURN
      END
