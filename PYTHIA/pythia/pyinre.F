 
C*********************************************************************
 
      SUBROUTINE PYINRE
 
C...Calculates full and effective widths of gauge bosons, stores
C...masses and widths, rescales coefficients to be used for
C...resonance production generation.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/LUDAT4/CHAF(500)
      CHARACTER CHAF*8
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/PYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT4/,/PYINT6/
      DIMENSION WDTP(0:40),WDTE(0:40,0:5),WDTPM(0:40),WDTEM(0:40,0:5)
      DIMENSION KCINP(16),KCORD(16),PMORD(16)
      DATA KCINP/23,24,25,6,7,8,17,18,32,34,35,36,37,38,39,40/
 
C...Born level couplings in MSSM Higgs doublet sector.
      XW=PARU(102)
      XWV=XW
      IF(MSTP(8).GE.2) XW=1.-(PMAS(24,1)/PMAS(23,1))**2
      XW1=1.-XW
      IF(MSTP(4).EQ.2) THEN
        TANBE=PARU(141)
        RATBE=((1.-TANBE**2)/(1.+TANBE**2))**2
        SQMZ=PMAS(23,1)**2
        SQMW=PMAS(24,1)**2
        SQMH=PMAS(25,1)**2
        SQMA=SQMH*(SQMZ-SQMH)/(SQMZ*RATBE-SQMH)
        SQMHP=0.5*(SQMA+SQMZ+SQRT((SQMA+SQMZ)**2-4.*SQMA*SQMZ*RATBE))
        SQMHC=SQMA+SQMW
        IF(SQMH.GE.SQMZ.OR.MIN(SQMA,SQMHP,SQMHC).LE.0.) THEN
          WRITE(MSTU(11),5000)
          STOP
        ENDIF
        PMAS(35,1)=SQRT(SQMHP)
        PMAS(36,1)=SQRT(SQMA)
        PMAS(37,1)=SQRT(SQMHC)
        ALSU=0.5*ATAN(2.*TANBE*(SQMA+SQMZ)/((1.-TANBE**2)*
     &  (SQMA-SQMZ)))
        BESU=ATAN(TANBE)
        PARU(142)=1.
        PARU(143)=1.
        PARU(161)=-SIN(ALSU)/COS(BESU)
        PARU(162)=COS(ALSU)/SIN(BESU)
        PARU(163)=PARU(161)
        PARU(164)=SIN(BESU-ALSU)
        PARU(165)=PARU(164)
        PARU(168)=SIN(BESU-ALSU)+0.5*COS(2.*BESU)*SIN(BESU+ALSU)/XW
        PARU(171)=COS(ALSU)/COS(BESU)
        PARU(172)=SIN(ALSU)/SIN(BESU)
        PARU(173)=PARU(171)
        PARU(174)=COS(BESU-ALSU)
        PARU(175)=PARU(174)
        PARU(176)=COS(2.*ALSU)*COS(BESU+ALSU)-2.*SIN(2.*ALSU)*
     &  SIN(BESU+ALSU)
        PARU(177)=COS(2.*BESU)*COS(BESU+ALSU)
        PARU(178)=COS(BESU-ALSU)-0.5*COS(2.*BESU)*COS(BESU+ALSU)/XW
        PARU(181)=TANBE
        PARU(182)=1./TANBE
        PARU(183)=PARU(181)
        PARU(184)=0.
        PARU(185)=PARU(184)
        PARU(186)=COS(BESU-ALSU)
        PARU(187)=SIN(BESU-ALSU)
        PARU(188)=PARU(186)
        PARU(189)=PARU(187)
        PARU(190)=0.
        PARU(195)=COS(BESU-ALSU)
      ENDIF
 
C...Change matrix element codes when top and 4th generation
C...decay before fragmentation.
      IF(MSTP(48).GE.1) THEN
        IOFF=MDCY(6,2)-1
        DO 100 I=4,7
        MDME(IOFF+I,2)=0
  100   CONTINUE
        MDME(IOFF+9,2)=0
      ENDIF
      IF(MSTP(6).EQ.1) THEN
        IOFF=MDCY(7,2)-1
        DO 110 I=1,4
        MDME(IOFF+I,2)=0
  110   CONTINUE
        IOFF=MDCY(8,2)-1
        DO 120 I=1,4
        MDME(IOFF+I,2)=0
  120   CONTINUE
        IOFF=MDCY(17,2)-1
        MDME(IOFF+2,2)=0
        MDME(IOFF+3,2)=0
        MDME(IOFF+4,2)=0
        IOFF=MDCY(18,2)-1
        MDME(IOFF+1,2)=0
        MDME(IOFF+2,2)=0
      ELSEIF(MSTP(49).GE.1) THEN
        IOFF=MDCY(7,2)-1
        DO 130 I=4,7
        MDME(IOFF+I,2)=0
  130   CONTINUE
        MDME(IOFF+9,2)=0
        MDME(IOFF+10,2)=0
        IOFF=MDCY(8,2)-1
        DO 140 I=4,7
        MDME(IOFF+I,2)=0
  140   CONTINUE
        MDME(IOFF+9,2)=0
        MDME(IOFF+10,2)=0
        IOFF=MDCY(17,2)-1
        MDME(IOFF+4,2)=0
        MDME(IOFF+6,2)=0
        IOFF=MDCY(18,2)-1
        MDME(IOFF+2,2)=0
        MDME(IOFF+3,2)=0
      ENDIF
 
C...Reset full and effective widths of gauge bosons.
      DO 160 I=21,40
      DO 150 J=0,40
      WIDP(I,J)=0.
      WIDE(I,J)=0.
  150 CONTINUE
      WIDS(I,1)=1.
      WIDS(I,2)=1.
      WIDS(I,3)=1.
  160 CONTINUE
 
C...Order resonances by increasing mass (except Z0 and W+/-).
      DO 170 I=1,3
      KCORD(I)=KCINP(I)
      PMORD(I)=PMAS(KCORD(I),1)
  170 CONTINUE
      DO 200 I=4,16
      KCIN=KCINP(I)
      PMIN=PMAS(KCIN,1)
      DO 180 I1=I-1,3,-1
      IF(PMIN.GE.PMORD(I1)) GOTO 190
      KCORD(I1+1)=KCORD(I1)
      PMORD(I1+1)=PMORD(I1)
  180 CONTINUE
  190 KCORD(I1+1)=KCIN
      PMORD(I1+1)=PMIN
  200 CONTINUE
 
C...Loop over possible resonances.
      DO 250 I=1,16
      KC=KCORD(I)
      IF(KC.EQ.6.AND.MSTP(48).LE.0) GOTO 250
      IF(KC.EQ.7.OR.KC.EQ.8.OR.KC.EQ.17.OR.KC.EQ.18) THEN
        IF(MSTP(6).NE.1.AND.(MSTP(49).LE.0.OR.MSTP(1).LE.3)) GOTO 250
        IF(KC.EQ.18.AND.PMORD(I).LT.1.) GOTO 250
      ENDIF
      KCL=KC
      IF(KC.GE.6.AND.KC.LE.8) KCL=KC+20
      IF(KC.EQ.17.OR.KC.EQ.18) KCL=KC+12
 
C...Change decay modes for q* and l*.
      IF(MSTP(6).EQ.1.AND.(KC.EQ.7.OR.KC.EQ.8.OR.KC.EQ.17.OR.
     &KC.EQ.18)) THEN
        DO 210 J=1,MDCY(KC,3)
        IDC=J+MDCY(KC,2)-1
        KF2=KFDP(IDC,2)
        IF(KF2.EQ.7.OR.KF2.EQ.8.OR.KF2.EQ.17.OR.KF2.EQ.18)
     &  KFDP(IDC,2)=KF2-6
  210   CONTINUE
      ENDIF
 
C...Check that no fourth generation channels on by mistake.
      IF(MSTP(1).LE.3) THEN
        DO 220 J=1,MDCY(KC,3)
        IDC=J+MDCY(KC,2)-1
        KFA1=IABS(KFDP(IDC,1))
        KFA2=IABS(KFDP(IDC,2))
        IF(KFA1.EQ.7.OR.KFA1.EQ.8.OR.KFA1.EQ.17.OR.KFA1.EQ.18.OR.KFA2
     &  .EQ.7.OR.KFA2.EQ.8.OR.KFA2.EQ.17.OR.KFA2.EQ.18) MDME(IDC,1)=-1
  220   CONTINUE
      ENDIF
 
C...Find mass and evaluate width.
      PMR=PMAS(KC,1)
      IF(KC.EQ.25.OR.KC.EQ.35.OR.KC.EQ.36) MINT(62)=1
      CALL PYWIDT(KC,PMR**2,WDTP,WDTE)
      IF(KC.EQ.6.OR.KC.EQ.7.OR.KC.EQ.8.OR.KC.EQ.17.OR.KC.EQ.18)
     &CALL PYWIDT(-KC,PMR**2,WDTPM,WDTEM)
      MINT(51)=0
 
C...Evaluate suppression factors due to non-simulated channels.
      IF(KCHG(KC,3).EQ.0) THEN
        WIDS(KCL,1)=((WDTE(0,1)+WDTE(0,2))**2+
     &  2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+
     &  2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        WIDS(KCL,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
        WIDS(KCL,3)=0.
      ELSEIF(KC.EQ.6.OR.KC.EQ.7.OR.KC.EQ.8.OR.KC.EQ.17.OR.KC.EQ.18) THEN
        WIDS(KCL,1)=((WDTE(0,1)+WDTE(0,2))*(WDTEM(0,1)+WDTEM(0,3))+
     &  (WDTE(0,1)+WDTE(0,2))*(WDTEM(0,4)+WDTEM(0,5))+
     &  (WDTE(0,4)+WDTE(0,5))*(WDTEM(0,1)+WDTEM(0,3))+
     &  WDTE(0,4)*WDTEM(0,5)+WDTE(0,5)*WDTEM(0,4))/WDTP(0)**2
        WIDS(KCL,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
        WIDS(KCL,3)=(WDTEM(0,1)+WDTEM(0,3)+WDTEM(0,4))/WDTP(0)
      ELSE
        WIDS(KCL,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+
     &  2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        WIDS(KCL,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
        WIDS(KCL,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)
        IF(KC.EQ.24) THEN
          VINT(91)=((WDTE(0,1)+WDTE(0,2))**2+2.*(WDTE(0,1)+WDTE(0,2))*
     &    (WDTE(0,4)+WDTE(0,5))+2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
          VINT(92)=((WDTE(0,1)+WDTE(0,3))**2+2.*(WDTE(0,1)+WDTE(0,3))*
     &    (WDTE(0,4)+WDTE(0,5))+2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        ENDIF
      ENDIF
 
C...Find factors to give widths in GeV.
      AEM=ULALEM(PMR**2)
      IF(MSTP(8).GE.1) AEM=SQRT(2.)*PARU(105)*PMAS(24,1)**2*XW/PARU(1)
      IF(KC.LE.20) THEN
        FAC=PMR
      ELSEIF(KC.EQ.23.OR.KC.EQ.32) THEN
        FAC=AEM/(48.*XW*XW1)*PMR
      ELSEIF(KC.EQ.24.OR.KC.EQ.34) THEN
        FAC=AEM/(24.*XW)*PMR
      ELSEIF(KC.EQ.25.OR.KC.EQ.35.OR.KC.EQ.36.OR.KC.EQ.37) THEN
        FAC=AEM/(8.*XW)*(PMR/PMAS(24,1))**2*PMR
      ELSEIF(KC.EQ.38) THEN
        FAC=PMR
      ELSEIF(KC.EQ.39) THEN
        FAC=AEM/4.*PMR
      ELSEIF(KC.EQ.40) THEN
        FAC=AEM/(12.*XW)*PMR
      ENDIF
 
C...Translate widths into GeV and save them.
      DO 230 J=0,40
      WIDP(KCL,J)=FAC*WDTP(J)
      WIDE(KCL,J)=FAC*WDTE(J,0)
  230 CONTINUE
 
C...Set resonance widths and branching ratios in JETSET;
C...also on/off switch for decays in PYTHIA/JETSET.
      PMAS(KC,2)=WIDP(KCL,0)
      PMAS(KC,3)=MIN(0.9*PMAS(KC,1),10.*PMAS(KC,2))
      MDCY(KC,1)=MSTP(41)
      DO 240 J=1,MDCY(KC,3)
      IDC=J+MDCY(KC,2)-1
      BRAT(IDC)=0.
      IF(WIDE(KCL,0).GT.0.) BRAT(IDC)=WIDE(KCL,J)/WIDE(KCL,0)
  240 CONTINUE
  250 CONTINUE
 
C...Flavours of leptoquark: redefine charge and name.
      KFLQQ=KFDP(MDCY(39,2),1)
      KFLQL=KFDP(MDCY(39,2),2)
      KCHG(39,1)=KCHG(IABS(KFLQQ),1)*ISIGN(1,KFLQQ)+
     &KCHG(IABS(KFLQL),1)*ISIGN(1,KFLQL)
      CHAF(39)(4:4)=CHAF(IABS(KFLQQ))(1:1)
      CHAF(39)(5:7)=CHAF(IABS(KFLQL))(1:3)
 
C...Scenario with q* and l*: redefine names.
      IF(MSTP(6).EQ.1) THEN
        CHAF(7)='d*'
        CHAF(8)='u*'
        CHAF(17)='e*'
        CHAF(18)='nu*_e'
      ENDIF
 
C...Special cases in treatment of gamma*/Z0: redefine process name.
      IF(MSTP(43).EQ.1) THEN
        PROC(1)='f + f~ -> gamma*'
        PROC(15)='f + f~ -> g + gamma*'
        PROC(19)='f + f~ -> gamma + gamma*'
        PROC(30)='f + g -> f + gamma*'
        PROC(35)='f + gamma -> f + gamma*'
      ELSEIF(MSTP(43).EQ.2) THEN
        PROC(1)='f + f~ -> Z0'
        PROC(15)='f + f~ -> g + Z0'
        PROC(19)='f + f~ -> gamma + Z0'
        PROC(30)='f + g -> f + Z0'
        PROC(35)='f + gamma -> f + Z0'
      ELSEIF(MSTP(43).EQ.3) THEN
        PROC(1)='f + f~ -> gamma*/Z0'
        PROC(15)='f + f~ -> g + gamma*/Z0'
        PROC(19)='f + f~ -> gamma + gamma*/Z0'
        PROC(30)='f + g -> f + gamma*/Z0'
        PROC(35)='f + gamma -> f + gamma*/Z0'
      ENDIF
 
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
      IF(MSTP(44).EQ.1) THEN
        PROC(141)='f + f~ -> gamma*'
      ELSEIF(MSTP(44).EQ.2) THEN
        PROC(141)='f + f~ -> Z0'
      ELSEIF(MSTP(44).EQ.3) THEN
        PROC(141)='f + f~ -> Z''0'
      ELSEIF(MSTP(44).EQ.4) THEN
        PROC(141)='f + f~ -> gamma*/Z0'
      ELSEIF(MSTP(44).EQ.5) THEN
        PROC(141)='f + f~ -> gamma*/Z''0'
      ELSEIF(MSTP(44).EQ.6) THEN
        PROC(141)='f + f~ -> Z0/Z''0'
      ELSEIF(MSTP(44).EQ.7) THEN
        PROC(141)='f + f~ -> gamma*/Z0/Z''0'
      ENDIF
 
C...Special cases in treatment of WW -> WW: redefine process name.
      IF(MSTP(45).EQ.1) THEN
        PROC(77)='W+ + W+ -> W+ + W+'
      ELSEIF(MSTP(45).EQ.2) THEN
        PROC(77)='W+ + W- -> W+ + W-'
      ELSEIF(MSTP(45).EQ.3) THEN
        PROC(77)='W+/- + W+/- -> W+/- + W+/-'
      ENDIF
 
C...Format for error information.
 5000 FORMAT(1X,'Error: unphysical input tan^2(beta) and m_H ',
     &'combination'/1X,'Execution stopped!')
 
      RETURN
      END
