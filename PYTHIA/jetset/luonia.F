 
C********************************************************************* 
 
      SUBROUTINE LUONIA(KFL,ECM) 
 
C...Purpose: to generate Upsilon and toponium decays into three 
C...gluons or two gluons and a photon. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Printout. Check input parameters. 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IF(KFL.LT.0.OR.KFL.GT.8) THEN 
        CALL LUERRM(16,'(LUONIA:) called with unknown flavour code') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(ECM.LT.PARJ(127)+2.02*PARF(101)) THEN 
        CALL LUERRM(16,'(LUONIA:) called with too small CM energy') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Initial e+e- and onium state (optional). 
      NC=0 
      IF(MSTJ(115).GE.2) THEN 
        NC=NC+2 
        CALL LU1ENT(NC-1,11,0.5*ECM,0.,0.) 
        K(NC-1,1)=21 
        CALL LU1ENT(NC,-11,0.5*ECM,PARU(1),0.) 
        K(NC,1)=21 
      ENDIF 
      KFLC=IABS(KFL) 
      IF(MSTJ(115).GE.3.AND.KFLC.GE.5) THEN 
        NC=NC+1 
        KF=110*KFLC+3 
        MSTU10=MSTU(10) 
        MSTU(10)=1 
        P(NC,5)=ECM 
        CALL LU1ENT(NC,KF,ECM,0.,0.) 
        K(NC,1)=21 
        K(NC,3)=1 
        MSTU(10)=MSTU10 
      ENDIF 
 
C...Choose x1 and x2 according to matrix element. 
      NTRY=0 
  100 X1=RLU(0) 
      X2=RLU(0) 
      X3=2.-X1-X2 
      IF(X3.GE.1..OR.((1.-X1)/(X2*X3))**2+((1.-X2)/(X1*X3))**2+ 
     &((1.-X3)/(X1*X2))**2.LE.2.*RLU(0)) GOTO 100 
      NTRY=NTRY+1 
      NJET=3 
      IF(MSTJ(101).LE.4) CALL LU3ENT(NC+1,21,21,21,ECM,X1,X3) 
      IF(MSTJ(101).GE.5) CALL LU3ENT(-(NC+1),21,21,21,ECM,X1,X3) 
 
C...Photon-gluon-gluon events. Small system modifications. Jet origin. 
      MSTU(111)=MSTJ(108) 
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1)) 
     &MSTU(111)=1 
      PARU(112)=PARJ(121) 
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122) 
      QF=0. 
      IF(KFLC.NE.0) QF=KCHG(KFLC,1)/3. 
      RGAM=7.2*QF**2*PARU(101)/ULALPS(ECM**2) 
      MK=0 
      ECMC=ECM 
      IF(RLU(0).GT.RGAM/(1.+RGAM)) THEN 
        IF(1.-MAX(X1,X2,X3).LE.MAX((PARJ(126)/ECM)**2,PARJ(125))) 
     &  NJET=2 
        IF(NJET.EQ.2.AND.MSTJ(101).LE.4) CALL LU2ENT(NC+1,21,21,ECM) 
        IF(NJET.EQ.2.AND.MSTJ(101).GE.5) CALL LU2ENT(-(NC+1),21,21,ECM) 
      ELSE 
        MK=1 
        ECMC=SQRT(1.-X1)*ECM 
        IF(ECMC.LT.2.*PARJ(127)) GOTO 100 
        K(NC+1,1)=1 
        K(NC+1,2)=22 
        K(NC+1,4)=0 
        K(NC+1,5)=0 
        IF(MSTJ(101).GE.5) K(NC+2,4)=MSTU(5)*(NC+3) 
        IF(MSTJ(101).GE.5) K(NC+2,5)=MSTU(5)*(NC+3) 
        IF(MSTJ(101).GE.5) K(NC+3,4)=MSTU(5)*(NC+2) 
        IF(MSTJ(101).GE.5) K(NC+3,5)=MSTU(5)*(NC+2) 
        NJET=2 
        IF(ECMC.LT.4.*PARJ(127)) THEN 
          MSTU10=MSTU(10) 
          MSTU(10)=1 
          P(NC+2,5)=ECMC 
          CALL LU1ENT(NC+2,83,0.5*(X2+X3)*ECM,PARU(1),0.) 
          MSTU(10)=MSTU10 
          NJET=0 
        ENDIF 
      ENDIF 
      DO 110 IP=NC+1,N 
      K(IP,3)=K(IP,3)+(MSTJ(115)/2)+(KFLC/5)*(MSTJ(115)/3)*(NC-1) 
  110 CONTINUE 
 
C...Differential cross-sections. Upper limit for cross-section. 
      IF(MSTJ(106).EQ.1) THEN 
        SQ2=SQRT(2.) 
        HF1=1.-PARJ(131)*PARJ(132) 
        HF3=PARJ(133)**2 
        CT13=(X1*X3-2.*X1-2.*X3+2.)/(X1*X3) 
        ST13=SQRT(1.-CT13**2) 
        SIGL=0.5*X3**2*((1.-X2)**2+(1.-X3)**2)*ST13**2 
        SIGU=(X1*(1.-X1))**2+(X2*(1.-X2))**2+(X3*(1.-X3))**2-SIGL 
        SIGT=0.5*SIGL 
        SIGI=(SIGL*CT13/ST13+0.5*X1*X3*(1.-X2)**2*ST13)/SQ2 
        SIGMAX=(2.*HF1+HF3)*ABS(SIGU)+2.*(HF1+HF3)*ABS(SIGL)+2.*(HF1+ 
     &  2.*HF3)*ABS(SIGT)+2.*SQ2*(HF1+2.*HF3)*ABS(SIGI) 
 
C...Angular orientation of event. 
  120   CHI=PARU(2)*RLU(0) 
        CTHE=2.*RLU(0)-1. 
        PHI=PARU(2)*RLU(0) 
        CCHI=COS(CHI) 
        SCHI=SIN(CHI) 
        C2CHI=COS(2.*CHI) 
        S2CHI=SIN(2.*CHI) 
        THE=ACOS(CTHE) 
        STHE=SIN(THE) 
        C2PHI=COS(2.*(PHI-PARJ(134))) 
        S2PHI=SIN(2.*(PHI-PARJ(134))) 
        SIG=((1.+CTHE**2)*HF1+STHE**2*C2PHI*HF3)*SIGU+2.*(STHE**2*HF1- 
     &  STHE**2*C2PHI*HF3)*SIGL+2.*(STHE**2*C2CHI*HF1+((1.+CTHE**2)* 
     &  C2CHI*C2PHI-2.*CTHE*S2CHI*S2PHI)*HF3)*SIGT-2.*SQ2*(2.*STHE*CTHE* 
     &  CCHI*HF1-2.*STHE*(CTHE*CCHI*C2PHI-SCHI*S2PHI)*HF3)*SIGI 
        IF(SIG.LT.SIGMAX*RLU(0)) GOTO 120 
        CALL LUDBRB(NC+1,N,0.,CHI,0D0,0D0,0D0) 
        CALL LUDBRB(NC+1,N,THE,PHI,0D0,0D0,0D0) 
      ENDIF 
 
C...Generate parton shower. Rearrange along strings and check. 
      IF(MSTJ(101).GE.5.AND.NJET.GE.2) THEN 
        CALL LUSHOW(NC+MK+1,-NJET,ECMC) 
        MSTJ14=MSTJ(14) 
        IF(MSTJ(105).EQ.-1) MSTJ(14)=-1 
        IF(MSTJ(105).GE.0) MSTU(28)=0 
        CALL LUPREP(0) 
        MSTJ(14)=MSTJ14 
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100 
      ENDIF 
 
C...Generate fragmentation. Information for LUTABU: 
      IF(MSTJ(105).EQ.1) CALL LUEXEC 
      MSTU(161)=110*KFLC+3 
      MSTU(162)=0 
 
      RETURN 
      END 
