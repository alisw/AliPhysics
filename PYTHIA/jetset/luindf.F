 
C********************************************************************* 
 
      SUBROUTINE LUINDF(IP) 
 
C...Purpose: to handle the fragmentation of a jet system (or a single 
C...jet) according to independent fragmentation models. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION DPS(5),PSI(4),NFI(3),NFL(3),IFET(3),KFLF(3), 
     &KFLO(2),PXO(2),PYO(2),WO(2) 
 
C...Reset counters. Identify parton system and take copy. Check flavour. 
      NSAV=N 
      MSTU90=MSTU(90) 
      NJET=0 
      KQSUM=0 
      DO 100 J=1,5 
      DPS(J)=0. 
  100 CONTINUE 
      I=IP-1 
  110 I=I+1 
      IF(I.GT.MIN(N,MSTU(4)-MSTU(32))) THEN 
        CALL LUERRM(12,'(LUINDF:) failed to reconstruct jet system') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(K(I,1).NE.1.AND.K(I,1).NE.2) GOTO 110 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 110 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 110 
      NJET=NJET+1 
      IF(KQ.NE.2) KQSUM=KQSUM+KQ 
      DO 120 J=1,5 
      K(NSAV+NJET,J)=K(I,J) 
      P(NSAV+NJET,J)=P(I,J) 
      DPS(J)=DPS(J)+P(I,J) 
  120 CONTINUE 
      K(NSAV+NJET,3)=I 
      IF(K(I,1).EQ.2.OR.(MSTJ(3).LE.5.AND.N.GT.I.AND. 
     &K(I+1,1).EQ.2)) GOTO 110 
      IF(NJET.NE.1.AND.KQSUM.NE.0) THEN 
        CALL LUERRM(12,'(LUINDF:) unphysical flavour combination') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Boost copied system to CM frame. Find CM energy and sum flavours. 
      IF(NJET.NE.1) THEN 
        MSTU(33)=1 
        CALL LUDBRB(NSAV+1,NSAV+NJET,0.,0.,-DPS(1)/DPS(4), 
     &  -DPS(2)/DPS(4),-DPS(3)/DPS(4)) 
      ENDIF 
      PECM=0. 
      DO 130 J=1,3 
      NFI(J)=0 
  130 CONTINUE 
      DO 140 I=NSAV+1,NSAV+NJET 
      PECM=PECM+P(I,4) 
      KFA=IABS(K(I,2)) 
      IF(KFA.LE.3) THEN 
        NFI(KFA)=NFI(KFA)+ISIGN(1,K(I,2)) 
      ELSEIF(KFA.GT.1000) THEN 
        KFLA=MOD(KFA/1000,10) 
        KFLB=MOD(KFA/100,10) 
        IF(KFLA.LE.3) NFI(KFLA)=NFI(KFLA)+ISIGN(1,K(I,2)) 
        IF(KFLB.LE.3) NFI(KFLB)=NFI(KFLB)+ISIGN(1,K(I,2)) 
      ENDIF 
  140 CONTINUE 
 
C...Loop over attempts made. Reset counters. 
      NTRY=0 
  150 NTRY=NTRY+1 
      IF(NTRY.GT.200) THEN 
        CALL LUERRM(14,'(LUINDF:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      N=NSAV+NJET 
      MSTU(90)=MSTU90 
      DO 160 J=1,3 
      NFL(J)=NFI(J) 
      IFET(J)=0 
      KFLF(J)=0 
  160 CONTINUE 
 
C...Loop over jets to be fragmented. 
      DO 230 IP1=NSAV+1,NSAV+NJET 
      MSTJ(91)=0 
      NSAV1=N 
      MSTU91=MSTU(90) 
 
C...Initial flavour and momentum values. Jet along +z axis. 
      KFLH=IABS(K(IP1,2)) 
      IF(KFLH.GT.10) KFLH=MOD(KFLH/1000,10) 
      KFLO(2)=0 
      WF=P(IP1,4)+SQRT(P(IP1,1)**2+P(IP1,2)**2+P(IP1,3)**2) 
 
C...Initial values for quark or diquark jet. 
  170 IF(IABS(K(IP1,2)).NE.21) THEN 
        NSTR=1 
        KFLO(1)=K(IP1,2) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        WO(1)=WF 
 
C...Initial values for gluon treated like random quark jet. 
      ELSEIF(MSTJ(2).LE.2) THEN 
        NSTR=1 
        IF(MSTJ(2).EQ.2) MSTJ(91)=1 
        KFLO(1)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        WO(1)=WF 
 
C...Initial values for gluon treated like quark-antiquark jet pair, 
C...sharing energy according to Altarelli-Parisi splitting function. 
      ELSE 
        NSTR=2 
        IF(MSTJ(2).EQ.4) MSTJ(91)=1 
        KFLO(1)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5) 
        KFLO(2)=-KFLO(1) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        PXO(2)=-PXO(1) 
        PYO(2)=-PYO(1) 
        WO(1)=WF*RLU(0)**(1./3.) 
        WO(2)=WF-WO(1) 
      ENDIF 
 
C...Initial values for rank, flavour, pT and W+. 
      DO 220 ISTR=1,NSTR 
  180 I=N 
      MSTU(90)=MSTU91 
      IRANK=0 
      KFL1=KFLO(ISTR) 
      PX1=PXO(ISTR) 
      PY1=PYO(ISTR) 
      W=WO(ISTR) 
 
C...New hadron. Generate flavour and hadron species. 
  190 I=I+1 
      IF(I.GE.MSTU(4)-MSTU(32)-NJET-5) THEN 
        CALL LUERRM(11,'(LUINDF:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IRANK=IRANK+1 
      K(I,1)=1 
      K(I,3)=IP1 
      K(I,4)=0 
      K(I,5)=0 
  200 CALL LUKFDI(KFL1,0,KFL2,K(I,2)) 
      IF(K(I,2).EQ.0) GOTO 180 
      IF(MSTJ(12).GE.3.AND.IRANK.EQ.1.AND.IABS(KFL1).LE.10.AND. 
     &IABS(KFL2).GT.10) THEN 
        IF(RLU(0).GT.PARJ(19)) GOTO 200 
      ENDIF 
 
C...Find hadron mass. Generate four-momentum. 
      P(I,5)=ULMASS(K(I,2)) 
      CALL LUPTDI(KFL1,PX2,PY2) 
      P(I,1)=PX1+PX2 
      P(I,2)=PY1+PY2 
      PR=P(I,5)**2+P(I,1)**2+P(I,2)**2 
      CALL LUZDIS(KFL1,KFL2,PR,Z) 
      MZSAV=0 
      IF(IABS(KFL1).GE.4.AND.IABS(KFL1).LE.8.AND.MSTU(90).LT.8) THEN 
        MZSAV=1 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I 
        PARU(90+MSTU(90))=Z 
      ENDIF 
      P(I,3)=0.5*(Z*W-PR/MAX(1E-4,Z*W)) 
      P(I,4)=0.5*(Z*W+PR/MAX(1E-4,Z*W)) 
      IF(MSTJ(3).GE.1.AND.IRANK.EQ.1.AND.KFLH.GE.4.AND. 
     &P(I,3).LE.0.001) THEN 
        IF(W.GE.P(I,5)+0.5*PARJ(32)) GOTO 180 
        P(I,3)=0.0001 
        P(I,4)=SQRT(PR) 
        Z=P(I,4)/W 
      ENDIF 
 
C...Remaining flavour and momentum. 
      KFL1=-KFL2 
      PX1=-PX2 
      PY1=-PY2 
      W=(1.-Z)*W 
      DO 210 J=1,5 
      V(I,J)=0. 
  210 CONTINUE 
 
C...Check if pL acceptable. Go back for new hadron if enough energy. 
      IF(MSTJ(3).GE.0.AND.P(I,3).LT.0.) THEN 
        I=I-1 
        IF(MZSAV.EQ.1) MSTU(90)=MSTU(90)-1 
      ENDIF 
      IF(W.GT.PARJ(31)) GOTO 190 
      N=I 
  220 CONTINUE 
      IF(MOD(MSTJ(3),5).EQ.4.AND.N.EQ.NSAV1) WF=WF+0.1*PARJ(32) 
      IF(MOD(MSTJ(3),5).EQ.4.AND.N.EQ.NSAV1) GOTO 170 
 
C...Rotate jet to new direction. 
      THE=ULANGL(P(IP1,3),SQRT(P(IP1,1)**2+P(IP1,2)**2)) 
      PHI=ULANGL(P(IP1,1),P(IP1,2)) 
      MSTU(33)=1 
      CALL LUDBRB(NSAV1+1,N,THE,PHI,0D0,0D0,0D0) 
      K(K(IP1,3),4)=NSAV1+1 
      K(K(IP1,3),5)=N 
 
C...End of jet generation loop. Skip conservation in some cases. 
  230 CONTINUE 
      IF(NJET.EQ.1.OR.MSTJ(3).LE.0) GOTO 490 
      IF(MOD(MSTJ(3),5).NE.0.AND.N-NSAV-NJET.LT.2) GOTO 150 
 
C...Subtract off produced hadron flavours, finished if zero. 
      DO 240 I=NSAV+NJET+1,N 
      KFA=IABS(K(I,2)) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      IF(KFLA.EQ.0) THEN 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)-ISIGN(1,K(I,2))*(-1)**KFLB 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)+ISIGN(1,K(I,2))*(-1)**KFLB 
      ELSE 
        IF(KFLA.LE.3) NFL(KFLA)=NFL(KFLA)-ISIGN(1,K(I,2)) 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)-ISIGN(1,K(I,2)) 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)-ISIGN(1,K(I,2)) 
      ENDIF 
  240 CONTINUE 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREQ.EQ.0) GOTO 320 
 
C...Take away flavour of low-momentum particles until enough freedom. 
      NREM=0 
  250 IREM=0 
      P2MIN=PECM**2 
      DO 260 I=NSAV+NJET+1,N 
      P2=P(I,1)**2+P(I,2)**2+P(I,3)**2 
      IF(K(I,1).EQ.1.AND.P2.LT.P2MIN) IREM=I 
      IF(K(I,1).EQ.1.AND.P2.LT.P2MIN) P2MIN=P2 
  260 CONTINUE 
      IF(IREM.EQ.0) GOTO 150 
      K(IREM,1)=7 
      KFA=IABS(K(IREM,2)) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      IF(KFLA.GE.4.OR.KFLB.GE.4) K(IREM,1)=8 
      IF(K(IREM,1).EQ.8) GOTO 250 
      IF(KFLA.EQ.0) THEN 
        ISGN=ISIGN(1,K(IREM,2))*(-1)**KFLB 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)+ISGN 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)-ISGN 
      ELSE 
        IF(KFLA.LE.3) NFL(KFLA)=NFL(KFLA)+ISIGN(1,K(IREM,2)) 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)+ISIGN(1,K(IREM,2)) 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)+ISIGN(1,K(IREM,2)) 
      ENDIF 
      NREM=NREM+1 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREQ.GT.NREM) GOTO 250 
      DO 270 I=NSAV+NJET+1,N 
      IF(K(I,1).EQ.8) K(I,1)=1 
  270 CONTINUE 
 
C...Find combination of existing and new flavours for hadron. 
  280 NFET=2 
      IF(NFL(1)+NFL(2)+NFL(3).NE.0) NFET=3 
      IF(NREQ.LT.NREM) NFET=1 
      IF(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3)).EQ.0) NFET=0 
      DO 290 J=1,NFET 
      IFET(J)=1+(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3)))*RLU(0) 
      KFLF(J)=ISIGN(1,NFL(1)) 
      IF(IFET(J).GT.IABS(NFL(1))) KFLF(J)=ISIGN(2,NFL(2)) 
      IF(IFET(J).GT.IABS(NFL(1))+IABS(NFL(2))) KFLF(J)=ISIGN(3,NFL(3)) 
  290 CONTINUE 
      IF(NFET.EQ.2.AND.(IFET(1).EQ.IFET(2).OR.KFLF(1)*KFLF(2).GT.0)) 
     &GOTO 280 
      IF(NFET.EQ.3.AND.(IFET(1).EQ.IFET(2).OR.IFET(1).EQ.IFET(3).OR. 
     &IFET(2).EQ.IFET(3).OR.KFLF(1)*KFLF(2).LT.0.OR.KFLF(1)*KFLF(3) 
     &.LT.0.OR.KFLF(1)*(NFL(1)+NFL(2)+NFL(3)).LT.0)) GOTO 280 
      IF(NFET.EQ.0) KFLF(1)=1+INT((2.+PARJ(2))*RLU(0)) 
      IF(NFET.EQ.0) KFLF(2)=-KFLF(1) 
      IF(NFET.EQ.1) KFLF(2)=ISIGN(1+INT((2.+PARJ(2))*RLU(0)),-KFLF(1)) 
      IF(NFET.LE.2) KFLF(3)=0 
      IF(KFLF(3).NE.0) THEN 
        KFLFC=ISIGN(1000*MAX(IABS(KFLF(1)),IABS(KFLF(3)))+ 
     &  100*MIN(IABS(KFLF(1)),IABS(KFLF(3)))+1,KFLF(1)) 
        IF(KFLF(1).EQ.KFLF(3).OR.(1.+3.*PARJ(4))*RLU(0).GT.1.) 
     &  KFLFC=KFLFC+ISIGN(2,KFLFC) 
      ELSE 
        KFLFC=KFLF(1) 
      ENDIF 
      CALL LUKFDI(KFLFC,KFLF(2),KFLDMP,KF) 
      IF(KF.EQ.0) GOTO 280 
      DO 300 J=1,MAX(2,NFET) 
      NFL(IABS(KFLF(J)))=NFL(IABS(KFLF(J)))-ISIGN(1,KFLF(J)) 
  300 CONTINUE 
 
C...Store hadron at random among free positions. 
      NPOS=MIN(1+INT(RLU(0)*NREM),NREM) 
      DO 310 I=NSAV+NJET+1,N 
      IF(K(I,1).EQ.7) NPOS=NPOS-1 
      IF(K(I,1).EQ.1.OR.NPOS.NE.0) GOTO 310 
      K(I,1)=1 
      K(I,2)=KF 
      P(I,5)=ULMASS(K(I,2)) 
      P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  310 CONTINUE 
      NREM=NREM-1 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREM.GT.0) GOTO 280 
 
C...Compensate for missing momentum in global scheme (3 options). 
  320 IF(MOD(MSTJ(3),5).NE.0.AND.MOD(MSTJ(3),5).NE.4) THEN 
        DO 340 J=1,3 
        PSI(J)=0. 
        DO 330 I=NSAV+NJET+1,N 
        PSI(J)=PSI(J)+P(I,J) 
  330   CONTINUE 
  340   CONTINUE 
        PSI(4)=PSI(1)**2+PSI(2)**2+PSI(3)**2 
        PWS=0. 
        DO 350 I=NSAV+NJET+1,N 
        IF(MOD(MSTJ(3),5).EQ.1) PWS=PWS+P(I,4) 
        IF(MOD(MSTJ(3),5).EQ.2) PWS=PWS+SQRT(P(I,5)**2+(PSI(1)*P(I,1)+ 
     &  PSI(2)*P(I,2)+PSI(3)*P(I,3))**2/PSI(4)) 
        IF(MOD(MSTJ(3),5).EQ.3) PWS=PWS+1. 
  350   CONTINUE 
        DO 370 I=NSAV+NJET+1,N 
        IF(MOD(MSTJ(3),5).EQ.1) PW=P(I,4) 
        IF(MOD(MSTJ(3),5).EQ.2) PW=SQRT(P(I,5)**2+(PSI(1)*P(I,1)+ 
     &  PSI(2)*P(I,2)+PSI(3)*P(I,3))**2/PSI(4)) 
        IF(MOD(MSTJ(3),5).EQ.3) PW=1. 
        DO 360 J=1,3 
        P(I,J)=P(I,J)-PSI(J)*PW/PWS 
  360   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  370   CONTINUE 
 
C...Compensate for missing momentum withing each jet separately. 
      ELSEIF(MOD(MSTJ(3),5).EQ.4) THEN 
        DO 390 I=N+1,N+NJET 
        K(I,1)=0 
        DO 380 J=1,5 
        P(I,J)=0. 
  380   CONTINUE 
  390   CONTINUE 
        DO 410 I=NSAV+NJET+1,N 
        IR1=K(I,3) 
        IR2=N+IR1-NSAV 
        K(IR2,1)=K(IR2,1)+1 
        PLS=(P(I,1)*P(IR1,1)+P(I,2)*P(IR1,2)+P(I,3)*P(IR1,3))/ 
     &  (P(IR1,1)**2+P(IR1,2)**2+P(IR1,3)**2) 
        DO 400 J=1,3 
        P(IR2,J)=P(IR2,J)+P(I,J)-PLS*P(IR1,J) 
  400   CONTINUE 
        P(IR2,4)=P(IR2,4)+P(I,4) 
        P(IR2,5)=P(IR2,5)+PLS 
  410   CONTINUE 
        PSS=0. 
        DO 420 I=N+1,N+NJET 
        IF(K(I,1).NE.0) PSS=PSS+P(I,4)/(PECM*(0.8*P(I,5)+0.2)) 
  420   CONTINUE 
        DO 440 I=NSAV+NJET+1,N 
        IR1=K(I,3) 
        IR2=N+IR1-NSAV 
        PLS=(P(I,1)*P(IR1,1)+P(I,2)*P(IR1,2)+P(I,3)*P(IR1,3))/ 
     &  (P(IR1,1)**2+P(IR1,2)**2+P(IR1,3)**2) 
        DO 430 J=1,3 
        P(I,J)=P(I,J)-P(IR2,J)/K(IR2,1)+(1./(P(IR2,5)*PSS)-1.)*PLS* 
     &  P(IR1,J) 
  430   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  440   CONTINUE 
      ENDIF 
 
C...Scale momenta for energy conservation. 
      IF(MOD(MSTJ(3),5).NE.0) THEN 
        PMS=0. 
        PES=0. 
        PQS=0. 
        DO 450 I=NSAV+NJET+1,N 
        PMS=PMS+P(I,5) 
        PES=PES+P(I,4) 
        PQS=PQS+P(I,5)**2/P(I,4) 
  450   CONTINUE 
        IF(PMS.GE.PECM) GOTO 150 
        NECO=0 
  460   NECO=NECO+1 
        PFAC=(PECM-PQS)/(PES-PQS) 
        PES=0. 
        PQS=0. 
        DO 480 I=NSAV+NJET+1,N 
        DO 470 J=1,3 
        P(I,J)=PFAC*P(I,J) 
  470   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
        PES=PES+P(I,4) 
        PQS=PQS+P(I,5)**2/P(I,4) 
  480   CONTINUE 
        IF(NECO.LT.10.AND.ABS(PECM-PES).GT.2E-6*PECM) GOTO 460 
      ENDIF 
 
C...Origin of produced particles and parton daughter pointers. 
  490 DO 500 I=NSAV+NJET+1,N 
      IF(MSTU(16).NE.2) K(I,3)=NSAV+1 
      IF(MSTU(16).EQ.2) K(I,3)=K(K(I,3),3) 
  500 CONTINUE 
      DO 510 I=NSAV+1,NSAV+NJET 
      I1=K(I,3) 
      K(I1,1)=K(I1,1)+10 
      IF(MSTU(16).NE.2) THEN 
        K(I1,4)=NSAV+1 
        K(I1,5)=NSAV+1 
      ELSE 
        K(I1,4)=K(I1,4)-NJET+1 
        K(I1,5)=K(I1,5)-NJET+1 
        IF(K(I1,5).LT.K(I1,4)) THEN 
          K(I1,4)=0 
          K(I1,5)=0 
        ENDIF 
      ENDIF 
  510 CONTINUE 
 
C...Document independent fragmentation system. Remove copy of jets. 
      NSAV=NSAV+1 
      K(NSAV,1)=11 
      K(NSAV,2)=93 
      K(NSAV,3)=IP 
      K(NSAV,4)=NSAV+1 
      K(NSAV,5)=N-NJET+1 
      DO 520 J=1,4 
      P(NSAV,J)=DPS(J) 
      V(NSAV,J)=V(IP,J) 
  520 CONTINUE 
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2)) 
      V(NSAV,5)=0. 
      DO 540 I=NSAV+NJET,N 
      DO 530 J=1,5 
      K(I-NJET+1,J)=K(I,J) 
      P(I-NJET+1,J)=P(I,J) 
      V(I-NJET+1,J)=V(I,J) 
  530 CONTINUE 
  540 CONTINUE 
      N=N-NJET+1 
      DO 550 IZ=MSTU90+1,MSTU(90) 
      MSTU(90+IZ)=MSTU(90+IZ)-NJET+1 
  550 CONTINUE 
 
C...Boost back particle system. Set production vertices. 
      IF(NJET.NE.1) CALL LUDBRB(NSAV+1,N,0.,0.,DPS(1)/DPS(4), 
     &DPS(2)/DPS(4),DPS(3)/DPS(4)) 
      DO 570 I=NSAV+1,N 
      DO 560 J=1,4 
      V(I,J)=V(IP,J) 
  560 CONTINUE 
  570 CONTINUE 
 
      RETURN 
      END 
