 
C********************************************************************* 
 
      SUBROUTINE LUSTRF(IP) 
C...Purpose: to handle the fragmentation of an arbitrary colour singlet 
C...jet system according to the Lund string fragmentation model. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION DPS(5),KFL(3),PMQ(3),PX(3),PY(3),GAM(3),IE(2),PR(2), 
     &IN(9),DHM(4),DHG(4),DP(5,5),IRANK(2),MJU(4),IJU(3),PJU(5,5), 
     &TJU(5),KFJH(2),NJS(2),KFJS(2),PJS(4,5),MSTU9T(8),PARU9T(8) 
 
C...Function: four-product of two vectors. 
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3) 
      DFOUR(I,J)=DP(I,4)*DP(J,4)-DP(I,1)*DP(J,1)-DP(I,2)*DP(J,2)- 
     &DP(I,3)*DP(J,3) 
 
C...Reset counters. Identify parton system. 
      MSTJ(91)=0 
      NSAV=N 
      MSTU90=MSTU(90) 
      NP=0 
      KQSUM=0 
      DO 100 J=1,5 
      DPS(J)=0D0 
  100 CONTINUE 
      MJU(1)=0 
      MJU(2)=0 
      I=IP-1 
  110 I=I+1 
      IF(I.GT.MIN(N,MSTU(4)-MSTU(32))) THEN 
        CALL LUERRM(12,'(LUSTRF:) failed to reconstruct jet system') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 110 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 110 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 110 
      IF(N+5*NP+11.GT.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Take copy of partons to be considered. Check flavour sum. 
      NP=NP+1 
      DO 120 J=1,5 
      K(N+NP,J)=K(I,J) 
      P(N+NP,J)=P(I,J) 
      IF(J.NE.4) DPS(J)=DPS(J)+P(I,J) 
  120 CONTINUE 
      DPS(4)=DPS(4)+SQRT(DBLE(P(I,1))**2+DBLE(P(I,2))**2+ 
     &DBLE(P(I,3))**2+DBLE(P(I,5))**2) 
      K(N+NP,3)=I 
      IF(KQ.NE.2) KQSUM=KQSUM+KQ 
      IF(K(I,1).EQ.41) THEN 
        KQSUM=KQSUM+2*KQ 
        IF(KQSUM.EQ.KQ) MJU(1)=N+NP 
        IF(KQSUM.NE.KQ) MJU(2)=N+NP 
      ENDIF 
      IF(K(I,1).EQ.2.OR.K(I,1).EQ.41) GOTO 110 
      IF(KQSUM.NE.0) THEN 
        CALL LUERRM(12,'(LUSTRF:) unphysical flavour combination') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Boost copied system to CM frame (for better numerical precision). 
      IF(ABS(DPS(3)).LT.0.99D0*DPS(4)) THEN 
        MBST=0 
        MSTU(33)=1 
        CALL LUDBRB(N+1,N+NP,0.,0.,-DPS(1)/DPS(4),-DPS(2)/DPS(4), 
     &  -DPS(3)/DPS(4)) 
      ELSE 
        MBST=1 
        HHBZ=SQRT(MAX(1D-6,DPS(4)+DPS(3))/MAX(1D-6,DPS(4)-DPS(3))) 
        DO 130 I=N+1,N+NP 
        HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2 
        IF(P(I,3).GT.0.) THEN 
          HHPEZ=(P(I,4)+P(I,3))/HHBZ 
          P(I,3)=0.5*(HHPEZ-HHPMT/HHPEZ) 
          P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ) 
        ELSE 
          HHPEZ=(P(I,4)-P(I,3))*HHBZ 
          P(I,3)=-0.5*(HHPEZ-HHPMT/HHPEZ) 
          P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ) 
        ENDIF 
  130   CONTINUE 
      ENDIF 
 
C...Search for very nearby partons that may be recombined. 
      NTRYR=0 
      PARU12=PARU(12) 
      PARU13=PARU(13) 
      MJU(3)=MJU(1) 
      MJU(4)=MJU(2) 
      NR=NP 
  140 IF(NR.GE.3) THEN 
        PDRMIN=2.*PARU12 
        DO 150 I=N+1,N+NR 
        IF(I.EQ.N+NR.AND.IABS(K(N+1,2)).NE.21) GOTO 150 
        I1=I+1 
        IF(I.EQ.N+NR) I1=N+1 
        IF(K(I,1).EQ.41.OR.K(I1,1).EQ.41) GOTO 150 
        IF(MJU(1).NE.0.AND.I1.LT.MJU(1).AND.IABS(K(I1,2)).NE.21) 
     &  GOTO 150 
        IF(MJU(2).NE.0.AND.I.GT.MJU(2).AND.IABS(K(I,2)).NE.21) GOTO 150 
        PAP=SQRT((P(I,1)**2+P(I,2)**2+P(I,3)**2)*(P(I1,1)**2+ 
     &  P(I1,2)**2+P(I1,3)**2)) 
        PVP=P(I,1)*P(I1,1)+P(I,2)*P(I1,2)+P(I,3)*P(I1,3) 
        PDR=4.*(PAP-PVP)**2/MAX(1E-6,PARU13**2*PAP+2.*(PAP-PVP)) 
        IF(PDR.LT.PDRMIN) THEN 
          IR=I 
          PDRMIN=PDR 
        ENDIF 
  150   CONTINUE 
 
C...Recombine very nearby partons to avoid machine precision problems. 
        IF(PDRMIN.LT.PARU12.AND.IR.EQ.N+NR) THEN 
          DO 160 J=1,4 
          P(N+1,J)=P(N+1,J)+P(N+NR,J) 
  160     CONTINUE 
          P(N+1,5)=SQRT(MAX(0.,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2- 
     &    P(N+1,3)**2)) 
          NR=NR-1 
          GOTO 140 
        ELSEIF(PDRMIN.LT.PARU12) THEN 
          DO 170 J=1,4 
          P(IR,J)=P(IR,J)+P(IR+1,J) 
  170     CONTINUE 
          P(IR,5)=SQRT(MAX(0.,P(IR,4)**2-P(IR,1)**2-P(IR,2)**2- 
     &    P(IR,3)**2)) 
          DO 190 I=IR+1,N+NR-1 
          K(I,2)=K(I+1,2) 
          DO 180 J=1,5 
          P(I,J)=P(I+1,J) 
  180     CONTINUE 
  190     CONTINUE 
          IF(IR.EQ.N+NR-1) K(IR,2)=K(N+NR,2) 
          NR=NR-1 
          IF(MJU(1).GT.IR) MJU(1)=MJU(1)-1 
          IF(MJU(2).GT.IR) MJU(2)=MJU(2)-1 
          GOTO 140 
        ENDIF 
      ENDIF 
      NTRYR=NTRYR+1 
 
C...Reset particle counter. Skip ahead if no junctions are present; 
C...this is usually the case! 
      NRS=MAX(5*NR+11,NP) 
      NTRY=0 
  200 NTRY=NTRY+1 
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN 
        PARU12=4.*PARU12 
        PARU13=2.*PARU13 
        GOTO 140 
      ELSEIF(NTRY.GT.100) THEN 
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      I=N+NRS 
      MSTU(90)=MSTU90 
      IF(MJU(1).EQ.0.AND.MJU(2).EQ.0) GOTO 580 
      DO 570 JT=1,2 
      NJS(JT)=0 
      IF(MJU(JT).EQ.0) GOTO 570 
      JS=3-2*JT 
 
C...Find and sum up momentum on three sides of junction. Check flavours. 
      DO 220 IU=1,3 
      IJU(IU)=0 
      DO 210 J=1,5 
      PJU(IU,J)=0. 
  210 CONTINUE 
  220 CONTINUE 
      IU=0 
      DO 240 I1=N+1+(JT-1)*(NR-1),N+NR+(JT-1)*(1-NR),JS 
      IF(K(I1,2).NE.21.AND.IU.LE.2) THEN 
        IU=IU+1 
        IJU(IU)=I1 
      ENDIF 
      DO 230 J=1,4 
      PJU(IU,J)=PJU(IU,J)+P(I1,J) 
  230 CONTINUE 
  240 CONTINUE 
      DO 250 IU=1,3 
      PJU(IU,5)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,3)**2) 
  250 CONTINUE 
      IF(K(IJU(3),2)/100.NE.10*K(IJU(1),2)+K(IJU(2),2).AND. 
     &K(IJU(3),2)/100.NE.10*K(IJU(2),2)+K(IJU(1),2)) THEN 
        CALL LUERRM(12,'(LUSTRF:) unphysical flavour combination') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Calculate (approximate) boost to rest frame of junction. 
      T12=(PJU(1,1)*PJU(2,1)+PJU(1,2)*PJU(2,2)+PJU(1,3)*PJU(2,3))/ 
     &(PJU(1,5)*PJU(2,5)) 
      T13=(PJU(1,1)*PJU(3,1)+PJU(1,2)*PJU(3,2)+PJU(1,3)*PJU(3,3))/ 
     &(PJU(1,5)*PJU(3,5)) 
      T23=(PJU(2,1)*PJU(3,1)+PJU(2,2)*PJU(3,2)+PJU(2,3)*PJU(3,3))/ 
     &(PJU(2,5)*PJU(3,5)) 
      T11=SQRT((2./3.)*(1.-T12)*(1.-T13)/(1.-T23)) 
      T22=SQRT((2./3.)*(1.-T12)*(1.-T23)/(1.-T13)) 
      TSQ=SQRT((2.*T11*T22+T12-1.)*(1.+T12)) 
      T1F=(TSQ-T22*(1.+T12))/(1.-T12**2) 
      T2F=(TSQ-T11*(1.+T12))/(1.-T12**2) 
      DO 260 J=1,3 
      TJU(J)=-(T1F*PJU(1,J)/PJU(1,5)+T2F*PJU(2,J)/PJU(2,5)) 
  260 CONTINUE 
      TJU(4)=SQRT(1.+TJU(1)**2+TJU(2)**2+TJU(3)**2) 
      DO 270 IU=1,3 
      PJU(IU,5)=TJU(4)*PJU(IU,4)-TJU(1)*PJU(IU,1)-TJU(2)*PJU(IU,2)- 
     &TJU(3)*PJU(IU,3) 
  270 CONTINUE 
 
C...Put junction at rest if motion could give inconsistencies. 
      IF(PJU(1,5)+PJU(2,5).GT.PJU(1,4)+PJU(2,4)) THEN 
        DO 280 J=1,3 
        TJU(J)=0. 
  280   CONTINUE 
        TJU(4)=1. 
        PJU(1,5)=PJU(1,4) 
        PJU(2,5)=PJU(2,4) 
        PJU(3,5)=PJU(3,4) 
      ENDIF 
 
C...Start preparing for fragmentation of two strings from junction. 
      ISTA=I 
      DO 550 IU=1,2 
      NS=IJU(IU+1)-IJU(IU) 
 
C...Junction strings: find longitudinal string directions. 
      DO 310 IS=1,NS 
      IS1=IJU(IU)+IS-1 
      IS2=IJU(IU)+IS 
      DO 290 J=1,5 
      DP(1,J)=0.5*P(IS1,J) 
      IF(IS.EQ.1) DP(1,J)=P(IS1,J) 
      DP(2,J)=0.5*P(IS2,J) 
      IF(IS.EQ.NS) DP(2,J)=-PJU(IU,J) 
  290 CONTINUE 
      IF(IS.EQ.NS) DP(2,4)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,3)**2) 
      IF(IS.EQ.NS) DP(2,5)=0. 
      DP(3,5)=DFOUR(1,1) 
      DP(4,5)=DFOUR(2,2) 
      DHKC=DFOUR(1,2) 
      IF(DP(3,5)+2.*DHKC+DP(4,5).LE.0.) THEN 
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2) 
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2) 
        DP(3,5)=0D0 
        DP(4,5)=0D0 
        DHKC=DFOUR(1,2) 
      ENDIF 
      DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5)) 
      DHK1=0.5*((DP(4,5)+DHKC)/DHKS-1.) 
      DHK2=0.5*((DP(3,5)+DHKC)/DHKS-1.) 
      IN1=N+NR+4*IS-3 
      P(IN1,5)=SQRT(DP(3,5)+2.*DHKC+DP(4,5)) 
      DO 300 J=1,4 
      P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J) 
      P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J) 
  300 CONTINUE 
  310 CONTINUE 
 
C...Junction strings: initialize flavour, momentum and starting pos. 
      ISAV=I 
      MSTU91=MSTU(90) 
  320 NTRY=NTRY+1 
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN 
        PARU12=4.*PARU12 
        PARU13=2.*PARU13 
        GOTO 140 
      ELSEIF(NTRY.GT.100) THEN 
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      I=ISAV 
      MSTU(90)=MSTU91 
      IRANKJ=0 
      IE(1)=K(N+1+(JT/2)*(NP-1),3) 
      IN(4)=N+NR+1 
      IN(5)=IN(4)+1 
      IN(6)=N+NR+4*NS+1 
      DO 340 JQ=1,2 
      DO 330 IN1=N+NR+2+JQ,N+NR+4*NS-2+JQ,4 
      P(IN1,1)=2-JQ 
      P(IN1,2)=JQ-1 
      P(IN1,3)=1. 
  330 CONTINUE 
  340 CONTINUE 
      KFL(1)=K(IJU(IU),2) 
      PX(1)=0. 
      PY(1)=0. 
      GAM(1)=0. 
      DO 350 J=1,5 
      PJU(IU+3,J)=0. 
  350 CONTINUE 
 
C...Junction strings: find initial transverse directions. 
      DO 360 J=1,4 
      DP(1,J)=P(IN(4),J) 
      DP(2,J)=P(IN(4)+1,J) 
      DP(3,J)=0. 
      DP(4,J)=0. 
  360 CONTINUE 
      DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2) 
      DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2) 
      DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4) 
      DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4) 
      DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4) 
      IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1. 
      IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1. 
      IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1. 
      IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1. 
      DHC12=DFOUR(1,2) 
      DHCX1=DFOUR(3,1)/DHC12 
      DHCX2=DFOUR(3,2)/DHC12 
      DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12) 
      DHCY1=DFOUR(4,1)/DHC12 
      DHCY2=DFOUR(4,2)/DHC12 
      DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12 
      DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2) 
      DO 370 J=1,4 
      DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J)) 
      P(IN(6),J)=DP(3,J) 
      P(IN(6)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)- 
     &DHCYX*DP(3,J)) 
  370 CONTINUE 
 
C...Junction strings: produce new particle, origin. 
  380 I=I+1 
      IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IRANKJ=IRANKJ+1 
      K(I,1)=1 
      K(I,3)=IE(1) 
      K(I,4)=0 
      K(I,5)=0 
 
C...Junction strings: generate flavour, hadron, pT, z and Gamma. 
  390 CALL LUKFDI(KFL(1),0,KFL(3),K(I,2)) 
      IF(K(I,2).EQ.0) GOTO 320 
      IF(MSTJ(12).GE.3.AND.IRANKJ.EQ.1.AND.IABS(KFL(1)).LE.10.AND. 
     &IABS(KFL(3)).GT.10) THEN 
        IF(RLU(0).GT.PARJ(19)) GOTO 390 
      ENDIF 
      P(I,5)=ULMASS(K(I,2)) 
      CALL LUPTDI(KFL(1),PX(3),PY(3)) 
      PR(1)=P(I,5)**2+(PX(1)+PX(3))**2+(PY(1)+PY(3))**2 
      CALL LUZDIS(KFL(1),KFL(3),PR(1),Z) 
      IF(IABS(KFL(1)).GE.4.AND.IABS(KFL(1)).LE.8.AND. 
     &MSTU(90).LT.8) THEN 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I 
        PARU(90+MSTU(90))=Z 
      ENDIF 
      GAM(3)=(1.-Z)*(GAM(1)+PR(1)/Z) 
      DO 400 J=1,3 
      IN(J)=IN(3+J) 
  400 CONTINUE 
 
C...Junction strings: stepping within or from 'low' string region easy. 
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)* 
     &P(IN(1),5)**2.GE.PR(1)) THEN 
        P(IN(1)+2,4)=Z*P(IN(1)+2,3) 
        P(IN(2)+2,4)=PR(1)/(P(IN(1)+2,4)*P(IN(1),5)**2) 
        DO 410 J=1,4 
        P(I,J)=(PX(1)+PX(3))*P(IN(3),J)+(PY(1)+PY(3))*P(IN(3)+1,J) 
  410   CONTINUE 
        GOTO 500 
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN 
        P(IN(2)+2,4)=P(IN(2)+2,3) 
        P(IN(2)+2,1)=1. 
        IN(2)=IN(2)+4 
        IF(IN(2).GT.N+NR+4*NS) GOTO 320 
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN 
          P(IN(1)+2,4)=P(IN(1)+2,3) 
          P(IN(1)+2,1)=0. 
          IN(1)=IN(1)+4 
        ENDIF 
      ENDIF 
 
C...Junction strings: find new transverse directions. 
  420 IF(IN(1).GT.N+NR+4*NS.OR.IN(2).GT.N+NR+4*NS.OR. 
     &IN(1).GT.IN(2)) GOTO 320 
      IF(IN(1).NE.IN(4).OR.IN(2).NE.IN(5)) THEN 
        DO 430 J=1,4 
        DP(1,J)=P(IN(1),J) 
        DP(2,J)=P(IN(2),J) 
        DP(3,J)=0. 
        DP(4,J)=0. 
  430   CONTINUE 
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2) 
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2) 
        DHC12=DFOUR(1,2) 
        IF(DHC12.LE.1E-2) THEN 
          P(IN(1)+2,4)=P(IN(1)+2,3) 
          P(IN(1)+2,1)=0. 
          IN(1)=IN(1)+4 
          GOTO 420 
        ENDIF 
        IN(3)=N+NR+4*NS+5 
        DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4) 
        DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4) 
        DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4) 
        IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1. 
        IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1. 
        IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1. 
        IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1. 
        DHCX1=DFOUR(3,1)/DHC12 
        DHCX2=DFOUR(3,2)/DHC12 
        DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12) 
        DHCY1=DFOUR(4,1)/DHC12 
        DHCY2=DFOUR(4,2)/DHC12 
        DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12 
        DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2) 
        DO 440 J=1,4 
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J)) 
        P(IN(3),J)=DP(3,J) 
        P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)- 
     &  DHCYX*DP(3,J)) 
  440   CONTINUE 
C...Express pT with respect to new axes, if sensible. 
        PXP=-(PX(3)*FOUR(IN(6),IN(3))+PY(3)*FOUR(IN(6)+1,IN(3))) 
        PYP=-(PX(3)*FOUR(IN(6),IN(3)+1)+PY(3)*FOUR(IN(6)+1,IN(3)+1)) 
        IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01) THEN 
          PX(3)=PXP 
          PY(3)=PYP 
        ENDIF 
      ENDIF 
 
C...Junction strings: sum up known four-momentum, coefficients for m2. 
      DO 470 J=1,4 
      DHG(J)=0. 
      P(I,J)=PX(1)*P(IN(6),J)+PY(1)*P(IN(6)+1,J)+PX(3)*P(IN(3),J)+ 
     &PY(3)*P(IN(3)+1,J) 
      DO 450 IN1=IN(4),IN(1)-4,4 
      P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J) 
  450 CONTINUE 
      DO 460 IN2=IN(5),IN(2)-4,4 
      P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J) 
  460 CONTINUE 
  470 CONTINUE 
      DHM(1)=FOUR(I,I) 
      DHM(2)=2.*FOUR(I,IN(1)) 
      DHM(3)=2.*FOUR(I,IN(2)) 
      DHM(4)=2.*FOUR(IN(1),IN(2)) 
 
C...Junction strings: find coefficients for Gamma expression. 
      DO 490 IN2=IN(1)+1,IN(2),4 
      DO 480 IN1=IN(1),IN2-1,4 
      DHC=2.*FOUR(IN1,IN2) 
      DHG(1)=DHG(1)+P(IN1+2,1)*P(IN2+2,1)*DHC 
      IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-P(IN2+2,1)*DHC 
      IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+P(IN1+2,1)*DHC 
      IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC 
  480 CONTINUE 
  490 CONTINUE 
 
C...Junction strings: solve (m2, Gamma) equation system for energies. 
      DHS1=DHM(3)*DHG(4)-DHM(4)*DHG(3) 
      IF(ABS(DHS1).LT.1E-4) GOTO 320 
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(2)*DHG(3)-DHG(4)* 
     &(P(I,5)**2-DHM(1))+DHG(2)*DHM(3) 
      DHS3=DHM(2)*(GAM(3)-DHG(1))-DHG(2)*(P(I,5)**2-DHM(1)) 
      P(IN(2)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/ABS(DHS1)- 
     &DHS2/DHS1) 
      IF(DHM(2)+DHM(4)*P(IN(2)+2,4).LE.0.) GOTO 320 
      P(IN(1)+2,4)=(P(I,5)**2-DHM(1)-DHM(3)*P(IN(2)+2,4))/ 
     &(DHM(2)+DHM(4)*P(IN(2)+2,4)) 
 
C...Junction strings: step to new region if necessary. 
      IF(P(IN(2)+2,4).GT.P(IN(2)+2,3)) THEN 
        P(IN(2)+2,4)=P(IN(2)+2,3) 
        P(IN(2)+2,1)=1. 
        IN(2)=IN(2)+4 
        IF(IN(2).GT.N+NR+4*NS) GOTO 320 
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN 
          P(IN(1)+2,4)=P(IN(1)+2,3) 
          P(IN(1)+2,1)=0. 
          IN(1)=IN(1)+4 
        ENDIF 
        GOTO 420 
      ELSEIF(P(IN(1)+2,4).GT.P(IN(1)+2,3)) THEN 
        P(IN(1)+2,4)=P(IN(1)+2,3) 
        P(IN(1)+2,1)=0. 
        IN(1)=IN(1)+JS 
        GOTO 820 
      ENDIF 
 
C...Junction strings: particle four-momentum, remainder, loop back. 
  500 DO 510 J=1,4 
      P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J) 
      PJU(IU+3,J)=PJU(IU+3,J)+P(I,J) 
  510 CONTINUE 
      IF(P(I,4).LT.P(I,5)) GOTO 320 
      PJU(IU+3,5)=TJU(4)*PJU(IU+3,4)-TJU(1)*PJU(IU+3,1)- 
     &TJU(2)*PJU(IU+3,2)-TJU(3)*PJU(IU+3,3) 
      IF(PJU(IU+3,5).LT.PJU(IU,5)) THEN 
        KFL(1)=-KFL(3) 
        PX(1)=-PX(3) 
        PY(1)=-PY(3) 
        GAM(1)=GAM(3) 
        IF(IN(3).NE.IN(6)) THEN 
          DO 520 J=1,4 
          P(IN(6),J)=P(IN(3),J) 
          P(IN(6)+1,J)=P(IN(3)+1,J) 
  520     CONTINUE 
        ENDIF 
        DO 530 JQ=1,2 
        IN(3+JQ)=IN(JQ) 
        P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4) 
        P(IN(JQ)+2,1)=P(IN(JQ)+2,1)-(3-2*JQ)*P(IN(JQ)+2,4) 
  530   CONTINUE 
        GOTO 380 
      ENDIF 
 
C...Junction strings: save quantities left after each string. 
      IF(IABS(KFL(1)).GT.10) GOTO 320 
      I=I-1 
      KFJH(IU)=KFL(1) 
      DO 540 J=1,4 
      PJU(IU+3,J)=PJU(IU+3,J)-P(I+1,J) 
  540 CONTINUE 
  550 CONTINUE 
 
C...Junction strings: put together to new effective string endpoint. 
      NJS(JT)=I-ISTA 
      KFJS(JT)=K(K(MJU(JT+2),3),2) 
      KFLS=2*INT(RLU(0)+3.*PARJ(4)/(1.+3.*PARJ(4)))+1 
      IF(KFJH(1).EQ.KFJH(2)) KFLS=3 
      IF(ISTA.NE.I) KFJS(JT)=ISIGN(1000*MAX(IABS(KFJH(1)), 
     &IABS(KFJH(2)))+100*MIN(IABS(KFJH(1)),IABS(KFJH(2)))+ 
     &KFLS,KFJH(1)) 
      DO 560 J=1,4 
      PJS(JT,J)=PJU(1,J)+PJU(2,J)+P(MJU(JT),J) 
      PJS(JT+2,J)=PJU(4,J)+PJU(5,J) 
  560 CONTINUE 
      PJS(JT,5)=SQRT(MAX(0.,PJS(JT,4)**2-PJS(JT,1)**2-PJS(JT,2)**2- 
     &PJS(JT,3)**2)) 
  570 CONTINUE 
 
C...Open versus closed strings. Choose breakup region for latter. 
  580 IF(MJU(1).NE.0.AND.MJU(2).NE.0) THEN 
        NS=MJU(2)-MJU(1) 
        NB=MJU(1)-N 
      ELSEIF(MJU(1).NE.0) THEN 
        NS=N+NR-MJU(1) 
        NB=MJU(1)-N 
      ELSEIF(MJU(2).NE.0) THEN 
        NS=MJU(2)-N 
        NB=1 
      ELSEIF(IABS(K(N+1,2)).NE.21) THEN 
        NS=NR-1 
        NB=1 
      ELSE 
        NS=NR+1 
        W2SUM=0. 
        DO 590 IS=1,NR 
        P(N+NR+IS,1)=0.5*FOUR(N+IS,N+IS+1-NR*(IS/NR)) 
        W2SUM=W2SUM+P(N+NR+IS,1) 
  590   CONTINUE 
        W2RAN=RLU(0)*W2SUM 
        NB=0 
  600   NB=NB+1 
        W2SUM=W2SUM-P(N+NR+NB,1) 
        IF(W2SUM.GT.W2RAN.AND.NB.LT.NR) GOTO 600 
      ENDIF 
 
C...Find longitudinal string directions (i.e. lightlike four-vectors). 
      DO 630 IS=1,NS 
      IS1=N+IS+NB-1-NR*((IS+NB-2)/NR) 
      IS2=N+IS+NB-NR*((IS+NB-1)/NR) 
      DO 610 J=1,5 
      DP(1,J)=P(IS1,J) 
      IF(IABS(K(IS1,2)).EQ.21) DP(1,J)=0.5*DP(1,J) 
      IF(IS1.EQ.MJU(1)) DP(1,J)=PJS(1,J)-PJS(3,J) 
      DP(2,J)=P(IS2,J) 
      IF(IABS(K(IS2,2)).EQ.21) DP(2,J)=0.5*DP(2,J) 
      IF(IS2.EQ.MJU(2)) DP(2,J)=PJS(2,J)-PJS(4,J) 
  610 CONTINUE 
      DP(3,5)=DFOUR(1,1) 
      DP(4,5)=DFOUR(2,2) 
      DHKC=DFOUR(1,2) 
      IF(DP(3,5)+2.*DHKC+DP(4,5).LE.0.) THEN 
        DP(3,5)=DP(1,5)**2 
        DP(4,5)=DP(2,5)**2 
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2+DP(1,5)**2) 
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2+DP(2,5)**2) 
        DHKC=DFOUR(1,2) 
      ENDIF 
      DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5)) 
      DHK1=0.5*((DP(4,5)+DHKC)/DHKS-1.) 
      DHK2=0.5*((DP(3,5)+DHKC)/DHKS-1.) 
      IN1=N+NR+4*IS-3 
      P(IN1,5)=SQRT(DP(3,5)+2.*DHKC+DP(4,5)) 
      DO 620 J=1,4 
      P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J) 
      P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J) 
  620 CONTINUE 
  630 CONTINUE 
 
C...Begin initialization: sum up energy, set starting position. 
      ISAV=I 
      MSTU91=MSTU(90) 
  640 NTRY=NTRY+1 
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN 
        PARU12=4.*PARU12 
        PARU13=2.*PARU13 
        GOTO 140 
      ELSEIF(NTRY.GT.100) THEN 
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      I=ISAV 
      MSTU(90)=MSTU91 
      DO 660 J=1,4 
      P(N+NRS,J)=0. 
      DO 650 IS=1,NR 
      P(N+NRS,J)=P(N+NRS,J)+P(N+IS,J) 
  650 CONTINUE 
  660 CONTINUE 
      DO 680 JT=1,2 
      IRANK(JT)=0 
      IF(MJU(JT).NE.0) IRANK(JT)=NJS(JT) 
      IF(NS.GT.NR) IRANK(JT)=1 
      IE(JT)=K(N+1+(JT/2)*(NP-1),3) 
      IN(3*JT+1)=N+NR+1+4*(JT/2)*(NS-1) 
      IN(3*JT+2)=IN(3*JT+1)+1 
      IN(3*JT+3)=N+NR+4*NS+2*JT-1 
      DO 670 IN1=N+NR+2+JT,N+NR+4*NS-2+JT,4 
      P(IN1,1)=2-JT 
      P(IN1,2)=JT-1 
      P(IN1,3)=1. 
  670 CONTINUE 
  680 CONTINUE 
 
C...Initialize flavour and pT variables for open string. 
      IF(NS.LT.NR) THEN 
        PX(1)=0. 
        PY(1)=0. 
        IF(NS.EQ.1.AND.MJU(1)+MJU(2).EQ.0) CALL LUPTDI(0,PX(1),PY(1)) 
        PX(2)=-PX(1) 
        PY(2)=-PY(1) 
        DO 690 JT=1,2 
        KFL(JT)=K(IE(JT),2) 
        IF(MJU(JT).NE.0) KFL(JT)=KFJS(JT) 
        MSTJ(93)=1 
        PMQ(JT)=ULMASS(KFL(JT)) 
        GAM(JT)=0. 
  690   CONTINUE 
 
C...Closed string: random initial breakup flavour, pT and vertex. 
      ELSE 
        KFL(3)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5) 
        CALL LUKFDI(KFL(3),0,KFL(1),KDUMP) 
        KFL(2)=-KFL(1) 
        IF(IABS(KFL(1)).GT.10.AND.RLU(0).GT.0.5) THEN 
          KFL(2)=-(KFL(1)+ISIGN(10000,KFL(1))) 
        ELSEIF(IABS(KFL(1)).GT.10) THEN 
          KFL(1)=-(KFL(2)+ISIGN(10000,KFL(2))) 
        ENDIF 
        CALL LUPTDI(KFL(1),PX(1),PY(1)) 
        PX(2)=-PX(1) 
        PY(2)=-PY(1) 
        PR3=MIN(25.,0.1*P(N+NR+1,5)**2) 
  700   CALL LUZDIS(KFL(1),KFL(2),PR3,Z) 
        ZR=PR3/(Z*P(N+NR+1,5)**2) 
        IF(ZR.GE.1.) GOTO 700 
        DO 710 JT=1,2 
        MSTJ(93)=1 
        PMQ(JT)=ULMASS(KFL(JT)) 
        GAM(JT)=PR3*(1.-Z)/Z 
        IN1=N+NR+3+4*(JT/2)*(NS-1) 
        P(IN1,JT)=1.-Z 
        P(IN1,3-JT)=JT-1 
        P(IN1,3)=(2-JT)*(1.-Z)+(JT-1)*Z 
        P(IN1+1,JT)=ZR 
        P(IN1+1,3-JT)=2-JT 
        P(IN1+1,3)=(2-JT)*(1.-ZR)+(JT-1)*ZR 
  710   CONTINUE 
      ENDIF 
 
C...Find initial transverse directions (i.e. spacelike four-vectors). 
      DO 750 JT=1,2 
      IF(JT.EQ.1.OR.NS.EQ.NR-1) THEN 
        IN1=IN(3*JT+1) 
        IN3=IN(3*JT+3) 
        DO 720 J=1,4 
        DP(1,J)=P(IN1,J) 
        DP(2,J)=P(IN1+1,J) 
        DP(3,J)=0. 
        DP(4,J)=0. 
  720   CONTINUE 
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2) 
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2) 
        DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4) 
        DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4) 
        DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4) 
        IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1. 
        IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1. 
        IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1. 
        IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1. 
        DHC12=DFOUR(1,2) 
        DHCX1=DFOUR(3,1)/DHC12 
        DHCX2=DFOUR(3,2)/DHC12 
        DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12) 
        DHCY1=DFOUR(4,1)/DHC12 
        DHCY2=DFOUR(4,2)/DHC12 
        DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12 
        DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2) 
        DO 730 J=1,4 
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J)) 
        P(IN3,J)=DP(3,J) 
        P(IN3+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)- 
     &  DHCYX*DP(3,J)) 
  730   CONTINUE 
      ELSE 
        DO 740 J=1,4 
        P(IN3+2,J)=P(IN3,J) 
        P(IN3+3,J)=P(IN3+1,J) 
  740   CONTINUE 
      ENDIF 
  750 CONTINUE 
 
C...Remove energy used up in junction string fragmentation. 
      IF(MJU(1)+MJU(2).GT.0) THEN 
        DO 770 JT=1,2 
        IF(NJS(JT).EQ.0) GOTO 770 
        DO 760 J=1,4 
        P(N+NRS,J)=P(N+NRS,J)-PJS(JT+2,J) 
  760   CONTINUE 
  770   CONTINUE 
      ENDIF 
 
C...Produce new particle: side, origin. 
  780 I=I+1 
      IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      JT=1.5+RLU(0) 
      IF(IABS(KFL(3-JT)).GT.10) JT=3-JT 
      IF(IABS(KFL(3-JT)).GE.4.AND.IABS(KFL(3-JT)).LE.8) JT=3-JT 
      JR=3-JT 
      JS=3-2*JT 
      IRANK(JT)=IRANK(JT)+1 
      K(I,1)=1 
      K(I,3)=IE(JT) 
      K(I,4)=0 
      K(I,5)=0 
 
C...Generate flavour, hadron and pT. 
  790 CALL LUKFDI(KFL(JT),0,KFL(3),K(I,2)) 
      IF(K(I,2).EQ.0) GOTO 640 
      IF(MSTJ(12).GE.3.AND.IRANK(JT).EQ.1.AND.IABS(KFL(JT)).LE.10.AND. 
     &IABS(KFL(3)).GT.10) THEN 
        IF(RLU(0).GT.PARJ(19)) GOTO 790 
      ENDIF 
      P(I,5)=ULMASS(K(I,2)) 
      CALL LUPTDI(KFL(JT),PX(3),PY(3)) 
      PR(JT)=P(I,5)**2+(PX(JT)+PX(3))**2+(PY(JT)+PY(3))**2 
 
C...Final hadrons for small invariant mass. 
      MSTJ(93)=1 
      PMQ(3)=ULMASS(KFL(3)) 
      PARJST=PARJ(33) 
      IF(MSTJ(11).EQ.2) PARJST=PARJ(34) 
      WMIN=PARJST+PMQ(1)+PMQ(2)+PARJ(36)*PMQ(3) 
      IF(IABS(KFL(JT)).GT.10.AND.IABS(KFL(3)).GT.10) WMIN= 
     &WMIN-0.5*PARJ(36)*PMQ(3) 
      WREM2=FOUR(N+NRS,N+NRS) 
      IF(WREM2.LT.0.10) GOTO 640 
      IF(WREM2.LT.MAX(WMIN*(1.+(2.*RLU(0)-1.)*PARJ(37)), 
     &PARJ(32)+PMQ(1)+PMQ(2))**2) GOTO 940 
 
C...Choose z, which gives Gamma. Shift z for heavy flavours. 
      CALL LUZDIS(KFL(JT),KFL(3),PR(JT),Z) 
      IF(IABS(KFL(JT)).GE.4.AND.IABS(KFL(JT)).LE.8.AND. 
     &MSTU(90).LT.8) THEN 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I 
        PARU(90+MSTU(90))=Z 
      ENDIF 
      KFL1A=IABS(KFL(1)) 
      KFL2A=IABS(KFL(2)) 
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10), 
     &MOD(KFL2A/1000,10)).GE.4) THEN 
        PR(JR)=(PMQ(JR)+PMQ(3))**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2 
        PW12=SQRT(MAX(0.,(WREM2-PR(1)-PR(2))**2-4.*PR(1)*PR(2))) 
        Z=(WREM2+PR(JT)-PR(JR)+PW12*(2.*Z-1.))/(2.*WREM2) 
        PR(JR)=(PMQ(JR)+PARJST)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2 
        IF((1.-Z)*(WREM2-PR(JT)/Z).LT.PR(JR)) GOTO 940 
      ENDIF 
      GAM(3)=(1.-Z)*(GAM(JT)+PR(JT)/Z) 
      DO 800 J=1,3 
      IN(J)=IN(3*JT+J) 
  800 CONTINUE 
 
C...Stepping within or from 'low' string region easy. 
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)* 
     &P(IN(1),5)**2.GE.PR(JT)) THEN 
        P(IN(JT)+2,4)=Z*P(IN(JT)+2,3) 
        P(IN(JR)+2,4)=PR(JT)/(P(IN(JT)+2,4)*P(IN(1),5)**2) 
        DO 810 J=1,4 
        P(I,J)=(PX(JT)+PX(3))*P(IN(3),J)+(PY(JT)+PY(3))*P(IN(3)+1,J) 
  810   CONTINUE 
        GOTO 900 
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN 
        P(IN(JR)+2,4)=P(IN(JR)+2,3) 
        P(IN(JR)+2,JT)=1. 
        IN(JR)=IN(JR)+4*JS 
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 640 
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN 
          P(IN(JT)+2,4)=P(IN(JT)+2,3) 
          P(IN(JT)+2,JT)=0. 
          IN(JT)=IN(JT)+4*JS 
        ENDIF 
      ENDIF 
 
C...Find new transverse directions (i.e. spacelike string vectors). 
  820 IF(JS*IN(1).GT.JS*IN(3*JR+1).OR.JS*IN(2).GT.JS*IN(3*JR+2).OR. 
     &IN(1).GT.IN(2)) GOTO 640 
      IF(IN(1).NE.IN(3*JT+1).OR.IN(2).NE.IN(3*JT+2)) THEN 
        DO 830 J=1,4 
        DP(1,J)=P(IN(1),J) 
        DP(2,J)=P(IN(2),J) 
        DP(3,J)=0. 
        DP(4,J)=0. 
  830   CONTINUE 
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2) 
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2) 
        DHC12=DFOUR(1,2) 
        IF(DHC12.LE.1E-2) THEN 
          P(IN(JT)+2,4)=P(IN(JT)+2,3) 
          P(IN(JT)+2,JT)=0. 
          IN(JT)=IN(JT)+4*JS 
          GOTO 820 
        ENDIF 
        IN(3)=N+NR+4*NS+5 
        DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4) 
        DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4) 
        DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4) 
        IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1. 
        IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1. 
        IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1. 
        IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1. 
        DHCX1=DFOUR(3,1)/DHC12 
        DHCX2=DFOUR(3,2)/DHC12 
        DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12) 
        DHCY1=DFOUR(4,1)/DHC12 
        DHCY2=DFOUR(4,2)/DHC12 
        DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12 
        DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2) 
        DO 840 J=1,4 
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J)) 
        P(IN(3),J)=DP(3,J) 
        P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)- 
     &  DHCYX*DP(3,J)) 
  840   CONTINUE 
C...Express pT with respect to new axes, if sensible. 
        PXP=-(PX(3)*FOUR(IN(3*JT+3),IN(3))+PY(3)* 
     &  FOUR(IN(3*JT+3)+1,IN(3))) 
        PYP=-(PX(3)*FOUR(IN(3*JT+3),IN(3)+1)+PY(3)* 
     &  FOUR(IN(3*JT+3)+1,IN(3)+1)) 
        IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01) THEN 
          PX(3)=PXP 
          PY(3)=PYP 
        ENDIF 
      ENDIF 
 
C...Sum up known four-momentum. Gives coefficients for m2 expression. 
      DO 870 J=1,4 
      DHG(J)=0. 
      P(I,J)=PX(JT)*P(IN(3*JT+3),J)+PY(JT)*P(IN(3*JT+3)+1,J)+ 
     &PX(3)*P(IN(3),J)+PY(3)*P(IN(3)+1,J) 
      DO 850 IN1=IN(3*JT+1),IN(1)-4*JS,4*JS 
      P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J) 
  850 CONTINUE 
      DO 860 IN2=IN(3*JT+2),IN(2)-4*JS,4*JS 
      P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J) 
  860 CONTINUE 
  870 CONTINUE 
      DHM(1)=FOUR(I,I) 
      DHM(2)=2.*FOUR(I,IN(1)) 
      DHM(3)=2.*FOUR(I,IN(2)) 
      DHM(4)=2.*FOUR(IN(1),IN(2)) 
 
C...Find coefficients for Gamma expression. 
      DO 890 IN2=IN(1)+1,IN(2),4 
      DO 880 IN1=IN(1),IN2-1,4 
      DHC=2.*FOUR(IN1,IN2) 
      DHG(1)=DHG(1)+P(IN1+2,JT)*P(IN2+2,JT)*DHC 
      IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-JS*P(IN2+2,JT)*DHC 
      IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+JS*P(IN1+2,JT)*DHC 
      IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC 
  880 CONTINUE 
  890 CONTINUE 
 
C...Solve (m2, Gamma) equation system for energies taken. 
      DHS1=DHM(JR+1)*DHG(4)-DHM(4)*DHG(JR+1) 
      IF(ABS(DHS1).LT.1E-4) GOTO 640 
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(JT+1)*DHG(JR+1)-DHG(4)* 
     &(P(I,5)**2-DHM(1))+DHG(JT+1)*DHM(JR+1) 
      DHS3=DHM(JT+1)*(GAM(3)-DHG(1))-DHG(JT+1)*(P(I,5)**2-DHM(1)) 
      P(IN(JR)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/ABS(DHS1)- 
     &DHS2/DHS1) 
      IF(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4).LE.0.) GOTO 640 
      P(IN(JT)+2,4)=(P(I,5)**2-DHM(1)-DHM(JR+1)*P(IN(JR)+2,4))/ 
     &(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4)) 
 
C...Step to new region if necessary. 
      IF(P(IN(JR)+2,4).GT.P(IN(JR)+2,3)) THEN 
        P(IN(JR)+2,4)=P(IN(JR)+2,3) 
        P(IN(JR)+2,JT)=1. 
        IN(JR)=IN(JR)+4*JS 
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 640 
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN 
          P(IN(JT)+2,4)=P(IN(JT)+2,3) 
          P(IN(JT)+2,JT)=0. 
          IN(JT)=IN(JT)+4*JS 
        ENDIF 
        GOTO 820 
      ELSEIF(P(IN(JT)+2,4).GT.P(IN(JT)+2,3)) THEN 
        P(IN(JT)+2,4)=P(IN(JT)+2,3) 
        P(IN(JT)+2,JT)=0. 
        IN(JT)=IN(JT)+4*JS 
        GOTO 820 
      ENDIF 
 
C...Four-momentum of particle. Remaining quantities. Loop back. 
  900 DO 910 J=1,4 
      P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J) 
      P(N+NRS,J)=P(N+NRS,J)-P(I,J) 
  910 CONTINUE 
      IF(P(I,4).LT.P(I,5)) GOTO 640 
      KFL(JT)=-KFL(3) 
      PMQ(JT)=PMQ(3) 
      PX(JT)=-PX(3) 
      PY(JT)=-PY(3) 
      GAM(JT)=GAM(3) 
      IF(IN(3).NE.IN(3*JT+3)) THEN 
        DO 920 J=1,4 
        P(IN(3*JT+3),J)=P(IN(3),J) 
        P(IN(3*JT+3)+1,J)=P(IN(3)+1,J) 
  920   CONTINUE 
      ENDIF 
      DO 930 JQ=1,2 
      IN(3*JT+JQ)=IN(JQ) 
      P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4) 
      P(IN(JQ)+2,JT)=P(IN(JQ)+2,JT)-JS*(3-2*JQ)*P(IN(JQ)+2,4) 
  930 CONTINUE 
      GOTO 780 
 
C...Final hadron: side, flavour, hadron, mass. 
  940 I=I+1 
      K(I,1)=1 
      K(I,3)=IE(JR) 
      K(I,4)=0 
      K(I,5)=0 
      CALL LUKFDI(KFL(JR),-KFL(3),KFLDMP,K(I,2)) 
      IF(K(I,2).EQ.0) GOTO 640 
      P(I,5)=ULMASS(K(I,2)) 
      PR(JR)=P(I,5)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2 
 
C...Final two hadrons: find common setup of four-vectors. 
      JQ=1 
      IF(P(IN(4)+2,3)*P(IN(5)+2,3)*FOUR(IN(4),IN(5)).LT.P(IN(7),3)* 
     &P(IN(8),3)*FOUR(IN(7),IN(8))) JQ=2 
      DHC12=FOUR(IN(3*JQ+1),IN(3*JQ+2)) 
      DHR1=FOUR(N+NRS,IN(3*JQ+2))/DHC12 
      DHR2=FOUR(N+NRS,IN(3*JQ+1))/DHC12 
      IF(IN(4).NE.IN(7).OR.IN(5).NE.IN(8)) THEN 
        PX(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3))-PX(JQ) 
        PY(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3)+1)-PY(JQ) 
        PR(3-JQ)=P(I+(JT+JQ-3)**2-1,5)**2+(PX(3-JQ)+(2*JQ-3)*JS* 
     &  PX(3))**2+(PY(3-JQ)+(2*JQ-3)*JS*PY(3))**2 
      ENDIF 
 
C...Solve kinematics for final two hadrons, if possible. 
      WREM2=WREM2+(PX(1)+PX(2))**2+(PY(1)+PY(2))**2 
      FD=(SQRT(PR(1))+SQRT(PR(2)))/SQRT(WREM2) 
      IF(MJU(1)+MJU(2).NE.0.AND.I.EQ.ISAV+2.AND.FD.GE.1.) GOTO 200 
      IF(FD.GE.1.) GOTO 640 
      FA=WREM2+PR(JT)-PR(JR) 
      IF(MSTJ(11).NE.2) PREV=0.5*EXP(MAX(-50.,LOG(FD)*PARJ(38)* 
     &(PR(1)+PR(2))**2)) 
      IF(MSTJ(11).EQ.2) PREV=0.5*FD**PARJ(39) 
      FB=SIGN(SQRT(MAX(0.,FA**2-4.*WREM2*PR(JT))),JS*(RLU(0)-PREV)) 
      KFL1A=IABS(KFL(1)) 
      KFL2A=IABS(KFL(2)) 
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10), 
     &MOD(KFL2A/1000,10)).GE.6) FB=SIGN(SQRT(MAX(0.,FA**2- 
     &4.*WREM2*PR(JT))),FLOAT(JS)) 
      DO 950 J=1,4 
      P(I-1,J)=(PX(JT)+PX(3))*P(IN(3*JQ+3),J)+(PY(JT)+PY(3))* 
     &P(IN(3*JQ+3)+1,J)+0.5*(DHR1*(FA+FB)*P(IN(3*JQ+1),J)+ 
     &DHR2*(FA-FB)*P(IN(3*JQ+2),J))/WREM2 
      P(I,J)=P(N+NRS,J)-P(I-1,J) 
  950 CONTINUE 
      IF(P(I-1,4).LT.P(I-1,5).OR.P(I,4).LT.P(I,5)) GOTO 640 
 
C...Mark jets as fragmented and give daughter pointers. 
      N=I-NRS+1 
      DO 960 I=NSAV+1,NSAV+NP 
      IM=K(I,3) 
      K(IM,1)=K(IM,1)+10 
      IF(MSTU(16).NE.2) THEN 
        K(IM,4)=NSAV+1 
        K(IM,5)=NSAV+1 
      ELSE 
        K(IM,4)=NSAV+2 
        K(IM,5)=N 
      ENDIF 
  960 CONTINUE 
 
C...Document string system. Move up particles. 
      NSAV=NSAV+1 
      K(NSAV,1)=11 
      K(NSAV,2)=92 
      K(NSAV,3)=IP 
      K(NSAV,4)=NSAV+1 
      K(NSAV,5)=N 
      DO 970 J=1,4 
      P(NSAV,J)=DPS(J) 
      V(NSAV,J)=V(IP,J) 
  970 CONTINUE 
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2)) 
      V(NSAV,5)=0. 
      DO 990 I=NSAV+1,N 
      DO 980 J=1,5 
      K(I,J)=K(I+NRS-1,J) 
      P(I,J)=P(I+NRS-1,J) 
      V(I,J)=0. 
  980 CONTINUE 
  990 CONTINUE 
      MSTU91=MSTU(90) 
      DO 1000 IZ=MSTU90+1,MSTU91 
      MSTU9T(IZ)=MSTU(90+IZ)-NRS+1-NSAV+N 
      PARU9T(IZ)=PARU(90+IZ) 
 1000 CONTINUE 
      MSTU(90)=MSTU90 
 
C...Order particles in rank along the chain. Update mother pointer. 
      DO 1020 I=NSAV+1,N 
      DO 1010 J=1,5 
      K(I-NSAV+N,J)=K(I,J) 
      P(I-NSAV+N,J)=P(I,J) 
 1010 CONTINUE 
 1020 CONTINUE 
      I1=NSAV 
      DO 1050 I=N+1,2*N-NSAV 
      IF(K(I,3).NE.IE(1)) GOTO 1050 
      I1=I1+1 
      DO 1030 J=1,5 
      K(I1,J)=K(I,J) 
      P(I1,J)=P(I,J) 
 1030 CONTINUE 
      IF(MSTU(16).NE.2) K(I1,3)=NSAV 
      DO 1040 IZ=MSTU90+1,MSTU91 
      IF(MSTU9T(IZ).EQ.I) THEN 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I1 
        PARU(90+MSTU(90))=PARU9T(IZ) 
      ENDIF 
 1040 CONTINUE 
 1050 CONTINUE 
      DO 1080 I=2*N-NSAV,N+1,-1 
      IF(K(I,3).EQ.IE(1)) GOTO 1080 
      I1=I1+1 
      DO 1060 J=1,5 
      K(I1,J)=K(I,J) 
      P(I1,J)=P(I,J) 
 1060 CONTINUE 
      IF(MSTU(16).NE.2) K(I1,3)=NSAV 
      DO 1070 IZ=MSTU90+1,MSTU91 
      IF(MSTU9T(IZ).EQ.I) THEN 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I1 
        PARU(90+MSTU(90))=PARU9T(IZ) 
      ENDIF 
 1070 CONTINUE 
 1080 CONTINUE 
 
C...Boost back particle system. Set production vertices. 
      IF(MBST.EQ.0) THEN 
        MSTU(33)=1 
        CALL LUDBRB(NSAV+1,N,0.,0.,DPS(1)/DPS(4),DPS(2)/DPS(4), 
     &  DPS(3)/DPS(4)) 
      ELSE 
        DO 1090 I=NSAV+1,N 
        HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2 
        IF(P(I,3).GT.0.) THEN 
          HHPEZ=(P(I,4)+P(I,3))*HHBZ 
          P(I,3)=0.5*(HHPEZ-HHPMT/HHPEZ) 
          P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ) 
        ELSE 
          HHPEZ=(P(I,4)-P(I,3))/HHBZ 
          P(I,3)=-0.5*(HHPEZ-HHPMT/HHPEZ) 
          P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ) 
        ENDIF 
 1090   CONTINUE 
      ENDIF 
      DO 1110 I=NSAV+1,N 
      DO 1100 J=1,4 
      V(I,J)=V(IP,J) 
 1100 CONTINUE 
 1110 CONTINUE 
 
      RETURN 
      END 
