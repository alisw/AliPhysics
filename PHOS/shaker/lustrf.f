*CMZ :          17/07/98  15.44.31  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUSTRF(IP)
C...Purpose: to handle the fragmentation of an arbitrary colour singlet
C...jet system according to the Lund string fragmentation model.
      IMPLICIT DOUBLE PRECISION(D)
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.
      DIMENSION DPS(5),KFL(3),PMQ(3),PX(3),PY(3),GAM(3),IE(2),PR(2),
     &IN(9),DHM(4),DHG(4),DP(5,5),IRANK(2),MJU(4),IJU(3),PJU(5,5),
     &TJU(5),KFJH(2),NJS(2),KFJS(2),PJS(4,5)

C...Function: four-product of two vectors.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
      DFOUR(I,J)=DP(I,4)*DP(J,4)-DP(I,1)*DP(J,1)-DP(I,2)*DP(J,2)-
     &DP(I,3)*DP(J,3)

C...Reset counters. Identify parton system.
      MSTJ(91)=0
      NSAV=N
      NP=0
      KQSUM=0
      DO 100 J=1,5
  100 DPS(J)=0.
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
  120 DPS(J)=DPS(J)+P(I,J)
      K(N+NP,3)=I
      IF(P(N+NP,4)**2.LT.P(N+NP,1)**2+P(N+NP,2)**2+P(N+NP,3)**2) THEN
        P(N+NP,4)=SQRT(P(N+NP,1)**2+P(N+NP,2)**2+P(N+NP,3)**2+
     &  P(N+NP,5)**2)
        DPS(4)=DPS(4)+MAX(0.,P(N+NP,4)-P(I,4))
      ENDIF
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
      MSTU(33)=1
      CALL LUDBRB(N+1,N+NP,0.,0.,-DPS(1)/DPS(4),-DPS(2)/DPS(4),
     &-DPS(3)/DPS(4))

C...Search for very nearby partons that may be recombined.
      NTRYR=0
      PARU12=PARU(12)
      PARU13=PARU(13)
      MJU(3)=MJU(1)
      MJU(4)=MJU(2)
      NR=NP
  130 IF(NR.GE.3) THEN
        PDRMIN=2.*PARU12
        DO 140 I=N+1,N+NR
        IF(I.EQ.N+NR.AND.IABS(K(N+1,2)).NE.21) GOTO 140
        I1=I+1
        IF(I.EQ.N+NR) I1=N+1
        IF(K(I,1).EQ.41.OR.K(I1,1).EQ.41) GOTO 140
        IF(MJU(1).NE.0.AND.I1.LT.MJU(1).AND.IABS(K(I1,2)).NE.21)
     &  GOTO 140
        IF(MJU(2).NE.0.AND.I.GT.MJU(2).AND.IABS(K(I,2)).NE.21) GOTO 140
        PAP=SQRT((P(I,1)**2+P(I,2)**2+P(I,3)**2)*(P(I1,1)**2+
     &  P(I1,2)**2+P(I1,3)**2))
        PVP=P(I,1)*P(I1,1)+P(I,2)*P(I1,2)+P(I,3)*P(I1,3)
        PDR=4.*(PAP-PVP)**2/(PARU13**2*PAP+2.*(PAP-PVP))
        IF(PDR.LT.PDRMIN) THEN
          IR=I
          PDRMIN=PDR
        ENDIF
  140   CONTINUE

C...Recombine very nearby partons to avoid machine precision problems.
        IF(PDRMIN.LT.PARU12.AND.IR.EQ.N+NR) THEN
          DO 150 J=1,4
  150     P(N+1,J)=P(N+1,J)+P(N+NR,J)
          P(N+1,5)=SQRT(MAX(0.,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &    P(N+1,3)**2))
          NR=NR-1
          GOTO 130
        ELSEIF(PDRMIN.LT.PARU12) THEN
          DO 160 J=1,4
  160     P(IR,J)=P(IR,J)+P(IR+1,J)
          P(IR,5)=SQRT(MAX(0.,P(IR,4)**2-P(IR,1)**2-P(IR,2)**2-
     &    P(IR,3)**2))
          DO 170 I=IR+1,N+NR-1
          K(I,2)=K(I+1,2)
          DO 170 J=1,5
  170     P(I,J)=P(I+1,J)
          IF(IR.EQ.N+NR-1) K(IR,2)=K(N+NR,2)
          NR=NR-1
          IF(MJU(1).GT.IR) MJU(1)=MJU(1)-1
          IF(MJU(2).GT.IR) MJU(2)=MJU(2)-1
          GOTO 130
        ENDIF
      ENDIF
      NTRYR=NTRYR+1

C...Reset particle counter. Skip ahead if no junctions are present;
C...this is usually the case!
      NRS=MAX(5*NR+11,NP)
      NTRY=0
  180 NTRY=NTRY+1
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN
        PARU12=4.*PARU12
        PARU13=2.*PARU13
        GOTO 130
      ELSEIF(NTRY.GT.100) THEN
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=N+NRS
      IF(MJU(1).EQ.0.AND.MJU(2).EQ.0) GOTO 500
      DO 490 JT=1,2
      NJS(JT)=0
      IF(MJU(JT).EQ.0) GOTO 490
      JS=3-2*JT

C...Find and sum up momentum on three sides of junction. Check flavours.
      DO 190 IU=1,3
      IJU(IU)=0
      DO 190 J=1,5
  190 PJU(IU,J)=0.
      IU=0
      DO 200 I1=N+1+(JT-1)*(NR-1),N+NR+(JT-1)*(1-NR),JS
      IF(K(I1,2).NE.21.AND.IU.LE.2) THEN
        IU=IU+1
        IJU(IU)=I1
      ENDIF
      DO 200 J=1,4
  200 PJU(IU,J)=PJU(IU,J)+P(I1,J)
      DO 210 IU=1,3
  210 PJU(IU,5)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,3)**2)
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
      DO 220 J=1,3
  220 TJU(J)=-(T1F*PJU(1,J)/PJU(1,5)+T2F*PJU(2,J)/PJU(2,5))
      TJU(4)=SQRT(1.+TJU(1)**2+TJU(2)**2+TJU(3)**2)
      DO 230 IU=1,3
  230 PJU(IU,5)=TJU(4)*PJU(IU,4)-TJU(1)*PJU(IU,1)-TJU(2)*PJU(IU,2)-
     &TJU(3)*PJU(IU,3)

C...Put junction at rest if motion could give inconsistencies.
      IF(PJU(1,5)+PJU(2,5).GT.PJU(1,4)+PJU(2,4)) THEN
        DO 240 J=1,3
  240   TJU(J)=0.
        TJU(4)=1.
        PJU(1,5)=PJU(1,4)
        PJU(2,5)=PJU(2,4)
        PJU(3,5)=PJU(3,4)
      ENDIF

C...Start preparing for fragmentation of two strings from junction.
      ISTA=I
      DO 470 IU=1,2
      NS=IJU(IU+1)-IJU(IU)

C...Junction strings: find longitudinal string directions.
      DO 260 IS=1,NS
      IS1=IJU(IU)+IS-1
      IS2=IJU(IU)+IS
      DO 250 J=1,5
      DP(1,J)=0.5*P(IS1,J)
      IF(IS.EQ.1) DP(1,J)=P(IS1,J)
      DP(2,J)=0.5*P(IS2,J)
  250 IF(IS.EQ.NS) DP(2,J)=-PJU(IU,J)
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
      DO 260 J=1,4
      P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J)
  260 P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J)

C...Junction strings: initialize flavour, momentum and starting pos.
      ISAV=I
  270 NTRY=NTRY+1
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN
        PARU12=4.*PARU12
        PARU13=2.*PARU13
        GOTO 130
      ELSEIF(NTRY.GT.100) THEN
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=ISAV
      IRANKJ=0
      IE(1)=K(N+1+(JT/2)*(NP-1),3)
      IN(4)=N+NR+1
      IN(5)=IN(4)+1
      IN(6)=N+NR+4*NS+1
      DO 280 JQ=1,2
      DO 280 IN1=N+NR+2+JQ,N+NR+4*NS-2+JQ,4
      P(IN1,1)=2-JQ
      P(IN1,2)=JQ-1
  280 P(IN1,3)=1.
      KFL(1)=K(IJU(IU),2)
      PX(1)=0.
      PY(1)=0.
      GAM(1)=0.
      DO 290 J=1,5
  290 PJU(IU+3,J)=0.

C...Junction strings: find initial transverse directions.
      DO 300 J=1,4
      DP(1,J)=P(IN(4),J)
      DP(2,J)=P(IN(4)+1,J)
      DP(3,J)=0.
  300 DP(4,J)=0.
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
      DO 310 J=1,4
      DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
      P(IN(6),J)=DP(3,J)
  310 P(IN(6)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &DHCYX*DP(3,J))

C...Junction strings: produce new particle, origin.
  320 I=I+1
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
  330 CALL LUKFDI(KFL(1),0,KFL(3),K(I,2))
      IF(K(I,2).EQ.0) GOTO 270
      IF(MSTJ(12).GE.3.AND.IRANKJ.EQ.1.AND.IABS(KFL(1)).LE.10.AND.
     &IABS(KFL(3)).GT.10) THEN
        IF(RLU(0).GT.PARJ(19)) GOTO 330
      ENDIF
      P(I,5)=ULMASS(K(I,2))
      CALL LUPTDI(KFL(1),PX(3),PY(3))
      PR(1)=P(I,5)**2+(PX(1)+PX(3))**2+(PY(1)+PY(3))**2
      CALL LUZDIS(KFL(1),KFL(3),PR(1),Z)
      GAM(3)=(1.-Z)*(GAM(1)+PR(1)/Z)
      DO 340 J=1,3
  340 IN(J)=IN(3+J)

C...Junction strings: stepping within or from 'low' string region easy.
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)*
     &P(IN(1),5)**2.GE.PR(1)) THEN
        P(IN(1)+2,4)=Z*P(IN(1)+2,3)
        P(IN(2)+2,4)=PR(1)/(P(IN(1)+2,4)*P(IN(1),5)**2)
        DO 350 J=1,4
  350   P(I,J)=(PX(1)+PX(3))*P(IN(3),J)+(PY(1)+PY(3))*P(IN(3)+1,J)
        GOTO 420
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN
        P(IN(2)+2,4)=P(IN(2)+2,3)
        P(IN(2)+2,1)=1.
        IN(2)=IN(2)+4
        IF(IN(2).GT.N+NR+4*NS) GOTO 270
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(1)+2,4)=P(IN(1)+2,3)
          P(IN(1)+2,1)=0.
          IN(1)=IN(1)+4
        ENDIF
      ENDIF

C...Junction strings: find new transverse directions.
  360 IF(IN(1).GT.N+NR+4*NS.OR.IN(2).GT.N+NR+4*NS.OR.
     &IN(1).GT.IN(2)) GOTO 270
      IF(IN(1).NE.IN(4).OR.IN(2).NE.IN(5)) THEN
        DO 370 J=1,4
        DP(1,J)=P(IN(1),J)
        DP(2,J)=P(IN(2),J)
        DP(3,J)=0.
  370   DP(4,J)=0.
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
        DHC12=DFOUR(1,2)
        IF(DHC12.LE.1E-2) THEN
          P(IN(1)+2,4)=P(IN(1)+2,3)
          P(IN(1)+2,1)=0.
          IN(1)=IN(1)+4
          GOTO 360
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
        DO 380 J=1,4
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
        P(IN(3),J)=DP(3,J)
  380   P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &  DHCYX*DP(3,J))
C...Express pT with respect to new axes, if sensible.
        PXP=-(PX(3)*FOUR(IN(6),IN(3))+PY(3)*FOUR(IN(6)+1,IN(3)))
        PYP=-(PX(3)*FOUR(IN(6),IN(3)+1)+PY(3)*FOUR(IN(6)+1,IN(3)+1))
        IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01) THEN
          PX(3)=PXP
          PY(3)=PYP
        ENDIF
      ENDIF

C...Junction strings: sum up known four-momentum, coefficients for m2.
      DO 400 J=1,4
      DHG(J)=0.
      P(I,J)=PX(1)*P(IN(6),J)+PY(1)*P(IN(6)+1,J)+PX(3)*P(IN(3),J)+
     &PY(3)*P(IN(3)+1,J)
      DO 390 IN1=IN(4),IN(1)-4,4
  390 P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
      DO 400 IN2=IN(5),IN(2)-4,4
  400 P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
      DHM(1)=FOUR(I,I)
      DHM(2)=2.*FOUR(I,IN(1))
      DHM(3)=2.*FOUR(I,IN(2))
      DHM(4)=2.*FOUR(IN(1),IN(2))

C...Junction strings: find coefficients for Gamma expression.
      DO 410 IN2=IN(1)+1,IN(2),4
      DO 410 IN1=IN(1),IN2-1,4
      DHC=2.*FOUR(IN1,IN2)
      DHG(1)=DHG(1)+P(IN1+2,1)*P(IN2+2,1)*DHC
      IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-P(IN2+2,1)*DHC
      IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+P(IN1+2,1)*DHC
  410 IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC

C...Junction strings: solve (m2, Gamma) equation system for energies.
      DHS1=DHM(3)*DHG(4)-DHM(4)*DHG(3)
      IF(ABS(DHS1).LT.1E-4) GOTO 270
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(2)*DHG(3)-DHG(4)*
     &(P(I,5)**2-DHM(1))+DHG(2)*DHM(3)
      DHS3=DHM(2)*(GAM(3)-DHG(1))-DHG(2)*(P(I,5)**2-DHM(1))
      P(IN(2)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/ABS(DHS1)-
     &DHS2/DHS1)
      IF(DHM(2)+DHM(4)*P(IN(2)+2,4).LE.0.) GOTO 270
      P(IN(1)+2,4)=(P(I,5)**2-DHM(1)-DHM(3)*P(IN(2)+2,4))/
     &(DHM(2)+DHM(4)*P(IN(2)+2,4))

C...Junction strings: step to new region if necessary.
      IF(P(IN(2)+2,4).GT.P(IN(2)+2,3)) THEN
        P(IN(2)+2,4)=P(IN(2)+2,3)
        P(IN(2)+2,1)=1.
        IN(2)=IN(2)+4
        IF(IN(2).GT.N+NR+4*NS) GOTO 270
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(1)+2,4)=P(IN(1)+2,3)
          P(IN(1)+2,1)=0.
          IN(1)=IN(1)+4
        ENDIF
        GOTO 360
      ELSEIF(P(IN(1)+2,4).GT.P(IN(1)+2,3)) THEN
        P(IN(1)+2,4)=P(IN(1)+2,3)
        P(IN(1)+2,1)=0.
        IN(1)=IN(1)+JS
        GOTO 710
      ENDIF

C...Junction strings: particle four-momentum, remainder, loop back.
  420 DO 430 J=1,4
      P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J)
  430 PJU(IU+3,J)=PJU(IU+3,J)+P(I,J)
      IF(P(I,4).LE.0.) GOTO 270
      PJU(IU+3,5)=TJU(4)*PJU(IU+3,4)-TJU(1)*PJU(IU+3,1)-
     &TJU(2)*PJU(IU+3,2)-TJU(3)*PJU(IU+3,3)
      IF(PJU(IU+3,5).LT.PJU(IU,5)) THEN
        KFL(1)=-KFL(3)
        PX(1)=-PX(3)
        PY(1)=-PY(3)
        GAM(1)=GAM(3)
        IF(IN(3).NE.IN(6)) THEN
          DO 440 J=1,4
          P(IN(6),J)=P(IN(3),J)
  440     P(IN(6)+1,J)=P(IN(3)+1,J)
        ENDIF
        DO 450 JQ=1,2
        IN(3+JQ)=IN(JQ)
        P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
  450   P(IN(JQ)+2,1)=P(IN(JQ)+2,1)-(3-2*JQ)*P(IN(JQ)+2,4)
        GOTO 320
      ENDIF

C...Junction strings: save quantities left after each string.
      IF(IABS(KFL(1)).GT.10) GOTO 270
      I=I-1
      KFJH(IU)=KFL(1)
      DO 460 J=1,4
  460 PJU(IU+3,J)=PJU(IU+3,J)-P(I+1,J)
  470 CONTINUE

C...Junction strings: put together to new effective string endpoint.
      NJS(JT)=I-ISTA
      KFJS(JT)=K(K(MJU(JT+2),3),2)
      KFLS=2*INT(RLU(0)+3.*PARJ(4)/(1.+3.*PARJ(4)))+1
      IF(KFJH(1).EQ.KFJH(2)) KFLS=3
      IF(ISTA.NE.I) KFJS(JT)=ISIGN(1000*MAX(IABS(KFJH(1)),
     &IABS(KFJH(2)))+100*MIN(IABS(KFJH(1)),IABS(KFJH(2)))+
     &KFLS,KFJH(1))
      DO 480 J=1,4
      PJS(JT,J)=PJU(1,J)+PJU(2,J)+P(MJU(JT),J)
  480 PJS(JT+2,J)=PJU(4,J)+PJU(5,J)
      PJS(JT,5)=SQRT(MAX(0.,PJS(JT,4)**2-PJS(JT,1)**2-PJS(JT,2)**2-
     &PJS(JT,3)**2))
  490 CONTINUE

C...Open versus closed strings. Choose breakup region for latter.
  500 IF(MJU(1).NE.0.AND.MJU(2).NE.0) THEN
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
        DO 510 IS=1,NR
        P(N+NR+IS,1)=0.5*FOUR(N+IS,N+IS+1-NR*(IS/NR))
  510   W2SUM=W2SUM+P(N+NR+IS,1)
        W2RAN=RLU(0)*W2SUM
        NB=0
  520   NB=NB+1
        W2SUM=W2SUM-P(N+NR+NB,1)
        IF(W2SUM.GT.W2RAN.AND.NB.LT.NR) GOTO 520
      ENDIF

C...Find longitudinal string directions (i.e. lightlike four-vectors).
      DO 540 IS=1,NS
      IS1=N+IS+NB-1-NR*((IS+NB-2)/NR)
      IS2=N+IS+NB-NR*((IS+NB-1)/NR)
      DO 530 J=1,5
      DP(1,J)=P(IS1,J)
      IF(IABS(K(IS1,2)).EQ.21) DP(1,J)=0.5*DP(1,J)
      IF(IS1.EQ.MJU(1)) DP(1,J)=PJS(1,J)-PJS(3,J)
      DP(2,J)=P(IS2,J)
      IF(IABS(K(IS2,2)).EQ.21) DP(2,J)=0.5*DP(2,J)
  530 IF(IS2.EQ.MJU(2)) DP(2,J)=PJS(2,J)-PJS(4,J)
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
      DO 540 J=1,4
      P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J)
  540 P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J)

C...Begin initialization: sum up energy, set starting position.
      ISAV=I
  550 NTRY=NTRY+1
      IF(NTRY.GT.100.AND.NTRYR.LE.4) THEN
        PARU12=4.*PARU12
        PARU13=2.*PARU13
        GOTO 130
      ELSEIF(NTRY.GT.100) THEN
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=ISAV
      DO 560 J=1,4
      P(N+NRS,J)=0.
      DO 560 IS=1,NR
  560 P(N+NRS,J)=P(N+NRS,J)+P(N+IS,J)
      DO 570 JT=1,2
      IRANK(JT)=0
      IF(MJU(JT).NE.0) IRANK(JT)=NJS(JT)
      IF(NS.GT.NR) IRANK(JT)=1
      IE(JT)=K(N+1+(JT/2)*(NP-1),3)
      IN(3*JT+1)=N+NR+1+4*(JT/2)*(NS-1)
      IN(3*JT+2)=IN(3*JT+1)+1
      IN(3*JT+3)=N+NR+4*NS+2*JT-1
      DO 570 IN1=N+NR+2+JT,N+NR+4*NS-2+JT,4
      P(IN1,1)=2-JT
      P(IN1,2)=JT-1
  570 P(IN1,3)=1.

C...Initialize flavour and pT variables for open string.
      IF(NS.LT.NR) THEN
        PX(1)=0.
        PY(1)=0.
        IF(NS.EQ.1.AND.MJU(1)+MJU(2).EQ.0) CALL LUPTDI(0,PX(1),PY(1))
        PX(2)=-PX(1)
        PY(2)=-PY(1)
        DO 580 JT=1,2
        KFL(JT)=K(IE(JT),2)
        IF(MJU(JT).NE.0) KFL(JT)=KFJS(JT)
        MSTJ(93)=1
        PMQ(JT)=ULMASS(KFL(JT))
  580   GAM(JT)=0.

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
  590   CALL LUZDIS(KFL(1),KFL(2),PR3,Z)
        ZR=PR3/(Z*P(N+NR+1,5)**2)
        IF(ZR.GE.1.) GOTO 590
        DO 600 JT=1,2
        MSTJ(93)=1
        PMQ(JT)=ULMASS(KFL(JT))
        GAM(JT)=PR3*(1.-Z)/Z
        IN1=N+NR+3+4*(JT/2)*(NS-1)
        P(IN1,JT)=1.-Z
        P(IN1,3-JT)=JT-1
        P(IN1,3)=(2-JT)*(1.-Z)+(JT-1)*Z
        P(IN1+1,JT)=ZR
        P(IN1+1,3-JT)=2-JT
  600   P(IN1+1,3)=(2-JT)*(1.-ZR)+(JT-1)*ZR
      ENDIF

C...Find initial transverse directions (i.e. spacelike four-vectors).
      DO 640 JT=1,2
      IF(JT.EQ.1.OR.NS.EQ.NR-1) THEN
        IN1=IN(3*JT+1)
        IN3=IN(3*JT+3)
        DO 610 J=1,4
        DP(1,J)=P(IN1,J)
        DP(2,J)=P(IN1+1,J)
        DP(3,J)=0.
  610   DP(4,J)=0.
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
        DO 620 J=1,4
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
        P(IN3,J)=DP(3,J)
  620   P(IN3+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &  DHCYX*DP(3,J))
      ELSE
        DO 630 J=1,4
        P(IN3+2,J)=P(IN3,J)
  630   P(IN3+3,J)=P(IN3+1,J)
      ENDIF
  640 CONTINUE

C...Remove energy used up in junction string fragmentation.
      IF(MJU(1)+MJU(2).GT.0) THEN
        DO 660 JT=1,2
        IF(NJS(JT).EQ.0) GOTO 660
        DO 650 J=1,4
  650   P(N+NRS,J)=P(N+NRS,J)-PJS(JT+2,J)
  660   CONTINUE
      ENDIF

C...Produce new particle: side, origin.
  670 I=I+1
      IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      JT=1.5+RLU(0)
      IF(IABS(KFL(3-JT)).GT.10) JT=3-JT
      JR=3-JT
      JS=3-2*JT
      IRANK(JT)=IRANK(JT)+1
      K(I,1)=1
      K(I,3)=IE(JT)
      K(I,4)=0
      K(I,5)=0

C...Generate flavour, hadron and pT.
  680 CALL LUKFDI(KFL(JT),0,KFL(3),K(I,2))
      IF(K(I,2).EQ.0) GOTO 550
      IF(MSTJ(12).GE.3.AND.IRANK(JT).EQ.1.AND.IABS(KFL(JT)).LE.10.AND.
     &IABS(KFL(3)).GT.10) THEN
        IF(RLU(0).GT.PARJ(19)) GOTO 680
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
      IF(WREM2.LT.0.10) GOTO 550
      IF(WREM2.LT.MAX(WMIN*(1.+(2.*RLU(0)-1.)*PARJ(37)),
     &PARJ(32)+PMQ(1)+PMQ(2))**2) GOTO 810

C...Choose z, which gives Gamma. Shift z for heavy flavours.
      CALL LUZDIS(KFL(JT),KFL(3),PR(JT),Z)
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     &MOD(KFL2A/1000,10)).GE.4) THEN
        PR(JR)=(PMQ(JR)+PMQ(3))**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        PW12=SQRT(MAX(0.,(WREM2-PR(1)-PR(2))**2-4.*PR(1)*PR(2)))
        Z=(WREM2+PR(JT)-PR(JR)+PW12*(2.*Z-1.))/(2.*WREM2)
        PR(JR)=(PMQ(JR)+PARJST)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        IF((1.-Z)*(WREM2-PR(JT)/Z).LT.PR(JR)) GOTO 810
      ENDIF
      GAM(3)=(1.-Z)*(GAM(JT)+PR(JT)/Z)
      DO 690 J=1,3
  690 IN(J)=IN(3*JT+J)

C...Stepping within or from 'low' string region easy.
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)*
     &P(IN(1),5)**2.GE.PR(JT)) THEN
        P(IN(JT)+2,4)=Z*P(IN(JT)+2,3)
        P(IN(JR)+2,4)=PR(JT)/(P(IN(JT)+2,4)*P(IN(1),5)**2)
        DO 700 J=1,4
  700   P(I,J)=(PX(JT)+PX(3))*P(IN(3),J)+(PY(JT)+PY(3))*P(IN(3)+1,J)
        GOTO 770
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1.
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 550
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
        ENDIF
      ENDIF

C...Find new transverse directions (i.e. spacelike string vectors).
  710 IF(JS*IN(1).GT.JS*IN(3*JR+1).OR.JS*IN(2).GT.JS*IN(3*JR+2).OR.
     &IN(1).GT.IN(2)) GOTO 550
      IF(IN(1).NE.IN(3*JT+1).OR.IN(2).NE.IN(3*JT+2)) THEN
        DO 720 J=1,4
        DP(1,J)=P(IN(1),J)
        DP(2,J)=P(IN(2),J)
        DP(3,J)=0.
  720   DP(4,J)=0.
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
        DHC12=DFOUR(1,2)
        IF(DHC12.LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
          GOTO 710
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
        DO 730 J=1,4
        DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
        P(IN(3),J)=DP(3,J)
  730   P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &  DHCYX*DP(3,J))
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
      DO 750 J=1,4
      DHG(J)=0.
      P(I,J)=PX(JT)*P(IN(3*JT+3),J)+PY(JT)*P(IN(3*JT+3)+1,J)+
     &PX(3)*P(IN(3),J)+PY(3)*P(IN(3)+1,J)
      DO 740 IN1=IN(3*JT+1),IN(1)-4*JS,4*JS
  740 P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
      DO 750 IN2=IN(3*JT+2),IN(2)-4*JS,4*JS
  750 P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
      DHM(1)=FOUR(I,I)
      DHM(2)=2.*FOUR(I,IN(1))
      DHM(3)=2.*FOUR(I,IN(2))
      DHM(4)=2.*FOUR(IN(1),IN(2))

C...Find coefficients for Gamma expression.
      DO 760 IN2=IN(1)+1,IN(2),4
      DO 760 IN1=IN(1),IN2-1,4
      DHC=2.*FOUR(IN1,IN2)
      DHG(1)=DHG(1)+P(IN1+2,JT)*P(IN2+2,JT)*DHC
      IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-JS*P(IN2+2,JT)*DHC
      IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+JS*P(IN1+2,JT)*DHC
  760 IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC

C...Solve (m2, Gamma) equation system for energies taken.
      DHS1=DHM(JR+1)*DHG(4)-DHM(4)*DHG(JR+1)
      IF(ABS(DHS1).LT.1E-4) GOTO 550
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(JT+1)*DHG(JR+1)-DHG(4)*
     &(P(I,5)**2-DHM(1))+DHG(JT+1)*DHM(JR+1)
      DHS3=DHM(JT+1)*(GAM(3)-DHG(1))-DHG(JT+1)*(P(I,5)**2-DHM(1))
      P(IN(JR)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/ABS(DHS1)-
     &DHS2/DHS1)
      IF(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4).LE.0.) GOTO 550
      P(IN(JT)+2,4)=(P(I,5)**2-DHM(1)-DHM(JR+1)*P(IN(JR)+2,4))/
     &(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4))

C...Step to new region if necessary.
      IF(P(IN(JR)+2,4).GT.P(IN(JR)+2,3)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1.
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 550
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
        ENDIF
        GOTO 710
      ELSEIF(P(IN(JT)+2,4).GT.P(IN(JT)+2,3)) THEN
        P(IN(JT)+2,4)=P(IN(JT)+2,3)
        P(IN(JT)+2,JT)=0.
        IN(JT)=IN(JT)+4*JS
        GOTO 710
      ENDIF

C...Four-momentum of particle. Remaining quantities. Loop back.
  770 DO 780 J=1,4
      P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J)
  780 P(N+NRS,J)=P(N+NRS,J)-P(I,J)
      IF(P(I,4).LE.0.) GOTO 550
      KFL(JT)=-KFL(3)
      PMQ(JT)=PMQ(3)
      PX(JT)=-PX(3)
      PY(JT)=-PY(3)
      GAM(JT)=GAM(3)
      IF(IN(3).NE.IN(3*JT+3)) THEN
        DO 790 J=1,4
        P(IN(3*JT+3),J)=P(IN(3),J)
  790   P(IN(3*JT+3)+1,J)=P(IN(3)+1,J)
      ENDIF
      DO 800 JQ=1,2
      IN(3*JT+JQ)=IN(JQ)
      P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
  800 P(IN(JQ)+2,JT)=P(IN(JQ)+2,JT)-JS*(3-2*JQ)*P(IN(JQ)+2,4)
      GOTO 670

C...Final hadron: side, flavour, hadron, mass.
  810 I=I+1
      K(I,1)=1
      K(I,3)=IE(JR)
      K(I,4)=0
      K(I,5)=0
      CALL LUKFDI(KFL(JR),-KFL(3),KFLDMP,K(I,2))
      IF(K(I,2).EQ.0) GOTO 550
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
      IF(MJU(1)+MJU(2).NE.0.AND.I.EQ.ISAV+2.AND.FD.GE.1.) GOTO 180
      IF(FD.GE.1.) GOTO 550
      FA=WREM2+PR(JT)-PR(JR)
      IF(MSTJ(11).NE.2) PREV=0.5*EXP(MAX(-80.,LOG(FD)*PARJ(38)*
     &(PR(1)+PR(2))**2))
      IF(MSTJ(11).EQ.2) PREV=0.5*FD**PARJ(39)
      FB=SIGN(SQRT(MAX(0.,FA**2-4.*WREM2*PR(JT))),JS*(RLU(0)-PREV))
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     &MOD(KFL2A/1000,10)).GE.6) FB=SIGN(SQRT(MAX(0.,FA**2-
     &4.*WREM2*PR(JT))),FLOAT(JS))
      DO 820 J=1,4
      P(I-1,J)=(PX(JT)+PX(3))*P(IN(3*JQ+3),J)+(PY(JT)+PY(3))*
     &P(IN(3*JQ+3)+1,J)+0.5*(DHR1*(FA+FB)*P(IN(3*JQ+1),J)+
     &DHR2*(FA-FB)*P(IN(3*JQ+2),J))/WREM2
  820 P(I,J)=P(N+NRS,J)-P(I-1,J)

C...Mark jets as fragmented and give daughter pointers.
      N=I-NRS+1
      DO 830 I=NSAV+1,NSAV+NP
      IM=K(I,3)
      K(IM,1)=K(IM,1)+10
      IF(MSTU(16).NE.2) THEN
        K(IM,4)=NSAV+1
        K(IM,5)=NSAV+1
      ELSE
        K(IM,4)=NSAV+2
        K(IM,5)=N
      ENDIF
  830 CONTINUE

C...Document string system. Move up particles.
      NSAV=NSAV+1
      K(NSAV,1)=11
      K(NSAV,2)=92
      K(NSAV,3)=IP
      K(NSAV,4)=NSAV+1
      K(NSAV,5)=N
      DO 840 J=1,4
      P(NSAV,J)=DPS(J)
  840 V(NSAV,J)=V(IP,J)
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))
      V(NSAV,5)=0.
      DO 850 I=NSAV+1,N
      DO 850 J=1,5
      K(I,J)=K(I+NRS-1,J)
      P(I,J)=P(I+NRS-1,J)
  850 V(I,J)=0.

C...Order particles in rank along the chain. Update mother pointer.
      DO 860 I=NSAV+1,N
      DO 860 J=1,5
      K(I-NSAV+N,J)=K(I,J)
  860 P(I-NSAV+N,J)=P(I,J)
      I1=NSAV
      DO 880 I=N+1,2*N-NSAV
      IF(K(I,3).NE.IE(1)) GOTO 880
      I1=I1+1
      DO 870 J=1,5
      K(I1,J)=K(I,J)
  870 P(I1,J)=P(I,J)
      IF(MSTU(16).NE.2) K(I1,3)=NSAV
  880 CONTINUE
      DO 900 I=2*N-NSAV,N+1,-1
      IF(K(I,3).EQ.IE(1)) GOTO 900
      I1=I1+1
      DO 890 J=1,5
      K(I1,J)=K(I,J)
  890 P(I1,J)=P(I,J)
      IF(MSTU(16).NE.2) K(I1,3)=NSAV
  900 CONTINUE

C...Boost back particle system. Set production vertices.
      MSTU(33)=1
      CALL LUDBRB(NSAV+1,N,0.,0.,DPS(1)/DPS(4),DPS(2)/DPS(4),
     &DPS(3)/DPS(4))
      DO 910 I=NSAV+1,N
      DO 910 J=1,4
  910 V(I,J)=V(IP,J)

      RETURN
      END
