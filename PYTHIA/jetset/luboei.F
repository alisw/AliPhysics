 
C********************************************************************* 
 
      SUBROUTINE LUBOEI(NSAV) 
 
C...Purpose: to modify event so as to approximately take into account 
C...Bose-Einstein effects according to a simple phenomenological 
C...parametrization. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUJETS/,/LUDAT1/ 
      DIMENSION DPS(4),KFBE(9),NBE(0:9),BEI(100) 
      DATA KFBE/211,-211,111,321,-321,130,310,221,331/ 
 
C...Boost event to overall CM frame. Calculate CM energy. 
      IF((MSTJ(51).NE.1.AND.MSTJ(51).NE.2).OR.N-NSAV.LE.1) RETURN 
      DO 100 J=1,4 
      DPS(J)=0. 
  100 CONTINUE 
      DO 120 I=1,N 
      KFA=IABS(K(I,2))
      IF(K(I,1).LE.10.AND.((KFA.GT.10.AND.KFA.LE.20).OR.KFA.EQ.22).AND.
     &K(I,3).GT.0) THEN
        KFMA=IABS(K(K(I,3),2))
        IF(KFMA.GT.10.AND.KFMA.LE.80) K(I,1)=-K(I,1)
      ENDIF
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 120 
      DO 110 J=1,4 
      DPS(J)=DPS(J)+P(I,J) 
  110 CONTINUE 
  120 CONTINUE 
      CALL LUDBRB(0,0,0.,0.,-DPS(1)/DPS(4),-DPS(2)/DPS(4), 
     &-DPS(3)/DPS(4)) 
      PECM=0. 
      DO 130 I=1,N 
      IF(K(I,1).GE.1.AND.K(I,1).LE.10) PECM=PECM+P(I,4) 
  130 CONTINUE 
 
C...Reserve copy of particles by species at end of record. 
      NBE(0)=N+MSTU(3) 
      DO 160 IBE=1,MIN(9,MSTJ(52)) 
      NBE(IBE)=NBE(IBE-1) 
      DO 150 I=NSAV+1,N 
      IF(K(I,2).NE.KFBE(IBE)) GOTO 150 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 150 
      IF(NBE(IBE).GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUBOEI:) no more memory left in LUJETS') 
        RETURN 
      ENDIF 
      NBE(IBE)=NBE(IBE)+1 
      K(NBE(IBE),1)=I 
      DO 140 J=1,3 
      P(NBE(IBE),J)=0. 
  140 CONTINUE 
  150 CONTINUE 
  160 CONTINUE 
      IF(NBE(MIN(9,MSTJ(52)))-NBE(0).LE.1) GOTO 280
 
C...Tabulate integral for subsequent momentum shift. 
      DO 220 IBE=1,MIN(9,MSTJ(52)) 
      IF(IBE.NE.1.AND.IBE.NE.4.AND.IBE.LE.7) GOTO 180 
      IF(IBE.EQ.1.AND.MAX(NBE(1)-NBE(0),NBE(2)-NBE(1),NBE(3)-NBE(2)) 
     &.LE.1) GOTO 180 
      IF(IBE.EQ.4.AND.MAX(NBE(4)-NBE(3),NBE(5)-NBE(4),NBE(6)-NBE(5), 
     &NBE(7)-NBE(6)).LE.1) GOTO 180 
      IF(IBE.GE.8.AND.NBE(IBE)-NBE(IBE-1).LE.1) GOTO 180 
      IF(IBE.EQ.1) PMHQ=2.*ULMASS(211) 
      IF(IBE.EQ.4) PMHQ=2.*ULMASS(321) 
      IF(IBE.EQ.8) PMHQ=2.*ULMASS(221) 
      IF(IBE.EQ.9) PMHQ=2.*ULMASS(331) 
      QDEL=0.1*MIN(PMHQ,PARJ(93)) 
      IF(MSTJ(51).EQ.1) THEN 
        NBIN=MIN(100,NINT(9.*PARJ(93)/QDEL)) 
        BEEX=EXP(0.5*QDEL/PARJ(93)) 
        BERT=EXP(-QDEL/PARJ(93)) 
      ELSE 
        NBIN=MIN(100,NINT(3.*PARJ(93)/QDEL)) 
      ENDIF 
      DO 170 IBIN=1,NBIN 
      QBIN=QDEL*(IBIN-0.5) 
      BEI(IBIN)=QDEL*(QBIN**2+QDEL**2/12.)/SQRT(QBIN**2+PMHQ**2) 
      IF(MSTJ(51).EQ.1) THEN 
        BEEX=BEEX*BERT 
        BEI(IBIN)=BEI(IBIN)*BEEX 
      ELSE 
        BEI(IBIN)=BEI(IBIN)*EXP(-(QBIN/PARJ(93))**2) 
      ENDIF 
      IF(IBIN.GE.2) BEI(IBIN)=BEI(IBIN)+BEI(IBIN-1) 
  170 CONTINUE 
 
C...Loop through particle pairs and find old relative momentum. 
  180 DO 210 I1M=NBE(IBE-1)+1,NBE(IBE)-1 
      I1=K(I1M,1) 
      DO 200 I2M=I1M+1,NBE(IBE) 
      I2=K(I2M,1) 
      Q2OLD=MAX(0.,(P(I1,4)+P(I2,4))**2-(P(I1,1)+P(I2,1))**2-(P(I1,2)+ 
     &P(I2,2))**2-(P(I1,3)+P(I2,3))**2-(P(I1,5)+P(I2,5))**2) 
      QOLD=SQRT(Q2OLD) 
 
C...Calculate new relative momentum. 
      IF(QOLD.LT.1E-3*QDEL) THEN 
        GOTO 200 
      ELSEIF(QOLD.LE.QDEL) THEN 
        QMOV=QOLD/3. 
      ELSEIF(QOLD.LT.(NBIN-0.1)*QDEL) THEN 
        RBIN=QOLD/QDEL 
        IBIN=RBIN 
        RINP=(RBIN**3-IBIN**3)/(3*IBIN*(IBIN+1)+1) 
        QMOV=(BEI(IBIN)+RINP*(BEI(IBIN+1)-BEI(IBIN)))* 
     &  SQRT(Q2OLD+PMHQ**2)/Q2OLD 
      ELSE 
        QMOV=BEI(NBIN)*SQRT(Q2OLD+PMHQ**2)/Q2OLD 
      ENDIF 
      Q2NEW=Q2OLD*(QOLD/(QOLD+3.*PARJ(92)*QMOV))**(2./3.) 
 
C...Calculate and save shift to be performed on three-momenta. 
      HC1=(P(I1,4)+P(I2,4))**2-(Q2OLD-Q2NEW) 
      HC2=(Q2OLD-Q2NEW)*(P(I1,4)-P(I2,4))**2 
      HA=0.5*(1.-SQRT(HC1*Q2NEW/(HC1*Q2OLD-HC2))) 
      DO 190 J=1,3 
      PD=HA*(P(I2,J)-P(I1,J)) 
      P(I1M,J)=P(I1M,J)+PD 
      P(I2M,J)=P(I2M,J)-PD 
  190 CONTINUE 
  200 CONTINUE 
  210 CONTINUE 
  220 CONTINUE 
 
C...Shift momenta and recalculate energies. 
      DO 240 IM=NBE(0)+1,NBE(MIN(9,MSTJ(52))) 
      I=K(IM,1) 
      DO 230 J=1,3 
      P(I,J)=P(I,J)+P(IM,J) 
  230 CONTINUE 
      P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  240 CONTINUE 
 
C...Rescale all momenta for energy conservation. 
      PES=0. 
      PQS=0. 
      DO 250 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 250 
      PES=PES+P(I,4) 
      PQS=PQS+P(I,5)**2/P(I,4) 
  250 CONTINUE 
      FAC=(PECM-PQS)/(PES-PQS) 
      DO 270 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 270 
      DO 260 J=1,3 
      P(I,J)=FAC*P(I,J) 
  260 CONTINUE 
      P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  270 CONTINUE 
 
C...Boost back to correct reference frame. 
  280 CALL LUDBRB(0,0,0.,0.,DPS(1)/DPS(4),DPS(2)/DPS(4),DPS(3)/DPS(4)) 
      DO 290 I=1,N
      IF(K(I,1).LT.0) K(I,1)=-K(I,1)
  290 CONTINUE
 
      RETURN 
      END 
