C $Id$
C
C
C
C
	SUBROUTINE QUENCH(JPJT,NTP)
	DIMENSION RDP(300),LQP(300),RDT(300),LQT(300)
#define BLANKET_SAVE
#include "hijcrdn.inc"
#include "hiparnt.inc"
C
#include "hijjet1.inc"
#include "hijjet2.inc"
#include "histrng.inc"
C
	SAVE
C
        BB=HINT1(19)                                            ! Uzhi
        PHI=HINT1(20)                                           ! Uzhi
        BBX=BB*COS(PHI)                                         ! Uzhi
        BBY=BB*SIN(PHI)                                         ! Uzhi
c
	IF(NTP.EQ.2) GO TO 400
	IF(NTP.EQ.3) GO TO 2000 
C*******************************************************
C Jet interaction for proj jet in the direction PHIP
C******************************************************
C
	IF(NFP(JPJT,7).NE.1) RETURN

	JP=JPJT
	DO 290 I=1,NPJ(JP)
	   PTJET0=SQRT(PJPX(JP,I)**2+PJPY(JP,I)**2)
	   IF(PTJET0.LE.HIPR1(11)) GO TO 290
	   PTOT=SQRT(PTJET0*PTJET0+PJPZ(JP,I)**2)
	   IF (IHPR2(50) .EQ. 1) THEN 
	      DEDX0 = HIPR1(14)*LOG10(PTJET0)/LOG10(5.)
	   ELSE
	      DEDX0 = HIPR1(14)
	   ENDIF
	   IF(PTOT.LT.HIPR1(8)) GO TO 290
	   PHIP=ULANGL_HIJING(PJPX(JP,I),PJPY(JP,I))
C******* find the wounded proj which can interact with jet***
	   KP=0
	   DO 100 I2=1,IHNT2(1)
	      IF(NFP(I2,5).NE.3 .OR. I2.EQ.JP) GO TO 100
	      DX=YP(1,I2)-YP(1,JP)
	      DY=YP(2,I2)-YP(2,JP)
	      PHI=ULANGL_HIJING(DX,DY)
	      DPHI=ABS(PHI-PHIP)
              IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI      ! Uzhi
	      IF(DPHI.GE.HIPR1(40)/2.0) GO TO 100
	      RD0=SQRT(DX*DX+DY*DY)
	      IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 100
	      KP=KP+1
	      LQP(KP)=I2
	      RDP(KP)=COS(DPHI)*RD0
 100	   CONTINUE
C*******	rearrange according decending rd************
	   DO 110 I2=1,KP-1
	      DO 110 J2=I2+1,KP
		 IF(RDP(I2).LT.RDP(J2)) GO TO 110
		 RD=RDP(I2)
		 LQ=LQP(I2)
		 RDP(I2)=RDP(J2)
		 LQP(I2)=LQP(J2)
		 RDP(J2)=RD
		 LQP(J2)=LQ
 110	      CONTINUE
C****** find wounded targ which can interact with jet********
	      KT=0
	      DO 120 I2=1,IHNT2(3)
		 IF(NFT(I2,5).NE.3) GO TO 120
		 DX=YT(1,I2)-YP(1,JP)-BBX
		 DY=YT(2,I2)-YP(2,JP)-BBY
		 PHI=ULANGL_HIJING(DX,DY)
		 DPHI=ABS(PHI-PHIP)
		 IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI ! Uzhi
		 IF(DPHI.GT.HIPR1(40)/2.0) GO TO 120
		 RD0=SQRT(DX*DX+DY*DY)
		 IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 120
		 KT=KT+1
		 LQT(KT)=I2
		 RDT(KT)=COS(DPHI)*RD0
 120	      CONTINUE
C*******	rearrange according decending rd************
	      DO 130 I2=1,KT-1
		 DO 130 J2=I2+1,KT
		    IF(RDT(I2).LT.RDT(J2)) GO TO 130
		    RD=RDT(I2)
		    LQ=LQT(I2)
		    RDT(I2)=RDT(J2)
		    LQT(I2)=LQT(J2)
		    RDT(J2)=RD
		    LQT(J2)=LQ
 130		 CONTINUE
		
		 MP=0
		 MT=0
		 R0=0.0
		 NQ=0
		 DP=0.0
		 PTOT=SQRT(PJPX(JP,I)**2+PJPY(JP,I)**2+PJPZ(JP,I)**2)
		 V1=PJPX(JP,I)/PTOT
		 V2=PJPY(JP,I)/PTOT
		 V3=PJPZ(JP,I)/PTOT

 200		 RN=RLU_HIJING(0)
 210		 IF(MT.GE.KT .AND. MP.GE.KP) GO TO 290
		 IF(MT.GE.KT) GO TO 220
		 IF(MP.GE.KP) GO TO 240
		 IF(RDP(MP+1).GT.RDT(MT+1)) GO TO 240
 220		 MP=MP+1
		 DRR=RDP(MP)-R0
		 IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 210
		 DP=DRR*DEDX0
		 IF(KFPJ(JP,I).NE.21) DP=0.5*DP
C	********string tension of quark jet is 0.5 of gluon's 
		 IF(DP.LE.0.2) GO TO 210
		 IF(PTOT.LE.0.4) GO TO 290
		 IF(PTOT.LE.DP) DP=PTOT-0.2
		 DE=DP

		 IF(KFPJ(JP,I).NE.21) THEN
		    PRSHU=PP(LQP(MP),1)**2+PP(LQP(MP),2)**2
     &                   +PP(LQP(MP),3)**2
		    DE=SQRT(PJPM(JP,I)**2+PTOT**2)
     &			-SQRT(PJPM(JP,I)**2+(PTOT-DP)**2)
		    ERSHU=(PP(LQP(MP),4)+DE-DP)**2
		    AMSHU=ERSHU-PRSHU
		    IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 210
		    PP(LQP(MP),4)=SQRT(ERSHU)
		    PP(LQP(MP),5)=SQRT(AMSHU)
		 ENDIF
C		********reshuffle the energy when jet has mass
		 R0=RDP(MP)
		 DP1=DP*V1
		 DP2=DP*V2
		 DP3=DP*V3
C		********momentum and energy transfer from jet
		 
		 NPJ(LQP(MP))=NPJ(LQP(MP))+1
		 KFPJ(LQP(MP),NPJ(LQP(MP)))=21
		 PJPX(LQP(MP),NPJ(LQP(MP)))=DP1
		 PJPY(LQP(MP),NPJ(LQP(MP)))=DP2
		 PJPZ(LQP(MP),NPJ(LQP(MP)))=DP3
		 PJPE(LQP(MP),NPJ(LQP(MP)))=DP
		 PJPM(LQP(MP),NPJ(LQP(MP)))=0.0
		 GO TO 260

 240		 MT=MT+1
		 DRR=RDT(MT)-R0
		 IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 210
		 DP=DRR*DEDX0
		 IF(DP.LE.0.2) GO TO 210
		 IF(PTOT.LE.0.4) GO TO 290
		 IF(PTOT.LE.DP) DP=PTOT-0.2
		 DE=DP

		 IF(KFPJ(JP,I).NE.21) THEN
		    PRSHU=PT(LQT(MT),1)**2+PT(LQT(MT),2)**2
     &                   +PT(LQT(MT),3)**2
		    DE=SQRT(PJPM(JP,I)**2+PTOT**2)
     &			-SQRT(PJPM(JP,I)**2+(PTOT-DP)**2)
		    ERSHU=(PT(LQT(MT),4)+DE-DP)**2
		    AMSHU=ERSHU-PRSHU
		    IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 210
		    PT(LQT(MT),4)=SQRT(ERSHU)
		    PT(LQT(MT),5)=SQRT(AMSHU)
		 ENDIF
C		********reshuffle the energy when jet has mass

		 R0=RDT(MT)
		 DP1=DP*V1
		 DP2=DP*V2
		 DP3=DP*V3
C		********momentum and energy transfer from jet
		 NTJ(LQT(MT))=NTJ(LQT(MT))+1
		 KFTJ(LQT(MT),NTJ(LQT(MT)))=21
		 PJTX(LQT(MT),NTJ(LQT(MT)))=DP1
		 PJTY(LQT(MT),NTJ(LQT(MT)))=DP2
		 PJTZ(LQT(MT),NTJ(LQT(MT)))=DP3
		 PJTE(LQT(MT),NTJ(LQT(MT)))=DP
		 PJTM(LQT(MT),NTJ(LQT(MT)))=0.0

 260		 PJPX(JP,I)=(PTOT-DP)*V1
		 PJPY(JP,I)=(PTOT-DP)*V2
		 PJPZ(JP,I)=(PTOT-DP)*V3
		 PJPE(JP,I)=PJPE(JP,I)-DE

		 PTOT=PTOT-DP
		 NQ=NQ+1
		 GO TO 200
 290	      CONTINUE

	      RETURN

C*******************************************************
C Jet interaction for target jet in the direction PHIT
C******************************************************
C
C******* find the wounded proj which can interact with jet***

 400	      IF(NFT(JPJT,7).NE.1) RETURN
	      JT=JPJT
	      DO 690 I=1,NTJ(JT)
		 PTJET0=SQRT(PJTX(JT,I)**2+PJTY(JT,I)**2)
		 IF(PTJET0.LE.HIPR1(11)) GO TO 690
		 PTOT=SQRT(PTJET0*PTJET0+PJTZ(JT,I)**2)
		 IF (IHPR2(50) .EQ. 1) THEN 
		    DEDX0 = HIPR1(14)*LOG10(PTJET0)/LOG10(5.)
		 ELSE
		    DEDX0 = HIPR1(14)
		 ENDIF
		 IF(PTOT.LT.HIPR1(8)) GO TO 690
		 PHIT=ULANGL_HIJING(PJTX(JT,I),PJTY(JT,I))
		 KP=0
		 DO 500 I2=1,IHNT2(1)
		    IF(NFP(I2,5).NE.3) GO TO 500
		    DX=YP(1,I2)+BBX-YT(1,JT)
		    DY=YP(2,I2)+BBY-YT(2,JT)
		    PHI=ULANGL_HIJING(DX,DY)
		    DPHI=ABS(PHI-PHIT)
		    IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI ! Uzhi
		    IF(DPHI.GT.HIPR1(40)/2.0) GO TO 500
		    RD0=SQRT(DX*DX+DY*DY)
		    IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 500
		    KP=KP+1
		    LQP(KP)=I2
		    RDP(KP)=COS(DPHI)*RD0
 500		 CONTINUE
C*******	rearrange according to decending rd************
		 DO 510 I2=1,KP-1
		    DO 510 J2=I2+1,KP
		       IF(RDP(I2).LT.RDP(J2)) GO TO 510
		       RD=RDP(I2)
		       LQ=LQP(I2)
		       RDP(I2)=RDP(J2)
		       LQP(I2)=LQP(J2)
		       RDP(J2)=RD
		       LQP(J2)=LQ
 510		    CONTINUE
C****** find wounded targ which can interact with jet********
		    KT=0
		    DO 520 I2=1,IHNT2(3)
		       IF(NFT(I2,5).NE.3 .OR. I2.EQ.JT) GO TO 520
		       DX=YT(1,I2)-YT(1,JT)
		       DY=YT(2,I2)-YT(2,JT)
		       PHI=ULANGL_HIJING(DX,DY)
		       DPHI=ABS(PHI-PHIT)
		       IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI ! Uzhi
		       IF(DPHI.GT.HIPR1(40)/2.0) GO TO 520
		       RD0=SQRT(DX*DX+DY*DY)
		       IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 520
		       KT=KT+1
		       LQT(KT)=I2
		       RDT(KT)=COS(DPHI)*RD0
 520		    CONTINUE
C*******	rearrange according to decending rd************
		    DO 530 I2=1,KT-1
		       DO 530 J2=I2+1,KT
			  IF(RDT(I2).LT.RDT(J2)) GO TO 530
			  RD=RDT(I2)
			  LQ=LQT(I2)
			  RDT(I2)=RDT(J2)
			  LQT(I2)=LQT(J2)
			  RDT(J2)=RD
			  LQT(J2)=LQ
 530		       CONTINUE
		       
		       MP=0
		       MT=0
		       NQ=0
		       DP=0.0
		       R0=0.0
		PTOT=SQRT(PJTX(JT,I)**2+PJTY(JT,I)**2+PJTZ(JT,I)**2)
		V1=PJTX(JT,I)/PTOT
		V2=PJTY(JT,I)/PTOT
		V3=PJTZ(JT,I)/PTOT

 600		RN=RLU_HIJING(0)
 610		IF(MT.GE.KT .AND. MP.GE.KP) GO TO 690
		IF(MT.GE.KT) GO TO 620
		IF(MP.GE.KP) GO TO 640
		IF(RDP(MP+1).GT.RDT(MT+1)) GO TO 640
620		MP=MP+1
		DRR=RDP(MP)-R0
		IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 610
		DP=DRR*DEDX0
		IF(KFTJ(JT,I).NE.21) DP=0.5*DP
C	********string tension of quark jet is 0.5 of gluon's 
		IF(DP.LE.0.2) GO TO 610
		IF(PTOT.LE.0.4) GO TO 690
		IF(PTOT.LE.DP) DP=PTOT-0.2
		DE=DP
C
		IF(KFTJ(JT,I).NE.21) THEN
		   PRSHU=PP(LQP(MP),1)**2+PP(LQP(MP),2)**2
     &                   +PP(LQP(MP),3)**2
		   DE=SQRT(PJTM(JT,I)**2+PTOT**2)
     &		     -SQRT(PJTM(JT,I)**2+(PTOT-DP)**2)
		   ERSHU=(PP(LQP(MP),4)+DE-DP)**2
		   AMSHU=ERSHU-PRSHU
		   IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 610
		   PP(LQP(MP),4)=SQRT(ERSHU)
		   PP(LQP(MP),5)=SQRT(AMSHU)
		ENDIF
C		********reshuffle the energy when jet has mass
C
		R0=RDP(MP)
		DP1=DP*V1
		DP2=DP*V2
		DP3=DP*V3
C		********momentum and energy transfer from jet
		NPJ(LQP(MP))=NPJ(LQP(MP))+1
		KFPJ(LQP(MP),NPJ(LQP(MP)))=21
		PJPX(LQP(MP),NPJ(LQP(MP)))=DP1
		PJPY(LQP(MP),NPJ(LQP(MP)))=DP2
		PJPZ(LQP(MP),NPJ(LQP(MP)))=DP3
		PJPE(LQP(MP),NPJ(LQP(MP)))=DP
		PJPM(LQP(MP),NPJ(LQP(MP)))=0.0

		GO TO 660

640		MT=MT+1
		DRR=RDT(MT)-R0
		IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 610
		DP=DRR*DEDX0
		IF(DP.LE.0.2) GO TO 610
		IF(PTOT.LE.0.4) GO TO 690
		IF(PTOT.LE.DP) DP=PTOT-0.2
		DE=DP

		IF(KFTJ(JT,I).NE.21) THEN
		   PRSHU=PT(LQT(MT),1)**2+PT(LQT(MT),2)**2
     &                   +PT(LQT(MT),3)**2
		   DE=SQRT(PJTM(JT,I)**2+PTOT**2)
     &		     -SQRT(PJTM(JT,I)**2+(PTOT-DP)**2)
		   ERSHU=(PT(LQT(MT),4)+DE-DP)**2
		   AMSHU=ERSHU-PRSHU
		   IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 610
		   PT(LQT(MT),4)=SQRT(ERSHU)
		   PT(LQT(MT),5)=SQRT(AMSHU)
		ENDIF
C		********reshuffle the energy when jet has mass

		R0=RDT(MT)
		DP1=DP*V1
		DP2=DP*V2
		DP3=DP*V3
C		********momentum and energy transfer from jet
		NTJ(LQT(MT))=NTJ(LQT(MT))+1
		KFTJ(LQT(MT),NTJ(LQT(MT)))=21
		PJTX(LQT(MT),NTJ(LQT(MT)))=DP1
		PJTY(LQT(MT),NTJ(LQT(MT)))=DP2
		PJTZ(LQT(MT),NTJ(LQT(MT)))=DP3
		PJTE(LQT(MT),NTJ(LQT(MT)))=DP
		PJTM(LQT(MT),NTJ(LQT(MT)))=0.0

660		PJTX(JT,I)=(PTOT-DP)*V1
		PJTY(JT,I)=(PTOT-DP)*V2
		PJTZ(JT,I)=(PTOT-DP)*V3
		PJTE(JT,I)=PJTE(JT,I)-DE

		PTOT=PTOT-DP
		NQ=NQ+1
		GO TO 600
690	CONTINUE
	RETURN
C********************************************************
C	Q-QBAR jet interaction
C********************************************************
2000	ISG=JPJT
	IF(IASG(ISG,3).NE.1) RETURN
C
	JP=IASG(ISG,1)
	JT=IASG(ISG,2)
	XJ=(YP(1,JP)+BBX+YT(1,JT))/2.0
	YJ=(YP(2,JP)+BBY+YT(2,JT))/2.0
	DO 2690 I=1,NJSG(ISG)
	   PTJET0=SQRT(PXSG(ISG,I)**2+PYSG(ISG,I)**2)
	   IF(PTJET0.LE.HIPR1(11).OR.PESG(ISG,I).LT.HIPR1(1))
     &            GO TO 2690
	   PTOT=SQRT(PTJET0*PTJET0+PZSG(ISG,I)**2)
	   IF (IHPR2(50) .EQ. 1) THEN 
	      DEDX0 = HIPR1(14)*LOG10(PTJET0)/LOG10(5.)
	   ELSE
	      DEDX0 = HIPR1(14)
	   ENDIF
	   IF(PTOT.LT.MAX(HIPR1(1),HIPR1(8))) GO TO 2690
	   PHIQ=ULANGL_HIJING(PXSG(ISG,I),PYSG(ISG,I))
	   KP=0
	   DO 2500 I2=1,IHNT2(1)
	      IF(NFP(I2,5).NE.3.OR.I2.EQ.JP) GO TO 2500
	      DX=YP(1,I2)+BBX-XJ
	      DY=YP(2,I2)+BBY-YJ
	      PHI=ULANGL_HIJING(DX,DY)
	      DPHI=ABS(PHI-PHIQ)
              IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI      ! Uzhi
	      IF(DPHI.GT.HIPR1(40)/2.0) GO TO 2500
	      RD0=SQRT(DX*DX+DY*DY)
	      IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 2500
	      KP=KP+1
	      LQP(KP)=I2
	      RDP(KP)=COS(DPHI)*RD0
 2500	   CONTINUE
C*******	rearrange according to decending rd************
	   DO 2510 I2=1,KP-1
	      DO 2510 J2=I2+1,KP
		 IF(RDP(I2).LT.RDP(J2)) GO TO 2510
		 RD=RDP(I2)
		 LQ=LQP(I2)
		 RDP(I2)=RDP(J2)
		 LQP(I2)=LQP(J2)
		 RDP(J2)=RD
		 LQP(J2)=LQ
 2510	      CONTINUE
C****** find wounded targ which can interact with jet********
	      KT=0
	      DO 2520 I2=1,IHNT2(3)
		 IF(NFT(I2,5).NE.3 .OR. I2.EQ.JT) GO TO 2520
		 DX=YT(1,I2)-XJ
		 DY=YT(2,I2)-YJ
		 PHI=ULANGL_HIJING(DX,DY)
		 DPHI=ABS(PHI-PHIQ)
		 IF(DPHI.GE.HIPR1(40)) DPHI=2.*HIPR1(40)-DPHI ! Uzhi
		 IF(DPHI.GT.HIPR1(40)/2.0) GO TO 2520
		 RD0=SQRT(DX*DX+DY*DY)
		 IF(RD0*SIN(DPHI).GT.HIPR1(12)) GO TO 2520
		 KT=KT+1
		 LQT(KT)=I2
		 RDT(KT)=COS(DPHI)*RD0
 2520	      CONTINUE
C*******	rearrange according to decending rd************
	      DO 2530 I2=1,KT-1
		 DO 2530 J2=I2+1,KT
		    IF(RDT(I2).LT.RDT(J2)) GO TO 2530
		    RD=RDT(I2)
		    LQ=LQT(I2)
		    RDT(I2)=RDT(J2)
		    LQT(I2)=LQT(J2)
		    RDT(J2)=RD
		    LQT(J2)=LQ
 2530		 CONTINUE
		
		 MP=0
		 MT=0
		 NQ=0
		 DP=0.0
		 R0=0.0
		 PTOT=SQRT(PXSG(ISG,I)**2+PYSG(ISG,I)**2
     &                +PZSG(ISG,I)**2)
		 V1=PXSG(ISG,I)/PTOT
		 V2=PYSG(ISG,I)/PTOT
		 V3=PZSG(ISG,I)/PTOT

 2600		 RN=RLU_HIJING(0)
 2610		 IF(MT.GE.KT .AND. MP.GE.KP) GO TO 2690
		 IF(MT.GE.KT) GO TO 2620
		 IF(MP.GE.KP) GO TO 2640
		 IF(RDP(MP+1).GT.RDT(MT+1)) GO TO 2640
 2620		 MP=MP+1
		 DRR=RDP(MP)-R0
		 IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 2610
		 DP=DRR*DEDX0/2.0
		 IF(DP.LE.0.2) GO TO 2610
		 IF(PTOT.LE.0.4) GO TO 2690
		 IF(PTOT.LE.DP) DP=PTOT-0.2
		 DE=DP
C
		 IF(K2SG(ISG,I).NE.21) THEN
		    IF(PTOT.LT.DP+HIPR1(1)) GO TO 2690
		    PRSHU=PP(LQP(MP),1)**2+PP(LQP(MP),2)**2
     &                    +PP(LQP(MP),3)**2
		    DE=SQRT(PMSG(ISG,I)**2+PTOT**2)
     &		       -SQRT(PMSG(ISG,I)**2+(PTOT-DP)**2)
		    ERSHU=(PP(LQP(MP),4)+DE-DP)**2
		    AMSHU=ERSHU-PRSHU
		    IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 2610
		    PP(LQP(MP),4)=SQRT(ERSHU)
		    PP(LQP(MP),5)=SQRT(AMSHU)
		 ENDIF
C		********reshuffle the energy when jet has mass
C
		 R0=RDP(MP)
		 DP1=DP*V1
		 DP2=DP*V2
		 DP3=DP*V3
C		********momentum and energy transfer from jet
		 NPJ(LQP(MP))=NPJ(LQP(MP))+1
		 KFPJ(LQP(MP),NPJ(LQP(MP)))=21
		 PJPX(LQP(MP),NPJ(LQP(MP)))=DP1
		 PJPY(LQP(MP),NPJ(LQP(MP)))=DP2
		 PJPZ(LQP(MP),NPJ(LQP(MP)))=DP3
		 PJPE(LQP(MP),NPJ(LQP(MP)))=DP
		 PJPM(LQP(MP),NPJ(LQP(MP)))=0.0

		 GO TO 2660

 2640		 MT=MT+1
		 DRR=RDT(MT)-R0
		 IF(RN.GE.1.0-EXP(-DRR/HIPR1(13))) GO TO 2610
		 DP=DRR*DEDX0
		 IF(DP.LE.0.2) GO TO 2610
		 IF(PTOT.LE.0.4) GO TO 2690
		 IF(PTOT.LE.DP) DP=PTOT-0.2
		 DE=DP

		 IF(K2SG(ISG,I).NE.21) THEN
		    IF(PTOT.LT.DP+HIPR1(1)) GO TO 2690
		    PRSHU=PT(LQT(MT),1)**2+PT(LQT(MT),2)**2
     &                    +PT(LQT(MT),3)**2
		    DE=SQRT(PMSG(ISG,I)**2+PTOT**2)
     &		       -SQRT(PMSG(ISG,I)**2+(PTOT-DP)**2)
		    ERSHU=(PT(LQT(MT),4)+DE-DP)**2
		    AMSHU=ERSHU-PRSHU
		    IF(AMSHU.LT.HIPR1(1)*HIPR1(1)) GO TO 2610
		    PT(LQT(MT),4)=SQRT(ERSHU)
		    PT(LQT(MT),5)=SQRT(AMSHU)
		 ENDIF
C               ********reshuffle the energy when jet has mass

		 R0=RDT(MT)
		 DP1=DP*V1
		 DP2=DP*V2
		 DP3=DP*V3
C		********momentum and energy transfer from jet
		 NTJ(LQT(MT))=NTJ(LQT(MT))+1
		 KFTJ(LQT(MT),NTJ(LQT(MT)))=21
		 PJTX(LQT(MT),NTJ(LQT(MT)))=DP1
		 PJTY(LQT(MT),NTJ(LQT(MT)))=DP2
		 PJTZ(LQT(MT),NTJ(LQT(MT)))=DP3
		 PJTE(LQT(MT),NTJ(LQT(MT)))=DP
		 PJTM(LQT(MT),NTJ(LQT(MT)))=0.0

 2660		 PXSG(ISG,I)=(PTOT-DP)*V1
		 PYSG(ISG,I)=(PTOT-DP)*V2
		 PZSG(ISG,I)=(PTOT-DP)*V3
		 PESG(ISG,I)=PESG(ISG,I)-DE

		 PTOT=PTOT-DP
		 NQ=NQ+1
		 GO TO 2600
 2690	CONTINUE
	RETURN
	END







