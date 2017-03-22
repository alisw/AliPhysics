* $Id$
C
C	
C*******************************************************************
C	   				                           *
C		Subroutine HIJSFT		                   *
C						                   *
C  Scatter two excited strings, JP from proj and JT from target    *
C*******************************************************************
	SUBROUTINE HIJSFT(JP,JT,JOUT,IERROR)
#define BLANKET_SAVE
#include "hijcrdn.inc"
#include "hiparnt.inc"
#include "hijdat.inc"
#include "hijjet1.inc"
#include "hijjet2.inc"
#include "histrng.inc"
#include "dpmcom1.inc"
#include "dpmcom2.inc"
	SAVE
C*******************************************************************
C	JOUT-> the number
C	of hard scatterings preceding this soft collision. 
C       IHNT2(13)-> 1=
C	double diffrac 2=single diffrac, 3=non-single diffrac.
C*******************************************************************
	IERROR=0
	JJP=JP
	JJT=JT
	NDPM=0
	IOPMAIN=0
	IF(JP.GT.IHNT2(1) .OR. JT.GT.IHNT2(3)) RETURN

	EPP=PP(JP,4)+PP(JP,3)
	EPM=PP(JP,4)-PP(JP,3)
	ETP=PT(JT,4)+PT(JT,3)
	ETM=PT(JT,4)-PT(JT,3)

	WP=EPP+ETP
	WM=EPM+ETM
	SW=WP*WM
C		********total W+,W- and center-of-mass energy

	IF(WP.LT.0.0 .OR. WM.LT.0.0) GO TO 1000

	IF(JOUT.EQ.0) THEN
		IF(EPP.LT.0.0) GO TO 1000
		IF(EPM.LT.0.0) GO TO 1000
		IF(ETP.LT.0.0) GO TO 1000
		IF(ETM.LT.0.0) GO TO 1000    
		IF(EPP/(EPM+0.01).LE.ETP/(ETM+0.01)) RETURN
	ENDIF
C		********For strings which does not follow a jet-prod,
C			scatter only if Ycm(JP)>Ycm(JT). When jets
C			are produced just before this collision
C			this requirement has already be enforced
C			(see SUBROUTINE HIJHRD)
	IHNT2(11)=JP
	IHNT2(12)=JT
C
C
C
	MISS=0
	PKC1=0.0
	PKC2=0.0
	PKC11=0.0
	PKC12=0.0
	PKC21=0.0
	PKC22=0.0
	DPKC11=0.0
	DPKC12=0.0
	DPKC21=0.0
	DPKC22=0.0
	IF(NFP(JP,10).EQ.1.OR.NFT(JT,10).EQ.1) THEN
	   IF(NFP(JP,10).EQ.1) THEN
	      PHI1=ULANGL_HIJING(PP(JP,10),PP(JP,11))
	      PPJET=SQRT(PP(JP,10)**2+PP(JP,11)**2)
	      PKC1=PPJET
	      PKC11=PP(JP,10)
	      PKC12=PP(JP,11)
	   ENDIF
	   IF(NFT(JT,10).EQ.1) THEN
	      PHI2=ULANGL_HIJING(PT(JT,10),PT(JT,11))
	      PTJET=SQRT(PT(JT,10)**2+PT(JT,11)**2)
	      PKC2=PTJET
	      PKC21=PT(JT,10)
	      PKC22=PT(JT,11)
	   ENDIF
C
C Quenching requested ?
C
	   IF(IHPR2(4).GT.0.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
              IF(NFP(JP,10).EQ.0) THEN
                 PHI=-PHI2
              ELSE IF(NFT(JT,10).EQ.0) THEN
                 PHI=PHI1
              ELSE
                 IF (PHI1.LT.0. .AND. PHI2.GT.0) THEN ! PH AM
                    PHI=(PHI1+PHI2-HIPR1(40))/2.0     ! PH AM
                 ELSE                                 ! PH AM
                    PHI=(PHI1+PHI2+HIPR1(40))/2.0     ! PH AM
                 ENDIF                                ! PH AM
              ENDIF
	      BX=HINT1(19)*COS(HINT1(20))
	      BY=HINT1(19)*SIN(HINT1(20))
	      XP0=YP(1,JP)
	      YP0=YP(2,JP)
	      XT0=YT(1,JT)+BX
	      YT0=YT(2,JT)+BY
	      R1=MAX(1.2*IHNT2(1)**0.3333333,
     &               SQRT(XP0**2+YP0**2))
	      R2=MAX(1.2*IHNT2(3)**0.3333333,
     &               SQRT((XT0-BX)**2+(YT0-BY)**2))
	      IF(ABS(COS(PHI)).LT.1.0E-5) THEN
C AM		 DD1=R1
C AM		 DD2=R1
                 DD1 = ABS(SQRT(R1**2-(XT0-BX)**2) - YT0)
                 DD2 = ABS(SQRT(R1**2-(XT0-BX)**2) + YT0)
C PH		 DD3=ABS(BY+SQRT(R2**2-(XP0-BX)**2)-YP0)
C PH		 DD4=ABS(BY-SQRT(R2**2-(XP0-BX)**2)-YP0)
		 DD3=ABS(BY+SQRT(R2**2-(XT0-BX)**2)-YP0)
		 DD4=ABS(BY-SQRT(R2**2-(XT0-BX)**2)-YP0)
		 GO TO 5
	      ENDIF
	      BB=2.0*SIN(PHI)*(COS(PHI)*YP0-SIN(PHI)*XP0)
	      CC=(YP0**2-R1**2)*COS(PHI)**2+XP0*SIN(PHI)*(
     &				XP0*SIN(PHI)-2.0*YP0*COS(PHI))
	      DD=BB**2-4.0*CC
	      IF(DD.LT.0.0) GO TO 10
	      XX1=(-BB+SQRT(DD))/2.0
	      XX2=(-BB-SQRT(DD))/2.0
	      DD1=ABS((XX1-XP0)/COS(PHI))
	      DD2=ABS((XX2-XP0)/COS(PHI))
C			
	      BB=2.0*SIN(PHI)*(COS(PHI)*(YT0-BY)-SIN(PHI)*XT0)-2.0*BX
	      CC=(BX**2+(YT0-BY)**2-R2**2)*COS(PHI)**2+XT0*SIN(PHI)
     &           *(XT0*SIN(PHI)-2.0*COS(PHI)*(YT0-BY))
     &		 -2.0*BX*SIN(PHI)*(COS(PHI)*(YT0-BY)-SIN(PHI)*XT0)
	      DD=BB**2-4.0*CC
	      IF(DD.LT.0.0) GO TO 10
	      XX1=(-BB+SQRT(DD))/2.0
	      XX2=(-BB-SQRT(DD))/2.0
	      DD3=ABS((XX1-XT0)/COS(PHI))
	      DD4=ABS((XX2-XT0)/COS(PHI))
C
 5	      DD1=MIN(DD1,DD3)
	      DD2=MIN(DD2,DD4)
	      IF(DD1.LT.HIPR1(13)) DD1=0.0
	      IF(DD2.LT.HIPR1(13)) DD2=0.0
	      IF(NFP(JP,10).EQ.1.AND.PPJET.GT.HIPR1(11)) THEN
		 IF (IHPR2(50) .EQ. 1) THEN 
		    DEDX0 = HIPR1(14)*LOG10(PPJET)/LOG10(5.)
		 ELSE
		    DEDX0 = HIPR1(14)
		 ENDIF
		 DP1=DD1*DEDX0/2.0
		 DP1=MIN(DP1,PPJET-HIPR1(11))
		 PKC1=PPJET-DP1
		 DPX1=COS(PHI1)*DP1
		 DPY1=SIN(PHI1)*DP1
		 PKC11=PP(JP,10)-DPX1
		 PKC12=PP(JP,11)-DPY1
		 IF(DP1.GT.0.0) THEN
		    CTHEP=PP(JP,12)/SQRT(PP(JP,12)**2+PPJET**2)
		    DPZ1=DP1*CTHEP/SQRT(1.0-CTHEP**2)
		    DPE1=SQRT(DPX1**2+DPY1**2+DPZ1**2)
		    EPPPRM=PP(JP,4)+PP(JP,3)-DPE1-DPZ1
		    EPMPRM=PP(JP,4)-PP(JP,3)-DPE1+DPZ1
		    IF(EPPPRM.LE.0.0.OR.EPMPRM.LE.0.0) GO TO 15
		    EPP=EPPPRM
		    EPM=EPMPRM
		    PP(JP,10)=PKC11
		    PP(JP,11)=PKC12
		    NPJ(JP)=NPJ(JP)+1
		    KFPJ(JP,NPJ(JP))=21
		    PJPX(JP,NPJ(JP))=DPX1
		    PJPY(JP,NPJ(JP))=DPY1
		    PJPZ(JP,NPJ(JP))=DPZ1
		    PJPE(JP,NPJ(JP))=DPE1
		    PJPM(JP,NPJ(JP))=0.0
		    PP(JP,3)=PP(JP,3)-DPZ1
		    PP(JP,4)=PP(JP,4)-DPE1
		 ENDIF
	      ENDIF
 15	      IF(NFT(JT,10).EQ.1.AND.PTJET.GT.HIPR1(11)) THEN
		 IF (IHPR2(50) .EQ. 1) THEN 
		    DEDX0 = HIPR1(14)*LOG10(PTJET)/LOG10(5.)
		 ELSE
		    DEDX0 = HIPR1(14)
		 ENDIF
		 DP2=DD2*DEDX0/2.0
		 DP2=MIN(DP2,PTJET-HIPR1(11))
		 PKC2=PTJET-DP2
		 DPX2=COS(PHI2)*DP2
		 DPY2=SIN(PHI2)*DP2
		 PKC21=PT(JT,10)-DPX2
		 PKC22=PT(JT,11)-DPY2
		 IF(DP2.GT.0.0) THEN
		    CTHET=PT(JT,12)/SQRT(PT(JT,12)**2+PTJET**2)
		    DPZ2=DP2*CTHET/SQRT(1.0-CTHET**2)
		    DPE2=SQRT(DPX2**2+DPY2**2+DPZ2**2)
		    ETPPRM=PT(JT,4)+PT(JT,3)-DPE2-DPZ2
		    ETMPRM=PT(JT,4)-PT(JT,3)-DPE2+DPZ2
		    IF(ETPPRM.LE.0.0.OR.ETMPRM.LE.0.0) GO TO 16
		    ETP=ETPPRM
		    ETM=ETMPRM
		    PT(JT,10)=PKC21
		    PT(JT,11)=PKC22
		    NTJ(JT)=NTJ(JT)+1
		    KFTJ(JT,NTJ(JT))=21
		    PJTX(JT,NTJ(JT))=DPX2
		    PJTY(JT,NTJ(JT))=DPY2
		    PJTZ(JT,NTJ(JT))=DPZ2
		    PJTE(JT,NTJ(JT))=DPE2
		    PJTM(JT,NTJ(JT))=0.0
		    PT(JT,3)=PT(JT,3)-DPZ2
		    PT(JT,4)=PT(JT,4)-DPE2
		 ENDIF
	      ENDIF
 16           CONTINUE
	      PKC11=0.
	      PKC12=0.
	      PKC21=0.
	      PKC22=0.
	      DPKC11=0.
	      DPKC12=0.
	      DPKC21=0.
	      DPKC22=0.

	      WP=EPP+ETP
	      WM=EPM+ETM
	      SW=WP*WM
	   ENDIF
	ENDIF

	IF(NFP(JP,10).EQ.1) THEN
	   PPJET=SQRT(PP(JP,10)**2+PP(JP,11)**2)
	   PKC1=PPJET
	ENDIF
	IF(NFT(JT,10).EQ.1) THEN
	   PTJET=SQRT(PT(JT,10)**2+PT(JT,11)**2)
	   PKC2=PTJET
	ENDIF
C		********If jet is quenched the pt from valence quark
C			hard scattering has to reduced by d*kapa
C
C   
10	PTP02=PP(JP,1)**2+PP(JP,2)**2
	PTT02=PT(JT,1)**2+PT(JT,2)**2
C	
	AMQ=MAX(PP(JP,14)+PP(JP,15),PT(JT,14)+PT(JT,15))
	AMX=HIPR1(1)+AMQ
C		********consider mass cut-off for strings which
C			must also include quark's mass
	AMP0=AMX
	DPM0=AMX
	NFDP=0
	IF(NFP(JP,5).LE.2.AND.NFP(JP,3).NE.0) THEN
		AMP0=ULMASS_HIJING(NFP(JP,3))
		NFDP=NFP(JP,3)+2*NFP(JP,3)/ABS(NFP(JP,3))
		DPM0=ULMASS_HIJING(NFDP)
		IF(DPM0.LE.0.0) THEN
			NFDP=NFDP-2*NFDP/ABS(NFDP)
			DPM0=ULMASS_HIJING(NFDP)
		ENDIF
	ENDIF
	AMT0=AMX
	DTM0=AMX
	NFDT=0
	IF(NFT(JT,5).LE.2.AND.NFT(JT,3).NE.0) THEN
		AMT0=ULMASS_HIJING(NFT(JT,3))
		NFDT=NFT(JT,3)+2*NFT(JT,3)/ABS(NFT(JT,3))
		DTM0=ULMASS_HIJING(NFDT)
		IF(DTM0.LE.0.0) THEN
			NFDT=NFDT-2*NFDT/ABS(NFDT)
			DTM0=ULMASS_HIJING(NFDT)
		ENDIF
	ENDIF
C	
	AMPN=SQRT(AMP0**2+PTP02)
	AMTN=SQRT(AMT0**2+PTT02)
	SNN=(AMPN+AMTN)**2+0.001
C
	IF(SW.LT.SNN+0.001) GO TO 4000
C		********Scatter only if SW>SNN
C*****give some PT kick to the two exited strings******************
20	SWPTN=4.0*(MAX(AMP0,AMT0)**2+MAX(PTP02,PTT02))
	SWPTD=4.0*(MAX(DPM0,DTM0)**2+MAX(PTP02,PTT02))
	SWPTX=4.0*(AMX**2+MAX(PTP02,PTT02))
	IF(SW.LE.SWPTN) THEN
		PKCMX=0.0
	ELSE IF(SW.GT.SWPTN .AND. SW.LE.SWPTD
     &		.AND.NPJ(JP).EQ.0.AND.NTJ(JT).EQ.0) THEN
	   PKCMX=SQRT(SW/4.0-MAX(AMP0,AMT0)**2)
     &           -SQRT(MAX(PTP02,PTT02))
	ELSE IF(SW.GT.SWPTD .AND. SW.LE.SWPTX
     &		.AND.NPJ(JP).EQ.0.AND.NTJ(JT).EQ.0) THEN
	   PKCMX=SQRT(SW/4.0-MAX(DPM0,DTM0)**2)
     &           -SQRT(MAX(PTP02,PTT02))
	ELSE IF(SW.GT.SWPTX) THEN
	   PKCMX=SQRT(SW/4.0-AMX**2)-SQRT(MAX(PTP02,PTT02))
	ENDIF
C		********maximun PT kick
C*********************************************************
	IF(NFP(JP,10).EQ.1.OR.NFT(JT,10).EQ.1) THEN
		IF(PKC1.GT.PKCMX) THEN
			PKC1=PKCMX
			PKC11=PKC1*COS(PHI1)
			PKC12=PKC1*SIN(PHI1)
			DPKC11=-(PP(JP,10)-PKC11)/2.0
			DPKC12=-(PP(JP,11)-PKC12)/2.0
		ENDIF
		IF(PKC2.GT.PKCMX) THEN
			PKC2=PKCMX
			PKC21=PKC2*COS(PHI2)
			PKC22=PKC2*SIN(PHI2)
			DPKC21=-(PT(JT,10)-PKC21)/2.0
			DPKC22=-(PT(JT,11)-PKC22)/2.0
		ENDIF
		DPKC1=DPKC11+DPKC21
		DPKC2=DPKC12+DPKC22
		NFP(JP,10)=-NFP(JP,10)
		NFT(JT,10)=-NFT(JT,10)
		GO TO 40
	ENDIF
C		********If the valence quarks had a hard-collision
C			the pt kick is the pt from hard-collision.
	I_SNG=0
	IF(IHPR2(13).NE.0 .AND. RLU_HIJING(0).LE.HIDAT(4)) I_SNG=1
	IF((NFP(JP,5).EQ.3 .OR.NFT(JT,5).EQ.3).OR.
     &		(NPJ(JP).NE.0.OR.NFP(JP,10).NE.0).OR.
     &		(NTJ(JT).NE.0.OR.NFT(JT,10).NE.0)) I_SNG=0
C
C               ********decite whether to have single-diffractive
	IF(IHPR2(5).EQ.0) THEN
		PKC=HIPR1(2)*SQRT(-ALOG(1.0-RLU_HIJING(0)
     &			*(1.0-EXP(-PKCMX**2/HIPR1(2)**2))))
		GO TO 30
	ENDIF
	PKC=HIRND2(3,0.0,PKCMX**2)
	PKC=SQRT(PKC)
	IF(PKC.GT.HIPR1(20)) 
     &	   PKC=HIPR1(2)*SQRT(-ALOG(EXP(-HIPR1(20)**2/HIPR1(2)**2)
     &	       -RLU_HIJING(0)*(EXP(-HIPR1(20)**2/HIPR1(2)**2)-
     &	       EXP(-PKCMX**2/HIPR1(2)**2))))
C
	IF(I_SNG.EQ.1) PKC=0.65*SQRT(
     &		-ALOG(1.0-RLU_HIJING(0)*(1.0-EXP(-PKCMX**2/0.65**2))))
C			********select PT kick
30	PHI0=2.0*HIPR1(40)*RLU_HIJING(0)
	PKC11=PKC*SIN(PHI0)
	PKC12=PKC*COS(PHI0)
	PKC21=-PKC11
	PKC22=-PKC12
	DPKC1=0.0
	DPKC2=0.0
40	PP11=PP(JP,1)+PKC11-DPKC1
	PP12=PP(JP,2)+PKC12-DPKC2
	PT11=PT(JT,1)+PKC21-DPKC1
	PT12=PT(JT,2)+PKC22-DPKC2
	PTP2=PP11**2+PP12**2
	PTT2=PT11**2+PT12**2
C
	AMPN=SQRT(AMP0**2+PTP2)
	AMTN=SQRT(AMT0**2+PTT2)
	SNN=(AMPN+AMTN)**2+0.001
C***************************************
	WP=EPP+ETP
	WM=EPM+ETM
	SW=WP*WM
C****************************************
	IF(SW.LT.SNN) THEN
	   MISS=MISS+1
	   IF(MISS.LE.100) then
	      PKC=0.0
	      GO TO 30
	   ENDIF
	   IF(IHPR2(10).NE.0) 
     &	     WRITE(6,*) 'Error occured in Pt kick section of HIJSFT'
	   GO TO 4000
	ENDIF
C******************************************************************
	AMPD=SQRT(DPM0**2+PTP2)
	AMTD=SQRT(DTM0**2+PTT2)

	AMPX=SQRT(AMX**2+PTP2)
	AMTX=SQRT(AMX**2+PTT2)

	DPN=AMPN**2/SW
	DTN=AMTN**2/SW
	DPD=AMPD**2/SW
	DTD=AMTD**2/SW
	DPX=AMPX**2/SW
	DTX=AMTX**2/SW
C
	SPNTD=(AMPN+AMTD)**2
	SPNTX=(AMPN+AMTX)**2
C			********CM energy if proj=N,targ=N*
	SPDTN=(AMPD+AMTN)**2
	SPXTN=(AMPX+AMTN)**2
C			********CM energy if proj=N*,targ=N
	SPDTX=(AMPD+AMTX)**2
	SPXTD=(AMPX+AMTD)**2
	SDD=(AMPD+AMTD)**2
	SXX=(AMPX+AMTX)**2

C
C	
C		********CM energy if proj=delta, targ=delta
C****************There are many different cases**********
c	IF(IHPR2(15).EQ.1) GO TO 500
C
C		********to have DPM type soft interactions
C
 45	CONTINUE
	IF(SW.GT.SXX+0.001) THEN
	   IF(I_SNG.EQ.0) THEN
 	      D1=DPX
	      D2=DTX
	      NFP3=0
	      NFT3=0
	      GO TO 400
	   ELSE
c**** 5/30/1998 this is identical to the above statement. Added to
c**** avoid questional branching to block.
	      IF((NFP(JP,5).EQ.3 .AND.NFT(JT,5).EQ.3).OR.
     &	         (NPJ(JP).NE.0.OR.NFP(JP,10).NE.0).OR.
     &		 (NTJ(JT).NE.0.OR.NFT(JT,10).NE.0)) THEN
                 D1=DPX
                 D2=DTX
                 NFP3=0
                 NFT3=0
                 GO TO 400
              ENDIF
C		********do not allow excited strings to have 
C			single-diffr 
	      IF(RLU_HIJING(0).GT.0.5.OR.(NFT(JT,5).GT.2.OR.
     &		      NTJ(JT).NE.0.OR.NFT(JT,10).NE.0)) THEN
		 D1=DPN
		 D2=DTX
		 NFP3=NFP(JP,3)
		 NFT3=0
		 GO TO 220
	      ELSE
		 D1=DPX
		 D2=DTN
		 NFP3=0
		 NFT3=NFT(JT,3)
		 GO TO 240
	      ENDIF
C		********have single diffractive collision
	   ENDIF
	ELSE IF(SW.GT.MAX(SPDTX,SPXTD)+0.001 .AND.
     &				SW.LE.SXX+0.001) THEN
	   IF(((NPJ(JP).EQ.0.AND.NTJ(JT).EQ.0.AND.
     &         RLU_HIJING(0).GT.0.5).OR.(NPJ(JP).EQ.0
     &         .AND.NTJ(JT).NE.0)).AND.NFP(JP,5).LE.2) THEN
	      D1=DPD
	      D2=DTX
	      NFP3=NFDP
	      NFT3=0
	      GO TO 220
	   ELSE IF(NTJ(JT).EQ.0.AND.NFT(JT,5).LE.2) THEN
	      D1=DPX
	      D2=DTD
	      NFP3=0
	      NFT3=NFDT
	      GO TO 240
	   ENDIF
	   GO TO 4000
	ELSE IF(SW.GT.MIN(SPDTX,SPXTD)+0.001.AND.
     &			SW.LE.MAX(SPDTX,SPXTD)+0.001) THEN
	   IF(SPDTX.LE.SPXTD.AND.NPJ(JP).EQ.0
     &                       .AND.NFP(JP,5).LE.2) THEN
	      D1=DPD
	      D2=DTX
	      NFP3=NFDP
	      NFT3=0
	      GO TO 220
	   ELSE IF(SPDTX.GT.SPXTD.AND.NTJ(JT).EQ.0
     &                       .AND.NFT(JT,5).LE.2) THEN
	      D1=DPX
	      D2=DTD
	      NFP3=0
	      NFT3=NFDT
	      GO TO 240
	   ENDIF
c*** 5/30/1998 added to avoid questional branching to another block
c*** this is identical to the statement following the next ELSE IF
	   IF(((NPJ(JP).EQ.0.AND.NTJ(JT).EQ.0
     &       .AND.RLU_HIJING(0).GT.0.5).OR.(NPJ(JP).EQ.0
     &        .AND.NTJ(JT).NE.0)).AND.NFP(JP,5).LE.2) THEN
	      D1=DPN
	      D2=DTX
	      NFP3=NFP(JP,3)
	      NFT3=0
	      GO TO 220
	   ELSE IF(NTJ(JT).EQ.0.AND.NFT(JT,5).LE.2) THEN
	      D1=DPX
	      D2=DTN
	      NFP3=0
	      NFT3=NFT(JT,3)
	      GO TO 240
	   ENDIF
	   GO TO 4000
	ELSE IF(SW.GT.MAX(SPNTX,SPXTN)+0.001 .AND.
     &			SW.LE.MIN(SPDTX,SPXTD)+0.001) THEN
	   IF(((NPJ(JP).EQ.0.AND.NTJ(JT).EQ.0
     &       .AND.RLU_HIJING(0).GT.0.5).OR.(NPJ(JP).EQ.0
     &        .AND.NTJ(JT).NE.0)).AND.NFP(JP,5).LE.2) THEN
	      D1=DPN
	      D2=DTX
	      NFP3=NFP(JP,3)
	      NFT3=0
	      GO TO 220
	   ELSE IF(NTJ(JT).EQ.0.AND.NFT(JT,5).LE.2) THEN
	      D1=DPX
	      D2=DTN
	      NFP3=0
	      NFT3=NFT(JT,3)
	      GO TO 240
	   ENDIF
	   GO TO 4000
	ELSE IF(SW.GT.MIN(SPNTX,SPXTN)+0.001 .AND.
     &			SW.LE.MAX(SPNTX,SPXTN)+0.001) THEN
	   IF(SPNTX.LE.SPXTN.AND.NPJ(JP).EQ.0
     &                           .AND.NFP(JP,5).LE.2) THEN
	      D1=DPN
	      D2=DTX
	      NFP3=NFP(JP,3)
	      NFT3=0
	      GO TO 220
	   ELSEIF(SPNTX.GT.SPXTN.AND.NTJ(JT).EQ.0
     &                           .AND.NFT(JT,5).LE.2) THEN
	      D1=DPX
	      D2=DTN
	      NFP3=0
	      NFT3=NFT(JT,3)
	      GO TO 240
	   ENDIF
	   GO TO 4000
	ELSE IF(SW.LE.MIN(SPNTX,SPXTN)+0.001 .AND.
     &			(NPJ(JP).NE.0 .OR.NTJ(JT).NE.0)) THEN
	   GO TO 4000
	ELSE IF(SW.LE.MIN(SPNTX,SPXTN)+0.001 .AND.
     &		NFP(JP,5).GT.2.AND.NFT(JT,5).GT.2) THEN
	   GO TO 4000
	ELSE IF(SW.GT.SDD+0.001.AND.SW.LE.
     &                     MIN(SPNTX,SPXTN)+0.001) THEN
	   D1=DPD
	   D2=DTD
	   NFP3=NFDP
	   NFT3=NFDT
	   GO TO 100
	ELSE IF(SW.GT.MAX(SPNTD,SPDTN)+0.001 
     &                      .AND. SW.LE.SDD+0.001) THEN
	   IF(RLU_HIJING(0).GT.0.5) THEN
	      D1=DPD
	      D2=DTN
	      NFP3=NFDP
	      NFT3=NFT(JT,3)
	      GO TO 100
	   ELSE
	      D1=DPN
	      D2=DTD
	      NFP3=NFP(JP,3)
	      NFT3=NFDT
	      GO TO 100
	   ENDIF
	ELSE IF(SW.GT.MIN(SPNTD,SPDTN)+0.001
     &		.AND. SW.LE.MAX(SPNTD,SPDTN)+0.001) THEN
	   IF(SPNTD.GT.SPDTN) THEN
	      D1=DPD
	      D2=DTN
	      NFP3=NFDP
	      NFT3=NFT(JT,3)
	      GO TO 100
	   ELSE
	      D1=DPN
	      D2=DTD
	      NFP3=NFP(JP,3)
	      NFT3=NFDT
	      GO TO 100
	   ENDIF
	ELSE IF(SW.LE.MIN(SPNTD,SPDTN)+0.001) THEN
	   D1=DPN
	   D2=DTN
	   NFP3=NFP(JP,3)
	   NFT3=NFT(JT,3)
	   GO TO 100
	ENDIF
	WRITE(6,*) ' Error in HIJSFT: There is no path to here'
	RETURN
C
C***************  elastic scattering ***************
C	this is like elastic, both proj and targ mass
C	must be fixed
C***************************************************
100	NFP5=MAX(2,NFP(JP,5))
	NFT5=MAX(2,NFT(JT,5))
	BB1=1.0+D1-D2
	BB2=1.0+D2-D1
	IF(BB1**2.LT.4.0*D1 .OR. BB2**2.LT.4.0*D2) THEN
		MISS=MISS+1
		IF(MISS.GT.100.OR.PKC.EQ.0.0) GO TO 3000
		PKC=PKC*0.5
		GO TO 30
	ENDIF
	IF(RLU_HIJING(0).LT.0.5) THEN
		X1=(BB1-SQRT(BB1**2-4.0*D1))/2.0
		X2=(BB2-SQRT(BB2**2-4.0*D2))/2.0
	ELSE
		X1=(BB1+SQRT(BB1**2-4.0*D1))/2.0
		X2=(BB2+SQRT(BB2**2-4.0*D2))/2.0
	ENDIF
	IHNT2(13)=2
	GO TO 600
C
C********** Single diffractive ***********************
C either proj or targ's mass is fixed
C*****************************************************
220	NFP5=MAX(2,NFP(JP,5))
	NFT5=3
	IF(NFP3.EQ.0) NFP5=3
	BB2=1.0+D2-D1
	IF(BB2**2.LT.4.0*D2) THEN
		MISS=MISS+1
		IF(MISS.GT.100.OR.PKC.EQ.0.0) GO TO 3000
		PKC=PKC*0.5
		GO TO 30
	ENDIF
	XMIN=(BB2-SQRT(BB2**2-4.0*D2))/2.0
	XMAX=(BB2+SQRT(BB2**2-4.0*D2))/2.0
	MISS4=0
222	X2=HIRND2(6,XMIN,XMAX)
	X1=D1/(1.0-X2)
	IF(X2*(1.0-X1).LT.(D2+1.E-4/SW)) THEN
		MISS4=MISS4+1
		IF(MISS4.LE.1000) GO TO 222
		GO TO 5000
	ENDIF
	IHNT2(13)=2
	GO TO 600
C			********Fix proj mass*********
240	NFP5=3
	NFT5=MAX(2,NFT(JT,5))
	IF(NFT3.EQ.0) NFT5=3
	BB1=1.0+D1-D2
	IF(BB1**2.LT.4.0*D1) THEN
		MISS=MISS+1
		IF(MISS.GT.100.OR.PKC.EQ.0.0) GO TO 3000
		PKC=PKC*0.5
		GO TO 30
	ENDIF
	XMIN=(BB1-SQRT(BB1**2-4.0*D1))/2.0
	XMAX=(BB1+SQRT(BB1**2-4.0*D1))/2.0
	MISS4=0
242	X1=HIRND2(6,XMIN,XMAX)
	X2=D2/(1.0-X1)
	IF(X1*(1.0-X2).LT.(D1+1.E-4/SW)) THEN
		MISS4=MISS4+1
		IF(MISS4.LE.1000) GO TO 242
		GO TO 5000
	ENDIF
	IHNT2(13)=2
	GO TO 600
C			********Fix targ mass*********
C
C*************non-single diffractive**********************
C	both proj and targ may not be fixed in mass 
C*********************************************************
C
400	NFP5=3
	NFT5=3
	BB1=1.0+D1-D2
	BB2=1.0+D2-D1
	IF(BB1**2.LT.4.0*D1 .OR. BB2**2.LT.4.0*D2) THEN
		MISS=MISS+1
		IF(MISS.GT.100.OR.PKC.EQ.0.0) GO TO 3000
		PKC=PKC*0.5
		GO TO 30
	ENDIF
	XMIN1=(BB1-SQRT(BB1**2-4.0*D1))/2.0
	XMAX1=(BB1+SQRT(BB1**2-4.0*D1))/2.0
	XMIN2=(BB2-SQRT(BB2**2-4.0*D2))/2.0
	XMAX2=(BB2+SQRT(BB2**2-4.0*D2))/2.0
	MISS4=0	
410	X1=HIRND2(4,XMIN1,XMAX1)
	X2=HIRND2(4,XMIN2,XMAX2)
	IF(NFP(JP,5).EQ.3.OR.NFT(JT,5).EQ.3) THEN
		X1=HIRND2(6,XMIN1,XMAX1)
		X2=HIRND2(6,XMIN2,XMAX2)
	ENDIF
C			********
	IF(ABS(NFP(JP,1)*NFP(JP,2)).GT.1000000.OR.
     &			ABS(NFP(JP,1)*NFP(JP,2)).LT.100) THEN
		X1=HIRND2(5,XMIN1,XMAX1)
	ENDIF
	IF(ABS(NFT(JT,1)*NFT(JT,2)).GT.1000000.OR.
     &			ABS(NFT(JT,1)*NFT(JT,2)).LT.100) THEN
		X2=HIRND2(5,XMIN2,XMAX2)
	ENDIF
c	IF(IOPMAIN.EQ.3) X1=HIRND2(6,XMIN1,XMAX1)
c	IF(IOPMAIN.EQ.2) X2=HIRND2(6,XMIN2,XMAX2) 
C	********For q-qbar or (qq)-(qq)bar system use symetric
C		distribution, for q-(qq) or qbar-(qq)bar use
C		unsymetrical distribution
C
	IF(ABS(NFP(JP,1)*NFP(JP,2)).GT.1000000) X1=1.0-X1
	XXP=X1*(1.0-X2)
	XXT=X2*(1.0-X1)
	IF(XXP.LT.(D1+1.E-4/SW) .OR. XXT.LT.(D2+1.E-4/SW)) THEN
		MISS4=MISS4+1
		IF(MISS4.LE.1000) GO TO 410
		GO TO 5000
	ENDIF
	IHNT2(13)=3
C***************************************************
C***************************************************
600	CONTINUE
	IF(X1*(1.0-X2).LT.(AMPN**2-1.E-4)/SW.OR.
     &			X2*(1.0-X1).LT.(AMTN**2-1.E-4)/SW) THEN
		MISS=MISS+1
		IF(MISS.GT.100.OR.PKC.EQ.0.0) GO TO 2000
		PKC=0.0
		GO TO 30
	ENDIF
C
	EPP=(1.0-X2)*WP
	EPM=X1*WM
	ETP=X2*WP
	ETM=(1.0-X1)*WM
	PP(JP,3)=(EPP-EPM)/2.0
	PP(JP,4)=(EPP+EPM)/2.0
	IF(EPP*EPM-PTP2.LT.0.0) GO TO 6000
	PP(JP,5)=SQRT(EPP*EPM-PTP2)
	NFP(JP,3)=NFP3
	NFP(JP,5)=NFP5

	PT(JT,3)=(ETP-ETM)/2.0
	PT(JT,4)=(ETP+ETM)/2.0
	IF(ETP*ETM-PTT2.LT.0.0) GO TO 6000
	PT(JT,5)=SQRT(ETP*ETM-PTT2)
	NFT(JT,3)=NFT3
	NFT(JT,5)=NFT5
C*****recoil PT from hard-inter is shared by two end-partons 
C       so that pt=p1+p2
	PP(JP,1)=PP11-PKC11
	PP(JP,2)=PP12-PKC12

	KICKDIP=1
	KICKDIT=1
	IF(ABS(NFP(JP,1)*NFP(JP,2)).GT.1000000.OR.
     &			ABS(NFP(JP,1)*NFP(JP,2)).LT.100) THEN
		KICKDIP=0
	ENDIF
	IF(ABS(NFT(JT,1)*NFT(JT,2)).GT.1000000.OR.
     &			ABS(NFT(JT,1)*NFT(JT,2)).LT.100) THEN
		KICKDIT=0
	ENDIF
	IF((KICKDIP.EQ.0.AND.RLU_HIJING(0).LT.0.5)
     &     .OR.(KICKDIP.NE.0.AND.RLU_HIJING(0)
     &     .LT.0.5/(1.0+(PKC11**2+PKC12**2)/HIPR1(22)**2))) THEN
	   PP(JP,6)=(PP(JP,1)-PP(JP,6)-PP(JP,8)-DPKC1)/2.0+PP(JP,6)
	   PP(JP,7)=(PP(JP,2)-PP(JP,7)-PP(JP,9)-DPKC2)/2.0+PP(JP,7)
	   PP(JP,8)=(PP(JP,1)-PP(JP,6)-PP(JP,8)-DPKC1)/2.0
     &              +PP(JP,8)+PKC11
	   PP(JP,9)=(PP(JP,2)-PP(JP,7)-PP(JP,9)-DPKC2)/2.0
     &              +PP(JP,9)+PKC12
	ELSE
	   PP(JP,8)=(PP(JP,1)-PP(JP,6)-PP(JP,8)-DPKC1)/2.0+PP(JP,8)
	   PP(JP,9)=(PP(JP,2)-PP(JP,7)-PP(JP,9)-DPKC2)/2.0+PP(JP,9)
	   PP(JP,6)=(PP(JP,1)-PP(JP,6)-PP(JP,8)-DPKC1)/2.0
     &              +PP(JP,6)+PKC11
	   PP(JP,7)=(PP(JP,2)-PP(JP,7)-PP(JP,9)-DPKC2)/2.0
     &              +PP(JP,7)+PKC12
	ENDIF
	PP(JP,1)=PP(JP,6)+PP(JP,8)
	PP(JP,2)=PP(JP,7)+PP(JP,9)
C				********pt kick for proj
	PT(JT,1)=PT11-PKC21
	PT(JT,2)=PT12-PKC22
	IF((KICKDIT.EQ.0.AND.RLU_HIJING(0).LT.0.5)
     &     .OR.(KICKDIT.NE.0.AND.RLU_HIJING(0)
     &     .LT.0.5/(1.0+(PKC21**2+PKC22**2)/HIPR1(22)**2))) THEN
	   PT(JT,6)=(PT(JT,1)-PT(JT,6)-PT(JT,8)-DPKC1)/2.0+PT(JT,6)
	   PT(JT,7)=(PT(JT,2)-PT(JT,7)-PT(JT,9)-DPKC2)/2.0+PT(JT,7)
	   PT(JT,8)=(PT(JT,1)-PT(JT,6)-PT(JT,8)-DPKC1)/2.0
     &              +PT(JT,8)+PKC21
	   PT(JT,9)=(PT(JT,2)-PT(JT,7)-PT(JT,9)-DPKC2)/2.0
     &              +PT(JT,9)+PKC22
	ELSE
	   PT(JT,8)=(PT(JT,1)-PT(JT,6)-PT(JT,8)-DPKC1)/2.0+PT(JT,8)
	   PT(JT,9)=(PT(JT,2)-PT(JT,7)-PT(JT,9)-DPKC2)/2.0+PT(JT,9)
	   PT(JT,6)=(PT(JT,1)-PT(JT,6)-PT(JT,8)-DPKC1)/2.0
     &              +PT(JT,6)+PKC21
	   PT(JT,7)=(PT(JT,2)-PT(JT,7)-PT(JT,9)-DPKC2)/2.0
     &              +PT(JT,7)+PKC22
	ENDIF
	PT(JT,1)=PT(JT,6)+PT(JT,8)
	PT(JT,2)=PT(JT,7)+PT(JT,9)
C			********pt kick for targ

	IF(NPJ(JP).NE.0) NFP(JP,5)=3
	IF(NTJ(JT).NE.0) NFT(JT,5)=3
C			********jets must be connected to string
	IF(EPP/(EPM+0.0001).LT.ETP/(ETM+0.0001).AND.
     &			ABS(NFP(JP,1)*NFP(JP,2)).LT.1000000)THEN
		DO 620 JSB=1,15
		PSB=PP(JP,JSB)
		PP(JP,JSB)=PT(JT,JSB)
		PT(JT,JSB)=PSB
		NSB=NFP(JP,JSB)
		NFP(JP,JSB)=NFT(JT,JSB)
		NFT(JT,JSB)=NSB
620		CONTINUE
C		********when Ycm(JP)<Ycm(JT) after the collision
C			exchange the positions of the two   
	ENDIF
C
	RETURN
C**************************************************
C**************************************************
1000	IERROR=1
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) '	Fatal HIJSFT start error,abandon this event'
	WRITE(6,*) '	PROJ E+,E-,W+',EPP,EPM,WP
	WRITE(6,*) '	TARG E+,E-,W-',ETP,ETM,WM
	WRITE(6,*) '	W+*W-, (APN+ATN)^2',SW,SNN
	RETURN
2000	IERROR=0
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) '	(2)energy partition fail,'
	WRITE(6,*) ' 	HIJSFT not performed, but continue'
	WRITE(6,*) '	MP1,MPN',X1*(1.0-X2)*SW,AMPN**2
	WRITE(6,*) '	MT2,MTN',X2*(1.0-X1)*SW,AMTN**2
	RETURN
3000	IERROR=0
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) '	(3)something is wrong with the pt kick, '
	WRITE(6,*) '	HIJSFT not performed, but continue'
	WRITE(6,*) '	D1=',D1,' D2=',D2,' SW=',SW
	WRITE(6,*) '	HISTORY NFP5=',NFP(JP,5),' NFT5=',NFT(JT,5)
	WRITE(6,*) '	THIS COLLISON NFP5=',NFP5, ' NFT5=',NFT5
	WRITE(6,*) '	# OF JET IN PROJ',NPJ(JP),' IN TARG',NTJ(JT)
	RETURN
4000	IERROR=0
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) '	(4)unable to choose process, but not harmful'
	WRITE(6,*) '	HIJSFT not performed, but continue'
	WRITE(6,*) '	PTP=',SQRT(PTP2),' PTT=',SQRT(PTT2),' SW=',SW
	WRITE(6,*) '	AMCUT=',AMX,' JP=',JP,' JT=',JT
	WRITE(6,*) '	HISTORY NFP5=',NFP(JP,5),' NFT5=',NFT(JT,5)
	RETURN
5000	IERROR=0
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) '	energy partition failed(5),for limited try'
	WRITE(6,*) '    HIJSFT not performed, but continue'
	WRITE(6,*) '	NFP5=',NFP5,' NFT5=',NFT5
	WRITE(6,*) '	D1',D1,' X1(1-X2)',X1*(1.0-X2)
	WRITE(6,*) '	D2',D2,' X2(1-X1)',X2*(1.0-X1)
	RETURN
6000	PKC=0.0
	MISS=MISS+1
	IF(MISS.LT.100) GO TO 30
	IERROR=1
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) ' ERROR OCCURED, HIJSFT NOT PERFORMED'
        WRITE(6,*) ' Abort this event'
	WRITE(6,*) 'MTP,PTP2',EPP*EPM,PTP2,'  MTT,PTT2',ETP*ETM,PTT2 
	RETURN
	END
