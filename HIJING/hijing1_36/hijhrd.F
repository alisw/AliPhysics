* $Id$
C
C
C
	SUBROUTINE HIJHRD(JP,JT,JOUT,JFLG,IOPJET0)
C
C	IOPTJET=1, ALL JET WILL FORM SINGLE STRING SYSTEM
C		0, ONLY Q-QBAR JET FORM SINGLE STRING SYSTEM
C*******Perform jets production and fragmentation when JP JT *******
C     scatter. JOUT-> number of hard scatterings precede this one  *
C     for the the same pair(JP,JT). JFLG->a flag to show whether   *
C     jets can be produced (with valence quark=1,gluon=2, q-qbar=3)*
C     or not(0). Information of jets are in  COMMON/ATTJET and     *
C     /MINJET. ABS(NFP(JP,6)) is the total number jets produced by *
C    JP. If NFP(JP,6)<0 JP can not produce jet anymore.	           *
C*******************************************************************
	DIMENSION IP(100,2),IPQ(50),IPB(50),IT(100,2),ITQ(50),ITB(50)
#define BLANKET_SAVE
#include "hijcrdn.inc"
#include "hiparnt.inc"
#include "hijdat.inc"
#include "histrng.inc"
#include "hijjet1.inc"
#include "hijjet2.inc"
#include "hijjet4.inc"
C************************************ HIJING common block
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint5_hijing.inc"
#include "hipyint.inc"
	SAVE
C*********************************** LU common block
	MXJT=500
C		SIZE OF COMMON BLOCK FOR # OF PARTON PER STRING
	MXSG=900
C		SIZE OF COMMON BLOCK FOR # OF SINGLE STRINGS
	MXSJ=100
C		SIZE OF COMMON BLOCK FOR # OF PARTON PER SINGLE
C		STRING
	JFLG=0
	IHNT2(11)=JP
	IHNT2(12)=JT
C
	IOPJET=IOPJET0
	IF(IOPJET.EQ.1.AND.(NFP(JP,6).NE.0.OR.NFT(JT,6).NE.0))
     &                   IOPJET=0
	IF(JP.GT.IHNT2(1) .OR. JT.GT.IHNT2(3)) RETURN
	IF(NFP(JP,6).LT.0 .OR. NFT(JT,6).LT.0) RETURN
C		******** JP or JT can not produce jet anymore
C
	IF(JOUT.EQ.0) THEN
		EPP=PP(JP,4)+PP(JP,3)
		EPM=PP(JP,4)-PP(JP,3)
		ETP=PT(JT,4)+PT(JT,3)
		ETM=PT(JT,4)-PT(JT,3)
		IF(EPP.LT.0.0) GO TO 1000
		IF(EPM.LT.0.0) GO TO 1000
		IF(ETP.LT.0.0) GO TO 1000
		IF(ETM.LT.0.0) GO TO 1000
		IF(EPP/(EPM+0.01).LE.ETP/(ETM+0.01)) RETURN
	ENDIF
C		********for the first hard scattering of (JP,JT)
C			have collision only when Ycm(JP)>Ycm(JT)

	ECUT1=HIPR1(1)+HIPR1(8)+PP(JP,14)+PP(JP,15)
	ECUT2=HIPR1(1)+HIPR1(8)+PT(JT,14)+PT(JT,15)
	IF(PP(JP,4).LE.ECUT1) THEN
		NFP(JP,6)=-ABS(NFP(JP,6))
		RETURN
	ENDIF
	IF(PT(JT,4).LE.ECUT2) THEN
		NFT(JT,6)=-ABS(NFT(JT,6))
		RETURN
	ENDIF
C		*********must have enough energy to produce jets

	MISS=0
	MISP=0
	MIST=0
C
	IF(NFP(JP,10).EQ.0 .AND. NFT(JT,10).EQ.0) THEN
		MINT(44)=MINT4
		MINT(45)=MINT5
		XSEC(0,1)=ATXS(0)
		XSEC(11,1)=ATXS(11)
		XSEC(12,1)=ATXS(12)
		XSEC(28,1)=ATXS(28)
		DO 120 I=1,20
		COEF(11,I)=ATCO(11,I)
		COEF(12,I)=ATCO(12,I)
		COEF(28,I)=ATCO(28,I)
120		CONTINUE
	ELSE
		ISUB11=0
		ISUB12=0
		ISUB28=0
		IF(XSEC(11,1).NE.0) ISUB11=1
		IF(XSEC(12,1).NE.0) ISUB12=1
		IF(XSEC(28,1).NE.0) ISUB28=1		
		MINT(44)=MINT4-ISUB11-ISUB12-ISUB28
		MINT(45)=MINT5-ISUB11-ISUB12-ISUB28
		XSEC(0,1)=ATXS(0)-ATXS(11)-ATXS(12)-ATXS(28)
		XSEC(11,1)=0.0
		XSEC(12,1)=0.0
		XSEC(28,1)=0.0	
		DO 110 I=1,20
		COEF(11,I)=0.0
		COEF(12,I)=0.0
		COEF(28,I)=0.0
110		CONTINUE
	ENDIF		
C	********Scatter the valence quarks only once per NN 
C       collision,
C		afterwards only gluon can have hard scattering.
 155	CALL PYTHIA_HIJING
	JJ=MINT(31)
	IF(JJ.NE.1) GO TO 155
C		*********one hard collision at a time
	IF(K(7,2).EQ.-K(8,2)) THEN
		QMASS2=(P(7,4)+P(8,4))**2-(P(7,1)+P(8,1))**2
     &			-(P(7,2)+P(8,2))**2-(P(7,3)+P(8,3))**2
		QM=ULMASS_HIJING(K(7,2))
		IF(QMASS2.LT.(2.0*QM+HIPR1(1))**2) GO TO 155
	ENDIF
C		********q-qbar jets must has minimum mass HIPR1(1)
	PXP=PP(JP,1)-P(3,1)
	PYP=PP(JP,2)-P(3,2)
	PZP=PP(JP,3)-P(3,3)
	PEP=PP(JP,4)-P(3,4)
	PXT=PT(JT,1)-P(4,1)
	PYT=PT(JT,2)-P(4,2)
	PZT=PT(JT,3)-P(4,3)
	PET=PT(JT,4)-P(4,4)

	IF(PEP.LE.ECUT1) THEN
		MISP=MISP+1
		IF(MISP.LT.50) GO TO 155
		NFP(JP,6)=-ABS(NFP(JP,6))
		RETURN
	ENDIF
	IF(PET.LE.ECUT2) THEN
		MIST=MIST+1
		IF(MIST.LT.50) GO TO 155
		NFT(JT,6)=-ABS(NFT(JT,6))
		RETURN
	ENDIF
C		******** if the remain energy<ECUT the proj or targ
C			 can not produce jet anymore

	WP=PEP+PZP+PET+PZT
	WM=PEP-PZP+PET-PZT
	IF(WP.LT.0.0 .OR. WM.LT.0.0) THEN
		MISS=MISS+1
		IF(MISS.LT.50) GO TO 155
		RETURN
	ENDIF
C		********the total W+, W- must be positive
	SW=WP*WM
	AMPX=SQRT((ECUT1-HIPR1(8))**2+PXP**2+PYP**2+0.01)
	AMTX=SQRT((ECUT2-HIPR1(8))**2+PXT**2+PYT**2+0.01)
	SXX=(AMPX+AMTX)**2
	IF(SW.LT.SXX.OR.VINT(43).LT.HIPR1(1)) THEN
		MISS=MISS+1
		IF(MISS.LT.50) GO TO 155
		RETURN
	ENDIF  
C		********the proj and targ remnants must have at least
C			a CM energy that can produce two strings
C			with minimum mass HIPR1(1)(see HIJSFT HIJFRG)
C
	HINT1(41)=P(7,1)
	HINT1(42)=P(7,2)
	HINT1(43)=P(7,3)
	HINT1(44)=P(7,4)
	HINT1(45)=P(7,5)
	HINT1(46)=SQRT(P(7,1)**2+P(7,2)**2)
	HINT1(51)=P(8,1)
	HINT1(52)=P(8,2)
	HINT1(53)=P(8,3)
	HINT1(54)=P(8,4)
	HINT1(55)=P(8,5)
	HINT1(56)=SQRT(P(8,1)**2+P(8,2)**2) 
	IHNT2(14)=K(7,2)
	IHNT2(15)=K(8,2)
C
	PINIRAD=(1.0-EXP(-2.0*(VINT(47)-HIDAT(1))))
     &		/(1.0+EXP(-2.0*(VINT(47)-HIDAT(1))))
	I_INIRAD=0
	IF(RLU_HIJING(0).LE.PINIRAD) I_INIRAD=1
	IF(K(7,2).EQ.-K(8,2)) GO TO 190
	IF(K(7,2).EQ.21.AND.K(8,2).EQ.21.AND.IOPJET.EQ.1) GO TO 190
C*******************************************************************
C	gluon  jets are going to be connectd with
C	the final leading string of quark-aintquark
C*******************************************************************
	JFLG=2
	JPP=0
	LPQ=0
	LPB=0
	JTT=0
	LTQ=0
	LTB=0
	IS7=0
	IS8=0
        HINT1(47)=0.0
        HINT1(48)=0.0
        HINT1(49)=0.0
        HINT1(50)=0.0
        HINT1(67)=0.0
        HINT1(68)=0.0
        HINT1(69)=0.0
        HINT1(70)=0.0
	DO 180 I=9,N
	   IF(K(I,3).EQ.1 .OR. K(I,3).EQ.2.OR.
     &                   ABS(K(I,2)).GT.30) GO TO 180
C************************************************************
           IF(K(I,3).EQ.7) THEN
              HINT1(47)=HINT1(47)+P(I,1)
              HINT1(48)=HINT1(48)+P(I,2)
              HINT1(49)=HINT1(49)+P(I,3)
              HINT1(50)=HINT1(50)+P(I,4)
           ENDIF
           IF(K(I,3).EQ.8) THEN
              HINT1(67)=HINT1(67)+P(I,1)
              HINT1(68)=HINT1(68)+P(I,2)
              HINT1(69)=HINT1(69)+P(I,3)
              HINT1(70)=HINT1(70)+P(I,4)
           ENDIF
C************************modifcation made on Apr 10. 1996*****
	   IF(K(I,2).GT.21.AND.K(I,2).LE.30) THEN
	      NDR=NDR+1
	      IADR(NDR,1)=JP
	      IADR(NDR,2)=JT
	      KFDR(NDR)=K(I,2)
	      PDR(NDR,1)=P(I,1)
	      PDR(NDR,2)=P(I,2)
	      PDR(NDR,3)=P(I,3)
	      PDR(NDR,4)=P(I,4)
	      PDR(NDR,5)=P(I,5)
C************************************************************
	      GO TO 180
C************************correction made on Oct. 14,1994*****
	   ENDIF
	   IF(K(I,3).EQ.7.OR.K(I,3).EQ.3) THEN
	      IF(K(I,3).EQ.7.AND.K(I,2).NE.21.AND.K(I,2).EQ.K(7,2)
     &	             .AND.IS7.EQ.0) THEN
		 PP(JP,10)=P(I,1)
		 PP(JP,11)=P(I,2)
		 PP(JP,12)=P(I,3)
		 PZP=PZP+P(I,3)
		 PEP=PEP+P(I,4)
		 NFP(JP,10)=1
		 IS7=1
		 GO TO 180
	      ENDIF
	      IF(K(I,3).EQ.3.AND.(K(I,2).NE.21.OR.
     &                               I_INIRAD.EQ.0)) THEN
		 PXP=PXP+P(I,1)
		 PYP=PYP+P(I,2)
		 PZP=PZP+P(I,3)
		 PEP=PEP+P(I,4)
		 GO TO 180 
	      ENDIF
	      JPP=JPP+1
	      IP(JPP,1)=I
	      IP(JPP,2)=0
	      IF(K(I,2).NE.21) THEN
		 IF(K(I,2).GT.0) THEN
		    LPQ=LPQ+1
		    IPQ(LPQ)=JPP
		    IP(JPP,2)=LPQ
		 ELSE IF(K(I,2).LT.0) THEN
		    LPB=LPB+1
		    IPB(LPB)=JPP
		    IP(JPP,2)=-LPB
		 ENDIF
	      ENDIF
	   ELSE IF(K(I,3).EQ.8.OR.K(I,3).EQ.4) THEN
	      IF(K(I,3).EQ.8.AND.K(I,2).NE.21.AND.K(I,2).EQ.K(8,2)
     &				.AND.IS8.EQ.0) THEN
		 PT(JT,10)=P(I,1)
		 PT(JT,11)=P(I,2)
		 PT(JT,12)=P(I,3)
		 PZT=PZT+P(I,3)
		 PET=PET+P(I,4)
		 NFT(JT,10)=1
		 IS8=1
		 GO TO 180
	      ENDIF			
	      IF(K(I,3).EQ.4.AND.(K(I,2).NE.21.OR.
     &                             I_INIRAD.EQ.0)) THEN
		 PXT=PXT+P(I,1)
		 PYT=PYT+P(I,2)
		 PZT=PZT+P(I,3)
		 PET=PET+P(I,4)
		 GO TO 180
	      ENDIF
	      JTT=JTT+1
	      IT(JTT,1)=I
	      IT(JTT,2)=0
	      IF(K(I,2).NE.21) THEN
		 IF(K(I,2).GT.0) THEN
		    LTQ=LTQ+1
		    ITQ(LTQ)=JTT
		    IT(JTT,2)=LTQ
		 ELSE IF(K(I,2).LT.0) THEN
		    LTB=LTB+1
		    ITB(LTB)=JTT
		    IT(JTT,2)=-LTB
		 ENDIF
	      ENDIF
	   ENDIF
 180	CONTINUE
c
c
	IF(LPQ.NE.LPB .OR. LTQ.NE.LTB) THEN
		MISS=MISS+1
		IF(MISS.LE.50) GO TO 155
		WRITE(6,*) ' Q -QBAR NOT MATCHED IN HIJHRD'
		JFLG=0
		RETURN
	ENDIF
C****The following will rearrange the partons so that a quark is***
C****allways followed by an anti-quark ****************************

	J=0
181	J=J+1
	IF(J.GT.JPP) GO TO 182
	IF(IP(J,2).EQ.0) THEN
		GO TO 181
	ELSE IF(IP(J,2).NE.0) THEN
		LP=ABS(IP(J,2))
		IP1=IP(J,1)
		IP2=IP(J,2)
		IP(J,1)=IP(IPQ(LP),1)
		IP(J,2)=IP(IPQ(LP),2)
		IP(IPQ(LP),1)=IP1
		IP(IPQ(LP),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPQ(LP)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPQ(LP)
		ENDIF
C		********replace J with a quark
		IP1=IP(J+1,1)
		IP2=IP(J+1,2)
		IP(J+1,1)=IP(IPB(LP),1)
		IP(J+1,2)=IP(IPB(LP),2)
		IP(IPB(LP),1)=IP1
		IP(IPB(LP),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPB(LP)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPB(LP)
		ENDIF
C		******** replace J+1 with anti-quark
		J=J+1
		GO TO 181
	ENDIF

182	J=0
183	J=J+1
	IF(J.GT.JTT) GO TO 184
	IF(IT(J,2).EQ.0) THEN
		GO TO 183
	ELSE IF(IT(J,2).NE.0) THEN
		LT=ABS(IT(J,2))
		IT1=IT(J,1)
		IT2=IT(J,2)
		IT(J,1)=IT(ITQ(LT),1)
		IT(J,2)=IT(ITQ(LT),2)
		IT(ITQ(LT),1)=IT1
		IT(ITQ(LT),2)=IT2
		IF(IT2.GT.0) THEN
			ITQ(IT2)=ITQ(LT)
		ELSE IF(IT2.LT.0) THEN
			ITB(-IT2)=ITQ(LT)
		ENDIF
C		********replace J with a quark
		IT1=IT(J+1,1)
		IT2=IT(J+1,2)
		IT(J+1,1)=IT(ITB(LT),1)
		IT(J+1,2)=IT(ITB(LT),2)
		IT(ITB(LT),1)=IT1
		IT(ITB(LT),2)=IT2
		IF(IT2.GT.0) THEN
			ITQ(IT2)=ITB(LT)
		ELSE IF(IT2.LT.0) THEN
			ITB(-IT2)=ITB(LT)
		ENDIF
C		******** replace J+1 with anti-quark
		J=J+1
		GO TO 183

	ENDIF

184	CONTINUE
	IF(NPJ(JP)+JPP.GT.MXJT.OR.NTJ(JT)+JTT.GT.MXJT) THEN
		JFLG=0
		WRITE(6,*) 'number of partons per string exceeds'
		WRITE(6,*) 'the common block size'
		RETURN
	ENDIF
C			********check the bounds of common blocks
	DO 186 J=1,JPP
		KFPJ(JP,NPJ(JP)+J)=K(IP(J,1),2)
		PJPX(JP,NPJ(JP)+J)=P(IP(J,1),1)
		PJPY(JP,NPJ(JP)+J)=P(IP(J,1),2)
		PJPZ(JP,NPJ(JP)+J)=P(IP(J,1),3)
		PJPE(JP,NPJ(JP)+J)=P(IP(J,1),4)
		PJPM(JP,NPJ(JP)+J)=P(IP(J,1),5)
186	CONTINUE
	NPJ(JP)=NPJ(JP)+JPP
	DO 188 J=1,JTT
		KFTJ(JT,NTJ(JT)+J)=K(IT(J,1),2)
		PJTX(JT,NTJ(JT)+J)=P(IT(J,1),1)
		PJTY(JT,NTJ(JT)+J)=P(IT(J,1),2)
		PJTZ(JT,NTJ(JT)+J)=P(IT(J,1),3)
		PJTE(JT,NTJ(JT)+J)=P(IT(J,1),4)
		PJTM(JT,NTJ(JT)+J)=P(IT(J,1),5)
188	CONTINUE
	NTJ(JT)=NTJ(JT)+JTT
	GO TO 900
C*****************************************************************
CThis is the case of a quark-antiquark jet it will fragment alone
C****************************************************************
190	JFLG=3
	IF(K(7,2).NE.21.AND.K(8,2).NE.21.AND.
     &                   K(7,2)*K(8,2).GT.0) GO TO 155
	JPP=0
	LPQ=0
	LPB=0
        DO 200 I=9,N
	   IF(K(I,3).EQ.1.OR.K(I,3).EQ.2.OR.
     &                  ABS(K(I,2)).GT.30) GO TO 200
		IF(K(I,2).GT.21.AND.K(I,2).LE.30) THEN
			NDR=NDR+1
			IADR(NDR,1)=JP
			IADR(NDR,2)=JT
			KFDR(NDR)=K(I,2)
			PDR(NDR,1)=P(I,1)
			PDR(NDR,2)=P(I,2)
			PDR(NDR,3)=P(I,3)
			PDR(NDR,4)=P(I,4)
			PDR(NDR,5)=P(I,5)
C************************************************************
			GO TO 200
C************************correction made on Oct. 14,1994*****
		ENDIF
		IF(K(I,3).EQ.3.AND.(K(I,2).NE.21.OR.
     &                              I_INIRAD.EQ.0)) THEN
			PXP=PXP+P(I,1)
			PYP=PYP+P(I,2)
			PZP=PZP+P(I,3)
			PEP=PEP+P(I,4)
			GO TO 200
		ENDIF
		IF(K(I,3).EQ.4.AND.(K(I,2).NE.21.OR.
     &                                I_INIRAD.EQ.0)) THEN
			PXT=PXT+P(I,1)
			PYT=PYT+P(I,2)
			PZT=PZT+P(I,3)
			PET=PET+P(I,4)
			GO TO 200
		ENDIF
		JPP=JPP+1
		IP(JPP,1)=I
		IP(JPP,2)=0
		IF(K(I,2).NE.21) THEN
			IF(K(I,2).GT.0) THEN
				LPQ=LPQ+1
				IPQ(LPQ)=JPP
				IP(JPP,2)=LPQ
			ELSE IF(K(I,2).LT.0) THEN
				LPB=LPB+1
				IPB(LPB)=JPP
				IP(JPP,2)=-LPB
			ENDIF
		ENDIF
200	CONTINUE
	IF(LPQ.NE.LPB) THEN
	   MISS=MISS+1
	   IF(MISS.LE.50) GO TO 155
	   WRITE(6,*) LPQ,LPB, 'Q-QBAR NOT CONSERVED OR NOT MATCHED'
	   JFLG=0
	   RETURN
	ENDIF

C**** The following will rearrange the partons so that a quark is***
C**** allways followed by an anti-quark ****************************
	J=0
220	J=J+1
	IF(J.GT.JPP) GO TO 222
	IF(IP(J,2).EQ.0) GO TO 220
		LP=ABS(IP(J,2))
		IP1=IP(J,1)
		IP2=IP(J,2)
		IP(J,1)=IP(IPQ(LP),1)
		IP(J,2)=IP(IPQ(LP),2)
		IP(IPQ(LP),1)=IP1
		IP(IPQ(LP),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPQ(LP)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPQ(LP)
		ENDIF
		IPQ(LP)=J
C		********replace J with a quark
		IP1=IP(J+1,1)
		IP2=IP(J+1,2)
		IP(J+1,1)=IP(IPB(LP),1)
		IP(J+1,2)=IP(IPB(LP),2)
		IP(IPB(LP),1)=IP1
		IP(IPB(LP),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPB(LP)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPB(LP)
		ENDIF
C		******** replace J+1 with an anti-quark
		IPB(LP)=J+1
		J=J+1
		GO TO 220

222	CONTINUE
	IF(LPQ.GE.1) THEN
		DO 240 L0=2,LPQ
			IP1=IP(2*L0-3,1)
			IP2=IP(2*L0-3,2)
			IP(2*L0-3,1)=IP(IPQ(L0),1)
			IP(2*L0-3,2)=IP(IPQ(L0),2)
			IP(IPQ(L0),1)=IP1
			IP(IPQ(L0),2)=IP2
			IF(IP2.GT.0) THEN
				IPQ(IP2)=IPQ(L0)
			ELSE IF(IP2.LT.0) THEN
				IPB(-IP2)=IPQ(L0)
			ENDIF
			IPQ(L0)=2*L0-3
C
			IP1=IP(2*L0-2,1)
			IP2=IP(2*L0-2,2)
			IP(2*L0-2,1)=IP(IPB(L0),1)
			IP(2*L0-2,2)=IP(IPB(L0),2)
			IP(IPB(L0),1)=IP1
			IP(IPB(L0),2)=IP2
			IF(IP2.GT.0) THEN
				IPQ(IP2)=IPB(L0)
			ELSE IF(IP2.LT.0) THEN
				IPB(-IP2)=IPB(L0)
			ENDIF
			IPB(L0)=2*L0-2
240		CONTINUE
C		********move all the qqbar pair to the front of 
C				the list, except the first pair
		IP1=IP(2*LPQ-1,1)
		IP2=IP(2*LPQ-1,2)
		IP(2*LPQ-1,1)=IP(IPQ(1),1)
		IP(2*LPQ-1,2)=IP(IPQ(1),2)
		IP(IPQ(1),1)=IP1
		IP(IPQ(1),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPQ(1)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPQ(1)
		ENDIF
		IPQ(1)=2*LPQ-1
C		********move the first quark to the beginning of
C				the last string system
		IP1=IP(JPP,1)
		IP2=IP(JPP,2)
		IP(JPP,1)=IP(IPB(1),1)
		IP(JPP,2)=IP(IPB(1),2)
		IP(IPB(1),1)=IP1
		IP(IPB(1),2)=IP2
		IF(IP2.GT.0) THEN
			IPQ(IP2)=IPB(1)
		ELSE IF(IP2.LT.0) THEN
			IPB(-IP2)=IPB(1)
		ENDIF
		IPB(1)=JPP
C		********move the first anti-quark to the end of the 
C			last string system
	ENDIF
	IF(NSG.GE.MXSG) THEN
	   JFLG=0
	   WRITE(6,*) 'number of jets forming single strings exceeds'
	   WRITE(6,*) 'the common block size'
	   RETURN
	ENDIF
	IF(JPP.GT.MXSJ) THEN
	   JFLG=0
	   WRITE(6,*) 'number of partons per single jet system'
	   WRITE(6,*) 'exceeds the common block size'
	   RETURN
	ENDIF
C		********check the bounds of common block size
	NSG=NSG+1
	NJSG(NSG)=JPP
	IASG(NSG,1)=JP
	IASG(NSG,2)=JT
	IASG(NSG,3)=0
	DO 300 I=1,JPP
		K1SG(NSG,I)=2
		K2SG(NSG,I)=K(IP(I,1),2)
		IF(K2SG(NSG,I).LT.0) K1SG(NSG,I)=1
		PXSG(NSG,I)=P(IP(I,1),1)
		PYSG(NSG,I)=P(IP(I,1),2)
		PZSG(NSG,I)=P(IP(I,1),3)
		PESG(NSG,I)=P(IP(I,1),4)
		PMSG(NSG,I)=P(IP(I,1),5)
300	CONTINUE
	K1SG(NSG,1)=2
	K1SG(NSG,JPP)=1
C******* reset the energy-momentum of incoming particles ********
900	PP(JP,1)=PXP
	PP(JP,2)=PYP
	PP(JP,3)=PZP
	PP(JP,4)=PEP
	PP(JP,5)=0.0
	PT(JT,1)=PXT
	PT(JT,2)=PYT
	PT(JT,3)=PZT
	PT(JT,4)=PET
	PT(JT,5)=0.0

	NFP(JP,6)=NFP(JP,6)+1
	NFT(JT,6)=NFT(JT,6)+1
	RETURN
C
1000	JFLG=-1
	IF(IHPR2(10).EQ.0) RETURN
	WRITE(6,*) 'Fatal HIJHRD error'
	WRITE(6,*) JP, ' proj E+,E-',EPP,EPM,' status',NFP(JP,5)
	WRITE(6,*) JT, ' targ E+,E_',ETP,ETM,' status',NFT(JT,5)
	RETURN
	END
