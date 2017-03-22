* $Id$
C
C
C
C********************************************************************
C	Sort the jets associated with a nucleon in order of their
C	rapdities
C********************************************************************
	SUBROUTINE HIJSRT(JPJT,NPT)
	DIMENSION KF(100),PX(100),PY(100),PZ(100),PE(100),PM(100)
	DIMENSION Y(100),IP(100,2)
#define BLANKET_SAVE
#include "hijjet1.inc"
	SAVE
	IF(NPT.EQ.2) GO TO 500
	JP=JPJT
	IQ=0
	I=1
100	KF(I)=KFPJ(JP,I)
	PX(I)=PJPX(JP,I)
	PY(I)=PJPY(JP,I)
	PZ(I)=PJPZ(JP,I)
	PE(I)=PJPE(JP,I)
	PM(I)=PJPM(JP,I)
	Y(I-IQ)=0.5*ALOG((ABS(PE(I)+PZ(I))+1.E-5)
     &          /(ABS(PE(I)-PZ(I))+1.E-5))
	IP(I-IQ,1)=I
	IP(I-IQ,2)=0
	IF(KF(I).NE.21) THEN
		IP(I-IQ,2)=1
		IQ=IQ+1
		I=I+1
		KF(I)=KFPJ(JP,I)
		PX(I)=PJPX(JP,I)
		PY(I)=PJPY(JP,I)
		PZ(I)=PJPZ(JP,I)
		PE(I)=PJPE(JP,I)
		PM(I)=PJPM(JP,I)
	ENDIF
	I=I+1
	IF(I.LE.NPJ(JP)) GO TO 100
			
	DO 200 I=1,NPJ(JP)-IQ
	DO 200 J=I+1,NPJ(JP)-IQ
		IF(Y(I).GT.Y(J)) GO TO 200
		IP1=IP(I,1)
		IP2=IP(I,2)
		IP(I,1)=IP(J,1)
		IP(I,2)=IP(J,2)
		IP(J,1)=IP1
		IP(J,2)=IP2
200	CONTINUE
C			********sort in decending y
	IQQ=0
	I=1
300	KFPJ(JP,I)=KF(IP(I-IQQ,1))
	PJPX(JP,I)=PX(IP(I-IQQ,1))
	PJPY(JP,I)=PY(IP(I-IQQ,1))
	PJPZ(JP,I)=PZ(IP(I-IQQ,1))
	PJPE(JP,I)=PE(IP(I-IQQ,1))
	PJPM(JP,I)=PM(IP(I-IQQ,1))
	IF(IP(I-IQQ,2).EQ.1) THEN
		KFPJ(JP,I+1)=KF(IP(I-IQQ,1)+1)
		PJPX(JP,I+1)=PX(IP(I-IQQ,1)+1)
		PJPY(JP,I+1)=PY(IP(I-IQQ,1)+1)
		PJPZ(JP,I+1)=PZ(IP(I-IQQ,1)+1)
		PJPE(JP,I+1)=PE(IP(I-IQQ,1)+1)
		PJPM(JP,I+1)=PM(IP(I-IQQ,1)+1)
		I=I+1
		IQQ=IQQ+1
	ENDIF
	I=I+1
	IF(I.LE.NPJ(JP)) GO TO 300

	RETURN

500	JT=JPJT
	IQ=0
	I=1
600	KF(I)=KFTJ(JT,I)
	PX(I)=PJTX(JT,I)
	PY(I)=PJTY(JT,I)
	PZ(I)=PJTZ(JT,I)
	PE(I)=PJTE(JT,I)
	PM(I)=PJTM(JT,I)
	Y(I-IQ)=0.5*ALOG((ABS(PE(I)+PZ(I))+1.E-5)
     &          /(ABS(PE(I)-PZ(I))+1.E-5))
	IP(I-IQ,1)=I
	IP(I-IQ,2)=0
	IF(KF(I).NE.21) THEN
		IP(I-IQ,2)=1
		IQ=IQ+1
		I=I+1
		KF(I)=KFTJ(JT,I)
		PX(I)=PJTX(JT,I)
		PY(I)=PJTY(JT,I)
		PZ(I)=PJTZ(JT,I)
		PE(I)=PJTE(JT,I)
		PM(I)=PJTM(JT,I)
	ENDIF
	I=I+1
	IF(I.LE.NTJ(JT)) GO TO 600
			
	DO 700 I=1,NTJ(JT)-IQ
	DO 700 J=I+1,NTJ(JT)-IQ
		IF(Y(I).LT.Y(J)) GO TO 700
		IP1=IP(I,1)
		IP2=IP(I,2)
		IP(I,1)=IP(J,1)
		IP(I,2)=IP(J,2)
		IP(J,1)=IP1
		IP(J,2)=IP2
700	CONTINUE
C			********sort in acending y
	IQQ=0
	I=1
800	KFTJ(JT,I)=KF(IP(I-IQQ,1))
	PJTX(JT,I)=PX(IP(I-IQQ,1))
	PJTY(JT,I)=PY(IP(I-IQQ,1))
	PJTZ(JT,I)=PZ(IP(I-IQQ,1))
	PJTE(JT,I)=PE(IP(I-IQQ,1))
	PJTM(JT,I)=PM(IP(I-IQQ,1))
	IF(IP(I-IQQ,2).EQ.1) THEN
		KFTJ(JT,I+1)=KF(IP(I-IQQ,1)+1)
		PJTX(JT,I+1)=PX(IP(I-IQQ,1)+1)
		PJTY(JT,I+1)=PY(IP(I-IQQ,1)+1)
		PJTZ(JT,I+1)=PZ(IP(I-IQQ,1)+1)
		PJTE(JT,I+1)=PE(IP(I-IQQ,1)+1)
		PJTM(JT,I+1)=PM(IP(I-IQQ,1)+1)
		I=I+1
		IQQ=IQQ+1
	ENDIF
	I=I+1
	IF(I.LE.NTJ(JT)) GO TO 800
	RETURN
	END	
