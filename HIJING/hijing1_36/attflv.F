* $Id$
C
C
C
	SUBROUTINE ATTFLV(ID,IDQ,IDQQ)
C
	SAVE
	IF(ABS(ID).LT.100) THEN
		NSIGN=1
		IDQ=ID/100
		IDQQ=-ID/10+IDQ*10
		IF(ABS(IDQ).EQ.3) NSIGN=-1
		IDQ=NSIGN*IDQ
		IDQQ=NSIGN*IDQQ
		IF(IDQ.LT.0) THEN
			ID0=IDQ
			IDQ=IDQQ
			IDQQ=ID0
		ENDIF
		RETURN
	ENDIF
C		********return ID of quark(IDQ) and anti-quark(IDQQ)
C			for pions and kaons
c
C	Return LU ID for quarks and diquarks for proton(ID=2212) 
C	anti-proton(ID=-2212) and nuetron(ID=2112)
C	LU ID for d=1,u=2, (ud)0=2101, (ud)1=2103, 
C       (dd)1=1103,(uu)1=2203.
C	Use SU(6)  weight  proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
C			  nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0
C 
	IDQ=2
	IF(ABS(ID).EQ.2112) IDQ=1
	IDQQ=2101
	X=RLU_HIJING(0)
	IF(X.LE.0.5) GO TO 30
	IF(X.GT.0.666667) GO TO 10
	IDQQ=2103
	GO TO 30
10	IDQ=1
	IDQQ=2203
	IF(ABS(ID).EQ.2112) THEN
		IDQ=2
		IDQQ=1103
	ENDIF
30	IF(ID.LT.0) THEN
		ID00=IDQQ
		IDQQ=-IDQ
		IDQ=-ID00
	ENDIF
	RETURN
	END	
