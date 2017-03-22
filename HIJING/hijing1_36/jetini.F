* $Id$
C
C
C
C
C
	SUBROUTINE JETINI(JP,JT,I_TRIG)
C*******Initialize PYTHIA for jet production**********************
C	I_TRIG=0: for normal processes
C	I_TRIG=1: for triggered processes
C       JP: sequence number of the projectile
C       JT: sequence number of the target
C     For A+A collisions, one has to initilize pythia
C     separately for each type of collisions, pp, pn,np and nn,
C     or hp and hn for hA collisions. In this subroutine we use the following
C     catalogue for different type of collisions:
C     h+h: h+h (I_TYPE=1)
C     h+A: h+p (I_TYPE=1), h+n (I_TYPE=2)
C     A+h: p+h (I_TYPE=1), n+h (I_TYPE=2)
C     A+A: p+p (I_TYPE=1), p+n (I_TYPE=2), n+p (I_TYPE=3), n+n (I_TYPE=4)
C*****************************************************************
	CHARACTER BEAM*16,TARG*16
	DIMENSION XSEC0(8,0:200),COEF0(8,200,20),INI(8),
     &		MINT44(8),MINT45(8)
#define BLANKET_SAVE
#include "hijcrdn.inc"
#include "hiparnt.inc"
#include "histrng.inc"
#include "hipyint.inc"
C
#include "ludat1_hijing.inc"
#include "ludat3_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint5_hijing.inc"
	DATA INI/8*0/I_LAST/-1/
	SAVE
C
        IHNT2(11)=JP
        IHNT2(12)=JT
        IF(IHNT2(5).NE.0 .AND. IHNT2(6).NE.0) THEN
           I_TYPE=1
        ELSE IF(IHNT2(5).NE.0 .AND. IHNT2(6).EQ.0) THEN
           I_TYPE=1
           IF(NFT(JT,4).EQ.2112) I_TYPE=2
        ELSE IF(IHNT2(5).EQ.0 .AND. IHNT2(6).NE.0) THEN
           I_TYPE=1
           IF(NFP(JP,4).EQ.2112) I_TYPE=2
        ELSE
           IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2212) THEN
              I_TYPE=1
           ELSE IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2112) THEN
              I_TYPE=2
           ELSE IF(NFP(JP,4).EQ.2112 .AND. NFT(JT,4).EQ.2212) THEN
              I_TYPE=3
           ELSE
              I_TYPE=4
           ENDIF
        ENDIF
c

	IF(I_TRIG.NE.0) GO TO 160
        IF(I_TRIG.EQ.I_LAST) GO TO 150
        MSTP(2)=2
c			********second order running alpha_strong
        MSTP(33)=1
        PARP(31)=HIPR1(17)
C			********inclusion of K factor
        MSTP(51)=3
C			********Duke-Owens set 1 structure functions
        MSTP(61)=1
C			********INITIAL STATE RADIATION
        MSTP(71)=1
C			********FINAL STATE RADIATION
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.2) MSTP(61)=0
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.1) MSTP(71)=0
	MSTP(91) = 1
c
        MSTP(81)=0
C			******** NO MULTIPLE INTERACTION
        MSTP(82)=1
C			*******STRUCTURE OF MUTLIPLE INTERACTION
        MSTP(111)=0
C		********frag off(have to be done by local call)
        IF(IHPR2(10).EQ.0) MSTP(122)=0
C		********No printout of initialization information
        PARP(81)=HIPR1(8)
        CKIN(5)=HIPR1(8)
        CKIN(3)=HIPR1(8)
        CKIN(4)=HIPR1(9)
        IF(HIPR1(9).LE.HIPR1(8)) CKIN(4)=-1.0
        CKIN(9)=-10.0
        CKIN(10)=10.0
        MSEL=0
        DO 100 ISUB=1,200
           MSUB(ISUB)=0
 100    CONTINUE
        MSUB(11)=1
        MSUB(12)=1
        MSUB(13)=1
        MSUB(28)=1
        MSUB(53)=1
        MSUB(68)=1
	MSUB(81)=1 
	MSUB(82)=1
        DO 110 J=1,MIN(8,MDCY(21,3))
 110    MDME(MDCY(21,2)+J-1,1)=0
        ISEL=4
        IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
        MDME(MDCY(21,2)+ISEL-1,1)=1
C			********QCD subprocesses
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
C                       ******* direct photon production
 150    IF(INI(I_TYPE).NE.0) GO TO 800
	GO TO 400
C
C	*****triggered subprocesses, jet, photon, heavy quark and DY
C
 160    I_TYPE=4+I_TYPE
        IF(I_TRIG.EQ.I_LAST) GO TO 260
        PARP(81)=ABS(HIPR1(10))-0.25
        CKIN(5)=ABS(HIPR1(10))-0.25
        CKIN(3)=ABS(HIPR1(10))-0.25
        CKIN(4)=ABS(HIPR1(10))+0.25
        IF(HIPR1(10).LT.HIPR1(8)) CKIN(4)=-1.0
C	WRITE(6,*) "TRIGGERED", CKIN(3), CKIN(4), CKIN(5), IHPR2(3)
C
C ALICE
C
	IF (HIPR1(50) .EQ. 1) THEN
	   MSTP(61) = 0
	   MSTP(71) = 0
	   MSTP(91) = 0
	ENDIF

c
        MSEL=0
        DO 101 ISUB=1,200
           MSUB(ISUB)=0
 101    CONTINUE
        IF(IHPR2(3).EQ.1) THEN
           MSUB(11)=1
           MSUB(12)=1
           MSUB(13)=1
           MSUB(28)=1
           MSUB(53)=1
           MSUB(68)=1
           MSUB(81)=1
           MSUB(82)=1
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
           DO 102 J=1,MIN(8,MDCY(21,3))
 102	   MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
C			********QCD subprocesses
        ELSE IF(IHPR2(3).EQ.2) THEN
C	   WRITE(6,*) "DIRECT PHOTON PRODUCTION"
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
C		********Direct photon production
c		q+qbar->g+gamma,q+qbar->gamma+gamma, q+g->q+gamma
        ELSE IF(IHPR2(3).EQ.3) THEN
           CKIN(3)=MAX(0.0,HIPR1(10))
           CKIN(5)=HIPR1(8)
           PARP(81)=HIPR1(8)
           MSUB(81)=1
           MSUB(82)=1
           DO 105 J=1,MIN(8,MDCY(21,3))
 105	   MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
C             **********Heavy quark production
        ENDIF
260	IF(INI(I_TYPE).NE.0) GO TO 800
C
C
400	INI(I_TYPE)=1
	IF(IHPR2(10).EQ.0) MSTP(122)=0
	IF(NFP(JP,4).EQ.2212) THEN
		BEAM='P'
	ELSE IF(NFP(JP,4).EQ.-2212) THEN
		BEAM='P~'
	ELSE IF(NFP(JP,4).EQ.2112) THEN
		BEAM='N'
	ELSE IF(NFP(JP,4).EQ.-2112) THEN
		BEAM='N~'
	ELSE IF(NFP(JP,4).EQ.211) THEN
		BEAM='PI+'
	ELSE IF(NFP(JP,4).EQ.-211) THEN
		BEAM='PI-'
	ELSE IF(NFP(JP,4).EQ.321) THEN
		BEAM='PI+'
	ELSE IF(NFP(JP,4).EQ.-321) THEN
		BEAM='PI-'
	ELSE
		WRITE(6,*) 'unavailable beam type', NFP(JP,4)
	ENDIF
	IF(NFT(JT,4).EQ.2212) THEN
		TARG='P'
	ELSE IF(NFT(JT,4).EQ.-2212) THEN
		TARG='P~'
	ELSE IF(NFT(JT,4).EQ.2112) THEN
		TARG='N'
	ELSE IF(NFT(JT,4).EQ.-2112) THEN
		TARG='N~'
	ELSE IF(NFT(JT,4).EQ.211) THEN
		TARG='PI+'
	ELSE IF(NFT(JT,4).EQ.-211) THEN
		TARG='PI-'
	ELSE IF(NFT(JT,4).EQ.321) THEN
		TARG='PI+'
	ELSE IF(NFT(JT,4).EQ.-321) THEN
		TARG='PI-'
	ELSE
		WRITE(6,*) 'unavailable target type', NFT(JT,4)
	ENDIF
C
	IHNT2(16)=1
C       ******************indicate for initialization use when
C                         structure functions are called in PYTHIA
C
C	WRITE(6,*) "CALL PYINIT WITH", CKIN(3), CKIN(4), MSTP(61), MSTP(71), MSTP(91)
C ALICE ->
	IF (IHPR2(49) .EQ. 1) THEN
C SWITCH OFF MASSIVE FINAL QUARKS 
	   MSEL     = 0
C QCD on
	   MSUB(11) = 1
	   MSUB(12) = 1
	   MSUB(13) = 1
	   MSUB(28) = 1
	   MSUB(53) = 1
	   MSUB(68) = 1
C HF OFF
	   MSUB(81) = 0 
	   MSUB(82) = 0
	   MDME(127, 1) = 1
	   MDME(128, 1) = 1
	   MDME(129, 1) = 1
C g->QQbar off
	   MDME(130, 1) = 0
	   MDME(131, 1) = 0
C highest flavor in final state shower
	   MSTJ(45) = 3
C highest flavor in initial state shower and pdf
	   MSTP(54) = 3
	ENDIF
C <- ALICE
	CALL PYINIT_HIJING('CMS',BEAM,TARG,HINT1(1))
	MINT4=MINT(44)
	MINT5=MINT(45)
	MINT44(I_TYPE)=MINT(44)
	MINT45(I_TYPE)=MINT(45)
	ATXS(0)=XSEC(0,1)
	XSEC0(I_TYPE,0)=XSEC(0,1)
	DO 500 I=1,200
		ATXS(I)=XSEC(I,1)
		XSEC0(I_TYPE,I)=XSEC(I,1)
		DO 500 J=1,20
			ATCO(I,J)=COEF(I,J)
			COEF0(I_TYPE,I,J)=COEF(I,J)
500	CONTINUE
C
	IHNT2(16)=0
C
        I_LAST=I_TRIG
	RETURN
C		********Store the initialization information for
C				late use
C
C
800	MINT(44)=MINT44(I_TYPE)
	MINT(45)=MINT45(I_TYPE)
	MINT4=MINT(44)
	MINT5=MINT(45)
	XSEC(0,1)=XSEC0(I_TYPE,0)
	ATXS(0)=XSEC(0,1)
	DO 900 I=1,200
		XSEC(I,1)=XSEC0(I_TYPE,I)
		ATXS(I)=XSEC(I,1)
	DO 900 J=1,20
		COEF(I,J)=COEF0(I_TYPE,I,J)
		ATCO(I,J)=COEF(I,J)
900	CONTINUE
        I_LAST=I_TRIG
        MINT(11)=NFP(JP,4)
        MINT(12)=NFT(JT,4)
	RETURN
	END









