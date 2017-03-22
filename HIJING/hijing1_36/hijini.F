* $Id$
C            
C
C
	SUBROUTINE HIJINI
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "histrng.inc"
#include "hijjet1.inc"
#include "hijjet2.inc"
#include "hijjet4.inc"
	SAVE
C****************Reset the momentum of initial particles************
C             and assign flavors to the proj and targ string       *
C*******************************************************************
	NSG=0
	NDR=0
	IPP=2212
	IPT=2212
	IF(IHNT2(5).NE.0) IPP=IHNT2(5)
	IF(IHNT2(6).NE.0) IPT=IHNT2(6)
C		********in case the proj or targ is a hadron.
C
	DO 100 I=1,IHNT2(1)
	PP(I,1)=0.0
	PP(I,2)=0.0
	PP(I,3)=SQRT(HINT1(1)**2/4.0-HINT1(8)**2)
	PP(I,4)=HINT1(1)/2
	PP(I,5)=HINT1(8)
	PP(I,6)=0.0
	PP(I,7)=0.0
	PP(I,8)=0.0
	PP(I,9)=0.0
	PP(I,10)=0.0
	NFP(I,3)=IPP
	NFP(I,4)=IPP
	NFP(I,5)=0
	NFP(I,6)=0
	NFP(I,7)=0
	NFP(I,8)=0
	NFP(I,9)=0
	NFP(I,10)=0
	NFP(I,11)=0
	NPJ(I)=0
	IF(I.GT.ABS(IHNT2(2))) NFP(I,3)=2112
	CALL ATTFLV(NFP(I,3),IDQ,IDQQ)
	NFP(I,1)=IDQ
	NFP(I,2)=IDQQ
	NFP(I,15)=-1
	IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &		RLU_HIJING(0).LT.0.5)) NFP(I,15)=1
	PP(I,14)=ULMASS_HIJING(IDQ)
	PP(I,15)=ULMASS_HIJING(IDQQ)
100	CONTINUE
C
	DO 200 I=1,IHNT2(3)
	PT(I,1)=0.0
	PT(I,2)=0.0
	PT(I,3)=-SQRT(HINT1(1)**2/4.0-HINT1(9)**2)
	PT(I,4)=HINT1(1)/2.0
	PT(I,5)=HINT1(9)
	PT(I,6)=0.0
	PT(I,7)=0.0
	PT(I,8)=0.0
	PT(I,9)=0.0
	PT(I,10)=0.0
	NFT(I,3)=IPT
	NFT(I,4)=IPT
	NFT(I,5)=0
	NFT(I,6)=0
	NFT(I,7)=0
	NFT(I,8)=0
	NFT(I,9)=0
	NFT(I,10)=0
	NFT(I,11)=0
	NTJ(I)=0
	IF(I.GT.ABS(IHNT2(4))) NFT(I,3)=2112
	CALL ATTFLV(NFT(I,3),IDQ,IDQQ)
	NFT(I,1)=IDQ
	NFT(I,2)=IDQQ
	NFT(I,15)=1
	IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &			RLU_HIJING(0).LT.0.5)) NFT(I,15)=-1
	PT(I,14)=ULMASS_HIJING(IDQ)
	PT(I,15)=ULMASS_HIJING(IDQQ)
200	CONTINUE
	RETURN
	END
