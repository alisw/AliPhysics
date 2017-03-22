* $Id$
C
C
C	THIS PROGRAM IS TO CALCULATE THE JET CROSS SECTION
C	THE INTEGRATION IS DONE BY USING VEGAS
C
	SUBROUTINE CRSJET
	IMPLICIT REAL*8(A-H,O-Z)
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "njet_hijing.inc"
#include "bveg1_hijing.inc"
#include "bveg2_hijing.inc"
#include "bveg3_hijing.inc"
	EXTERNAL FJET,FJETRIG
	SAVE
C
c************************
c	NCALL give the number of inner-iteration, ITMX 
C       gives the limit of out-iteration. Nprn is an option
C       ( 1: print the integration process. 0: do not print)
C
	NDIM=3
	IP_CRS=0
	CALL VEGAS(FJET,AVGI,SD,CHI2A)
	HINT1(14)=AVGI/2.5682
	IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
		IP_CRS=1
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(15)=AVGI/2.5682
	ENDIF
	IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
		IP_CRS=2
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(16)=AVGI/2.5682
	ENDIF
	IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
		IP_CRS=3
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(17)=AVGI/2.5682
	ENDIF
C		********Total inclusive jet cross section(Pt>P0) 
C
	IF(IHPR2(3).NE.0) THEN
	   IP_CRS=0
	   CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	   HINT1(61)=AVGI/2.5682
	   IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
	      IP_CRS=1
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(62)=AVGI/2.5682
	   ENDIF
	   IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
	      IP_CRS=2
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(63)=AVGI/2.5682
	   ENDIF
	   IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
	      IP_CRS=3
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(64)=AVGI/2.5682
	   ENDIF
	ENDIF
C			********cross section of trigger jet
C
	RETURN
	END
