* $Id$
C
C
C
	SUBROUTINE PARTON(F,X1,X2,QQ)
	IMPLICIT REAL*8(A-H,O-Z)
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "njet_hijing.inc"
	SAVE
	DIMENSION F(2,7) 
	DLAM=HIPR1(15)
	Q0=HIPR1(16)
	S=DLOG(DLOG(QQ/DLAM**2)/DLOG(Q0**2/DLAM**2))
	IF(IHPR2(7).EQ.2) GO TO 200
C*******************************************************
	AT1=0.419+0.004*S-0.007*S**2
	AT2=3.460+0.724*S-0.066*S**2
	GMUD=4.40-4.86*S+1.33*S**2
	AT3=0.763-0.237*S+0.026*S**2
	AT4=4.00+0.627*S-0.019*S**2
	GMD=-0.421*S+0.033*S**2
C*******************************************************
	CAS=1.265-1.132*S+0.293*S**2
	AS=-0.372*S-0.029*S**2
	BS=8.05+1.59*S-0.153*S**2
	APHS=6.31*S-0.273*S**2
	BTAS=-10.5*S-3.17*S**2
	GMS=14.7*S+9.80*S**2
C********************************************************
C	CAC=0.135*S-0.075*S**2
C	AC=-0.036-0.222*S-0.058*S**2
C	BC=6.35+3.26*S-0.909*S**2
C	APHC=-3.03*S+1.50*S**2
C	BTAC=17.4*S-11.3*S**2
C	GMC=-17.9*S+15.6*S**2
C***********************************************************
	CAG=1.56-1.71*S+0.638*S**2
	AG=-0.949*S+0.325*S**2
	BG=6.0+1.44*S-1.05*S**2
	APHG=9.0-7.19*S+0.255*S**2
	BTAG=-16.5*S+10.9*S**2
	GMG=15.3*S-10.1*S**2
	GO TO 300
C********************************************************
200	AT1=0.374+0.014*S
	AT2=3.33+0.753*S-0.076*S**2
	GMUD=6.03-6.22*S+1.56*S**2
	AT3=0.761-0.232*S+0.023*S**2
	AT4=3.83+0.627*S-0.019*S**2
	GMD=-0.418*S+0.036*S**2
C************************************
	CAS=1.67-1.92*S+0.582*S**2
	AS=-0.273*S-0.164*S**2
	BS=9.15+0.530*S-0.763*S**2
	APHS=15.7*S-2.83*S**2
	BTAS=-101.0*S+44.7*S**2
	GMS=223.0*S-117.0*S**2
C*********************************
C	CAC=0.067*S-0.031*S**2
C	AC=-0.120-0.233*S-0.023*S**2
C	BC=3.51+3.66*S-0.453*S**2
C	APHC=-0.474*S+0.358*S**2
C	BTAC=9.50*S-5.43*S**2
C	GMC=-16.6*S+15.5*S**2
C**********************************
	CAG=0.879-0.971*S+0.434*S**2
	AG=-1.16*S+0.476*S**2
	BG=4.0+1.23*S-0.254*S**2
	APHG=9.0-5.64*S-0.817*S**2
	BTAG=-7.54*S+5.50*S**2
	GMG=-0.596*S+1.26*S**2
C*********************************
300	B12=DEXP(GMRE(AT1)+GMRE(AT2+1.D0)-GMRE(AT1+AT2+1.D0))
	B34=DEXP(GMRE(AT3)+GMRE(AT4+1.D0)-GMRE(AT3+AT4+1.D0))
	CNUD=3.D0/B12/(1.D0+GMUD*AT1/(AT1+AT2+1.D0))
	CND=1.D0/B34/(1.D0+GMD*AT3/(AT3+AT4+1.D0))
C********************************************************
C	FUD=X*(U+D)
C	FS=X*2(UBAR+DBAR+SBAR)  AND UBAR=DBAR=SBAR
C*******************************************************
C PH    Make sure X1 belongs to [0,1]
	IF (X1.LT.0.D0) X1=0.D0
	IF (X1.GT.1.D0) X1=1.D0
	FUD1=CNUD*X1**AT1*(1.D0-X1)**AT2*(1.D0+GMUD*X1)
	FS1=CAS*X1**AS*(1.D0-X1)**BS*(1.D0+APHS*X1
     &      +BTAS*X1**2+GMS*X1**3)
	F(1,3)=CND*X1**AT3*(1.D0-X1)**AT4*(1.D0+GMD*X1)+FS1/6.D0
	F(1,1)=FUD1-F(1,3)+FS1/3.D0
	F(1,2)=FS1/6.D0
	F(1,4)=FS1/6.D0
	F(1,5)=FS1/6.D0
	F(1,6)=FS1/6.D0
	F(1,7)=CAG*X1**AG*(1.D0-X1)**BG*(1.D0+APHG*X1
     &         +BTAG*X1**2+GMG*X1**3)
C
C PH    Make sure X2 belongs to [0,1]
	IF (X2.LT.0.D0) X2=0.D0
	IF (X2.GT.1.D0) X2=1.D0
	FUD2=CNUD*X2**AT1*(1.D0-X2)**AT2*(1.D0+GMUD*X2)
	FS2=CAS*X2**AS*(1.D0-X2)**BS*(1.D0+APHS*X2
     &      +BTAS*X2**2+GMS*X2**3)
	F(2,3)=CND*X2**AT3*(1.D0-X2)**AT4*(1.D0+GMD*X2)+FS2/6.D0
	F(2,1)=FUD2-F(2,3)+FS2/3.D0
	F(2,2)=FS2/6.D0
	F(2,4)=FS2/6.D0
	F(2,5)=FS2/6.D0
	F(2,6)=FS2/6.D0
	F(2,7)=CAG*X2**AG*(1.D0-X2)**BG*(1.D0+APHG*X2
     &         +BTAG*X2**2+GMG*X2**3)
C***********Nuclear effect on the structure function****************
C
        IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
	   AAX=1.193*ALOG(FLOAT(IHNT2(1)))**0.16666666
	   RRX=AAX*(X1**3-1.2*X1**2+0.21*X1)+1.0
     &	       +1.079*(FLOAT(IHNT2(1))**0.33333333-1.0)
     &         /DLOG(IHNT2(1)+1.0D0)*DSQRT(X1)*DEXP(-X1**2/0.01)
	   IF(IP_CRS.EQ.1 .OR.IP_CRS.EQ.3) RRX=DEXP(-X1**2/0.01)
	   DO 400 I=1,7
	      F(1,I)=RRX*F(1,I)
 400	   CONTINUE
        ENDIF
        IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
	   AAX=1.193*ALOG(FLOAT(IHNT2(3)))**0.16666666
	   RRX=AAX*(X2**3-1.2*X2**2+0.21*X2)+1.0
     &	       +1.079*(FLOAT(IHNT2(3))**0.33333-1.0)
     &         /DLOG(IHNT2(3)+1.0D0)*DSQRT(X2)*DEXP(-X2**2/0.01)
	   IF(IP_CRS.EQ.2 .OR. IP_CRS.EQ.3) RRX=DEXP(-X2**2/0.01)
	   DO 500 I=1,7
	      F(2,I)=RRX*F(2,I)
 500	   CONTINUE
        ENDIF
c
	RETURN
	END
