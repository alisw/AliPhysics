* $Id$
C*************************
C
C
C
C
C ********************************************************
C ************************              WOOD-SAX
        SUBROUTINE HIJWDS(IA,IDH,XHIGH)
C     SETS UP HISTOGRAM IDH WITH RADII FOR
C     NUCLEUS IA DISTRIBUTED ACCORDING TO THREE PARAM WOOD SAXON
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "wood_hijing.inc"
        DIMENSION IAA(20),RR(20),DD(20),WW(20),RMS(20)
        EXTERNAL RWDSAX,WDSAX

C
C   PARAMETERS OF SPECIAL NUCLEI FROM ATOMIC DATA AND NUC DATA TABLES
C     VOL 14, 5-6 1974
        DATA IAA/2,4,12,16,27,32,40,56,63,93,184,197,208,7*0./
        DATA RR/0.01,.964,2.355,2.608,2.84,3.458,3.766,3.971,4.214,
     1        4.87,6.51,6.38,6.624,7*0./
        DATA DD/0.5882,.322,.522,.513,.569,.61,.586,.5935,.586,.573,
     1        .535,.535,.549,7*0./
        DATA WW/0.0,.517,-0.149,-0.051,0.,-0.208,-0.161,13*0./
        DATA RMS/2.11,1.71,2.46,2.73,3.05,3.247,3.482,3.737,3.925,4.31,
     1        5.42,5.33,5.521,7*0./
        SAVE
C
      	A=IA
C
C 		********SET WOOD-SAX PARAMS FIRST  AS IN DATE ET AL
      	D=0.54
C			********D IS WOOD SAX DIFFUSE PARAM IN FM
	R=1.19*A**(1./3.) - 1.61*A**(-1./3.)
C 			********R IS RADIUS PARAM
	W=0.
C 		********W IS The third of three WOOD-SAX PARAM
C
C      		********CHECK TABLE FOR SPECIAL CASES
	DO 10 I=1,13
		IF (IA.EQ.IAA(I)) THEN
			R=RR(I)
     			D=DD(I)
      			W=WW(I)
      			RS=RMS(I)
      		END IF
10    	CONTINUE
C     			********FNORM is the normalize factor
      	FNORM=1.0
      	XLOW=0.
      	XHIGH=R+ 12.*D
      	IF (W.LT.-0.01)  THEN
      		IF (XHIGH.GT.R/SQRT(ABS(W))) XHIGH=R/SQRT(ABS(W))
      	END IF
      	FGAUS=GAUSS1(RWDSAX,XLOW,XHIGH,0.001)
      	FNORM=1./FGAUS
C
        IF (IDH.EQ.1) THEN
           HINT1(72)=R
           HINT1(73)=D
           HINT1(74)=W
           HINT1(75)=FNORM/4.0/HIPR1(40)
        ELSE IF (IDH.EQ.2) THEN
           HINT1(76)=R
           HINT1(77)=D
           HINT1(78)=W
           HINT1(79)=FNORM/4.0/HIPR1(40)
        ENDIF
C
C     	NOW SET UP HBOOK FUNCTIONS IDH FOR  R**2*RHO(R)
C     	THESE HISTOGRAMS ARE USED TO GENERATE RANDOM RADII
      	CALL HIFUN(IDH,XLOW,XHIGH,RWDSAX)
      	RETURN
      	END
