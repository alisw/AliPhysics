*CMZ :          22/02/99  17.57.04  by  Federico Carminati
*CMZ :  2.03/01 10/09/98  11.48.45  by  Roberto Barbera
*CMZ :  2.00/05 25/05/98  14.39.00  by  Federico Carminati
*CMZ :  2.00/04 22/05/98  09.37.53  by  Roberto Barbera (Catania)
*CMZ :  2.00/02 24/04/98  13.35.44  by  Federico Carminati
*CMZ :  2.00/01 22/04/98  19.18.39  by  Federico Carminati
*CMZ :  1.05/06 24/10/95  14.00.47  by  Nick van Eijndhoven (RUU/CERN)
*CMZ :  1.05/00 10/12/92  16.45.55  by  Nick van Eijndhoven (RUU/CERN)
*-- Author :    Nick van Eijndhoven (CERN)   24/09/90
      SUBROUTINE ITS_GEO4(IVERS,IDTMED)
C
C *** DEFINITION OF THE GEOMETRY OF THE ITS ***
C *** RB APR-1990 CATANIA ***
C
C CALLED BY : SXGEOM
C ORIGIN    : ROBERTO BARBERA
C
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
*KEND.
      DIMENSION IDTMED(*)
C
      DIMENSION IDROTM(5250)
C
      CHARACTER*4 NATRA(8),NATRA1(16),NATRA2(6),NATRA3(6),NATRA4(16)
C
      REAL*8 BIGA,BIGB,BIGA1,BIGB1
      REAL*8 COEFFA,COEFFB,COEFFC
      REAL*8 XCC1,YCC1,XCC2,YCC2
C
      DIMENSION DITS(3),DBUS(3),DCHI(3)
      DIMENSION DSUP(3),DTUB(3),DWAT(3)
      DIMENSION DCEI(3),DFRA(10)
      DIMENSION DAL1(3),DAL2(3),DKAP(3)
      DIMENSION DPCB(3),DCOP(3),DCER(3),DSIL(3)
      DIMENSION DPLA(3),DEPX(3)
      DIMENSION DGH(15)
      DIMENSION XX(13),YY(13),XBEG(12),YBEG(12),XEND(12),YEND(12)	
      DIMENSION DARC(5),XARC(12),YARC(12),RARC(12)
      DIMENSION DTRA(3),DTRA1(3),DTRA2(3),DTRA3(3),DTRA4(3)
      DIMENSION XTRA(8),YTRA(8),ZTRA(8)
      DIMENSION XTRA1(6),YTRA1(6),ZTRA1(6)
      DIMENSION DCONE(5),DTUBE(3),DPGON(10)
      DIMENSION DBOX1(3),DBOX2(3),DSRV(3),DELA(3)
C
      DATA NATRA/'TR01','TR02','TR03','TR04','TR05','TR06','TR07',
     $           'TR08'/
      DATA NATRA1/'TR11','TR12','TR13','TR14','TR15','TR16','TR17',
     $            'TR18','TR19','TR20','TR21','TR22','TR23','TR24',
     $            'TR25','TR26'/
      DATA NATRA2/'TR31','TR32','TR33','TR34','TR35','TR36'/
      DATA NATRA3/'TR41','TR42','TR43','TR44','TR45','TR46'/
      DATA NATRA4/'TR51','TR52','TR53','TR54','TR55','TR56','TR57',
     $            'TR58','TR59','TR60','TR61','TR62','TR63','TR64',
     $            'TR65','TR66'/
C
	DATA XX/0.,0.,-4.824,-4.833,-22.167,-22.585,-28.070,-27.626,
     $  -38.139,-19.749,-13.449,-14.726,0./
        DATA YY/0.,34.028,34.064,32.594,29.984,30.914,29.180,27.777,
     $  22.522,-7.918,-4.769,-2.216,0./
	DATA XBEG/0.,-0.497,-4.827,-5.425,-22.395,-23.009,-27.919,
     $  -28.161,-37.859,-19.337,-13.673,-13.773/
        DATA YBEG/0.430,34.032,33.564,32.505,30.493,30.781,28.702,
     $  27.509,22.058,-7.711,-4.322,-2.073/
        DATA XEND/0.,-4.323,-4.830,-21.643,-22.395,-27.593,-27.807,
     $  -37.655,-19.988,-13.897,-14.294,-0.425/
        DATA YEND/33.531,34.060,33.192,30.064,30.493,29.330,28.347,
     $  22.764,-7.523,-4.992,-3.077,-0.064/
	DATA XARC/-19.560,-14.120,-13.669,-0.500,-0.500,-4.327,-5.529,
     $  -21.746,-22.858,-27.442,-28.474,-37.431/
        DATA YARC/-7.264,-4.545,-2.765,0.431,33.531,33.560,33.197,
     $  30.756,30.304,28.854,28.136,22.316/
	DATA RARC/0.5,0.5,0.7,0.5,0.5,0.5,0.7,0.7,0.5,0.5,0.7,0.5/
C
	DATA RR,TTETA,PPHI,GTETA/4.08332,63.,-27.,89.4/  ! 89.4 = 90.-0.6
C
C --- Define ghost volume containing the whole ITS and fill it with air
C     or vacuum
C
	DGH(1)=0.
        DGH(2)=360.
        DGH(3)=4.
        DGH(4)=-70.
        DGH(5)=49.999
        DGH(6)=49.999
        DGH(7)=-25.
        DGH(8)=3.
        DGH(9)=49.999
        DGH(10)=25.
        DGH(11)=3.
        DGH(12)=49.999
        DGH(13)=70.
        DGH(14)=49.999
        DGH(15)=49.999	
        CALL GSVOLU('ITSV','PCON',IDTMED(276),DGH,15,IOUT)
C
C --- Place the ghost volume in its mother volume (ALIC) and make it
C     invisible
C
        CALL GSPOS('ITSV',1,'ALIC',0.,0.,0.,0,'ONLY')	
 	CALL GSATT('ITSV','SEEN',0)
C	
C ************************************************************************
C *                                                                      *
C *                               P I X E L S                            *
C *                               ===========                            *
C *                                                                      *
C ************************************************************************
C
C      GOTO 2345             ! skip ITS layer no. 1 and 2
C
C --- Define ghost volume containing the Pixel Detectors and fill it with air
C     or vacuum
C
        XXM=(49.999-3.)/(70.-25.)
        DGH(1)=0.
        DGH(2)=360.
        DGH(3)=4.
        DGH(4)=-25.-(9.-3.01)/XXM
        DGH(5)=9.
        DGH(6)=9.
        DGH(7)=-25.
        DGH(8)=3.01
        DGH(9)=9.
        DGH(10)=25.
        DGH(11)=3.01
        DGH(12)=9.
        DGH(13)=25.+(9.-3.01)/XXM
        DGH(14)=9.
        DGH(15)=9.
        CALL GSVOLU('IT12','PCON',IDTMED(276),DGH,15,IOUT)
C
C --- Place the ghost volume in its mother volume (ITSV) and make it
C     invisible
C
        CALL GSPOS('IT12',1,'ITSV',0.,0.,0.,0,'ONLY')
        CALL GSATT('IT12','SEEN',0)
C
C --- Define a ghost volume containing a single element of layer #1
C     and fill it with air or vacuum
C
        DBOX1(1)=0.005+0.01+0.0075
	DBOX1(2)=0.79
	DBOX1(3)=12.67
	CALL GSVOLU('IPV1','BOX ',IDTMED(204),DBOX1,3,IOUT)
C
C --- Divide each element of layer #1 in three ladders along the beam direction
C
	CALL GSDVN('IPB1','IPV1',3,3)
C
C --- Make the ghost volumes invisible
C
        CALL GSATT('IPV1','SEEN',0)
        CALL GSATT('IPB1','SEEN',0)  	
C
C --- Define a volume containing the chip of pixels (silicon, layer #1)
C
        DCHI(1)=0.005	
        DCHI(2)=0.79
        DCHI(3)=DBOX1(3)/3.
	CALL GSVOLU('ICH1','BOX ',IDTMED(201),DCHI,3,IOUT)
C
C --- Define a volume containing the bus of pixels (silicon, layer #1)
C
        DBUS(1)=0.01
        DBUS(2)=0.64
        DBUS(3)=4.19
	CALL GSVOLU('IBU1','BOX ',IDTMED(202),DBUS,3,IOUT)
C
C --- Define a volume containing the sensitive part of pixels
C     (silicon, layer #1)
C	
        DITS(1)=0.0075
        DITS(2)=0.64
        DITS(3)=4.19
	CALL GSVOLU('ITS1','BOX ',IDTMED(200),DITS,3,IOUT)
C
C --- Place the chip into its mother (IPB1)
C
        XPOS=DBOX1(1)-DCHI(1)
	YPOS=0.
	ZPOS=0.
	CALL GSPOS('ICH1',1,'IPB1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive volume into its mother (IPB1)
C
        XPOS=DBOX1(1)-2.*DCHI(1)-DITS(1)
	YPOS=DCHI(2)-DITS(2)
	ZPOS=-(DCHI(3)-DITS(3))
	CALL GSPOS('ITS1',1,'IPB1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the bus into its mother (IPB1)
C
        XPOS=DBOX1(1)-2.*DCHI(1)-2.*DITS(1)-DBUS(1)
	YPOS=DCHI(2)-DBUS(2)
	ZPOS=-(DCHI(3)-DBUS(3))
	CALL GSPOS('IBU1',1,'IPB1',XPOS,YPOS,ZPOS,0,'ONLY')	
C
C --- Define a ghost volume containing a single element of layer #2
C     and fill it with air or vacuum
C
        DBOX2(1)=0.005+0.01+0.0075
	DBOX2(2)=0.79
	DBOX2(3)=16.91
	CALL GSVOLU('IPV2','BOX ',IDTMED(204),DBOX2,3,IOUT)
C
C --- Divide each element of layer #2 in four ladders along the beam direction
C
	CALL GSDVN('IPB2','IPV2',4,3)
C
C --- Make the ghost volumes invisible
C
        CALL GSATT('IPV2','SEEN',0)
        CALL GSATT('IPB2','SEEN',0) 	
C
C --- Define a volume containing the chip of pixels (silicon, layer #2)
C
        DCHI(1)=0.005	
        DCHI(2)=0.79
        DCHI(3)=DBOX2(3)/4.
	CALL GSVOLU('ICH2','BOX ',IDTMED(201),DCHI,3,IOUT)
C
C --- Define a volume containing the bus of pixels (silicon, layer #2)
C
        DBUS(1)=0.01
        DBUS(2)=0.64
        DBUS(3)=4.19
	CALL GSVOLU('IBU2','BOX ',IDTMED(202),DBUS,3,IOUT)
C
C --- Define a volume containing the sensitive part of pixels
C     (silicon, layer #2)
C
        DITS(1)=0.0075
        DITS(2)=0.64
        DITS(3)=4.19
	CALL GSVOLU('ITS2','BOX ',IDTMED(200),DITS,3,IOUT)
C
C --- Place the chip into its mother (IPB2)
C
        XPOS=DBOX1(1)-2.*DBUS(1)-2.*DITS(1)-DCHI(1)
	YPOS=0.
	ZPOS=0.
	CALL GSPOS('ICH2',1,'IPB2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive volume into its mother (IPB2)
C
        XPOS=DBOX1(1)-2.*DBUS(1)-DITS(1)
	YPOS=-(DCHI(2)-DITS(2))
	ZPOS=-(DCHI(3)-DITS(3))
	CALL GSPOS('ITS2',1,'IPB2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the bus into its mother (IPB2)
C
        XPOS=DBOX1(1)-DBUS(1)
	YPOS=-(DCHI(2)-DBUS(2))
	ZPOS=-(DCHI(3)-DBUS(3))
	CALL GSPOS('IBU2',1,'IPB2',XPOS,YPOS,ZPOS,0,'ONLY')		
C
C --- Define a generic segment of an element of the mechanical support
C
        DSUP(1)=0.
	DSUP(2)=0.
	DSUP(3)=0.
        CALL GSVOLU('SPIX','BOX ',IDTMED(203),DSUP,0,IOUT)
C
C --- Define a generic arc of an element of the mechanical support
C
        DARC(1)=0.
	DARC(2)=0.
	DARC(3)=0.
	CALL GSVOLU('SARC','TUBS',IDTMED(203),DARC,0,IOUT)
C
C --- Define the mechanical supports of layers #1 and #2 and place the
C     elements of the layers in it
C
      JBOX1=0        ! counter over the number of elements of layer #1 (1-20)
      JBOX2=0        ! counter over the number of elements of layer #2 (1-40)
C
      DO I=1,10	      ! number of carbon fiber supports (see sketch)
C
C --- Place part # 1-2 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(1)-XBEG(1))*(XEND(1)-XBEG(1))+
     $   (YEND(1)-YBEG(1))*(YEND(1)-YBEG(1)))/20.
         DSUP(3)=25.
         XCC=(XX(1)+XX(2))/20.
         YCC=(YY(1)+YY(2))/20.
         XCCC=(XBEG(1)+XEND(1))/20.
         YCCC=(YBEG(1)+YEND(1))/20.
         IF(XX(1).EQ.XX(2)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(2)-YY(1),XX(2)-XX(1))*RADDEG-90.
	 ENDIF 		
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA12=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+1),90.,ATHETA12,90.,
     $   90.+ATHETA12,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+1,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+1),'ONLY',DSUP,3)
C
C --- Place part # 2-3 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(2)-XBEG(2))*(XEND(2)-XBEG(2))+
     $   (YEND(2)-YBEG(2))*(YEND(2)-YBEG(2)))/20.
         DSUP(3)=25.
         XCC=(XX(2)+XX(3))/20.
         YCC=(YY(2)+YY(3))/20.
         XCCC=(XBEG(2)+XEND(2))/20.
         YCCC=(YBEG(2)+YEND(2))/20.
         IF(XX(2).EQ.XX(3)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(3)-YY(2),XX(3)-XX(2))*RADDEG-90.
	 ENDIF 		
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA23=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+2),90.,ATHETA23,90.,
     $   90.+ATHETA23,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+2),'ONLY',DSUP,3)
C
C --- Place an element of layer #2
C
         BIGA=(YY(3)-YY(2))/(XX(3)-XX(2))
         BIGB=(XX(3)*YY(2)-XX(2)*YY(3))/(XX(3)-XX(2))/10.
         COEFFA=(BIGA*BIGA+1.)
         COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
         COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.04713*0.04713)
         XCC1=(-COEFFB+SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC1=BIGA*XCC1+BIGB
         BIGA1=-1./BIGA
         BIGB1=(XCC1/BIGA+YCC1)
         COEFFA=(BIGA1*BIGA1+1.)
         COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
         COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX2(1))*(DSUP(1)+DBOX2(1)))
         XCC2=(-COEFFB+SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC2=BIGA1*XCC2+BIGB1
         XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
         YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
	 JBOX2=JBOX2+1
	 CALL GSPOS('IPV2',JBOX2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+2),'ONLY')	
C
C --- Place part # 3-4 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(3)-XBEG(3))*(XEND(3)-XBEG(3))+
     $   (YEND(3)-YBEG(3))*(YEND(3)-YBEG(3)))/20.
         DSUP(3)=25.
         XCC=(XX(2)+XX(3))/20.
         YCC=(YY(2)+YY(3))/20.
         XCCC=(XBEG(3)+XEND(3))/20.
         YCCC=(YBEG(3)+YEND(3))/20.
         IF(XX(3).EQ.XX(4)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(4)-YY(3),XX(4)-XX(3))*RADDEG-90.
	 ENDIF 		
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA34=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+3),90.,ATHETA34,90.,
     $   90.+ATHETA34,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+3,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+3),'ONLY',DSUP,3)
C
C --- Place part # 4-5 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(4)-XBEG(4))*(XEND(4)-XBEG(4))+
     $   (YEND(4)-YBEG(4))*(YEND(4)-YBEG(4)))/20.
         DSUP(3)=25.
         XCC=(XX(4)+XX(5))/20.
         YCC=(YY(4)+YY(5))/20.
         XCCC=(XBEG(4)+XEND(4))/20.
         YCCC=(YBEG(4)+YEND(4))/20.
         IF(XX(4).EQ.XX(5)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(5)-YY(4),XX(5)-XX(4))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA45=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+4),90.,ATHETA45,90.,
     $   90.+ATHETA45,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+4,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+4),'ONLY',DSUP,3)
C
C --- Place an element of layer #2
C
	 BIGA=(YY(5)-YY(4))/(XX(5)-XX(4))
	 BIGB=(XX(5)*YY(4)-XX(4)*YY(5))/(XX(5)-XX(4))/10.
	 COEFFA=(BIGA*BIGA+1.)
	 COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
	 COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.*0.)
	 XCC1=XCCC
	 YCC1=YCCC
	 BIGA1=-1./BIGA
	 BIGB1=(XCC1/BIGA+YCC1)
	 COEFFA=(BIGA1*BIGA1+1.)
	 COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
	 COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX2(1))*(DSUP(1)+DBOX2(1)))
	 XCC2=(-COEFFB-SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
	 YCC2=BIGA1*XCC2+BIGB1
	 XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
	 YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
	 XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
	 YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
	 ZPOS=0.
	 JBOX2=JBOX2+1
	 CALL GSPOS('IPV2',JBOX2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+4),'ONLY')	
C
C --- Place part # 5-6 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(5)-XBEG(5))*(XEND(5)-XBEG(5))+
     $   (YEND(5)-YBEG(5))*(YEND(5)-YBEG(5)))/20.
         DSUP(3)=25.
         XCC=(XX(5)+XX(6))/20.
         YCC=(YY(5)+YY(6))/20.
         XCCC=(XBEG(5)+XEND(5))/20.
         YCCC=(YBEG(5)+YEND(5))/20.
         IF(XX(5).EQ.XX(6)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(6)-YY(5),XX(6)-XX(5))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA56=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+5),90.,ATHETA56,90.,
     $   90.+ATHETA56,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+5,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+5),'ONLY',DSUP,3)
C
C --- Place part # 6-7 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(6)-XBEG(6))*(XEND(6)-XBEG(6))+
     $   (YEND(6)-YBEG(6))*(YEND(6)-YBEG(6)))/20.
         DSUP(3)=25.
         XCC=(XX(6)+XX(7))/20.
         YCC=(YY(6)+YY(7))/20.
         XCCC=(XBEG(6)+XEND(6))/20.
         YCCC=(YBEG(6)+YEND(6))/20.
         IF(XX(6).EQ.XX(7)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(7)-YY(6),XX(7)-XX(6))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA67=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+6),90.,ATHETA67,90.,
     $   90.+ATHETA67,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+6,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+6),'ONLY',DSUP,3)
C
C --- Place an element of layer #2
C
         BIGA=(YY(7)-YY(6))/(XX(7)-XX(6))
         BIGB=(XX(7)*YY(6)-XX(6)*YY(7))/(XX(7)-XX(6))/10.
         COEFFA=(BIGA*BIGA+1.)
         COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
         COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.*0.)
         XCC1=XCCC
         YCC1=YCCC
         BIGA1=-1./BIGA
         BIGB1=(XCC1/BIGA+YCC1)
         COEFFA=(BIGA1*BIGA1+1.)
         COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
         COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX2(1))*(DSUP(1)+DBOX2(1)))
         XCC2=(-COEFFB-DSQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC2=BIGA1*XCC2+BIGB1
         XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
         YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
	 JBOX2=JBOX2+1
	 CALL GSPOS('IPV2',JBOX2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+6),'ONLY')	
C
C --- Place part # 7-8 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(7)-XBEG(7))*(XEND(7)-XBEG(7))+
     $        (YEND(7)-YBEG(7))*(YEND(7)-YBEG(7)))/20.
         DSUP(3)=25.
         XCC=(XX(7)+XX(8))/20.
         YCC=(YY(7)+YY(8))/20.
         XCCC=(XBEG(7)+XEND(7))/20.
         YCCC=(YBEG(7)+YEND(7))/20.
         IF(XX(7).EQ.XX(8)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(8)-YY(7),XX(8)-XX(7))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA78=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+7),90.,ATHETA78,90.,
     $   90.+ATHETA78,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+7,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+7),'ONLY',DSUP,3)
C
C --- Place part # 8-9 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(8)-XBEG(8))*(XEND(8)-XBEG(8))+
     $        (YEND(8)-YBEG(8))*(YEND(8)-YBEG(8)))/20.
         DSUP(3)=25.
         XCC=(XX(8)+XX(9))/20.
         YCC=(YY(8)+YY(9))/20.
         XCCC=(XBEG(8)+XEND(8))/20.
         YCCC=(YBEG(8)+YEND(8))/20.
         IF(XX(2).EQ.XX(3)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(9)-YY(8),XX(9)-XX(8))*RADDEG-90.
	 ENDIF	 	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA89=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+8),90.,ATHETA89,90.,
     $        90.+ATHETA89,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+8,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+8),'ONLY',DSUP,3)
C
C --- Place an element of layer #2
C
         BIGA=(YY(9)-YY(8))/(XX(9)-XX(8))
         BIGB=(XX(9)*YY(8)-XX(8)*YY(9))/(XX(9)-XX(8))/10.
         COEFFA=(BIGA*BIGA+1.)
         COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
         COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.23564*0.23564)
         XCC1=(-COEFFB-SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC1=BIGA*XCC1+BIGB
         BIGA1=-1./BIGA
         BIGB1=(XCC1/BIGA+YCC1)
         COEFFA=(BIGA1*BIGA1+1.)
         COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
         COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX2(1))*(DSUP(1)+DBOX2(1)))
         XCC2=(-COEFFB-SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC2=BIGA1*XCC2+BIGB1
         XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
         YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
	 JBOX2=JBOX2+1
	 CALL GSPOS('IPV2',JBOX2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+8),'ONLY')
C
C --- Place part # 9-10 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(9)-XBEG(9))*(XEND(9)-XBEG(9))+
     $   (YEND(9)-YBEG(9))*(YEND(9)-YBEG(9)))/20.
         DSUP(3)=25.
         XCC=(XX(9)+XX(10))/20.
         YCC=(YY(9)+YY(10))/20.
         XCCC=(XBEG(9)+XEND(9))/20.
         YCCC=(YBEG(9)+YEND(9))/20.
         IF(XX(9).EQ.XX(10)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(10)-YY(9),XX(10)-XX(9))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA910=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+9),90.,ATHETA910,90.,
     $   90.+ATHETA910,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+9,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+9),'ONLY',DSUP,3)
C
C --- Place part # 12-13 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(12)-XBEG(12))*(XEND(12)-XBEG(12))+
     $   (YEND(12)-YBEG(12))*(YEND(12)-YBEG(12)))/20.
         DSUP(3)=25.
         XCC=(XX(12)+XX(13))/20.
         YCC=(YY(12)+YY(13))/20.
         XCCC=(XBEG(12)+XEND(12))/20.
         YCCC=(YBEG(12)+YEND(12))/20.
         IF(XX(12).EQ.XX(13)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(13)-YY(12),XX(13)-XX(12))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA1213=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+12),270.,ATHETA1213,90.,
     $   270.+ATHETA1213,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+12,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+12),'ONLY',DSUP,3)
C
C --- Place an element of layer #1
C
	 BIGA=(YY(13)-YY(12))/(XX(13)-XX(12))
	 BIGB=(XX(13)*YY(12)-XX(12)*YY(13))/(XX(13)-XX(12))/10.
	 COEFFA=(BIGA*BIGA+1.)
	 COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
	 COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.12567*0.12567)
	 XCC1=(-COEFFB+SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
	 YCC1=BIGA*XCC1+BIGB
	 BIGA1=-1./BIGA
	 BIGB1=(XCC1/BIGA+YCC1)
	 COEFFA=(BIGA1*BIGA1+1.)
	 COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
	 COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX1(1))*(DSUP(1)+DBOX1(1)))
	 XCC2=(-COEFFB+SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
	 YCC2=BIGA1*XCC2+BIGB1
	 XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
	 YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
	 XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
	 YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
	 ZPOS=0.
	 JBOX1=JBOX1+1
	 CALL GSPOS('IPV1',JBOX1,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+12),'ONLY')
C
C --- Place part # 11-12 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(11)-XBEG(11))*(XEND(11)-XBEG(11))+
     $   (YEND(11)-YBEG(11))*(YEND(11)-YBEG(11)))/20.
         DSUP(3)=25.
         XCC=(XX(11)+XX(12))/20.
         YCC=(YY(11)+YY(12))/20.
         XCCC=(XBEG(11)+XEND(11))/20.
         YCCC=(YBEG(11)+YEND(11))/20.
         IF(XX(11).EQ.XX(12)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(12)-YY(11),XX(12)-XX(11))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA1112=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+11),90.,ATHETA1112,90.,
     $   90.+ATHETA1112,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+11,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+11),'ONLY',DSUP,3)
C
C --- Place part # 10-11 (see sketch)
C
 	 OFFSET1=-27.
         DSUP(1)=0.01
         DSUP(2)=SQRT((XEND(10)-XBEG(10))*(XEND(10)-XBEG(10))+
     $   (YEND(10)-YBEG(10))*(YEND(10)-YBEG(10)))/20.
         DSUP(3)=25.
         XCC=(XX(10)+XX(11))/20.
         YCC=(YY(10)+YY(11))/20.
         XCCC=(XBEG(10)+XEND(10))/20.
         YCCC=(YBEG(10)+YEND(10))/20.
         IF(XX(10).EQ.XX(11)) THEN
            OFFSET2=0.
         ELSE
	    OFFSET2=ATG(YY(11)-YY(10),XX(11)-XX(10))*RADDEG-90.
	 ENDIF	
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCCC*COS(APHI)-YCCC*SIN(APHI)+XZERO
         YPOS1=XCCC*SIN(APHI)+YCCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         ATHETA1011=FLOAT(I-1)*36.+OFFSET1+OFFSET2-GTETA
         CALL SXSROT(IDROTM(1100+(I-1)*13+10),270.,ATHETA1011,90.,
     $   270.+ATHETA1011,0.,0.)
         CALL GSPOSP('SPIX',(I-1)*13+10,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+10),'ONLY',DSUP,3)
C
C --- Place an element of layer #1
C
         BIGA=(YY(11)-YY(10))/(XX(11)-XX(10))
         BIGB=(XX(11)*YY(10)-XX(10)*YY(11))/(XX(11)-XX(10))/10.
         COEFFA=(BIGA*BIGA+1.)
         COEFFB=(BIGA*BIGB-BIGA*YCC-XCC)
         COEFFC=(XCC*XCC+YCC*YCC-2.*YCC*BIGB+BIGB*BIGB-0.001*0.001)
         XCC1=(-COEFFB-SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC1=BIGA*XCC1+BIGB
         BIGA1=-1./BIGA
         BIGB1=(XCC1/BIGA+YCC1)
         COEFFA=(BIGA1*BIGA1+1.)
         COEFFB=(BIGA1*BIGB1-BIGA1*YCC1-XCC1)
         COEFFC=(XCC1*XCC1+YCC1*YCC1-2.*YCC1*BIGB1+BIGB1*BIGB1-
     $   (DSUP(1)+DBOX1(1))*(DSUP(1)+DBOX1(1)))
         XCC2=(-COEFFB+SQRT(COEFFB*COEFFB-COEFFA*COEFFC))/COEFFA
         YCC2=BIGA1*XCC2+BIGB1
         XPOS1=XCC2*COS(APHI)-YCC2*SIN(APHI)+XZERO
         YPOS1=XCC2*SIN(APHI)+YCC2*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
	 JBOX1=JBOX1+1
	 CALL GSPOS('IPV1',JBOX1,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+10),'ONLY')
C
C --- Place arc # 5 (between part 1-2 and part 2-3) (see sketch)
C
	 DARC(1)=RARC(5)/10.-0.02
         DARC(2)=RARC(5)/10.
         DARC(3)=25.
         DARC(4)=ATHETA12+90.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA23+90.-FLOAT(I-1)*(360./10.)
         XCC=XARC(5)/10.
         YCC=YARC(5)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+5,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+5),'ONLY',DARC,5)
C
C --- Place arc # 6 (between part 2-3 and part 3-4) (see sketch)
C
	 DARC(1)=RARC(6)/10.-0.02
         DARC(2)=RARC(6)/10.
         DARC(3)=25.
         DARC(4)=ATHETA23-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA34-FLOAT(I-1)*(360./10.)
         XCC=XARC(6)/10.
         YCC=YARC(6)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+6,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+6),'ONLY',DARC,5)
C
C --- Place arc # 7 (between part 3-4 and part 4-5) (see sketch)
C
	 DARC(1)=RARC(7)/10.-0.02
         DARC(2)=RARC(7)/10.
         DARC(3)=25.
         DARC(4)=ATHETA45+90.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA34+90.-FLOAT(I-1)*(360./10.)
         XCC=XARC(7)/10.
         YCC=YARC(7)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+7,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+7),'ONLY',DARC,5)
C
C --- Place arc # 8 (between part 4-5 and part 5-6) (see sketch)
C
	 DARC(1)=RARC(8)/10.-0.02
         DARC(2)=RARC(8)/10.
         DARC(3)=25.
         DARC(4)=ATHETA56+180.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA45+180.-FLOAT(I-1)*(360./10.)
         XCC=XARC(8)/10.
         YCC=YARC(8)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+8,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+8),'ONLY',DARC,5)
C
C --- Place arc # 9 (between part 5-6 and part 6-7) (see sketch)
C
	 DARC(1)=RARC(9)/10.-0.02
         DARC(2)=RARC(9)/10.
         DARC(3)=25.
         DARC(4)=ATHETA56-90.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA67-90.-FLOAT(I-1)*(360./10.)
         XCC=XARC(9)/10.
         YCC=YARC(9)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+9,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+9),'ONLY',DARC,5)
C
C --- Place arc # 10 (between part 6-7 and part 7-8) (see sketch)
C
	 DARC(1)=RARC(10)/10.-0.02
         DARC(2)=RARC(10)/10.
         DARC(3)=25.
         DARC(4)=ATHETA67-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA78-FLOAT(I-1)*(360./10.)
         XCC=XARC(10)/10.
         YCC=YARC(10)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+10,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+10),'ONLY',DARC,5)
C
C --- Place arc # 11 (between part 7-8 and part 8-9) (see sketch)
C
	 DARC(1)=RARC(11)/10.-0.02
         DARC(2)=RARC(11)/10.
         DARC(3)=25.
         DARC(4)=ATHETA89-90.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA78-90.-FLOAT(I-1)*(360./10.)
         XCC=XARC(11)/10.
         YCC=YARC(11)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+11,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+11),'ONLY',DARC,5)	
C
C --- Place arc # 12 (between part 8-9 and part 9-10) (see sketch)
C
	 DARC(1)=RARC(12)/10.-0.02
         DARC(2)=RARC(12)/10.
         DARC(3)=25.
         DARC(4)=ATHETA89-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA910-FLOAT(I-1)*(360./10.)
         XCC=XARC(12)/10.
         YCC=YARC(12)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+12,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+12),'ONLY',DARC,5)
C
C --- Place arc # 1 (between part 9-10 and part 10-11) (see sketch)
C
	 DARC(1)=RARC(1)/10.-0.02
         DARC(2)=RARC(1)/10.
         DARC(3)=25.
         DARC(4)=ATHETA1011+10.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA910-120.-FLOAT(I-1)*(360./10.)
         XCC=XARC(1)/10.
         YCC=YARC(1)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+1,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+1),'ONLY',DARC,5)
C
C --- Place arc # 2 (between part 10-11 and part 11-12) (see sketch)
C
	 DARC(1)=RARC(2)/10.-0.02
         DARC(2)=RARC(2)/10.
         DARC(3)=25.
         DARC(4)=ATHETA1011-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA1112+45.-FLOAT(I-1)*(360./10.)
         XCC=XARC(2)/10.
         YCC=YARC(2)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+2,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+2),'ONLY',DARC,5)
C
C --- Place arc # 3 (between part 11-12 and part 12-13) (see sketch)
C
	 DARC(1)=RARC(3)/10.-0.02
         DARC(2)=RARC(3)/10.
         DARC(3)=25.
         DARC(4)=ATHETA1213+115.-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA1112+115.-FLOAT(I-1)*(360./10.)
         XCC=XARC(3)/10.
         YCC=YARC(3)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+3,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+3),'ONLY',DARC,5)
C
C --- Place arc # 4 (between part 12-13 and part 1-2) (see sketch)
C
	 DARC(1)=RARC(4)/10.-0.02
         DARC(2)=RARC(4)/10.
         DARC(3)=25.
         DARC(4)=ATHETA1213-FLOAT(I-1)*(360./10.)
         DARC(5)=ATHETA12-FLOAT(I-1)*(360./10.)
         XCC=XARC(4)/10.
         YCC=YARC(4)/10.
         APHI=(PPHI+FLOAT(I-1)*36.)*DEGRAD
	 XZERO=RR*COS((TTETA+FLOAT(I-1)*36.)*DEGRAD)
	 YZERO=RR*SIN((TTETA+FLOAT(I-1)*36.)*DEGRAD)
         XPOS1=XCC*COS(APHI)-YCC*SIN(APHI)+XZERO
         YPOS1=XCC*SIN(APHI)+YCC*COS(APHI)+YZERO
         XPOS=XPOS1*COS(GTETA*DEGRAD)+YPOS1*SIN(GTETA*DEGRAD)
         YPOS=-XPOS1*SIN(GTETA*DEGRAD)+YPOS1*COS(GTETA*DEGRAD)
         ZPOS=0.
         CALL GSPOSP('SARC',(I-1)*13+4,'IT12',XPOS,YPOS,ZPOS,IDROTM(
     $   1100+(I-1)*13+4),'ONLY',DARC,5)
C       	
      END DO
2345  CONTINUE
C	
C ************************************************************************
C *                                                                      *
C *                               D R I F T S                            *
C *                               ===========                            *
C *                                                                      *
C ************************************************************************
C
C --- Define a ghost volume containing the Silicon Drift Detectors
C     (layer #3 and #4) and fill it with air or vacuum
C
        XXM=(49.999-3.)/(70.-25.)
        DGH(1)=0.
        DGH(2)=360.
        DGH(3)=4.
        DGH(4)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM-(27.-9.01)/XXM
        DGH(5)=27.
        DGH(6)=27.
        DGH(7)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM
        DGH(8)=9.01
        DGH(9)=27.
        DGH(10)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM
        DGH(11)=9.01
        DGH(12)=27.
        DGH(13)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM+(27.-9.01)/XXM
        DGH(14)=27.
        DGH(15)=27.
        CALL GSVOLU('IT34','PCON',IDTMED(276),DGH,15,IOUT)
C
C --- Place the ghost volume in its mother volume (ITSV) and make it
C     invisible
C
        CALL GSPOS('IT34',1,'ITSV',0.,0.,0.,0,'ONLY')
        CALL GSATT('IT34','SEEN',0)
C
C --- Layer #3
C
C        GOTO 3456           ! skip ITS layer no. 3
C
C --- Define a ghost volume containing a single ladder of layer #3 (with the
C     smaller lenght of ribs) and fill it with air or vacuum
C
      DBOX1(1)=0.5+(0.0172+0.03+0.0252+0.04+0.003)
      DBOX1(2)=3.85      ! the widest element is the sensitive element
      DBOX1(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5	  ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IDV1','BOX ',IDTMED(229),DBOX1,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('IDV1','SEEN',0)
C
C --- Define a volume containing the sensitive part of drifts
C     (silicon, layer #3)
C
      DITS(1)=0.0172	  ! see material budget report by G. Feofilov
      DITS(2)=3.85
      DITS(3)=4.35
      CALL GSVOLU('ITS3','BOX ',IDTMED(225),DITS,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts made of
C     carbon (layer #3)
C
      DSUP(1)=0.5-DITS(1)
      DSUP(2)=0.01
      DSUP(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR11','BOX ',IDTMED(228),DSUP,3,IOUT)
C
C --- Define the first part of the (smaller) rib between two sensitive parts
C     made of aluminum (layer #3)
C
      DAL1(1)=0.5-DITS(1)
      DAL1(2)=0.00096/2.
      DAL1(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR12','BOX ',IDTMED(231),DAL1,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts made of
C     kapton (layer #3)
C
      DKAP(1)=0.5-DITS(1)
      DKAP(2)=0.0317/2.
      DKAP(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR13','BOX ',IDTMED(237),DKAP,3,IOUT)
C
C --- Define the second part of the (smaller) rib between two sensitive parts
C     made of aluminum (layer #3)
C
      DAL2(1)=0.5-DITS(1)
      DAL2(2)=0.0027/2.
      DAL2(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR14','BOX ',IDTMED(231),DAL2,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts
C     made of silicon (the electronics) (layer #3)
C
      DCHI(1)=0.5-DITS(1)
      DCHI(2)=0.0071/2.
      DCHI(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR15','BOX ',IDTMED(226),DAL2,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts
C     made of water (the cooler) (layer #3)
C
      DWAT(1)=0.5-DITS(1)
      DWAT(2)=0.0093/2.
      DWAT(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5    ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR16','BOX ',IDTMED(232),DWAT,3,IOUT)
C
C --- Define the third part of the (smaller) rib between two sensitive parts
C     made of aluminum (the cooling tubes) (layer #3)
C
      DTUB(1)=0.5-DITS(1)
      DTUB(2)=0.00134/2.
      DTUB(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR17','BOX ',IDTMED(231),DTUB,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of PCB (layer #3)
C
      DPCB(1)=0.03            ! twice the foreseen thickness
      DPCB(2)=3.5
      DPCB(3)=7.5
      CALL GSVOLU('IEL1','BOX ',IDTMED(234),DPCB,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #3)
C
      DCOP(1)=0.0252            ! twice the foreseen thickness
      DCOP(2)=3.5
      DCOP(3)=7.5
      CALL GSVOLU('IEL2','BOX ',IDTMED(235),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of ceramics (layer #3)
C
      DCER(1)=0.04            ! twice the foreseen thickness
      DCER(2)=3.5
      DCER(3)=7.5
      CALL GSVOLU('IEL3','BOX ',IDTMED(236),DCER,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (layer #3)
C
      DSIL(1)=0.003            ! twice the foreseen thickness
      DSIL(2)=3.5
      DSIL(3)=7.5
      CALL GSVOLU('IEL4','BOX ',IDTMED(227),DSIL,3,IOUT)
C
C --- Place the sensitive part of the drifts (smaller ribs) into its mother
C     (IDV1)
C
      YPOS=0.
      DO J=1,5    ! odd elements are up and even elements are down
	 IF(J.EQ.1) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)
	    ZPOS=0.-DITS(3)+1.-2.*DITS(3)-0.1-DITS(3)
	 ELSE IF(J.EQ.2) THEN
	    XPOS=-DBOX1(1)+DITS(1)	
	    ZPOS=0.-DITS(3)+1.-DITS(3)
	 ELSE IF(J.EQ.3) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.
	 ELSE IF(J.EQ.4) THEN
	    XPOS=-DBOX1(1)+DITS(1)	
	    ZPOS=0.+DITS(3)-1.+DITS(3)
	 ELSE IF(J.EQ.5) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-1.+2.*DITS(3)+0.1+DITS(3)
	 ENDIF
	 CALL GSPOS('ITS3',J,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')	  	
      END DO
C
C --- Place the smaller ribs into their mother (IDV1)
C
C --- Right ribs (just a matter of convention)
C
      XPOS=0.5-DBOX1(1)+DITS(1)
      ZPOS=0.
C
C --- Carbon
C
      YPOS=2.81
      CALL GSPOS('IR11',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=2.81+DSUP(2)+DAL1(2)
      CALL GSPOS('IR12',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+DKAP(2)
      CALL GSPOS('IR13',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2)
      CALL GSPOS('IR14',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2)
      CALL GSPOS('IR15',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2)
      CALL GSPOS('IR16',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2)
      CALL GSPOS('IR17',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right ribs (just a matter of convention)
C
C --- Carbon
C
      YPOS=-2.81
      CALL GSPOS('IR11',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=-(2.81+DSUP(2)+DAL1(2))
      CALL GSPOS('IR12',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+DKAP(2))
      CALL GSPOS('IR13',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2))
      CALL GSPOS('IR14',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2))
      CALL GSPOS('IR15',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2))
      CALL GSPOS('IR16',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2))
      CALL GSPOS('IR17',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the end-ladder stuff into its mother (IDV1)
C
C
C --- Negative-Z end-ladder
C
      YPOS=0.
      ZPOS=-(8.7*5.-2.*1.+2.*0.1)/2.-7.5
C
C --- PCB
C
      XPOS=DBOX1(1)-DPCB(1)
      CALL GSPOS('IEL1',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX1(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL2',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL3',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL4',1,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Positive-Z end-ladder
C
      YPOS=0.
      ZPOS=(8.7*5.-2.*1.+2.*0.1)/2.+7.5
C
C --- PCB
C
      XPOS=DBOX1(1)-DPCB(1)
      CALL GSPOS('IEL1',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX1(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL2',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL3',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL4',2,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Define a ghost volume containing a single ladder of layer #3 (with the
C     larger lenght of ribs) and fill it with air or vacuum
C
      DBOX2(1)=0.65+(0.0172+0.03+0.0252+0.04+0.003)
      DBOX2(2)=3.85     ! the widest element is the sensitive element
      DBOX2(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5	  ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IDV2','BOX ',IDTMED(229),DBOX2,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('IDV2','SEEN',0)
C
C --- Define the part of the (larger) rib between two sensitive parts made of
C     carbon (layer #3)
C
      DSUP(1)=0.65-DITS(1)
      DSUP(2)=0.01
      DSUP(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR21','BOX ',IDTMED(228),DSUP,3,IOUT)
C
C --- Define the first part of the (larger) rib between two sensitive parts
C     made of aluminum (layer #3)
C
      DAL1(1)=0.65-DITS(1)
      DAL1(2)=0.00096/2.
      DAL1(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR22','BOX ',IDTMED(231),DAL1,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts made of
C     kapton (layer #3)
C
      DKAP(1)=0.65-DITS(1)
      DKAP(2)=0.0317/2.
      DKAP(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR23','BOX ',IDTMED(237),DKAP,3,IOUT)
C
C --- Define the second part of the (larger) rib between two sensitive parts
C     made of aluminum (layer #3)
C
      DAL2(1)=0.65-DITS(1)
      DAL2(2)=0.0027/2.
      DAL2(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR24','BOX ',IDTMED(231),DAL2,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts
C     made of silicon (the electronics) (layer #3)
C
      DCHI(1)=0.65-DITS(1)
      DCHI(2)=0.0071/2.
      DCHI(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR25','BOX ',IDTMED(226),DAL2,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts
C     made of water (the cooler) (layer #3)
C
      DWAT(1)=0.65-DITS(1)
      DWAT(2)=0.0093/2.
      DWAT(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR26','BOX ',IDTMED(232),DWAT,3,IOUT)
C
C --- Define the third part of the (larger) rib between two sensitive parts
C     made of aluminum (the cooling tubes) (layer #3)
C
      DTUB(1)=0.65-DITS(1)
      DTUB(2)=0.00134/2.
      DTUB(3)=(8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR27','BOX ',IDTMED(231),DTUB,3,IOUT)
C
C --- Place the sensitive part of the drifts (smaller ribs) into its mother
C     (IDV2)
C
      YPOS=0.
      DO J=1,5      ! odd element are up and even elements are down
	 IF(J.EQ.1) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)
	    ZPOS=0.-DITS(3)+1.-2.*DITS(3)-0.1-DITS(3)
	 ELSE IF(J.EQ.2) THEN
	    XPOS=-DBOX2(1)+DITS(1)	
	    ZPOS=0.-DITS(3)+1.-DITS(3)
	 ELSE IF(J.EQ.3) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.
	 ELSE IF(J.EQ.4) THEN
	    XPOS=-DBOX2(1)+DITS(1)	
	    ZPOS=0.+DITS(3)-1.+DITS(3)
	 ELSE IF(J.EQ.5) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-1.+2.*DITS(3)+0.1+DITS(3)
	 ENDIF
	 CALL GSPOS('ITS3',J,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')	  	
      END DO
C
C --- Place the larger ribs into their mother (IDV2)
C
C
C --- Right ribs (just a matter of convention)
C
      XPOS=0.65-DBOX2(1)+DITS(1)
      ZPOS=0.
C
C --- Carbon
C
      YPOS=2.81
      CALL GSPOS('IR21',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=2.81+DSUP(2)+DAL1(2)
      CALL GSPOS('IR22',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+DKAP(2)
      CALL GSPOS('IR23',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2)
      CALL GSPOS('IR24',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2)
      CALL GSPOS('IR25',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2)
      CALL GSPOS('IR26',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2)
      CALL GSPOS('IR27',1,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right ribs (just a matter of convention)
C
C --- Carbon
C
      YPOS=-2.81
      CALL GSPOS('IR21',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=-(2.81+DSUP(2)+DAL1(2))
      CALL GSPOS('IR22',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+DKAP(2))
      CALL GSPOS('IR23',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2))
      CALL GSPOS('IR24',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2))
      CALL GSPOS('IR25',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2))
      CALL GSPOS('IR26',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2))
      CALL GSPOS('IR27',2,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the end-ladder stuff into its mother (IDV1)
C
C
C --- Negative-Z end-ladder
C
      YPOS=0.
      ZPOS=-(8.7*5.-2.*1.+2.*0.1)/2.-7.5
C
C --- PCB
C
      XPOS=DBOX2(1)-DPCB(1)
      CALL GSPOS('IEL1',3,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX2(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL2',3,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL3',3,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL4',3,'IDV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Positive-Z end-ladder
C
      YOS=0.
      ZPOS=(8.7*5.-2.*1.+2.*0.1)/2.+7.5
C
C --- PCB
C
      XPOS=DBOX2(1)-DPCB(1)
      CALL GSPOS('IEL1',4,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX2(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL2',4,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL3',4,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL4',4,'IDV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the drift ladders of layer #3 in their
C     mother volume (IT34)
C     Odd elements have large ribs and even elements have small ribs
C
      DO I=1,12
        ATHETA=FLOAT(I-1)*(360./12.)
        CALL SXSROT(IDROTM(1300+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
        IF(MOD(I,2).EQ.0) THEN
           RZERO=(14.5+13.5)/2.
	   XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/12.)
	   YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/12.)
	   ZPOS=0.
	   CALL GSPOS('IDV1',I,'IT34',XPOS,YPOS,ZPOS,IDROTM(1300+I),
     $     'ONLY')
        ELSE
           RZERO=(14.5+13.2)/2.	
	   XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/12.)
	   YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/12.)
	   ZPOS=0.
	   CALL GSPOS('IDV2',I,'IT34',XPOS,YPOS,ZPOS,IDROTM(1300+I),
     $     'ONLY')		   	
	ENDIF
      ENDDO
C
3456  CONTINUE
C
C --- Layer #4
C
C        GOTO 4567           ! skip ITS layer no. 4
C
C --- Define a ghost volume containing a single ladder of layer #4 (with the
C     smaller lenght of ribs) and fill it with air or vacuum
C
      DBOX1(1)=0.5+(0.0172+0.03+0.0252+0.04+0.003)
      DBOX1(2)=3.5       ! the widest element is the end-ladder stuff
      DBOX1(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5	  ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IDV3','BOX ',IDTMED(229),DBOX1,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('IDV3','SEEN',0)
C
C --- Define a volume containing the sensitive part of drifts
C     (silicon, layer #4)
C
      DITS(1)=0.0172	  ! see material budget report by G. Feofilov
      DITS(2)=3.125
      DITS(3)=4.35
      CALL GSVOLU('ITS4','BOX ',IDTMED(225),DITS,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts made of
C     carbon (layer #4)
C
      DSUP(1)=0.5-DITS(1)
      DSUP(2)=0.01
      DSUP(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR31','BOX ',IDTMED(228),DSUP,3,IOUT)
C
C --- Define the first part of the (smaller) rib between two sensitive parts
C     made of aluminum (layer #4)
C
      DAL1(1)=0.5-DITS(1)
      DAL1(2)=0.00096/2.
      DAL1(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR32','BOX ',IDTMED(231),DAL1,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts made of
C     kapton (layer #4)
C
      DKAP(1)=0.5-DITS(1)
      DKAP(2)=0.0317/2.
      DKAP(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR33','BOX ',IDTMED(237),DKAP,3,IOUT)
C
C --- Define the second part of the (smaller) rib between two sensitive parts
C     made of aluminum (layer #4)
C
      DAL2(1)=0.5-DITS(1)
      DAL2(2)=0.0027/2.
      DAL2(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR34','BOX ',IDTMED(231),DAL2,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts
C     made of silicon (the electronics) (layer #4)
C
      DCHI(1)=0.5-DITS(1)
      DCHI(2)=0.0071/2.
      DCHI(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR35','BOX ',IDTMED(226),DAL2,3,IOUT)
C
C --- Define the part of the (smaller) rib between two sensitive parts
C     made of water (the cooler) (layer #4)
C
      DWAT(1)=0.5-DITS(1)
      DWAT(2)=0.0093/2.
      DWAT(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5    ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR36','BOX ',IDTMED(232),DWAT,3,IOUT)
C
C --- Define the third part of the (smaller) rib between two sensitive parts
C     made of aluminum (the cooling tubes) (layer #4)
C
      DTUB(1)=0.5-DITS(1)
      DTUB(2)=0.00134/2.
      DTUB(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR37','BOX ',IDTMED(231),DTUB,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of PCB (layer #4)
C
      DPCB(1)=0.03            ! twice the foreseen thickness
      DPCB(2)=3.5
      DPCB(3)=7.5
      CALL GSVOLU('IEL5','BOX ',IDTMED(234),DPCB,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #4)
C
      DCOP(1)=0.0252            ! twice the foreseen thickness
      DCOP(2)=3.5
      DCOP(3)=7.5
      CALL GSVOLU('IEL6','BOX ',IDTMED(235),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of ceramics (layer #4)
C
      DCER(1)=0.04            ! twice the foreseen thickness
      DCER(2)=3.5
      DCER(3)=7.5
      CALL GSVOLU('IEL7','BOX ',IDTMED(236),DCER,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (layer #4)
C
      DSIL(1)=0.003            ! twice the foreseen thickness
      DSIL(2)=3.5
      DSIL(3)=7.5
      CALL GSVOLU('IEL8','BOX ',IDTMED(227),DSIL,3,IOUT)
C
C --- Place the sensitive part of the drifts (smaller ribs) into its mother
C     (IDV3)
C
      YPOS=0.
      DO J=1,7    ! odd elements are down and even elements are up
	 IF(J.EQ.1) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)
	    ZPOS=0.-DITS(3)+0.7-2.*DITS(3)-0.-2.*DITS(3)+1.3-DITS(3)
	 ELSE IF(J.EQ.2) THEN
	    XPOS=-DBOX1(1)+DITS(1)	
	    ZPOS=0.-DITS(3)+0.7-2.*DITS(3)-0.-DITS(3)
	 ELSE IF(J.EQ.3) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.-DITS(3)+0.7-DITS(3)
	 ELSE IF(J.EQ.4) THEN
	    XPOS=-DBOX1(1)+DITS(1)	
	    ZPOS=0.
	 ELSE IF(J.EQ.5) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-0.7+DITS(3)	    	
	 ELSE IF(J.EQ.6) THEN
	    XPOS=-DBOX1(1)+DITS(1)	    	
	    ZPOS=0.+DITS(3)-0.7+2.*DITS(3)+0.+DITS(3)
	 ELSE IF(J.EQ.7) THEN
	    XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-0.7+2.*DITS(3)+0.+2.*DITS(3)-1.3+DITS(3)
	 ENDIF
	 CALL GSPOS('ITS4',J,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')	  	
      END DO
C
C --- Place the smaller ribs into their mother (IDV3)
C
C --- Right ribs (just a matter of convention)
C
      XPOS=0.5-DBOX1(1)+DITS(1)
      ZPOS=0.
C
C --- Carbon
C
      YPOS=2.81
      CALL GSPOS('IR31',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=2.81+DSUP(2)+DAL1(2)
      CALL GSPOS('IR32',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+DKAP(2)
      CALL GSPOS('IR33',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2)
      CALL GSPOS('IR34',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2)
      CALL GSPOS('IR35',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2)
      CALL GSPOS('IR36',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2)
      CALL GSPOS('IR37',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right ribs (just a matter of convention)
C
C --- Carbon
C
      YPOS=-2.81
      CALL GSPOS('IR31',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=-(2.81+DSUP(2)+DAL1(2))
      CALL GSPOS('IR32',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+DKAP(2))
      CALL GSPOS('IR33',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2))
      CALL GSPOS('IR34',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2))
      CALL GSPOS('IR35',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2))
      CALL GSPOS('IR36',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2))
      CALL GSPOS('IR37',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the end-ladder stuff into its mother (IDV1)
C
C
C --- Negative-Z end-ladder
C
      YPOS=0.
      ZPOS=-(8.7*7.-2.*0.7-2.*1.3)/2.-7.5
C
C --- PCB
C
      XPOS=DBOX1(1)-DPCB(1)
      CALL GSPOS('IEL5',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX1(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL6',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL7',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL8',1,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Positive-Z end-ladder
C
      YPOS=0.
      ZPOS=(8.7*7.-2.*0.7-2.*1.3)/2.+7.5
C
C --- PCB
C
      XPOS=DBOX1(1)-DPCB(1)
      CALL GSPOS('IEL5',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX1(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL6',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL7',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX1(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL8',2,'IDV3',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Define a ghost volume containing a single ladder of layer #4 (with the
C     larger lenght of ribs) and fill it with air or vacuum
C
      DBOX2(1)=0.65+(0.0172+0.03+0.0252+0.04+0.003)
      DBOX2(2)=3.5       ! the widest element is the end-ladder stuff
      DBOX2(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5	  ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IDV4','BOX ',IDTMED(229),DBOX2,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('IDV4','SEEN',0)
C
C --- Define the part of the (larger) rib between two sensitive parts made of
C     carbon (layer #4)
C
      DSUP(1)=0.65-DITS(1)
      DSUP(2)=0.01
      DSUP(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR41','BOX ',IDTMED(228),DSUP,3,IOUT)
C
C --- Define the first part of the (larger) rib between two sensitive parts
C     made of aluminum (layer #4)
C
      DAL1(1)=0.65-DITS(1)
      DAL1(2)=0.00096/2.
      DAL1(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR42','BOX ',IDTMED(231),DAL1,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts made of
C     kapton (layer #4)
C
      DKAP(1)=0.65-DITS(1)
      DKAP(2)=0.0317/2.
      DKAP(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR43','BOX ',IDTMED(237),DKAP,3,IOUT)
C
C --- Define the second part of the (larger) rib between two sensitive parts
C     made of aluminum (layer #4)
C
      DAL2(1)=0.65-DITS(1)
      DAL2(2)=0.0027/2.
      DAL2(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR44','BOX ',IDTMED(231),DAL2,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts
C     made of silicon (the electronics) (layer #4)
C
      DCHI(1)=0.65-DITS(1)
      DCHI(2)=0.0071/2.
      DCHI(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR45','BOX ',IDTMED(226),DAL2,3,IOUT)
C
C --- Define the part of the (larger) rib between two sensitive parts
C     made of water (the cooler) (layer #4)
C
      DWAT(1)=0.65-DITS(1)
      DWAT(2)=0.0093/2.
      DWAT(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR46','BOX ',IDTMED(232),DWAT,3,IOUT)
C
C --- Define the third part of the (larger) rib between two sensitive parts
C     made of aluminum (the cooling tubes) (layer #4)
C
      DTUB(1)=0.65-DITS(1)
      DTUB(2)=0.00134/2.
      DTUB(3)=(8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5     ! 7.5 cm is the lenght of the end-ladder module
      CALL GSVOLU('IR47','BOX ',IDTMED(231),DTUB,3,IOUT)
C
C --- Place the sensitive part of the drifts (smaller ribs) into its mother
C     (IDV4)
C
      YPOS=0.
      DO J=1,7    ! odd elements are down and even elements are up
	 IF(J.EQ.1) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)
	    ZPOS=0.-DITS(3)+0.7-2.*DITS(3)-0.-2.*DITS(3)+1.3-DITS(3)
	 ELSE IF(J.EQ.2) THEN
	    XPOS=-DBOX2(1)+DITS(1)	
	    ZPOS=0.-DITS(3)+0.7-2.*DITS(3)-0.-DITS(3)
	 ELSE IF(J.EQ.3) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.-DITS(3)+0.7-DITS(3)
	 ELSE IF(J.EQ.4) THEN
	    XPOS=-DBOX2(1)+DITS(1)	
	    ZPOS=0.
	 ELSE IF(J.EQ.5) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-0.7+DITS(3)	    	
	 ELSE IF(J.EQ.6) THEN
	    XPOS=-DBOX2(1)+DITS(1)	
	    ZPOS=0.+DITS(3)-0.7+2.*DITS(3)+0.+DITS(3)
	 ELSE IF(J.EQ.7) THEN
	    XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-2.*DSIL(1)-
     $      DITS(1)	
	    ZPOS=0.+DITS(3)-0.7+2.*DITS(3)+0.+2.*DITS(3)-1.3+DITS(3)
	 ENDIF
	 CALL GSPOS('ITS4',J,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')	  	
      END DO
C
C --- Place the larger ribs into their mother (IDV4)
C
C
C --- Right ribs (just a matter of convention)
C
      XPOS=0.65-DBOX2(1)+DITS(1)
      ZPOS=0.
C
C --- Carbon
C
      YPOS=2.81
      CALL GSPOS('IR41',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=2.81+DSUP(2)+DAL1(2)
      CALL GSPOS('IR42',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+DKAP(2)
      CALL GSPOS('IR43',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2)
      CALL GSPOS('IR44',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2)
      CALL GSPOS('IR45',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2)
      CALL GSPOS('IR46',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2)
      CALL GSPOS('IR47',1,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right ribs (just a matter of convention)
C
C --- Carbon
C
      YPOS=-2.81
      CALL GSPOS('IR41',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #1
C
      YPOS=-(2.81+DSUP(2)+DAL1(2))
      CALL GSPOS('IR42',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Kapton
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+DKAP(2))
      CALL GSPOS('IR43',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #2
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+DAL2(2))
      CALL GSPOS('IR44',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (chip)
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+DCHI(2))
      CALL GSPOS('IR45',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Water
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $DWAT(2))
      CALL GSPOS('IR46',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Aluminum #3
C
      YPOS=-(2.81+DSUP(2)+2.*DAL1(2)+2.*DKAP(2)+2.*DAL2(2)+2.*DCHI(2)+
     $2.*DWAT(2)+DTUB(2))
      CALL GSPOS('IR47',2,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the end-ladder stuff into its mother (IDV1)
C
C
C --- Negative-Z end-ladder
C
      YPOS=0.
      ZPOS=-(8.7*7.-2.*0.7-2.*1.3)/2.-7.5
C
C --- PCB
C
      XPOS=DBOX2(1)-DPCB(1)
      CALL GSPOS('IEL5',3,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX2(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL6',3,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL7',3,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL8',3,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Positive-Z end-ladder
C
      YOS=0.
      ZPOS=(8.7*7.-2.*0.7-2.*1.3)/2.+7.5
C
C --- PCB
C
      XPOS=DBOX2(1)-DPCB(1)
      CALL GSPOS('IEL5',4,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
      XPOS=DBOX2(1)-2.*DPCB(1)-DCOP(1)
      CALL GSPOS('IEL6',4,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Ceramics
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-DCER(1)
      CALL GSPOS('IEL7',4,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
      XPOS=DBOX2(1)-2.*DPCB(1)-2.*DCOP(1)-2.*DCER(1)-DSIL(1)
      CALL GSPOS('IEL8',4,'IDV4',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the drift ladders of layer #4 in their
C     mother volume (IT34)
C     Odd elements have large ribs and even elements have small ribs
C
      DO I=1,24
        ATHETA=FLOAT(I-1)*(360./24.)
        CALL SXSROT(IDROTM(1400+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
        IF(MOD(I,2).EQ.0) THEN
           RZERO=(24.0+23.0)/2.
	   XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/24.)
	   YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/24.)
	   ZPOS=0.
	   CALL GSPOS('IDV3',I,'IT34',XPOS,YPOS,ZPOS,IDROTM(1400+I),
     $     'ONLY')
        ELSE
           RZERO=(24.0+22.8)/2.	
	   XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/24.)
	   YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/24.)
	   ZPOS=0.
	   CALL GSPOS('IDV4',I,'IT34',XPOS,YPOS,ZPOS,IDROTM(1400+I),
     $     'ONLY')		   	
	ENDIF
      ENDDO
C
4567  CONTINUE
C	
C ************************************************************************
C *                                                                      *
C *                               S T R I P S                            *
C *                               ===========                            *
C *                                                                      *
C ************************************************************************
C
C
C --- Define SSD with the 35+39 lay-out
C
      IF(IVERS.EQ.40 .OR. IVERS.EQ.41 .OR. IVERS.EQ.42) THEN 	
C
C --- Define ghost volume containing the Strip Detectors and fill it with air
C     or vacuum
C
        XXM=(49.999-3.)/(70.-25.)
        DGH(1)=0.
        DGH(2)=360.
        DGH(3)=4.
        DGH(4)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM-(27.-9.01)/XXM-
     $  (37.-27)/XXM-(49.998-37.)/XXM
        DGH(5)=49.998
        DGH(6)=49.998
        DGH(7)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM-(27.-9.01)/XXM-
     $  (37.-27)/XXM
        DGH(8)=37.
        DGH(9)=49.998
        DGH(10)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM+(27.-9.01)/XXM+
     $  (37.-27)/XXM
        DGH(11)=37.
        DGH(12)=49.998
        DGH(13)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM+(27.-9.01)/XXM+
     $  (37.-27)/XXM+(49.998-37.)/XXM
        DGH(14)=49.998
        DGH(15)=49.998
        CALL GSVOLU('IT56','PCON',IDTMED(276),DGH,15,IOUT)
        CALL GSPOS('IT56',1,'ITSV',0.,0.,0.,0,'ONLY')
        CALL GSATT('IT56','SEEN',0)
C
C --- Layer #5
C
C        GOTO 5678           ! skip ITS layer no. 5
C
C --- Define a ghost volume containing a single ladder of layer #5 and fill
C     it with air or vacuum
C
        DBOX1(1)=(0.0600+2.*0.0150)/2.
        DBOX1(2)=3.75
        DBOX1(3)=90.22/2.
        CALL GSVOLU('ISV1','BOX ',IDTMED(254),DBOX1,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ISV1','SEEN',0)
C
C --- Define a ghost volume containing the electronics and cooling of
C     a single ladder of layer #5 and fill it with air or vacuum
C
        DSRV(1)=(0.47+(SQRT(3.)/2.)*4.2+0.05)/2.
        DSRV(2)=3.75
        DSRV(3)=90.22/2.
        CALL GSVOLU('SSV1','BOX ',IDTMED(254),DSRV,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('SSV1','SEEN',0)
C
C --- Define a ghost volume containing the end-ladder stuff of
C     a single ladder of layer #5 and fill it with air or vacuum
C
        DELA(1)=2.
        DELA(2)=3.5
        DELA(3)=4.0
        CALL GSVOLU('ELL5','BOX ',IDTMED(254),DELA,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ELL5','SEEN',0)
C
C --- Define a volume containing the sensitive part of the strips
C     (silicon, layer #5)
C
        DITS(1)=0.0150
        DITS(2)=3.75
        DITS(3)=2.1
        CALL GSVOLU('ITS5','BOX ',IDTMED(250),DITS,3,IOUT)
C
C --- Define a volume containing the electronics of the strips
C     (silicon, layer #5)
C
        DCHI(1)=0.02
        DCHI(2)=3.4
        DCHI(3)=0.525
        CALL GSVOLU('SCH5','BOX ',IDTMED(251),DCHI,3,IOUT)
C
C --- Define the cooling tubes (aluminum, layer #5)
C
        DTUB(1)=0.09
        DTUB(2)=DTUB(1)+0.01
        DTUB(3)=90.22/2.
        CALL GSVOLU('STB5','TUBE',IDTMED(256),DTUB,3,IOUT)
C
C --- Define the cooling fluid (water or freon, layer #5)
C
        DWAT(1)=0.
        DWAT(2)=0.09
        DWAT(3)=90.22/2.
        CALL GSVOLU('SWT5','TUBE',IDTMED(257),DWAT,3,IOUT)   ! water
C        CALL GSVOLU('SWT5','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
C
C --- Define the (triangular) element of the heat bridge (carbon, layer #5)
C
        DFRA(1)=120.
        DFRA(2)=360.
        DFRA(3)=3.
        DFRA(4)=2.
        DFRA(5)=-0.015
        DFRA(6)=4.2*SQRT(3.)/6.
        DFRA(7)=DFRA(6)+0.03
        DFRA(8)=0.015
        DFRA(9)=DFRA(6)
        DFRA(10)=DFRA(7)
        CALL GSVOLU('SFR5','PGON',IDTMED(253),DFRA,10,IOUT)
C
C --- Define the element connecting the triangles of the heat bridge
C     (carbon, layer #5)
C
        DCEI(1)=0.
        DCEI(2)=0.03
        DCEI(3)=90.22/2.
        CALL GSVOLU('SCE5','TUBE',IDTMED(253),DCEI,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of plastic (G10FR4)
C     (layer #5)
C
	DPLA(1)=(10./(8.*7.))/2.
	DPLA(2)=3.5
	DPLA(3)=4.
	CALL GSVOLU('EPL5','BOX ',IDTMED(263),DPLA,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #5)
C
	DCOP(1)=(2./(8.*7.))/2.
	DCOP(2)=3.5
	DCOP(3)=4.
	CALL GSVOLU('ECU5','BOX ',IDTMED(260),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of epoxy (layer #5)
C
	DEPX(1)=(30./(8.*7.))/2.
	DEPX(2)=3.5
	DEPX(3)=4.
	CALL GSVOLU('EPX5','BOX ',IDTMED(263),DEPX,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (bus)
C     (layer #5)
C
	DSIL(1)=(20./(8.*7.))/2.
	DSIL(2)=3.5
	DSIL(3)=4.
	CALL GSVOLU('ESI5','BOX ',IDTMED(252),DSIL,3,IOUT)
C
C --- Place the end-ladder stuff into its mother (ELL5)
C
        SEP=(4.-2.*(DPLA(1)+DCOP(1)+DEPX(1)+DSIL(1)))/3.
	YPOS=0.
	ZPOS=0.
C
C --- Plastic
C
	XPOS=-DELA(1)+DPLA(1)
	CALL GSPOS('EPL5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+DCOP(1)
	CALL GSPOS('ECU5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Epoxy
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+DEPX(1)
	CALL GSPOS('EPX5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+2.*DEPX(1)+
     $  SEP+DSIL(1)
	CALL GSPOS('ESI5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive part of the strips into its mother (ISV1)
C
        YPOS=0.
        DO J=1,23
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=DBOX1(1)-DITS(1)
	  ELSE
	    XPOS=-DBOX1(1)+DITS(1)
          ENDIF
	  ZPOS=(-11.+FLOAT(J-1))*3.91
	  CALL GSPOS('ITS5',J,'ISV1',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO
C
C --- Place the electronics of the strips into its mother (SSV1)
C
        YPOS=0.
        DO J=1,23
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=-DSRV(1)+0.28
	  ELSE
	    XPOS=-DSRV(1)+0.28-2.*DITS(1)-0.03
          ENDIF
	  ZPOS=(-11.+FLOAT(J-1))*3.91+0.85
	  CALL GSPOS('SCH5',J,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO      	
C
C --- Place the cooling tubes and the cooling fluid into their mother (SSV1)
C
      XPOS=-DSRV(1)+0.41
      ZPOS=0.
C
C --- Left tube (just a matter of convention)
C
      YPOS=-2.25-0.1
      CALL GSPOS('STB5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right tube (just a matter of convention)
C
      YPOS=2.25+0.1
      CALL GSPOS('STB5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the heat bridge elements into their mother (SSV1)
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/6.)*4.2
      YPOS=0.
      DO J=1,24
         ZPOS=(-11.+FLOAT(J-1))*3.91-4.2/2.
	 CALL GSPOS('SFR5',J,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      END DO
C
C --- Place the elements connecting the triangles of the heat bridge
C     into their mother (SSV1)
C
      ZPOS=0.
C
C --- Left element (just a matter of convention)
C
      XPOS=-DSRV(1)+0.47
      YPOS=-(2.1+0.015)
      CALL GSPOS('SCE5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right element
C
      XPOS=-DSRV(1)+0.47
      YPOS=2.1+0.015
      CALL GSPOS('SCE5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Top element
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/2.)*4.2+0.015
      YPOS=0.
      CALL GSPOS('SCE5',3,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the strip ladders (ISV1),
C     electronics/cooling (SSV1) and end-ladder stuff (ELL5) of layer #5 in
C     their mother volume (IT56)
C
        OFFSET1=ATG(0.9,40.)
        OFFSET2=5.2
        RZERO=40.+DBOX1(1)
	RUNO=40.+2.*DBOX1(1)+DSRV(1)
        RTWO=40.+2.*DBOX1(1)+DELA(1)	
	DO I=1,35
	  ATHETA=TWOPI*FLOAT(I-1)*RADDEG/35.+OFFSET2
	  CALL SXSROT(IDROTM(1500+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
C
C --- Strip ladders
C
	  XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/35.+OFFSET1)
	  YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/35.+OFFSET1)
          ZPOS=0.
	  CALL GSPOS('ISV1',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
C
C --- Electronics/cooling
C
	  XPOS=RUNO*COS(TWOPI*FLOAT(I-1)/35.+OFFSET1)
	  YPOS=RUNO*SIN(TWOPI*FLOAT(I-1)/35.+OFFSET1)
	  ZPOS=0.
	  CALL GSPOS('SSV1',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
C
C --- End-ladders (nagative-Z and positive-Z)
C
	  XPOS=RTWO*COS(TWOPI*FLOAT(I-1)/35.+OFFSET1)
	  YPOS=RTWO*SIN(TWOPI*FLOAT(I-1)/35.+OFFSET1)
	  ZPOS=-(DBOX1(3)+DELA(3)+6.)
	  CALL GSPOS('ELL5',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
	  ZPOS=DBOX1(3)+DELA(3)+6.
	  CALL GSPOS('ELL5',I+35,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
	END DO
C	
5678    CONTINUE
C
C --- Layer #6
C
C        GOTO 5778           ! skip ITS layer no. 6
C
C --- Define a ghost volume containing a single ladder of layer #6 and fill
C     it with air or vacuum
C
        DBOX2(1)=(0.0600+2.*0.0150)/2.
        DBOX2(2)=3.75
        DBOX2(3)=101.95/2.
        CALL GSVOLU('ISV2','BOX ',IDTMED(254),DBOX2,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ISV2','SEEN',0)
C
C --- Define a ghost volume containing the electronics and cooling of
C     a single ladder of layer #6 and fill it with air or vacuum
C
        DSRV(1)=(0.47+(SQRT(3.)/2.)*4.2+0.05)/2.
        DSRV(2)=3.75
        DSRV(3)=101.95/2.
        CALL GSVOLU('SSV2','BOX ',IDTMED(254), DSRV,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('SSV2','SEEN',0)
C
C --- Define a ghost volume containing the end-ladder stuff of
C     a single ladder of layer #6 and fill it with air or vacuum
C
        DELA(1)=2.
        DELA(2)=3.5
        DELA(3)=4.0
        CALL GSVOLU('ELL6','BOX ',IDTMED(254),DELA,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ELL6','SEEN',0)
C
C --- Define a volume containing the sensitive part of the strips
C     (silicon, layer #6)
C
        DITS(1)=0.0150
        DITS(2)=3.75
        DITS(3)=2.1
        CALL GSVOLU('ITS6','BOX ',IDTMED(250),DITS,3,IOUT)
C
C --- Define a volume containing the electronics of the strips
C     (silicon, layer #6)
C
        DCHI(1)=0.02
        DCHI(2)=3.4
        DCHI(3)=0.525
        CALL GSVOLU('SCH6','BOX ',IDTMED(251),DCHI,3,IOUT)
C
C --- Define the cooling tubes (aluminum, layer #6)
C
        DTUB(1)=0.09
        DTUB(2)=DTUB(1)+0.01
        DTUB(3)=101.95/2.
        CALL GSVOLU('STB6','TUBE',IDTMED(256),DTUB,3,IOUT)
C
C --- Define the cooling fluid (water or freon, layer #6)
C
        DWAT(1)=0.
        DWAT(2)=0.09
        DWAT(3)=101.95/2.
        CALL GSVOLU('SWT6','TUBE',IDTMED(257),DWAT,3,IOUT)   ! water
C        CALL GSVOLU('SWT6','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
C
C --- Define the (triangular) element of the heat bridge (carbon, layer #6)
C
        DFRA(1)=120.
        DFRA(2)=360.
        DFRA(3)=3.
        DFRA(4)=2.
        DFRA(5)=-0.015
        DFRA(6)=4.2*SQRT(3.)/6.
        DFRA(7)=DFRA(6)+0.03
        DFRA(8)=0.015
        DFRA(9)=DFRA(6)
        DFRA(10)=DFRA(7)
        CALL GSVOLU('SFR6','PGON',IDTMED(253),DFRA,10,IOUT)
C
C --- Define the element connecting the triangles of the heat bridge
C     (carbon, layer #6)
C
        DCEI(1)=0.
        DCEI(2)=0.03
        DCEI(3)=101.95/2.
        CALL GSVOLU('SCE6','TUBE',IDTMED(253),DCEI,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of plastic (G10FR4)
C     (layer #6)
C
	DPLA(1)=(10./(8.*7.))/2.
	DPLA(2)=3.5
	DPLA(3)=4.
	CALL GSVOLU('EPL6','BOX ',IDTMED(263),DPLA,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #6)
C
	DCOP(1)=(2./(8.*7.))/2.
	DCOP(2)=3.5
	DCOP(3)=4.
	CALL GSVOLU('ECU6','BOX ',IDTMED(260),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of epoxy (layer #6)
C
	DEPX(1)=(30./(8.*7.))/2.
	DEPX(2)=3.5
	DEPX(3)=4.
	CALL GSVOLU('EPX6','BOX ',IDTMED(263),DEPX,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (bus)
C     (layer #6)
C
	DSIL(1)=(20./(8.*7.))/2.
	DSIL(2)=3.5
	DSIL(3)=4.
	CALL GSVOLU('ESI6','BOX ',IDTMED(252),DSIL,3,IOUT)
C
C --- Place the end-ladder stuff into its mother (ELL5)
C
        SEP=(4.-2.*(DPLA(1)+DCOP(1)+DEPX(1)+DSIL(1)))/3.
	YPOS=0.
	ZPOS=0.
C
C --- Plastic
C
	XPOS=-DELA(1)+DPLA(1)
	CALL GSPOS('EPL6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+DCOP(1)
	CALL GSPOS('ECU6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Epoxy
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+DEPX(1)
	CALL GSPOS('EPX6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+2.*DEPX(1)+
     $  SEP+DSIL(1)
	CALL GSPOS('ESI6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive part of the strips into its mother (ISV2)
C
        YPOS=0.
        DO J=1,26
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=DBOX2(1)-DITS(1)
	  ELSE
	    XPOS=-DBOX2(1)+DITS(1)
          ENDIF
	  ZPOS=(-12.+FLOAT(J-1))*3.91-1.96
	  CALL GSPOS('ITS6',J,'ISV2',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO
C
C --- Place the electronics of the strips into its mother (SSV2)
C
        YPOS=0.
        DO J=1,26
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=-DSRV(1)+0.28
	  ELSE
	    XPOS=-DSRV(1)+0.28-2.*DITS(1)-0.03
          ENDIF
	  ZPOS=(-12.+FLOAT(J-1))*3.91-1.96+0.85
	  CALL GSPOS('SCH5',J,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO      	
C
C --- Place the cooling tubes and the cooling fluid into their mother (SSV1)
C
      XPOS=-DSRV(1)+0.41
      ZPOS=0.
C
C --- Left tube (just a matter of convention)
C
      YPOS=-2.25-0.1
      CALL GSPOS('STB6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right tube (just a matter of convention)
C
      YPOS=2.25+0.1
      CALL GSPOS('STB6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the heat bridge elements into their mother (SSV2)
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/6.)*4.2
      YPOS=0.
      DO J=1,27
         ZPOS=(-12.+FLOAT(J-1))*3.91-1.96-4.2/2.
	 CALL GSPOS('SFR6',J,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      END DO
C
C --- Place the elements connecting the triangles of the heat bridge
C     into their mother (SSV2)
C
      ZPOS=0.
C
C --- Left element (just a matter of convention)
C
      XPOS=-DSRV(1)+0.47
      YPOS=-(2.1+0.015)
      CALL GSPOS('SCE6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right element
C
      XPOS=-DSRV(1)+0.47
      YPOS=2.1+0.015
      CALL GSPOS('SCE6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Top element
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/2.)*4.2+0.015
      YPOS=0.
      CALL GSPOS('SCE6',3,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the strip ladders (ISV2),
C     electronics/cooling (SSV2) and end-ladder stuff (ELL6) of layer #6 in
C     their mother volume (IT56)
C
        OFFSET1=ATG(1.0,45.)
        OFFSET2=5.2
        RZERO=45.+DBOX2(1)
	RUNO=45.+2.*DBOX2(1)+DSRV(1)
        RTWO=45.+2.*DBOX2(1)+DELA(1)	
	DO I=1,39
	  ATHETA=TWOPI*FLOAT(I-1)*RADDEG/39.+OFFSET2
	  CALL SXSROT(IDROTM(1600+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
C
C --- Strip ladders
C
	  XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/39.+OFFSET1)
	  YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/39.+OFFSET1)
          ZPOS=0.
	  CALL GSPOS('ISV2',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
C
C --- Electronics/cooling
C
	  XPOS=RUNO*COS(TWOPI*FLOAT(I-1)/39.+OFFSET1)
	  YPOS=RUNO*SIN(TWOPI*FLOAT(I-1)/39.+OFFSET1)
	  ZPOS=0.
	  CALL GSPOS('SSV2',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
C
C --- End-ladders (nagative-Z and positive-Z)
C
	  XPOS=RTWO*COS(TWOPI*FLOAT(I-1)/39.+OFFSET1)
	  YPOS=RTWO*SIN(TWOPI*FLOAT(I-1)/39.+OFFSET1)
	  ZPOS=-(DBOX2(3)+DELA(3)+6.)
	  CALL GSPOS('ELL6',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
	  ZPOS=DBOX2(3)+DELA(3)+6.
	  CALL GSPOS('ELL6',I+39,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
	END DO
C	
5778    CONTINUE
C
      ENDIF
C
C --- Define SSD with the 32+36 lay-out
C
      IF(IVERS.EQ.43 .OR. IVERS.EQ.44 .OR. IVERS.EQ.45) THEN 	
C
C --- Define ghost volume containing the Strip Detectors and fill it with air
C     or vacuum
C
        XXM=(49.999-3.)/(70.-25.)
        DGH(1)=0.
        DGH(2)=360.
        DGH(3)=4.
        DGH(4)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM-(27.-9.01)/XXM-
     $  (36.-27)/XXM-(49.998-36.)/XXM
        DGH(5)=49.998
        DGH(6)=49.998
        DGH(7)=-25.-(9.-3.01)/XXM-(9.01-9.)/XXM-(27.-9.01)/XXM-
     $  (36.-27)/XXM
        DGH(8)=36.
        DGH(9)=49.998
        DGH(10)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM+(27.-9.01)/XXM+
     $  (36.-27)/XXM
        DGH(11)=36.
        DGH(12)=49.998
        DGH(13)=25.+(9.-3.01)/XXM+(9.01-9.)/XXM+(27.-9.01)/XXM+
     $  (36.-27)/XXM+(49.998-36.)/XXM
        DGH(14)=49.998
        DGH(15)=49.998
        CALL GSVOLU('IT56','PCON',IDTMED(276),DGH,15,IOUT)
        CALL GSPOS('IT56',1,'ITSV',0.,0.,0.,0,'ONLY')
        CALL GSATT('IT56','SEEN',0)
C
C --- Layer #5
C
C        GOTO 6678           ! skip ITS layer no. 5
C
C --- Define a ghost volume containing a single ladder of layer #5 and fill
C     it with air or vacuum
C
        DBOX1(1)=(0.0600+2.*0.0150)/2.
        DBOX1(2)=3.75
        DBOX1(3)=86.31/2.
        CALL GSVOLU('ISV1','BOX ',IDTMED(254),DBOX1,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ISV1','SEEN',0)
C
C --- Define a ghost volume containing the electronics and cooling of
C     a single ladder of layer #5 and fill it with air or vacuum
C
        DSRV(1)=(0.47+(SQRT(3.)/2.)*4.2+0.05)/2.
        DSRV(2)=3.75
        DSRV(3)=86.31/2.
        CALL GSVOLU('SSV1','BOX ',IDTMED(254),DSRV,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('SSV1','SEEN',0)
C
C --- Define a ghost volume containing the end-ladder stuff of
C     a single ladder of layer #5 and fill it with air or vacuum
C
        DELA(1)=2.
        DELA(2)=3.5
        DELA(3)=4.0
        CALL GSVOLU('ELL5','BOX ',IDTMED(254),DELA,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ELL5','SEEN',0)
C
C --- Define a volume containing the sensitive part of the strips
C     (silicon, layer #5)
C
        DITS(1)=0.0150
        DITS(2)=3.75
        DITS(3)=2.1
        CALL GSVOLU('ITS5','BOX ',IDTMED(250),DITS,3,IOUT)
C
C --- Define a volume containing the electronics of the strips
C     (silicon, layer #5)
C
        DCHI(1)=0.02
        DCHI(2)=3.4
        DCHI(3)=0.525
        CALL GSVOLU('SCH5','BOX ',IDTMED(251),DCHI,3,IOUT)
C
C --- Define the cooling tubes (aluminum, layer #5)
C
        DTUB(1)=0.09
        DTUB(2)=DTUB(1)+0.01
        DTUB(3)=86.31/2.
        CALL GSVOLU('STB5','TUBE',IDTMED(256),DTUB,3,IOUT)
C
C --- Define the cooling fluid (water or freon, layer #5)
C
        DWAT(1)=0.
        DWAT(2)=0.09
        DWAT(3)=86.31/2.
        CALL GSVOLU('SWT5','TUBE',IDTMED(257),DWAT,3,IOUT)   ! water
C        CALL GSVOLU('SWT5','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
C
C --- Define the (triangular) element of the heat bridge (carbon, layer #5)
C
        DFRA(1)=120.
        DFRA(2)=360.
        DFRA(3)=3.
        DFRA(4)=2.
        DFRA(5)=-0.015
        DFRA(6)=4.2*SQRT(3.)/6.
        DFRA(7)=DFRA(6)+0.03
        DFRA(8)=0.015
        DFRA(9)=DFRA(6)
        DFRA(10)=DFRA(7)
        CALL GSVOLU('SFR5','PGON',IDTMED(253),DFRA,10,IOUT)
C
C --- Define the element connecting the triangles of the heat bridge
C     (carbon, layer #5)
C
        DCEI(1)=0.
        DCEI(2)=0.03
        DCEI(3)=86.31/2.
        CALL GSVOLU('SCE5','TUBE',IDTMED(253),DCEI,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of plastic (G10FR4)
C     (layer #5)
C
	DPLA(1)=(10./(8.*7.))/2.
	DPLA(2)=3.5
	DPLA(3)=4.
	CALL GSVOLU('EPL5','BOX ',IDTMED(263),DPLA,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #5)
C
	DCOP(1)=(2./(8.*7.))/2.
	DCOP(2)=3.5
	DCOP(3)=4.
	CALL GSVOLU('ECU5','BOX ',IDTMED(260),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of epoxy (layer #5)
C
	DEPX(1)=(30./(8.*7.))/2.
	DEPX(2)=3.5
	DEPX(3)=4.
	CALL GSVOLU('EPX5','BOX ',IDTMED(263),DEPX,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (bus)
C     (layer #5)
C
	DSIL(1)=(20./(8.*7.))/2.
	DSIL(2)=3.5
	DSIL(3)=4.
	CALL GSVOLU('ESI5','BOX ',IDTMED(252),DSIL,3,IOUT)
C
C --- Place the end-ladder stuff into its mother (ELL5)
C
        SEP=(4.-2.*(DPLA(1)+DCOP(1)+DEPX(1)+DSIL(1)))/3.
	YPOS=0.
	ZPOS=0.
C
C --- Plastic
C
	XPOS=-DELA(1)+DPLA(1)
	CALL GSPOS('EPL5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+DCOP(1)
	CALL GSPOS('ECU5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Epoxy
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+DEPX(1)
	CALL GSPOS('EPX5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+2.*DEPX(1)+
     $  SEP+DSIL(1)
	CALL GSPOS('ESI5',1,'ELL5',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive part of the strips into its mother (ISV1)
C
        YPOS=0.
        DO J=1,22
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=DBOX1(1)-DITS(1)
	  ELSE
	    XPOS=-DBOX1(1)+DITS(1)
          ENDIF
	  ZPOS=(-10.+FLOAT(J-1))*3.91-1.96
	  CALL GSPOS('ITS5',J,'ISV1',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO
C
C --- Place the electronics of the strips into its mother (SSV1)
C
        YPOS=0.
        DO J=1,22
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=-DSRV(1)+0.28
	  ELSE
	    XPOS=-DSRV(1)+0.28-2.*DITS(1)-0.03
          ENDIF
	  ZPOS=(-10.+FLOAT(J-1))*3.91-1.96+0.85
	  CALL GSPOS('SCH5',J,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO      	
C
C --- Place the cooling tubes and the cooling fluid into their mother (SSV1)
C
      XPOS=-DSRV(1)+0.41
      ZPOS=0.
C
C --- Left tube (just a matter of convention)
C
      YPOS=-2.25-0.1
      CALL GSPOS('STB5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right tube (just a matter of convention)
C
      YPOS=2.25+0.1
      CALL GSPOS('STB5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the heat bridge elements into their mother (SSV1)
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/6.)*4.2
      YPOS=0.
      DO J=1,23
         ZPOS=(-10.+FLOAT(J-1))*3.91-1.96-4.2/2.
	 CALL GSPOS('SFR5',J,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
      END DO
C
C --- Place the elements connecting the triangles of the heat bridge
C     into their mother (SSV1)
C
      ZPOS=0.
C
C --- Left element (just a matter of convention)
C
      XPOS=-DSRV(1)+0.47
      YPOS=-(2.1+0.015)
      CALL GSPOS('SCE5',1,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right element
C
      XPOS=-DSRV(1)+0.47
      YPOS=2.1+0.015
      CALL GSPOS('SCE5',2,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Top element
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/2.)*4.2+0.015
      YPOS=0.
      CALL GSPOS('SCE5',3,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the strip ladders (ISV1),
C     electronics/cooling (SSV1) and end-ladder stuff (ELL5) of layer #5 in
C     their mother volume (IT56)
C
        OFFSET1=ATG(0.8,36.6)
        OFFSET2=5.2
        RZERO=36.6+DBOX1(1)
	RUNO=36.6+2.*DBOX1(1)+DSRV(1)
        RTWO=36.6+2.*DBOX1(1)+DELA(1)	
	DO I=1,32
	  ATHETA=TWOPI*FLOAT(I-1)*RADDEG/32.+OFFSET2
	  CALL SXSROT(IDROTM(1500+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
C
C --- Strip ladders
C
	  XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/32.+OFFSET1)
	  YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/32.+OFFSET1)
          ZPOS=0.
	  CALL GSPOS('ISV1',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
C
C --- Electronics/cooling
C
	  XPOS=RUNO*COS(TWOPI*FLOAT(I-1)/32.+OFFSET1)
	  YPOS=RUNO*SIN(TWOPI*FLOAT(I-1)/32.+OFFSET1)
	  ZPOS=0.
	  CALL GSPOS('SSV1',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
C
C --- End-ladders (nagative-Z and positive-Z)
C
	  XPOS=RTWO*COS(TWOPI*FLOAT(I-1)/32.+OFFSET1)
	  YPOS=RTWO*SIN(TWOPI*FLOAT(I-1)/32.+OFFSET1)
	  ZPOS=-(DBOX1(3)+DELA(3)+6.)
	  CALL GSPOS('ELL5',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
	  ZPOS=DBOX1(3)+DELA(3)+6.
	  CALL GSPOS('ELL5',I+35,'IT56',XPOS,YPOS,ZPOS,IDROTM(1500+I),
     $    'ONLY')
	END DO
C	
6678    CONTINUE
C
C --- Layer #6
C
C        GOTO 6778           ! skip ITS layer no. 6
C
C --- Define a ghost volume containing a single ladder of layer #6 and fill
C     it with air or vacuum
C
        DBOX2(1)=(0.0600+2.*0.0150)/2.
        DBOX2(2)=3.75
        DBOX2(3)=94.13/2.
        CALL GSVOLU('ISV2','BOX ',IDTMED(254),DBOX2,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ISV2','SEEN',0)
C
C --- Define a ghost volume containing the electronics and cooling of
C     a single ladder of layer #6 and fill it with air or vacuum
C
        DSRV(1)=(0.47+(SQRT(3.)/2.)*4.2+0.05)/2.
        DSRV(2)=3.75
        DSRV(3)=94.13/2.
        CALL GSVOLU('SSV2','BOX ',IDTMED(254),DSRV,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('SSV2','SEEN',0)
C
C --- Define a ghost volume containing the end-ladder stuff of
C     a single ladder of layer #6 and fill it with air or vacuum
C
        DELA(1)=2.
        DELA(2)=3.5
        DELA(3)=4.0
        CALL GSVOLU('ELL6','BOX ',IDTMED(254),DELA,3,IOUT)
C
C --- Make the ghost volume invisible
C
        CALL GSATT('ELL6','SEEN',0)
C
C --- Define a volume containing the sensitive part of the strips
C     (silicon, layer #6)
C
        DITS(1)=0.0150
        DITS(2)=3.75
        DITS(3)=2.1
        CALL GSVOLU('ITS6','BOX ',IDTMED(250),DITS,3,IOUT)
C
C --- Define a volume containing the electronics of the strips
C     (silicon, layer #6)
C
        DCHI(1)=0.02
        DCHI(2)=3.4
        DCHI(3)=0.525
        CALL GSVOLU('SCH6','BOX ',IDTMED(251),DCHI,3,IOUT)
C
C --- Define the cooling tubes (aluminum, layer #6)
C
        DTUB(1)=0.09
        DTUB(2)=DTUB(1)+0.01
        DTUB(3)=94.13/2.
        CALL GSVOLU('STB6','TUBE',IDTMED(256),DTUB,3,IOUT)
C
C --- Define the cooling fluid (water or freon, layer #6)
C
        DWAT(1)=0.
        DWAT(2)=0.09
        DWAT(3)=94.13/2.
        CALL GSVOLU('SWT6','TUBE',IDTMED(257),DWAT,3,IOUT)   ! water
C        CALL GSVOLU('SWT6','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
C
C --- Define the (triangular) element of the heat bridge (carbon, layer #6)
C
        DFRA(1)=120.
        DFRA(2)=360.
        DFRA(3)=3.
        DFRA(4)=2.
        DFRA(5)=-0.015
        DFRA(6)=4.2*SQRT(3.)/6.
        DFRA(7)=DFRA(6)+0.03
        DFRA(8)=0.015
        DFRA(9)=DFRA(6)
        DFRA(10)=DFRA(7)
        CALL GSVOLU('SFR6','PGON',IDTMED(253),DFRA,10,IOUT)
C
C --- Define the element connecting the triangles of the heat bridge
C     (carbon, layer #6)
C
        DCEI(1)=0.
        DCEI(2)=0.03
        DCEI(3)=94.13/2.
        CALL GSVOLU('SCE6','TUBE',IDTMED(253),DCEI,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of plastic (G10FR4)
C     (layer #6)
C
	DPLA(1)=(10./(8.*7.))/2.
	DPLA(2)=3.5
	DPLA(3)=4.
	CALL GSVOLU('EPL6','BOX ',IDTMED(263),DPLA,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of copper (layer #6)
C
	DCOP(1)=(2./(8.*7.))/2.
	DCOP(2)=3.5
	DCOP(3)=4.
	CALL GSVOLU('ECU6','BOX ',IDTMED(260),DCOP,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of epoxy (layer #6)
C
	DEPX(1)=(30./(8.*7.))/2.
	DEPX(2)=3.5
	DEPX(3)=4.
	CALL GSVOLU('EPX6','BOX ',IDTMED(263),DEPX,3,IOUT)
C
C --- Define the part of the end-ladder stuff made of silicon (bus)
C     (layer #6)
C
	DSIL(1)=(20./(8.*7.))/2.
	DSIL(2)=3.5
	DSIL(3)=4.
	CALL GSVOLU('ESI6','BOX ',IDTMED(252),DSIL,3,IOUT)
C
C --- Place the end-ladder stuff into its mother (ELL5)
C
        SEP=(4.-2.*(DPLA(1)+DCOP(1)+DEPX(1)+DSIL(1)))/3.
	YPOS=0.
	ZPOS=0.
C
C --- Plastic
C
	XPOS=-DELA(1)+DPLA(1)
	CALL GSPOS('EPL6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Copper
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+DCOP(1)
	CALL GSPOS('ECU6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Epoxy
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+DEPX(1)
	CALL GSPOS('EPX6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Silicon (bus)
C
	XPOS=-DELA(1)+2.*DPLA(1)+SEP+2.*DCOP(1)+SEP+2.*DEPX(1)+
     $  SEP+DSIL(1)
	CALL GSPOS('ESI6',1,'ELL6',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the sensitive part of the strips into its mother (ISV2)
C
        YPOS=0.
        DO J=1,24
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=-DBOX2(1)+DITS(1)
	  ELSE
	    XPOS=DBOX2(1)-DITS(1)
          ENDIF
	  ZPOS=(-11.+FLOAT(J-1))*3.91-1.96
	  CALL GSPOS('ITS6',J,'ISV2',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO
C
C --- Place the electronics of the strips into its mother (SSV2)
C
        YPOS=0.
        DO J=1,24
          IF(MOD(J,2).EQ.0) THEN
	    XPOS=-DSRV(1)+0.28-2.*DITS(1)-0.03
	  ELSE
	    XPOS=-DSRV(1)+0.28
          ENDIF
	  ZPOS=(-11.+FLOAT(J-1))*3.91-1.96+0.85
	  CALL GSPOS('SCH5',J,'SSV1',XPOS,YPOS,ZPOS,0,'ONLY')
        END DO      	
C
C --- Place the cooling tubes and the cooling fluid into their mother (SSV2)
C
      XPOS=-DSRV(1)+0.41
      ZPOS=0.
C
C --- Left tube (just a matter of convention)
C
      YPOS=-2.25-0.1
      CALL GSPOS('STB6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right tube (just a matter of convention)
C
      YPOS=2.25+0.1
      CALL GSPOS('STB6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      CALL GSPOS('SWT6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the heat bridge elements into their mother (SSV2)
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/6.)*4.2
      YPOS=0.
      DO J=1,25
         ZPOS=(-11.+FLOAT(J-1))*3.91-1.96-4.2/2.
	 CALL GSPOS('SFR6',J,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
      END DO
C
C --- Place the elements connecting the triangles of the heat bridge
C     into their mother (SSV2)
C
      ZPOS=0.
C
C --- Left element (just a matter of convention)
C
      XPOS=-DSRV(1)+0.47
      YPOS=-(2.1+0.015)
      CALL GSPOS('SCE6',1,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Right element
C
      XPOS=-DSRV(1)+0.47
      YPOS=2.1+0.015
      CALL GSPOS('SCE6',2,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Top element
C
      XPOS=-DSRV(1)+0.47+(SQRT(3.)/2.)*4.2+0.015
      YPOS=0.
      CALL GSPOS('SCE6',3,'SSV2',XPOS,YPOS,ZPOS,0,'ONLY')
C
C --- Place the ghost volumes containing the strip ladders (ISV2),
C     electronics/cooling (SSV2) and end-ladder stuff (ELL6) of layer #6 in
C     their mother volume (IT56)
C
        OFFSET1=ATG(0.9,41.2)
        OFFSET2=5.2
        RZERO=41.2+DBOX2(1)
	RUNO=41.2+2.*DBOX2(1)+DSRV(1)
        RTWO=41.2+2.*DBOX2(1)+DELA(1)	
	DO I=1,36
	  ATHETA=TWOPI*FLOAT(I-1)*RADDEG/36.+OFFSET2
	  CALL SXSROT(IDROTM(1600+I),90.,ATHETA,90.,90.+ATHETA,0.,0.)
C
C --- Strip ladders
C
	  XPOS=RZERO*COS(TWOPI*FLOAT(I-1)/36.+OFFSET1)
	  YPOS=RZERO*SIN(TWOPI*FLOAT(I-1)/36.+OFFSET1)
          ZPOS=0.
	  CALL GSPOS('ISV2',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
C
C --- Electronics/cooling
C
	  XPOS=RUNO*COS(TWOPI*FLOAT(I-1)/36.+OFFSET1)
	  YPOS=RUNO*SIN(TWOPI*FLOAT(I-1)/36.+OFFSET1)
	  ZPOS=0.
	  CALL GSPOS('SSV2',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
C
C --- End-ladders (nagative-Z and positive-Z)
C
	  XPOS=RTWO*COS(TWOPI*FLOAT(I-1)/36.+OFFSET1)
	  YPOS=RTWO*SIN(TWOPI*FLOAT(I-1)/36.+OFFSET1)
	  ZPOS=-(DBOX2(3)+DELA(3)+6.)
	  CALL GSPOS('ELL6',I,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
	  ZPOS=DBOX2(3)+DELA(3)+6.
	  CALL GSPOS('ELL6',I+39,'IT56',XPOS,YPOS,ZPOS,IDROTM(1600+I),
     $    'ONLY')
	END DO
C	
6778    CONTINUE
C
      ENDIF
C	
C ************************************************************************
C *                                                                      *
C *             E N D - C A P S     A N D     F R A M E S                *
C *             =========================================                *
C *                                                                      *
C ************************************************************************
C
C --- Define a dummy cylinder for multiple scattering tests
C
C      GOTO 7890             ! skip dummy cylinder for multiple scattering tests
C
C      DITS(1)=49.
C      DITS(2)=DITS(1)+0.1
C      DITS(3)=60.3
C      CALL GSVOLU('ITST','TUBE',IDTMED(255),DITS,3,IOUT)
C      CALL GSPOS('ITST',1,'ITSV',0.,0.,0.,0,'ONLY')
C7890  CONTINUE
C
C --- The 0.74% X0 outer wall (C) of the gas vessel at r=50cm ---
C
      IF(IVERS.EQ.40 .OR. IVERS.EQ.43) THEN
C
C        GOTO 8901                                  ! skip outer wall
C	
        DITS(1)=49.9
        DITS(2)=DITS(1)+0.06926
        DITS(3)=62.7+2.*DPCB(3)-10.5         ! old value 60.3
        CALL GSVOLU('ITSG','TUBE',IDTMED(275),DITS,3,I)
        CALL GSPOS('ITSG',1,'ITSV',0.,0.,0.,0,'ONLY')
      ELSE
        GOTO 8901
      ENDIF
8901  CONTINUE
C
C --- The frame between the end-caps (octagonal lay-out) ---
C
C      GOTO 9012                                    ! skip octagonal frame
C
      IF(IVERS.EQ.41) THEN
C
        RZERO=34.0
        DTRA(1)=0.92
        DTRA(2)=1.0
        DTRA(3)=50.5+2.*DPCB(3)-10.5
        DTRA1(1)=0.92
        DTRA1(2)=1.0
	DTRA1(3)=SQRT(DTRA(3)*DTRA(3)+(55.4*55.4-50.5*50.5))/2.
        ANGLE=360./8.
        OFFSET=ANGLE/2.
        DO I=1,8
           XTRA(I)=RZERO*COS(FLOAT(I-1)*ANGLE*DEGRAD)
           YTRA(I)=RZERO*SIN(FLOAT(I-1)*ANGLE*DEGRAD)
           ZTRA(I)=0.
           CALL GSVOLU(NATRA(I),'TUBE',IDTMED(275),DTRA,3,IOUT)
           CALL GSPOS(NATRA(I),1,'ITSV',XTRA(I),YTRA(I),ZTRA(I),0,
     $     'ONLY')
        END DO
C
        ATHETA=22.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(1)+XTRA(2))/2.
        YPOS=(YTRA(1)+YTRA(2))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(1),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5101),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5101),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(2),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5102),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5102),'ONLY'
     $       )
C
        ATHETA=67.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(2)+XTRA(3))/2.
        YPOS=(YTRA(2)+YTRA(3))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(3),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5103),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5103),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(4),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5104),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5104),'ONLY'
     $       )
C
        ATHETA=112.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(3)+XTRA(4))/2.
        YPOS=(YTRA(3)+YTRA(4))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(5),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5105),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5105),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(6),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5106),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5106),'ONLY'
     $       )
C
        ATHETA=157.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(4)+XTRA(5))/2.
        YPOS=(YTRA(4)+YTRA(5))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(7),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5107),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(7),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5107),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(8),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5108),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(8),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5108),'ONLY'
     $       )
C
        ATHETA=22.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(5)+XTRA(6))/2.
        YPOS=(YTRA(5)+YTRA(6))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(9),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5109),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(9),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5109),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(10),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5110),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(10),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5110)
     $       ,'ONLY')
C
        ATHETA=67.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(6)+XTRA(7))/2.
        YPOS=(YTRA(6)+YTRA(7))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(11),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5111),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(11),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5111)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(12),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5112),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(12),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5112)
     $       ,'ONLY')
C
        ATHETA=112.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(7)+XTRA(8))/2.
        YPOS=(YTRA(7)+YTRA(8))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(13),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5113),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(13),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5113)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(14),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5114),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(14),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5114)
     $       ,'ONLY')
C
        ATHETA=157.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(8)+XTRA(1))/2.
        YPOS=(YTRA(8)+YTRA(1))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(15),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5115),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(15),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5115)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(16),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5116),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(16),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5116)
     $       ,'ONLY')
C
C
      ELSEIF(IVERS.EQ.44) THEN
C
C
        RZERO=34.0
        DTRA(1)=0.92
        DTRA(2)=1.0
        DTRA(3)=50.5+2.*DPCB(3)-10.5
        DTRA1(1)=0.92
        DTRA1(2)=1.0
	DTRA1(3)=SQRT(DTRA(3)*DTRA(3)+(55.4*55.4-50.5*50.5))/2.
        ANGLE=360./8.
        OFFSET=ANGLE/2.
        DO I=1,8
           XTRA(I)=RZERO*COS(FLOAT(I-1)*ANGLE*DEGRAD)
           YTRA(I)=RZERO*SIN(FLOAT(I-1)*ANGLE*DEGRAD)
           ZTRA(I)=0.
           CALL GSVOLU(NATRA(I),'TUBE',IDTMED(275),DTRA,3,IOUT)
           CALL GSPOS(NATRA(I),1,'ITSV',XTRA(I),YTRA(I),ZTRA(I),0,
     $     'ONLY')
        END DO
C
        ATHETA=22.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(1)+XTRA(2))/2.
        YPOS=(YTRA(1)+YTRA(2))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(1),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5101),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5101),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(2),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5102),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5102),'ONLY'
     $       )
C
        ATHETA=67.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(2)+XTRA(3))/2.
        YPOS=(YTRA(2)+YTRA(3))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(3),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5103),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5103),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(4),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5104),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5104),'ONLY'
     $       )
C
        ATHETA=112.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(3)+XTRA(4))/2.
        YPOS=(YTRA(3)+YTRA(4))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(5),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5105),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5105),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(6),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5106),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5106),'ONLY'
     $       )
C
        ATHETA=157.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(4)+XTRA(5))/2.
        YPOS=(YTRA(4)+YTRA(5))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(7),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5107),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(7),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5107),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(8),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5108),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(8),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5108),'ONLY'
     $       )
C
        ATHETA=22.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(5)+XTRA(6))/2.
        YPOS=(YTRA(5)+YTRA(6))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(9),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5109),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(9),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5109),'ONLY'
     $       )
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(10),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5110),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(10),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5110)
     $       ,'ONLY')
C
        ATHETA=67.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(6)+XTRA(7))/2.
        YPOS=(YTRA(6)+YTRA(7))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(11),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5111),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(11),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5111)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(12),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5112),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(12),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5112)
     $       ,'ONLY')
C
        ATHETA=112.5
        APHI2=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA(7)+XTRA(8))/2.
        YPOS=(YTRA(7)+YTRA(8))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(13),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5113),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(13),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5113)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(14),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5114),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(14),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5114)
     $       ,'ONLY')
C
        ATHETA=157.5
        APHI1=ACOS(DTRA(3)/SQRT(DTRA(3)*DTRA(3)+
     $  ((50.5/COS(28.*DEGRAD))*(50.5/COS(28.*DEGRAD))-50.5*50.5)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA(8)+XTRA(1))/2.
        YPOS=(YTRA(8)+YTRA(1))/2.
        ZPOS=DTRA(3)/2.
        CALL GSVOLU(NATRA1(15),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5115),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(15),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5115)
     $       ,'ONLY')
        ZPOS=-DTRA(3)/2.
        CALL GSVOLU(NATRA1(16),'TUBE',IDTMED(275),DTRA1,3,IOUT)
        CALL SXSROT(IDROTM(5116),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA1(16),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5116)
     $       ,'ONLY')
      ELSE
        GOTO 9012
      ENDIF
C
9012  CONTINUE
C
C
C --- The frame between the end-caps (hexagonal lay-out) ---
C
C      GOTO 9123                                    ! skip hexagonal frame
C
      IF(IVERS.EQ.42) THEN
C
        RZERO=33.5
        DTRA2(1)=0.92
        DTRA2(2)=1.0
        DTRA2(3)=50.+2.*DPCB(3)-10.5
        DTRA3(1)=0.92
        DTRA3(2)=1.0
        DTRA3(3)=16.75
        DTRA4(1)=0.92
        DTRA4(2)=1.0
        DTRA4(3)=SQRT(DTRA2(3)*DTRA2(3)+(59.9*59.9-50.*50.))/2.
        ANGLE=360./6.
        OFFSET=ANGLE/2.
        DO I=1,6
           XTRA1(I)=RZERO*COS((FLOAT(I-1)*ANGLE+OFFSET)*DEGRAD)
           YTRA1(I)=RZERO*SIN((FLOAT(I-1)*ANGLE+OFFSET)*DEGRAD)
           ZTRA1(I)=0.
           CALL GSVOLU(NATRA2(I),'TUBE',IDTMED(275),DTRA2,3,IOUT)
           CALL GSPOS(NATRA2(I),1,'ITSV',XTRA1(I),YTRA1(I),ZTRA1(I)
     $          ,0,'ONLY')
        END DO
C
        ATHETA=60.
        APHI=90.
        XPOS=(XTRA1(1)+XTRA1(2))/2.
        YPOS=(YTRA1(1)+YTRA1(2))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(1),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5201),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5201),'ONLY'
     $       )
C
        ATHETA=120.
        APHI=90.
        XPOS=(XTRA1(2)+XTRA1(3))/2.
        YPOS=(YTRA1(2)+YTRA1(3))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(2),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5202),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5202),'ONLY'
     $       )
C
        ATHETA=180.
        APHI=90.
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(3),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5203),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5203),'ONLY'
     $       )
C
        ATHETA=60.
        APHI=90.
        XPOS=(XTRA1(4)+XTRA1(5))/2.
        YPOS=(YTRA1(4)+YTRA1(5))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(4),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5204),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5204),'ONLY'
     $       )
C
        ATHETA=120.
        APHI=90.
        XPOS=(XTRA1(5)+XTRA1(6))/2.
        YPOS=(YTRA1(5)+YTRA1(6))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(5),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5205),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5205),'ONLY'
     $       )
C
        ATHETA=180.
        APHI=90.
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(6),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5206),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5206),'ONLY'
     $       )
C
        ATHETA=60.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(1)+XTRA1(2))/2.
        YPOS=(YTRA1(1)+YTRA1(2))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(1),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5211),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5211),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(2),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5212),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5212),'ONLY'
     $       )
C
        ATHETA=120.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(2)+XTRA1(3))/2.
        YPOS=(YTRA1(2)+YTRA1(3))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(3),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5213),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5213),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(4),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5214),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5214),'ONLY'
     $       )
C
        ATHETA=180.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(5),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5215),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5215),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(6),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5216),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5216),'ONLY'
     $       )
        ATHETA=180.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(7),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5217),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(7),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5217),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(8),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5218),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(8),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5218),'ONLY'
     $       )
C
        ATHETA=60.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(4)+XTRA1(5))/2.
        YPOS=(YTRA1(4)+YTRA1(5))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(9),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5219),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(9),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5219),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(10),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5220),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(10),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5220)
     $       ,'ONLY')
C
        ATHETA=120.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(5)+XTRA1(6))/2.
        YPOS=(YTRA1(5)+YTRA1(6))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(11),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5221),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(11),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5221)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(12),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5222),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(12),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5222)
     $       ,'ONLY')
C
        ATHETA=180.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(13),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5223),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(13),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5223)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(14),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5224),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(14),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5224)
     $       ,'ONLY')
        ATHETA=180.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(15),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5225),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(15),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5225)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(16),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5226),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(16),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5226)
     $       ,'ONLY')

C
C
      ELSEIF(IVERS.EQ.45) THEN
C
C
        RZERO=33.5
        DTRA2(1)=0.92
        DTRA2(2)=1.0
        DTRA2(3)=50.+2.*DPCB(3)-10.5
        DTRA3(1)=0.92
        DTRA3(2)=1.0
        DTRA3(3)=16.75
        DTRA4(1)=0.92
        DTRA4(2)=1.0
        DTRA4(3)=SQRT(DTRA2(3)*DTRA2(3)+(59.9*59.9-50.*50.))/2.
        ANGLE=360./6.
        OFFSET=ANGLE/2.
        DO I=1,6
           XTRA1(I)=RZERO*COS((FLOAT(I-1)*ANGLE+OFFSET)*DEGRAD)
           YTRA1(I)=RZERO*SIN((FLOAT(I-1)*ANGLE+OFFSET)*DEGRAD)
           ZTRA1(I)=0.
           CALL GSVOLU(NATRA2(I),'TUBE',IDTMED(275),DTRA2,3,IOUT)
           CALL GSPOS(NATRA2(I),1,'ITSV',XTRA1(I),YTRA1(I),ZTRA1(I)
     $          ,0,'ONLY')
        END DO
C
        ATHETA=60.
        APHI=90.
        XPOS=(XTRA1(1)+XTRA1(2))/2.
        YPOS=(YTRA1(1)+YTRA1(2))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(1),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5201),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5201),'ONLY'
     $       )
C
        ATHETA=120.
        APHI=90.
        XPOS=(XTRA1(2)+XTRA1(3))/2.
        YPOS=(YTRA1(2)+YTRA1(3))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(2),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5202),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5202),'ONLY'
     $       )
C
        ATHETA=180.
        APHI=90.
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(3),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5203),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5203),'ONLY'
     $       )
C
        ATHETA=60.
        APHI=90.
        XPOS=(XTRA1(4)+XTRA1(5))/2.
        YPOS=(YTRA1(4)+YTRA1(5))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(4),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5204),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5204),'ONLY'
     $       )
C
        ATHETA=120.
        APHI=90.
        XPOS=(XTRA1(5)+XTRA1(6))/2.
        YPOS=(YTRA1(5)+YTRA1(6))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(5),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5205),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5205),'ONLY'
     $       )
C
        ATHETA=180.
        APHI=90.
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=0.
        CALL GSVOLU(NATRA3(6),'TUBE',IDTMED(275),DTRA3,3,IOUT)
        CALL SXSROT(IDROTM(5206),90.,ATHETA,90.+APHI,90.+ATHETA,APHI,90.
     $       +ATHETA)
        CALL GSPOS(NATRA3(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5206),'ONLY'
     $       )
C
        ATHETA=60.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(1)+XTRA1(2))/2.
        YPOS=(YTRA1(1)+YTRA1(2))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(1),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5211),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(1),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5211),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(2),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5212),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(2),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5212),'ONLY'
     $       )
C
        ATHETA=120.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(2)+XTRA1(3))/2.
        YPOS=(YTRA1(2)+YTRA1(3))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(3),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5213),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(3),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5213),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(4),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5214),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(4),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5214),'ONLY'
     $       )
C
        ATHETA=180.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(5),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5215),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(5),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5215),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(6),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5216),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(6),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5216),'ONLY'
     $       )
        ATHETA=180.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(3)+XTRA1(4))/2.
        YPOS=(YTRA1(3)+YTRA1(4))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(7),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5217),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(7),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5217),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(8),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5218),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(8),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5218),'ONLY'
     $       )
C
        ATHETA=60.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(4)+XTRA1(5))/2.
        YPOS=(YTRA1(4)+YTRA1(5))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(9),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5219),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(9),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5219),'ONLY'
     $       )
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(10),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5220),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(10),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5220)
     $       ,'ONLY')
C
        ATHETA=120.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(5)+XTRA1(6))/2.
        YPOS=(YTRA1(5)+YTRA1(6))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(11),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5221),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(11),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5221)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(12),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5222),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(12),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5222)
     $       ,'ONLY')
C
        ATHETA=180.
        APHI2=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI1=180.-APHI2
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(13),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5223),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(13),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5223)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(14),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5224),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(14),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5224)
     $       ,'ONLY')
        ATHETA=180.
        APHI1=ACOS(DTRA2(3)/SQRT(DTRA2(3)*DTRA2(3)+
     $  ((50./COS(34.*DEGRAD))*(50./COS(34.*DEGRAD))-50.*50.)))*
     $  RADDEG
        APHI2=180.-APHI1
        XPOS=(XTRA1(6)+XTRA1(1))/2.
        YPOS=(YTRA1(6)+YTRA1(1))/2.
        ZPOS=DTRA2(3)/2.
        CALL GSVOLU(NATRA4(15),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5225),90.,ATHETA,90.+APHI1,90.+ATHETA,APHI1
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(15),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5225)
     $       ,'ONLY')
        ZPOS=-DTRA2(3)/2.
        CALL GSVOLU(NATRA4(16),'TUBE',IDTMED(275),DTRA4,3,IOUT)
        CALL SXSROT(IDROTM(5226),90.,ATHETA,90.+APHI2,90.+ATHETA,APHI2
     $       ,90.+ATHETA)
        CALL GSPOS(NATRA4(16),1,'ITSV',XPOS,YPOS,ZPOS,IDROTM(5226)
     $       ,'ONLY')
      ELSE
        GOTO 9123
      ENDIF
C
 9123 CONTINUE
C
C --- Define the end-caps
C
C      GOTO 9234             ! skip both end-caps
C
C --- Define the Z>0 end-cap
C
C      GOTO 9345             ! skip the Z>0 end-cap
C
      DCONE(1)=(338.-3.)/2./10.
      DCONE(2)=12.
      DCONE(3)=12.02
      DCONE(4)=(338.-3.)*455./(338.-3.-10.)/10.
      DCONE(5)=(338.-3.)*455./(338.-3.-10.)/10.+0.02/COS(45.*DEGRAD)
      XPOS=0.
      YPOS=0.
      ZPOS=(583.+(338.-3.))/2./10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      CALL GSVOLU('RCON','CONE',IDTMED(275),DCONE,5,IOUT)
      CALL GSPOS('RCON',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DTUBE(1)=(338.-3.)*455./(338.-3.-10.)/10.+0.02/COS(45.*DEGRAD)
      DTUBE(2)=49.9   ! In the Simonetti's drawings 52. In the TP 50.
      DTUBE(3)=0.15
      XPOS=0.
      YPOS=0.
      ZPOS=(583./2.+(338-1.5))/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      CALL GSVOLU('RTB1','TUBE',IDTMED(275),DTUBE,3,IOUT)
      CALL GSPOS('RTB1',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DTUBE(1)=10.5
      DTUBE(2)=12.0
      DTUBE(3)=26.8/2./10.
      XPOS=0.
      YPOS=0.
      ZPOS=(583./2.-89.+26.8/2.)/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      CALL GSVOLU('RTB2','TUBE',IDTMED(275),DTUBE,3,IOUT)
      CALL GSPOS('RTB2',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DPGON(1)=15.
      DPGON(2)=360.
      DPGON(3)=12.
      DPGON(4)=2
      DPGON(5)=(583./2.-62.2)/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      DPGON(6)=12.0
      DPGON(7)=13.5
      DPGON(8)=583./2./10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      DPGON(9)=12.0
      DPGON(10)=13.5
      XPOS=0.
      YPOS=0.
      ZPOS=0.
      CALL GSVOLU('RP03','PGON',IDTMED(275),DPGON,10,IOUT)
      CALL GSPOS('RP03',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DPGON(1)=7.5
      DPGON(2)=360.
      DPGON(3)=24.
      DPGON(4)=2
      DPGON(5)=(583./2.+(338.-273.))/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
      DPGON(6)=21.
      DPGON(7)=23.
      DPGON(8)=(583./2.+(338.-273.+15.))/10.+2.*DPCB(3)-10.5 ! end-ladder electronics
      DPGON(9)=21.
      DPGON(10)=23.
      XPOS=0.
      YPOS=0.
      ZPOS=0.
      CALL GSVOLU('RP04','PGON',IDTMED(275),DPGON,10,IOUT)
      CALL GSPOS('RP04',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      IF(IVERS.EQ.40 .OR. IVERS.EQ.41 .OR. IVERS.EQ.42) THEN
         OFFSET2=5.2
         DPGON(1)=360./(2.*35.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=35.
         DPGON(4)=2
         DPGON(5)=(583./2.+(338.-106.))/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(6)=37.7
         DPGON(7)=40.
         DPGON(8)=(583./2.+(338.-106.+15.))/10.+2.*DPCB(3)-10.5 ! end-ladder electr.
         DPGON(9)=37.7
         DPGON(10)=40.
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('RP05','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('RP05',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
         DPGON(1)=360./(2.*39.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=39.
         DPGON(4)=2
         DPGON(5)=(583./2.+(338.-56.))/10.+2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(6)=42.7
         DPGON(7)=45.
         DPGON(8)=(583./2.+(338.-56.+15.))/10.+2.*DPCB(3)-10.5  ! end-ladder electr.
         DPGON(9)=42.7
         DPGON(10)=45.
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('RP06','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('RP06',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
      ENDIF	
      IF(IVERS.EQ.43 .OR. IVERS.EQ.44 .OR. IVERS.EQ.45) THEN
         OFFSET2=5.2
         DPGON(1)=360./(2.*32.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=32.
         DPGON(4)=2
         DPGON(5)=(583./2.+(338.-106.))/10.-(40.-36.6)/TAN(45.*DEGRAD)+
     $   2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(6)=34.3
         DPGON(7)=36.6
         DPGON(8)=(583./2.+(338.-106.+15.))/10.-
     $   (40.-36.6)/TAN(45.*DEGRAD)+2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(9)=34.3
         DPGON(10)=36.6
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('RP05','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('RP05',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
         DPGON(1)=360./(2.*36.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=36.
         DPGON(4)=2
         DPGON(5)=(583./2.+(338.-56.))/10.-(45.-41.2)/TAN(45.*DEGRAD)+
     $   2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(6)=38.9
         DPGON(7)=41.2
         DPGON(8)=(583./2.+(338.-56.+15.))/10.-
     $   (45.-41.2)/TAN(45.*DEGRAD)+2.*DPCB(3)-10.5  ! end-ladder electronics
         DPGON(9)=38.9
         DPGON(10)=41.2
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('RP06','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('RP06',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
      ENDIF	
C
C 9345 CONTINUE
C
C --- Define the Z<0 end-cap
C
C      GOTO 9456             ! skip the Z<0 end-cap
C
      DCONE(1)=(338.-3.)/2./10.
      DCONE(2)=(338.-3.)*455./(338.-3.-10.)/10.
      DCONE(3)=(338.-3.)*455./(338.-3.-10.)/10.+0.02/COS(45.*DEGRAD)
      DCONE(4)=12.
      DCONE(5)=12.02
      XPOS=0.
      YPOS=0.
      ZPOS=-(583.+(338.-3.))/2./10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      CALL GSVOLU('LCON','CONE',IDTMED(275),DCONE,5,IOUT)
      CALL GSPOS('LCON',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DTUBE(1)=(338.-3.)*455./(338.-3.-10.)/10.+0.02/COS(45.*DEGRAD)
      DTUBE(2)=49.9	! In the Simonetti's drawings 52. In the TP 50.
      DTUBE(3)=0.15
      XPOS=0.
      YPOS=0.
      ZPOS=-(583./2.+(338-1.5))/10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      CALL GSVOLU('LTB1','TUBE',IDTMED(275),DTUBE,3,IOUT)
      CALL GSPOS('LTB1',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DTUBE(1)=10.5
      DTUBE(2)=12.0
      DTUBE(3)=26.8/2./10.
      XPOS=0.
      YPOS=0.
      ZPOS=-(583./2.-89.+26.8/2.)/10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      CALL GSVOLU('LTB2','TUBE',IDTMED(275),DTUBE,3,IOUT)
      CALL GSPOS('LTB2',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DPGON(1)=15.
      DPGON(2)=360.
      DPGON(3)=12.
      DPGON(4)=2
      DPGON(5)=-583./2./10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      DPGON(6)=12.0
      DPGON(7)=13.5
      DPGON(8)=-(583./2.-62.2)/10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      DPGON(9)=12.0
      DPGON(10)=13.5
      XPOS=0.
      YPOS=0.
      ZPOS=0.
      CALL GSVOLU('LP03','PGON',IDTMED(275),DPGON,10,IOUT)
      CALL GSPOS('LP03',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      DPGON(1)=7.5
      DPGON(2)=360.
      DPGON(3)=24.
      DPGON(4)=2
      DPGON(5)=-(583./2.+(338.-273.+15.))/10.-2.*DPCB(3)+10.5 ! end-ladder electr.
      DPGON(6)=21.
      DPGON(7)=23.
      DPGON(8)=-(583./2.+(338.-273.))/10.-2.*DPCB(3)+10.5  ! end-ladder electronics
      DPGON(9)=21.
      DPGON(10)=23.
      XPOS=0.
      YPOS=0.
      ZPOS=0.
      CALL GSVOLU('LP04','PGON',IDTMED(275),DPGON,10,IOUT)
      CALL GSPOS('LP04',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
      IF(IVERS.EQ.40 .OR. IVERS.EQ.41 .OR. IVERS.EQ.42) THEN
         OFFSET2=5.2
         DPGON(1)=360./(2.*35.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=35.
         DPGON(4)=2
         DPGON(5)=-(583./2.+(338.-106.))/10.-2.*DPCB(3)+10.5  ! end-ladder electr.
         DPGON(6)=37.7
         DPGON(7)=40.
         DPGON(8)=-(583./2.+(338.-106.+15.))/10.-2.*DPCB(3)+10.5 ! end-ladder electr.
         DPGON(9)=37.7
         DPGON(10)=40.
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('LP05','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('LP05',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
         DPGON(1)=360./(2.*39.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=39.
         DPGON(4)=2
         DPGON(5)=-(583./2.+(338.-56.))/10.-2.*DPCB(3)+10.5  ! end-ladder electronics
         DPGON(6)=42.7
         DPGON(7)=45.
         DPGON(8)=-(583./2.+(338.-56.+15.))/10.-2.*DPCB(3)+10.5  ! end-ladder electr.
         DPGON(9)=42.7
         DPGON(10)=45.
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('LP06','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('LP06',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
      ENDIF	
      IF(IVERS.EQ.43 .OR. IVERS.EQ.44 .OR. IVERS.EQ.45) THEN
         OFFSET2=5.2
         DPGON(1)=360./(2.*32.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=32.
         DPGON(4)=2
         DPGON(5)=-(583./2.+(338.-106.))/10.+(40.-36.6)/TAN(45.*DEGRAD)-
     $   2.*DPCB(3)+10.5  ! end-ladder electronics
         DPGON(6)=34.3
         DPGON(7)=36.6
         DPGON(8)=-(583./2.+(338.-106.+15.))/10.+
     $   (40.-36.6)/TAN(45.*DEGRAD)-2.*DPCB(3)+10.5  ! end-ladder electronics
         DPGON(9)=34.3
         DPGON(10)=36.6
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('LP05','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('LP05',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
C
         DPGON(1)=360./(2.*36.)+OFFSET2
         DPGON(2)=360.
         DPGON(3)=36.
         DPGON(4)=2
         DPGON(5)=-(583./2.+(338.-56.))/10.+(45.-41.2)/TAN(45.*DEGRAD)-
     $   2.*DPCB(3)+10.5  ! end-ladder electronics
         DPGON(6)=38.9
         DPGON(7)=41.2
         DPGON(8)=-(583./2.+(338.-56.+15.))/10.+
     $   (45.-41.2)/TAN(45.*DEGRAD)-2.*DPCB(3)+10.5  ! end-ladder electronics
         DPGON(9)=38.9
         DPGON(10)=41.2
         XPOS=0.
         YPOS=0.
         ZPOS=0.
         CALL GSVOLU('LP06','PGON',IDTMED(275),DPGON,10,IOUT)
         CALL GSPOS('LP06',1,'ITSV',XPOS,YPOS,ZPOS,0,'ONLY')
      ENDIF	
C
C 9456 CONTINUE
C
 9234 CONTINUE
C
 9999 CONTINUE
      RETURN
      END
