*
* $Id$
*
* $Log$
* Revision 1.3  1996/10/01 14:12:05  ravndal
* right units for SPHE specs.
*
* Revision 1.2  1996/09/30 14:54:13  ravndal
* Right units for the spec of PARA
*
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 20/07/94  18.08.26  by  S.Ravndal
*-- Author :
      SUBROUTINE GDSPEC(NAME)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws specifications of volume NAME         *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDFSPC                      *
C.    *       Authors : P.Zanarini   *********                         *
C.    *                 A.McPherson  *****                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcshno.inc"
      CHARACTER*4 ICTUB(11)
      CHARACTER*4 NAME,NAMSEE,ISON
      CHARACTER*4 IBOX(5),ITRD1(5),ITUBE(5),ITUBS(5)
      CHARACTER*4 ITRD2(6),ICON(5),ICONS(7),ISPH(6),ITRAP(11),IPARA(6),
     +IPGON(7),IPCON(6),IGTRA(12),IHYPE(4),IELTU(3)
      DIMENSION PAR(50),IPAR(12),IPA(3),ISPAR(3)
      DIMENSION U0(3),V0(3),THE(3),PHI(3),ISHT(2)
      DIMENSION U01(3),V01(3)
      DIMENSION NNDM(100),INDM(5,100),ATT(10)
      CHARACTER*12 CHTEXT
      SAVE IBOX,ITRD1,ITRD2,ITRAP,ITUBE,ITUBS,ICON,ICONS,ISPH,IPARA
      SAVE IPGON,IPCON,IGTRA,ICTUB,IHYPE,IELTU
      SAVE NNDM,INDM,U01,V01,THE,PHI,XMAN1,YMAN1
C
      DATA IBOX /'DX  ','DY  ','DZ  ','    ','    '/
      DATA ITRD1/'DX1 ','DX2 ','DY  ','DZ  ','    '/
      DATA ITRD2/'DX1 ','DX2 ','DY1 ','DY2 ','DZ  ','    '/
      DATA ITRAP/'DZ  ','THET','PHI ','H1  ','BL1 ','TL1 ','ALP1',
     +'H2  ','BL2 ','TL2 ','ALP2'/
      DATA ITUBE/'RMIN','RMAX','DZ  ','    ','    '/
      DATA ITUBS/'RMIN','RMAX','DZ  ','PHI1','PHI2'/
      DATA ICON /'DZ  ','RMN1','RMX1','RMN2','RMX2'/
      DATA ICONS/'DZ  ','RMN1','RMX1','RMN2','RMX2','PHI1','PHI2'/
      DATA ISPH /'RMIN','RMAX','THE1','THE2','PHI1','PHI2'/
      DATA IPARA/'DX  ','DY  ','DZ  ','ALPH','THET','PHI '/
      DATA IPGON/'PHI1','DPHI','NPDV','NZ  ','Z   ','RMIN','RMAX'/
      DATA IPCON/'PHI1','DPHI','NZ  ','Z   ','RMIN','RMAX'/
      DATA IHYPE/'RMIN','RMAX','DZ  ','TWST'/
      DATA IGTRA/'DZ  ','THET','PHI ','TWIS','H1  ','BL1 ','TL1 ',
     +'ALP1','H2  ','BL2 ','TL2 ','ALP2'/
      DATA ICTUB/'RMIN','RMAX','DZ  ','PHI1','PHI2','LXL ','LYL ',
     +           'LZL ','LXH ','LYH ','LZH '/
      DATA IELTU /'A   ','B   ','DZ  '/
C
      DATA NNDM/0,0,0,4,0,2,0,2,4,3,4,3,0,0,13*0,5,72*0/
      DATA INDM/0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 2,3,7,11,0, 0,0,0,0,0,
     +          4,5,0,0,0, 0,0,0,0,0, 6,7,0,0,0, 3,4,5,6,0, 4,5,6,0,0,
     +          1,2,3,4,0, 1,2,3,0,0, 0,0,0,0,0, 70*0, 2,3,4,8,12,
     +          360*0/
C
      DATA U01/14.5,5.5,14.5/
      DATA V01/14.,5.,5./
      DATA THE/45.,0.,90./
      DATA PHI/135.,0.,180./
      DATA XMAN1/8.8/
      DATA YMAN1/11.6/
C.
C.    ------------------------------------------------------------------
C.
C             Is NAME an existing volume ?
C
***      CALL IGRNG(20.,20.)
      CALL HPLFRA(0.,20.,0.,20.,'AB')
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF (IVO.LE.0) GO TO 999
C
C             Normalize to PLTRNX,PLTRNY
C
      DO 10 I=1,3
         U0(I)=U01(I)*PLTRNX/20.
         V0(I)=V01(I)*PLTRNY/20.
   10 CONTINUE
      XMAN=XMAN1*PLTRNX/20.
      YMAN=YMAN1*PLTRNY/20.
C
C             Save GDRAW calling parameters
C             and ZOOM internal parameters
C
      SAVTHE=GTHETA
      SAVPHI=GPHI
      SAVPSI=GPSI
      SAVU0=GU0
      SAVV0=GV0
      SAVSCU=GSCU
      SAVSCV=GSCV
      SVGZUA=GZUA
      SVGZVA=GZVA
      SVGZUB=GZUB
      SVGZVB=GZVB
      SVGZUC=GZUC
      SVGZVC=GZVC
      GZUA=1
      GZVA=1
      GZUB=0
      GZVB=0
      GZUC=0
      GZVC=0
C
C             Get shape type
C
      JVO=LQ(JVOLUM-IVO)
      ISHAPE=Q(JVO+2)
C
C             Get user parameters
C
      CALL GFPARA(NAME,1,0,NPAR,NATT,PAR,ATT)
      IF(NPAR.LE.0) GO TO 250
C
C
C             Check parameter sizes
C
      PARMAX=-1.
      DO 40 I=1,NPAR
         IF(NNDM(ISHAPE).LE.0) GO TO 30
         NDM=NNDM(ISHAPE)
         DO 20 IDM=1,NDM
            IF(I.EQ.INDM(IDM,ISHAPE)) GO TO 40
   20    CONTINUE
   30    ABSPAR=ABS(PAR(I))
         PARMAX=MAX(PARMAX,ABSPAR)
   40 CONTINUE
C
      GSCU=MIN(PLTRNX,PLTRNY)/(7.*PARMAX)
      GSCV=GSCU
      AXSIZ=PARMAX*0.35
C
C             Draw header
C
      CALL GDHEAD(-1,NAME,0.)
C
C             Draw parameters list
C
      DY=0.4
      IF(NPAR.GT.20) NPAR=20
      IF(NPAR.GT.10) DY=5.0/NPAR
      H=DY*0.7
C
      DO 210 I=1,NPAR
         CALL UCTOH(' = <',ISPAR(2),4,4)
         CALL UCTOH('CM $',ISPAR(3),4,4)
         IF (ISHAPE.NE.1) GO TO 50
         CALL UCTOH('B<OX',ISHT(1),4,4)
         CALL UCTOH('   $',ISHT(2),4,4)
         CALL UCTOH(IBOX(I),ISPAR(1),4,4)
         GO TO 200
   50    IF (ISHAPE.NE.2) GO TO 60
         CALL UCTOH('T<RD',ISHT(1),4,4)
         CALL UCTOH('>1 $',ISHT(2),4,4)
         CALL UCTOH(ITRD1(I),ISPAR(1),4,4)
         GO TO 200
   60    IF(ISHAPE.NE.3) GO TO 70
         CALL UCTOH('T<RD',ISHT(1),4,4)
         CALL UCTOH('>2 $',ISHT(2),4,4)
         CALL UCTOH(ITRD2(I),ISPAR(1),4,4)
         GO TO 200
   70    IF(ISHAPE.NE.4) GO TO 80
         CALL UCTOH('T<RA',ISHT(1),4,4)
         CALL UCTOH('P  $',ISHT(2),4,4)
         CALL UCTOH(ITRAP(I),ISPAR(1),4,4)
         IF(I.EQ.2.OR.I.EQ.3.OR.I.EQ.7.OR.I.EQ.11)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
   80    IF (ISHAPE.NE.5) GO TO 90
         CALL UCTOH('T<UB',ISHT(1),4,4)
         CALL UCTOH('E  $',ISHT(2),4,4)
         CALL UCTOH(ITUBE(I),ISPAR(1),4,4)
         GO TO 200
   90    IF (ISHAPE.NE.6) GO TO 100
         CALL UCTOH('T<UB',ISHT(1),4,4)
         CALL UCTOH('S  $',ISHT(2),4,4)
         CALL UCTOH(ITUBS(I),ISPAR(1),4,4)
         IF(I.GT.3)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
  100    IF(ISHAPE.NE.7) GO TO 110
         CALL UCTOH('C<ON',ISHT(1),4,4)
         CALL UCTOH('E  $',ISHT(2),4,4)
         CALL UCTOH(ICON(I),ISPAR(1),4,4)
         GO TO 200
  110    IF(ISHAPE.NE.8) GO TO 120
         CALL UCTOH('C<ON',ISHT(1),4,4)
         CALL UCTOH('S  $',ISHT(2),4,4)
         CALL UCTOH(ICONS(I),ISPAR(1),4,4)
         IF(I.GT.5)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
  120    IF(ISHAPE.NE.9) GO TO 130
         CALL UCTOH('S<PH',ISHT(1),4,4)
         CALL UCTOH('E  $',ISHT(2),4,4)
         CALL UCTOH(ISPH(I),ISPAR(1),4,4)
         IF(I.GT.2)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
  130    IF(ISHAPE.NE.10) GO TO 140
         CALL UCTOH('P<AR',ISHT(1),4,4)
         CALL UCTOH('A  $',ISHT(2),4,4)
         CALL UCTOH(IPARA(I),ISPAR(1),4,4)
         IF(I.GT.3)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
  140    IF(ISHAPE.NE.11.AND.ISHAPE.NE.12) GO TO 170
         CALL UCTOH('P<GO',ISHT(1),4,4)
         CALL UCTOH('N  $',ISHT(2),4,4)
         IU=I
         IF(IU.LT.8) GO TO 150
         I2=I-5
         I3=I2/3
         I4=I2-I3*3
         IU=I4+5
  150    CALL UCTOH(IPGON(IU),ISPAR(1),4,4)
         IF(I.LT.3)THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         IF(I.EQ.3.OR.I.EQ.4)CALL UCTOH('   $', ISPAR(3),4,4)
         IF(ISHAPE.EQ.11) GO TO 200
         CALL UCTOH('P<CO',ISHT(1),4,4)
         IU=I
         IF(IU.LT.7) GO TO 160
         I2=I-4
         I3=I2/3
         I4=I2-I3*3
         IU=I4+4
  160    CALL UCTOH(IPCON(IU),ISPAR(1),4,4)
         IF(I.EQ.4) CALL UCTOH('CM $',ISPAR(3),4,4)
C
         GO TO 200
  170    CONTINUE
         IF(ISHAPE.NE.13) GO TO 180
         CALL UCTOH('E<LT',ISHT(1),4,4)
         CALL UCTOH('U  $',ISHT(2),4,4)
         CALL UCTOH(IELTU(I),ISPAR(1),4,4)
         GO TO 200
  180    CONTINUE
         IF(ISHAPE .NE. 14) GO TO 190
         CALL UCTOH('H<YP',ISHT(1),4,4)
         CALL UCTOH('E  $',ISHT(2),4,4)
         CALL UCTOH(IHYPE(I),ISPAR(1),4,4)
         IF(I .GT. 3) THEN
            CALL UCTOH(' = <',ISPAR(2),4,4)
            CALL UCTOH('DEG$',ISPAR(3),4,4)
         ENDIF
         GO TO 200
C
  190    CONTINUE
C
C             Cut tube
C
         IF (ISHAPE.EQ.NSCTUB) THEN
            CALL UCTOH('CT<U',ISHT(1),4,4)
            CALL UCTOH('B  $',ISHT(2),4,4)
            CALL UCTOH(ICTUB(I),ISPAR(1),4,4)
            IF(I.GE.4.AND.I.LE.5)THEN
               CALL UCTOH(' = <',ISPAR(2),4,4)
               CALL UCTOH('DEG$',ISPAR(3),4,4)
            ELSE IF(I.GE.6)THEN
               CALL UCTOH(' =  ',ISPAR(2),4,4)
               CALL UCTOH('  $ ',ISPAR(3),4,4)
            END IF
            GO TO 200
         ENDIF
C
C             General twisted trapezoid.
C
         IF(ISHAPE.NE.28) GO TO 230
         CALL UCTOH('G>TR',ISHT(1),4,4)
         CALL UCTOH('A  $',ISHT(2),4,4)
         CALL UCTOH(IGTRA(I),ISPAR(1),4,4)
         IF(I.EQ.2.OR.I.EQ.3.OR.I.EQ.4.OR.I.EQ.8.OR. I.EQ.12) CALL
     +   UCTOH('DEG$',ISPAR(3),4,4)
C
  200    CONTINUE
         IF (I.EQ.1) THEN
            XTEXT=4.*PLTRNX/20.
            YTEXT=16.5*PLTRNY/20.
            CSIZE=DY*MIN(PLTRNX,PLTRNY)/20.
            CALL UHTOC(ISHT,4,CHTEXT,12)
            CALL GDRAWT(XTEXT,YTEXT,CHTEXT,CSIZE,0.,1,-1)
         ENDIF
         Y=16.5-(I+0.5)*DY
         XTEXT=3.*PLTRNX/20.
         YTEXT=Y*PLTRNY/20.
         CSIZE=H*MIN(PLTRNX,PLTRNY)/20.
         CALL UHTOC(ISPAR,4,CHTEXT,12)
         CALL GDRAWT(XTEXT,YTEXT,CHTEXT,CSIZE,0.,1,-1)
         CALL HBCDF(PAR(I),8,IPAR)
         IF (PAR(I).EQ.0.) CALL UCTOH('0',IPAR(1),1,1)
         CALL UCTOH('$',IPAR(9),1,1)
         CALL UBUNCH(IPAR,IPA,9)
         XTEXT=(H*10.+3.)*PLTRNX/20.
         YTEXT=Y*PLTRNY/20.
         CSIZE=H*MIN(PLTRNX,PLTRNY)/20.
         CALL UHTOC(IPA,4,CHTEXT,12)
         CALL GDRAWT(XTEXT,YTEXT,CHTEXT,CSIZE,0.,1,-1)
  210 CONTINUE
C
C             Draw views
C
      CALL GFATT(NAME,'SEEN',KSEEN)
C
C             Add local value SEEN 1 to starting node of tree
C
      KSEEN=KSEEN+110
      ISEEN=KSEEN
      CALL GSATT(NAME,'SEEN',ISEEN)
C
      CALL GSATT(NAME,'COLO',2)
      CALL GDNSON(NAME,NSON,IDIV)
      DO 220 N=1,NSON
         CALL GDSON(N,NAME,ISON)
         CALL GFATT(ISON,'SEEN',KSEEN)
C
C             ISON is a volume with multiplicity;
C             first occurrence has already been set
C
         IF (KSEEN.GT.50) GO TO 220
C
C             Add local value SEEN -2 to each one-level-down node
C
         KSEEN=KSEEN+80
         ISEEN=KSEEN
         CALL GSATT(ISON,'SEEN',ISEEN)
C
         CALL GSATT(ISON,'COLO',4)
  220 CONTINUE
C
      CALL GDCOL(3)
      XSCAL=PLTRNX/4.
      YSCAL=PLTRNY/2.
      CALL GDSCAL(XSCAL,YSCAL)
**    IF (GSCU.LE.0.05) CALL GDMAN(XMAN,YMAN)
      IF (GSCU.LE.0.05) CALL GDWMN1(XMAN,YMAN)
C
      CALL GDRAW(NAME,THE(1),PHI(1),0.,U0(1),V0(1),GSCU,GSCV)
      CALL GDAXIS(0.,0.,0.,AXSIZ)
      CALL GDRAWC(NAME,3,0.005,U0(2),V0(2),GSCU,GSCV)
      CALL GDAXIS(0.,0.,0.,AXSIZ)
      CALL GDRAWC(NAME,1,0.005,U0(3),V0(3),GSCU,GSCV)
      CALL GDAXIS(0.,0.,0.,AXSIZ)
C
  230 CALL GDCOL(0)
C
C             Reset global SEEN values
C
      DO 240 IVO=1,NVOLUM
         CALL UHTOC(IQ(JVOLUM+IVO),4,NAMSEE,4)
         CALL GFATT(NAMSEE,'SEEN',KSEEN)
         IF (KSEEN.LT.50) GO TO 240
         ISEENL=KSEEN/10.+0.5
         ISEENG=KSEEN-ISEENL*10.
         CALL GSATT(NAMSEE,'SEEN',ISEENG)
  240 CONTINUE
C
  250 CONTINUE
C
C             Restore GDRAW calling parameters
C             and ZOOM internal parameters
C
      GTHETA=SAVTHE
      GPHI=SAVPHI
      GPSI=SAVPSI
      GU0=SAVU0
      GV0=SAVV0
      GSCU=SAVSCU
      GSCV=SAVSCV
      NGVIEW=0
      GZUA=SVGZUA
      GZVA=SVGZVA
      GZUB=SVGZUB
      GZVB=SVGZVB
      GZUC=SVGZUC
      GZVC=SVGZVC
      CALL ISELNT(1)
  999 END
