*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:04  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GSPE3D(NAME,TETA3,PHI3,PSI3,U03,V03,ZM3)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws specifications of volume NAME         *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDFSPC                      *
C.    *                                                                *
C.    *       Author :    S.Giani 29-10-92    ************             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcshno.inc"
      COMMON/SP3D/ISPFLA
*
      CHARACTER*4 NAME,NAMSEE,ISON
      DIMENSION PAR(100)
      DIMENSION U0(3),V0(3),THE(3),PHI(3)
      DIMENSION U01(3),V01(3)
      DIMENSION NNDM(100),INDM(5,100),ATT(10)
      SAVE NNDM,INDM,U01,V01,THE,PHI,XMAN1,YMAN1
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
      DATA XMAN1/3./
      DATA YMAN1/3./
C.
C.    ------------------------------------------------------------------
C.
      CALL ICLRWK(1,1)
      ISPFLA=0
C             Is NAME an existing volume ?
C
      CALL IGRNG(20.,20.)
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
      GTHETA=TETA3
      GPHI=PHI3
      GPSI=PSI3
*
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
      IF(NPAR.LE.0) GO TO 70
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
      AXSIZ=PARMAX*.5
      PARMAX=PARMAX/ZM3
      GSCU=MIN(PLTRNX,PLTRNY)/(7.*PARMAX)
      GSCV=GSCU
C
C             Draw header
C
      CALL GDHEAD(-1,NAME,0.)
C
C             Draw parameters list
C
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
      DO 50  N=1,NSON
         CALL GDSON(N,NAME,ISON)
         CALL GFATT(ISON,'SEEN',KSEEN)
C
C             ISON is a volume with multiplicity;
C             first occurrence has already been set
C
         IF (KSEEN.GT.50) GO TO 50
C
C             Add local value SEEN -2 to each one-level-down node
C
         KSEEN=KSEEN+80
         ISEEN=KSEEN
         CALL GSATT(ISON,'SEEN',ISEEN)
C
         CALL GSATT(ISON,'COLO',4)
   50 CONTINUE
C
      ISPFLA=1
      CALL GDRAW(NAME,TETA3,PHI3,PSI3,U03,V03,GSCU,GSCV)
      CALL GDCOL(3)
*      XSCAL=PLTRNX/4.
*      YSCAL=PLTRNY/2.
      XSCAL=17.
      YSCAL=3.
      CALL GDSCAL(XSCAL,YSCAL)
**    IF (GSCU.LE.0.05) CALL GDMAN(XMAN,YMAN)
      IF (GSCU.LE.0.05) CALL GDWMN1(XMAN,YMAN)
      CALL GDAXIS(0.,0.,0.,AXSIZ)
      ISPFLA=0
C
      CALL GPVOLX(IVO)
*
      CALL GDCOL(0)
C
C             Reset global SEEN values
C
      DO 60  IVO=1,NVOLUM
         CALL UHTOC(IQ(JVOLUM+IVO),4,NAMSEE,4)
         CALL GFATT(NAMSEE,'SEEN',KSEEN)
         IF (KSEEN.LT.50) GO TO 60
         ISEENL=KSEEN/10.+0.5
         ISEENG=KSEEN-ISEENL*10.
         CALL GSATT(NAMSEE,'SEEN',ISEENG)
   60 CONTINUE
C
   70 CONTINUE
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
  999 END
