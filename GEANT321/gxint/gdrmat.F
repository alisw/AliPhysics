*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 10/10/94  15.45.50  by  S.Ravndal
*-- Author :
      SUBROUTINE GDRMAT(IMATES,IPART,CHMECA,MMEC)
C.    ******************************************************************
C.    *                                                                *
C.    *       DRAW cross sections and energy loss tables               *
C.    *       material IMATE, particle IPART, mecanism name CHMECA,    *
C.    *                                                                *
C.    *      The CHMECAnism name can be :                              *
C.    *      'HADF'   'INEF'   'ELAF'   'FISF'   'CAPF'                *
C.    *      'HADG'   'INEG'   'ELAG'   'FISG'   'CAPG'                *
C.    *      'LOSS'   'PHOT'   'ANNI'   'COMP'   'BREM'                *
C.    *      'PAIR'   'DRAY'   'PFIS'   'RAYL'   'HADG'                *
C.    *      'MUNU'   'RANG'   'STEP'                                  *
C.    *                                                                *
C.    *       For Hadronic particles it also computes the              *
C.    *       hadronic cross section from FLUKA ( '***F' ) or          *
C.    *       GHEISHA ( '***G' ) programs:                             *
C.    *       HADF or HADG -- total                                    *
C.    *       INEF or INEG -- inelastic                                *
C.    *       ELAF or ELAG -- elastic                                  *
C.    *       FISF or FISG -- fission (0.0 for FLUKA)                  *
C.    *       CAPF or CAPG -- neutron capture (0.0 for FLUKA)          *
C.    *                                                                *
C.    *             Input parameters                                   *
C.    *  IMATE   Geant material number                                 *
C.    *  IPART   Geant particle number                                 *
C.    *  CHMECA   mechanism name of the bank to be fetched             *
C.    *                                                                *
C.    *    ==>Called by : <USER>  GXCONT                               *
C.    *       Authors    R.Brun, M.Maire    *********                  *
C.    *                                                                *
C.    ******************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
*
      CHARACTER*4 CHMECA(MMEC)
      CHARACTER*20 MESSGE,NAMEM,NAMEP
      CHARACTER*64 TITLE
      LOGICAL HEXIST,LXBARN
#include "geant321/gcnmec.inc"
      DIMENSION CRY(201),IYES(NMECA+2),XC(2)
*
      IF(IMATES.LT.0) THEN
         LXBARN=.TRUE.
      ELSE
         LXBARN=.FALSE.
      ENDIF
      IMATE=ABS(IMATES)
      IDM=-1
      CALL IGSA(0)
      DO 10 JMEC=1,MMEC
         CALL GPLMAT(IMATES,IPART,CHMECA(JMEC),NEKBIN,ELOW,IDM)
   10 CONTINUE
      YMIN=BIG
      YMAX=-BIG
      IOK=0
      IBASE = 10000*IMATE+100*IPART
      DO 30 IDB=1,NMECA+2
         IYES(IDB)=0
         ID=IBASE+IDB
         IF(HEXIST(ID)) THEN
            CALL HNOENT(ID,NOENT)
            IF(NOENT.LE.0)THEN
               CALL HDELET(ID)
            ELSE
               IYES(IDB)=1
               IOK=1
               DO 20 I=1,NEKBIN
                  CRS=HI(ID,I)
                  IF(CRS.LE.0.)GO TO 20
                  IF(CRS.GT.1.E10)GO TO 20
                  IF(CRS.LT.YMIN)YMIN=CRS
                  IF(CRS.GT.YMAX)YMAX=CRS
   20          CONTINUE
            ENDIF
         ENDIF
   30 CONTINUE
      IF(IOK.EQ.0)GO TO 999
*
      CALL HPLOPT('LOGX',1)
      CALL HPLOPT('LOGY',1)
      CALL IGSET('MSCF',0.8)
      CALL IGSET('TXAL',10.)
      YLI=LOG10(YMIN)
      YLA=LOG10(YMAX)
      DELTA=0.05*(YLA-YLI)
      YMIN=10.**(YLI-DELTA)
      YMAX=10.**(YLA+DELTA)
      CALL HPLFRA(ELOW(1),ELOW(NEKBIN+1),YMIN,YMAX,' ')
      XC(1)=ELOW(8)
      XC(2)=ELOW(NEKBIN/2)
      JMA=LQ(JMATE-IMATE)
      CALL UHTOC(IQ(JMA+1),4,NAMEM,20)
      NCHM=LNBLNK(NAMEM)
      IF(NAMEM(NCHM:NCHM).EQ.'$')THEN
         NCHM=NCHM-1
      ENDIF
      JPA=LQ(JPART-IPART)
      CALL UHTOC(IQ(JPA+1),4,NAMEP,20)
      NCHP=LNBLNK(NAMEP)
      IF(NAMEP(NCHP:NCHP).EQ.'$')THEN
         NCHP=NCHP-1
      ENDIF
      TITLE='Tables for '//NAMEP(1:NCHP)//' in '//NAMEM(1:NCHM)
      XCT=XC(1)
      CALL HPLTOC(XCT,YMAX,XCTT,YCTT,NT)
      YCTT=YCTT+0.3
      CALL ISELNT(1)
      CALL IGSET('CHHE',0.28)
      CALL ITX(XCTT,YCTT,TITLE)
      CALL ISELNT(NT)
*
      CALL IGSA(0)
      NSHIF=0
      NPLOT=0
      HSYM = 20
      HCOL = 0
      DO 50 IDB=1,NMECA+2
         IF(IYES(IDB).NE.0)THEN
            IF(HSYM.EQ.31) HSYM = 20
            HSYM = HSYM+1
            IF(HCOL.EQ.7) HCOL = 0
            HCOL = HCOL+1
            ID=10000*IMATE+100*IPART+IDB
            CALL HUNPAK(ID,CRY,'HIST',1)
            CALL HDELET(ID)
            CALL IGSET('MTYP',HSYM)
            CALL IGSET('PMCI',HCOL)
            CALL IGSET('PLCI',HCOL)
            KMIN = 0
            KMAX = 0
            DO 40 IBIN=1,NEKBIN
               IF(CRY(IBIN).GT.0.) THEN
                  IF(KMIN.EQ.0) THEN
                     KMIN = IBIN
                  ENDIF
                  KMAX = IBIN
               ENDIF
   40       CONTINUE
            IF(KMIN.NE.0) THEN
               CALL IGRAPH(KMAX-KMIN+1,ELOW(KMIN),CRY(KMIN),'PLGXY')
            ENDIF
            NPLOT = NPLOT+1
            IND=2-MOD(NPLOT,2)
            IF(IND.EQ.1) THEN
               NSHIF = NSHIF+1
            ENDIF
            YC=10.**(YLA-(NSHIF-1)*0.8*DELTA)
            CALL HPLTOC(XC(IND),YC,XCM,YCM,NT)
            CALL ISELNT(1)
            CALL IPM(1,XCM,YCM)
            IF(IDB.LE.NMECA) THEN
               IF (CHNMEC(IDB).EQ.'LOSS') THEN
                  MESSGE=CHNMEC(IDB)//' (MeV/cm)'
               ELSEIF (CHNMEC(IDB).EQ.'RANG') THEN
                  MESSGE=CHNMEC(IDB)//' (cm)'
               ELSEIF (CHNMEC(IDB).EQ.'STEP') THEN
                  MESSGE=CHNMEC(IDB)//' (cm)'
               ELSE
                  MESSGE=CHNMEC(IDB)//' X-sec'
                  IF(LXBARN) THEN
                    MESSGE(LNBLNK(MESSGE)+1:) = ' (barn)'
                  ELSE
                    MESSGE(LNBLNK(MESSGE)+1:) = ' (1/cm)'
                  ENDIF
               ENDIF
            ELSEIF(IDB.EQ.NMECA+1) THEN
               MESSGE='Tot X-sec'
                  IF(LXBARN) THEN
                    MESSGE(LNBLNK(MESSGE)+1:) = ' (barn)'
                  ELSE
                    MESSGE(LNBLNK(MESSGE)+1:) = ' (1/cm)'
                  ENDIF
            ELSEIF(IDB.EQ.NMECA+2) THEN
               MESSGE='Mean free path (cm)'
            ENDIF
            CALL IGSET('CHHE',0.28)
            CALL ITX(XCM+0.5,YCM,MESSGE)
            CALL ISELNT(NT)
         ENDIF
   50 CONTINUE
      CALL HPLOPT('LINX',1)
      CALL HPLOPT('LINY',1)
      CALL IGSET('PMCI',1.)
      CALL IGSET('PLCI',1.)
  999 END
