*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:29:30  mclareni
* Import version 8.04
*
* Revision 1.5  2000/04/15 14:10:47  plothow
* Version 8.03: Nuclear PDFs introduced
*
* Revision 1.4  2000/04/09 14:47:31  plothow
* version 8.02 - add CTEQ4 ref2 PDFs
*
* Revision 1.3  2000/04/09 10:06:54  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
*
* Revision 1.2  1996/10/28 15:13:59  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:06  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STRUCTM(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C
C  *********************************************************************
C  *                                                                   *
C  *   Main steering routine for all sets of structure functions       *
C  *                                                                   *
C  *                                                                   *
C  *   Input:    X     = x value of parton                             *
C  *             SCALE = QCD scale in GeV                              *
C  *                                                                   *
C  *   Output:   UPV   = up valence quark                              *
C  *             DNV   = down valence quark                            *
C  *             USEA  = sea (up_bar)                                  *
C  *             DSEA  = sea (down_bar)                                *
C  *             STR   = strange quark                                 *
C  *             CHM   = charm quark                                   *
C  *             BOT   = bottom quark                                  *
C  *             TOP   = top quark                                     *
C  *             GL    = gluon                                         *
C  *                                                                   *
C  *                                                                   *
C  *   The variables  NPTYPE, NGROUP and NSET  should be,              *
C  *   the variables  NFL, LO, TMAS  and                               *
C  *   QCDL4, QCDL5, XMIN, XMAX, Q2MIN, Q2MAX                          *
C  *   could be provided by the user via a call to the                 *
C  *   subroutine PDFSET at the initialization phase, where            *
C  *                                                                   *
C  *             NPTYPE = Particle type                                *
C  *                      (number or character string,                 *
C  *                       1,2,3 or 'NU','PI','PH')                    *
C  *                      of desired structure functions set           *
C  *                      (Default: NPTYPE = 1  or  'NU'               *
C  *             NGROUP = author group                                 *
C  *                      (number or character string,                 *
C  *                       1 to 7 or i.e. 'DFLM','MRS','MT','GRV',etc) *
C  *                      of desired structure functions set           *
C  *                      (Default: NGROUP = 5  or  'GRV')             *
C  *             NSET   = number of desired structure functions set    *
C  *                      (Default: NSET = 3)                          *
C  *             NFL    = desired number of flavours for alpha(s)      *
C  *                      (Default: NFL = 5)                           *
C  *             LO     = order of alpha(s) calculation                *
C  *                      (Default: LO = 2)                            *
C  *             TMAS   = top quark mass in GeV/c**2  (optional)       *
C  *                      (Default: TMAS = 100.0D0)                    *
C  *             QCDL4  = QCD scale in GeV for four flavours           *
C  *             QCDL5  = QCD scale in GeV for five flavours           *
C  *                      corresponding to QCDL4                       *
C  *             XMIN   = minimal allowed x value                      *
C  *             XMAX   = maximal allowed x value                      *
C  *             Q2MIN  = minimal allowed Q**2 value                   *
C  *             Q2MAX  = maximal allowed Q**2 value                   *
C  *                                                                   *
C  *   for each set of structure fuction.                              *
C  *                                                                   *
C  *   The internal COMMON blocks                                      *
C  *                                                                   *
C  *   COMMON/W50511/ NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS              *
C  *   COMMON/W50512/ QCDL4,QCDL5                                      *
C  *   COMMON/W50513/ XMIN,XMAX,Q2MIN,Q2MAX                            *
C  *                                                                   *
C  *   are filled then by the subroutine PDFSET.                       *
C  *                                                                   *
C  *                                                                   *
C  *     Note: STRUCTM returns X * parton distribution function !      *
C  *                                                                   *
C  *                                                                   *
C  *     Author:   H. Plothow-Besch                                    *
C  *               CERN-PPE, CH - 1211 Geneva 23, Switzerland          *
C  *                                                                   *
C  *     Please return any problems, questions, suggestions            *
C  *     to the author                                                 *
C  *                                                                   *
C  *********************************************************************
C
C
#include "pdf/impdp.inc"
C
#include "pdf/w5051p1.inc"
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50514.inc"
#include "pdf/w50514w.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w50519.inc"
#include "pdf/w505120.inc"
#include "pdf/w505121.inc"
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
      DATA ZEROD/0.D0/,ONED/1.D0/,TWOD/2.D0/
      SAVE /W50514/, /W50514W/, /W50516/
C.
#include "pdf/w50511c.inc"
C.
C User wants new version (4.0 or bigger) of PDFLIB format
      IF (IFLSET.NE.1) THEN
         IF(FIRST) THEN
         WRITE(N6,*) ' Warning : NO initialisation via PDFSET made !!'
         WRITE(N6,*)
     +   '           ALL Parameters set to default (Nucleon PDFs) !!'
         ENDIF
         NPTYPE = LPTYPE
         NGROUP = LGROUP
         NSET   = LNSET
C
         PARM(1) = 'Nptype'
         VAL(1)  = NPTYPE
         PARM(2) = 'Ngroup'
         VAL(2)  = NGROUP
         PARM(3) = 'Nset'
         VAL(3)  = NSET
C
         CALL PDFSET(PARM,VAL)
      ENDIF
C.
C User wants old version (3.0 or less) of PDFLIB format
      IF(.NOT.NEWVER) THEN
        IF (MODE .GE.0 .AND. MODE .LE.MODEMX) THEN
           IF(MODE.EQ.0) THEN
             NPTYPE = LPTYPO
             NGROUP = LGROPO
             NSET   = LNSETO
           ELSE
             NPTYPE = NPTYCR(MODE)
             NGROUP = NGROCR(MODE)
             NSET   = NSETCR(MODE)
C   Check on validity of parameter values
             IF (NPTYPE.LT.0 .OR. NGROUP.LT.0 .OR. NSET.LT.0) THEN
                IF(FIRST) THEN
                WRITE(N6,*) ' PDFLIB  : MODE value INCORRECT, MODE = ',
     +          MODE
               WRITE(N6,*) ' Warning : MODE value set to OLD default !!'
                ENDIF
                NPTYPE = LPTYPO
                NGROUP = LGROPO
                NSET   = LNSETO
             ENDIF
           ENDIF
        ENDIF
      ENDIF
C.
C... Define printer
           IF(N6.LE.0) N6 = L6
C
        DUPV = ZEROD
        DDNV = ZEROD
        DUSEA = ZEROD
        DDSEA = ZEROD
        DSTR = ZEROD
        DCHM = ZEROD
        DBOT = ZEROD
        DTOP = ZEROD
        DGL  = ZEROD
C
        ZUPV = ZEROD
        ZDNV = ZEROD
        ZUSEA = ZEROD
        ZDSEA = ZEROD
        ZSTR = ZEROD
        ZCHM = ZEROD
        ZBOT = ZEROD
        ZTOP = ZEROD
        ZGL  = ZEROD
C
        DX  = X
        DQ  = SCALE
        DQ2 = SCALE*SCALE
C.
        IF (DX.LE.ZEROD .OR. DX.GE.ONED) THEN
           WRITE(N6,*) ' PDFLIB  : Nptype = ',NPTYPE,',   Ngroup = ',
     +                 NGROUP,',   Nset = ',NSET
           WRITE(N6,*)
     +     ' Error   : X value outside physical range , X = ',DX
           GOTO 100
        ENDIF
        WTXMIN = WTXMIN + PDFWGT
        IF (DX.LT.XMIN) THEN
           WXMIN = WXMIN + PDFWGT
           IF(IFLPRT.GE.3) THEN
              WRITE(N6,*) ' PDFLIB  : Nptype = ',NPTYPE,',   Ngroup = ',
     +                 NGROUP,',   Nset = ',NSET
              WRITE(N6,*)
     +        ' Error   : X value smaller allowed range , X = ',DX
           ENDIF
        ENDIF
        WTXMAX = WTXMAX + PDFWGT
        IF (DX.GT.XMAX) THEN
           WXMAX = WXMAX + PDFWGT
           IF(IFLPRT.GE.3) THEN
              WRITE(N6,*) ' PDFLIB  : Nptype = ',NPTYPE,',   Ngroup = ',
     +                 NGROUP,',   Nset = ',NSET
              WRITE(N6,*)
     +        ' Error   : X value bigger  allowed range , X = ',DX
           ENDIF
           GOTO 100
        ENDIF
        WTQ2MIN = WTQ2MIN + PDFWGT
        IF (DQ2.LT.Q2MIN) THEN
           WQ2MIN = WQ2MIN + PDFWGT
           IF(IFLPRT.GE.3) THEN
              WRITE(N6,*) ' PDFLIB  : Nptype = ',NPTYPE,',   Ngroup = ',
     +                 NGROUP,',   Nset = ',NSET
              WRITE(N6,*)
     +      ' Error   : SCALE value smaller allowed range , Q**2 = ',DQ2
           ENDIF
        ENDIF
        WTQ2MAX = WTQ2MAX + PDFWGT
        IF (DQ2.GT.Q2MAX) THEN
           WQ2MAX = WQ2MAX + PDFWGT
           IF(IFLPRT.GE.3) THEN
              WRITE(N6,*) ' PDFLIB  : Nptype = ',NPTYPE,',   Ngroup = ',
     +                 NGROUP,',   Nset = ',NSET
              WRITE(N6,*)
     +       ' Error   : SCALE value bigger allowed range , Q**2 = ',DQ2
           ENDIF
        ENDIF
C
      IF (FIRST) WRITE(N6,*) ' '
      IF (FIRST) WRITE(N6,*) ' '
C Call default structure functions set
      IF    (NPTYPE.EQ.0 .OR. NGROUP.EQ.0 .OR. NSET.EQ.0) THEN
C       IF (FIRST) WRITE(N6,*) 
C     + ' Nucleon PDFs :  GRV Set HO Structure Functions'
C       IF (FIRST) WRITE(N6,*) 
C     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
C       CALL GRVHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C       DDSEA = DUSEA
C       GOTO 100
C          IF (FIRST) WRITE(N6,*) 
C     + ' Nucleon PDFs :  MRS Set (G) (L255-MSb) Structure Functions'
C          IF (FIRST) WRITE(N6,*) 
C     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
C       CALL SFMRSG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
C          DTOP=ZEROD
C          GOTO 100
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  MRS central-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFMRS_99A1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
      ENDIF
C Start with NUCLEON structure functions
      IF(NPTYPE.EQ.1) THEN
      IF(NGROUP.EQ.1) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  PRIVATE Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL STRPRIV(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  BEBC Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRBEBC(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  OR Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRUCOR(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  BEP Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRBEP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GHR Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRGHR(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  DO Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRDO1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  DO Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRDO2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  EHLQ Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFEHLQ1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*)
     + ' Nucleon PDFs :  EHLQ Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFEHLQ2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  New DO Set 1.1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL NEWDO1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.11) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',   Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.2) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMsoft (Valence) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM4(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMhard (Valence) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM5(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMsoft (Gluon) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM6(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMhard (Gluon) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM7(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMaverage (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM8(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLMaverage (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM9(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLM160 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLM260 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  DFLM360 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DFLM3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.GE.10) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',   Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.3) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set 1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRMRS1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set 2 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRMRS2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set 3 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRMRS3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set E Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRUCE(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set B Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRUCB(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set EP Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRCEP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set BP Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRCBP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS1 Set E Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL STRH1E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
C       ELSEIF(NSET.EQ.9) THEN
C         IF (FIRST) WRITE(N6,*) 
C     + ' Nucleon PDFs :  HMRS1 Set B Structure Functions'
C          IF (FIRST) WRITE(N6,*) 
C     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
C         IF (FIRST) WRITE(N6,*)
C    +' Warning !!   :  in principle, set has been retracted by the auth
C    +ors !'
C         FIRST = .FALSE.
C         CALL STRH1B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
C         DDSEA = DUSEA
C         DTOP=ZEROD
C         GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  KMRS Set B0 (190) Structure Functions (old)'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB0O(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS2 Set E Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH2E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS2 Set B Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH2B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.12) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set E+ Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3EP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.13) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set E- Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3EM(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.14) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set E Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.15) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set B (190) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.16) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set B (100) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3B1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.17) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  HMRS Set B (300) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRH3B3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.18) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  KMRS Set B- Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKBM(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.19) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  KMRS Set B--R2-SH Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFKBMR2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.20) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  KMRS Set B--R5-SH Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFKBMR5(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.21) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  KMRS Set B0 (190) updated Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.22) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set B0 (L135) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB02(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.23) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set B0 (L160) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB03(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.24) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set B0 (L200) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB04(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.25) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set B0 (L235) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STRKB05(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.26) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set S0 (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.27) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D0 (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.28) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D- (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.29) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set S0p (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSBS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.30) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D0p (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSBD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.31) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D-p (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSBDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.32) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set S0p (L230-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFDISS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.33) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D0p (L230-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFDISD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.34) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set D-p (L230-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFDISDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.35) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (H) (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSBH(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.36) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (H) (L230-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFDISH(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.37) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (A) (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSBA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.38) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Fit (A) (L230-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSFA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.39) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L231-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSAP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.40) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Fit (Ap) (L231-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMFAP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.41) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (G) (L255-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.42) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Fit (G) (L255-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSFG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.43) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (A) (L230-MSb) and low Q**2 Structure F
     +unctions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF(DQ2.GT.5.D0)
     +    CALL SFMSBA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          IF(DQ2.LE.5.D0)
     +    CALL SFMSQA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.44) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (A) (L230-DIS) and low Q**2 Structure F
     +unctions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF(DQ2.GT.5.D0)
     +    CALL SFDISA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          IF(DQ2.LE.5.D0)
     +    CALL SFDSQA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.45) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L150-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS105(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.46) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L201-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS110(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.47) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L266-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS115(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.48) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L344-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS120(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.49) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L435-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS125(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.50) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Ap) (L542-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMS130(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.51) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (J) (L344-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSJ(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.52) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (Jp) (L507-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSJP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.53) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (R1) (L241-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSR1(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.54) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (R2) (L344-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSR2(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.55) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (R3) (L241-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSR3(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.56) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (R4) (L344-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMSR4(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.57) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Fit (R1) (L241-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMFR1(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.58) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Fit (R2) (L344-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMFR2(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.59) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRRS charm1 (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSC1(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.60) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRRS charm2 (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSC2(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.61) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRRS charm3 (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSC3(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.62) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST central-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98D1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.63) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98D2
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.64) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98D3
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.65) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-a_s (L229-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98D4
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.66) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-a_s (L383-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98D5
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.67) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST central-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98A1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.68) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98A2
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.69) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98A3
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.70) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-a_s (L229-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98A4
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.71) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-a_s (L383-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98A5
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.72) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST central-g (L174-lo) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98LO1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.73) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-g (L174-lo) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98LO2
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.74) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-g (L174-lo) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98LO3
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.75) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-a_s (L136-lo) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98LO4
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.76) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-a_s (L216-lo) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98LO5
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.77) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST ht (L292.5-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_98ht1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.78) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST central-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.79) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D2
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.80) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-g (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D3
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.81) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-a_s (L229-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D4
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.82) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-a_s (L383-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D5
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.83) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST q-up (L303.3-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D6
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.84) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST q-down (L290.3-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D7
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.85) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST s-up (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D8
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.86) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST s-down (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D9
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.87) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST c-up (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D10
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.88) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST c-down (L300-DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99D11
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.89) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST central-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A1
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.90) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A2
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.91) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-g (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A3
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.92) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST lower-a_s (L229-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A4
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.93) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST higher-a_s (L383-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A5
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.94) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST q-up (L303.3-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A6
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.95) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST q-down (L290.3-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A7
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.96) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST s-up (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A8
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.97) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST s-down (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A9
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.98) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST c-up (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A10
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.99) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST c-down (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A11
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.100) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRST larger-d/u (L300-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRS_99A12
     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.101) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.4) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set S1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set B1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set B2 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set E1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU4(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set 6 (1/2s) (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU5(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set S1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU6(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set B1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU7(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set B2 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU8(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set E1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU9(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set 6 (1/2s) (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMT10(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MT Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMT11(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.12) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 1L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ11(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.13) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 1M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ12(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.14) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 1MS (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ13(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.15) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 1ML (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ14(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.16) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 1D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ15(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.17) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ21(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.18) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ22(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.19) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2MS (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ23(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.20) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2MF (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ24(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.21) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2ML (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ25(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.22) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ26(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.23) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pL (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ21(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.24) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pM (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ22(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.25) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pMS (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ23(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.26) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pMF (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ24(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.27) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pML (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ25(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.28) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 2pD (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ26(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.29) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 3L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ31(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.30) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 3M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ32(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.31) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 3D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ33(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.32) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ41(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.33) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ42(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.34) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ43(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.35) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4A1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ44(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.36) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4A2 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ45(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.37) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4A3 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ43(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.38) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4A4 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ46(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.39) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4A5 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ47(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.40) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4HJ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ48(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.41) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4LQ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ49(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
        ELSEIF(NSET.EQ.42) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4HQ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ45(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
        ELSEIF(NSET.EQ.43) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4HQ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ46(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.44) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4F3  Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ47(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
        ELSEIF(NSET.EQ.45) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 4F4  Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SRCTQ48(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.46) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ53(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.47) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ52(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.48) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ51(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.49) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5HJ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ54(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.50) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5HQ (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ55(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.51) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5F3  Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ56(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.52) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5F4  Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ57(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.53) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5M1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ58(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.54) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5HQ1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFCTQ59(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.55) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5L (LO) parametrized Structure Functio
     +ns'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SPCTQ53(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.56) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  CTEQ Set 5HM1 (MSb) parametrized Structure Func
     +tions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SPCTQ51(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.GE.57) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.5) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  old GRV Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFGRVH(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  old GRV Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFGRVL(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVLO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV-94 Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVLO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94LO(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV-94 Set HO (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVHO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94HO(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV-94 Set HO (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVHO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94DI(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRSV Set LO polarized Structure Functions (stan
     +dard)'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRSVL1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRSV Set LO polarized Structure Functions (vale
     +nce)'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRSVL2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRSV Set NLO polarized Structure Functions (sta
     +ndard)'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRSVH1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRSV Set NLO polarized Structure Functions (val
     +ence)'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRSVH2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.12) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV98 Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRV98_lo
     +         (DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.13) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV98 Set nlm (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRV98_nlm
     +         (DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.14) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  GRV98 Set nld (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRV98_nld
     +         (DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.GE.15) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.6) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  ABFOW Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFABFOW(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.2) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.7) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  BM Set A Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFBMA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  BM Set B Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFBMB(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.3) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSE
        WRITE(N6,*)' NUCLEON Structure function call for unknown GROUP '
     +  ,NGROUP,',   and unknown NSET ',NSET
        STOP
      ENDIF
      ENDIF
C Now start PION structure functions
      IF(NPTYPE.EQ.2) THEN
      IF(NGROUP.EQ.1) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  OW-P Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STROWP1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  OW-P Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL STROWP2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.3) THEN
          WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.2) THEN
        WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
        STOP
      ELSEIF(NGROUP.EQ.3) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  SMRS-P Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SMRSP31(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  SMRS-P Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SMRSP32(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  SMRS-P Set 3 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SMRSP33(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.4) THEN
          WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.4) THEN
        WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
        STOP
      ELSEIF(NGROUP.EQ.5) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  GRV-P Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVPIHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  GRV-P Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVPILO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.GE.3) THEN
          WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.6) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  ABFKW-P Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL ABFKW1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  ABFKW-P Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL ABFKW2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Pion PDFs :  ABFKW-P Set 3 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '              Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL ABFKW3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.4) THEN
          WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSE
        WRITE(N6,*) ' PION Structure function call for unknown GROUP ',
     +  NGROUP,',   and unknown NSET ',NSET
        STOP
      ENDIF
      ENDIF
C Now start PHOTON structure functions
      IF(NPTYPE.EQ.3) THEN
      IF(NGROUP.EQ.1) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DO-G Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DOPHO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DO-G Set NLL Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DOPHO2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.3) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.2) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DG-G Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DGPHO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DG-G Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DGPHO2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DG-G Set 3 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DGPHO3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  DG-G Set 4 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL DGPHO4(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.5) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.3) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  LAC-G Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFLACG1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  LAC-G Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFLACG2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  LAC-G Set 3 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFLACG3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GAL-G Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
       CALL SFGAL(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
      ELSEIF(NSET.GE.5) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.4) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GS-G Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFGSHO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GS-G LO Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFGSLO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GS-G LO Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFGSLO2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GS-96-G HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GS96HO(DX,DQ,DUPV,DDNV,DSTR,DCHM,DBOT,DGL)
          DUSEA=DUPV
          DDSEA=DDNV
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GS-96-G LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GS96LO(DX,DQ,DUPV,DDNV,DSTR,DCHM,DBOT,DGL)
          DUSEA=DUPV
          DDSEA=DDNV
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.6) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.5) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GRV-G Set L-HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVGAH0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GRV-G Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVGAHO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GRV-G Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVGALO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) THEN
              WRITE(N6,*) 
     + ' Photon PDFs :  GRS-G Set LO Photon Structure Functions'
              WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
              WRITE(N6,*) 
     + '             !! User MUST call STRUCTP with the CORRECT INPUT pa
     +rameters !!'
          ENDIF
          GOTO 100
        ELSEIF(NSET.GE.5) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.6) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  ACFGP-G Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL ACFGP1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  ACFGP-G Set HO-mc Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL ACFGP2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  AFG-G Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFAFG1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.4) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.7) THEN
        WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
        STOP
      ELSEIF(NGROUP.EQ.8) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 1 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 2 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 3 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 4 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI4(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 5 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI5(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  WHIT-G Set 6 Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFWHI6(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.7) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.9) THEN
        IF(NSET.LE.8) THEN
          IF (FIRST) THEN
              WRITE(N6,*) 
     + ' Photon PDFs :  SaS Photon Structure Functions'
              WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
              WRITE(N6,*) 
     + '             !! User MUST call STRUCTP with the CORRECT INPUT pa
     +rameters !!'
          ENDIF
          GOTO 100
        ELSEIF(NSET.GE.9) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSE
        WRITE(N6,*)' PHOTON Structure function call for unknown GROUP ',
     +  NGROUP,',   and unknown NSET ',NSET
        STOP
      ENDIF
      ENDIF
C Now start NPDF's with Nuclear Corrections
      IF(NPTYPE.EQ.4) THEN
          IF (FIRST) THEN
              WRITE(N6,*) 
     + ' PDFLIB :  User wants NPDFs,  BUT  NPtype = ',NPTYPE,
     + '  is not valid !!'
              WRITE(N6,*) 
     + '        !! User MUST call STRUCTA for NPDFs !!' 
             WRITE(N6,*)
     + ' PDFLIB :  Please check ALL INPUT parameters for PDFSET and STRU
     +CTA !!'
          ENDIF
          GOTO 100
      ENDIF
  100 IF(FIRST) WRITE(N6,*)
     +' ----------------------------------------------------------------
     +----------------------'
      IF (FIRST) FIRST = .FALSE.
C
      UPV = MAX( ZEROD,DUPV)
      DNV = MAX( ZEROD,DDNV)
      USEA = MAX( ZEROD,DUSEA)
      DSEA = MAX( ZEROD,DDSEA)
      STR = MAX( ZEROD,DSTR)
      CHM = MAX( ZEROD,DCHM)
      BOT = MAX( ZEROD,DBOT)
      TOP = MAX( ZEROD,DTOP)
      GL  = MAX( ZEROD,DGL)
C
      RETURN
C
  200 IF(FIRST) WRITE(N6,*)
     +' ----------------------------------------------------------------
     +----------------------'
      IF (FIRST) FIRST = .FALSE.
C
      UPV = DUPV
      DNV = DDNV
      USEA = DUSEA
      DSEA = DDSEA
      STR = DSTR
      CHM = DCHM
      BOT = DBOT
      TOP = DTOP
      GL  = DGL
C
      RETURN
      END
