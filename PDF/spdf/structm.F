*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:56  cernlib
* Version 7.04
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
C       IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
C    +             ',   Nset = ',NSET,
C    +             ',  for GRV Set HO Structure Functions'
C       CALL GRVHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C       DDSEA = DUSEA
C       GOTO 100
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (G) (L255-MSb) Structure Functions'
          CALL SFMRSG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
      ENDIF
C Start with NUCLEON structure functions
      IF(NPTYPE.EQ.1) THEN
      IF(NGROUP.EQ.1) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for PRIVATE Structure Functions'
       CALL STRPRIV(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for BEBC Structure Functions'
          CALL STRBEBC(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for OR Structure Functions'
          CALL STRUCOR(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for BEP Structure Functions'
          CALL STRBEP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DGL)
          DDSEA = DUSEA
          DCHM=ZEROD
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GHR Structure Functions'
          CALL STRGHR(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DO Set 1 Structure Functions'
          CALL STRDO1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DO Set 2 Structure Functions'
          CALL STRDO2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for EHLQ Set 1 Structure Functions'
          CALL SFEHLQ1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for EHLQ Set 2 Structure Functions'
          CALL SFEHLQ2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for New DO Set 1.1 Structure Functions'
          CALL NEWDO1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GT.11) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',   Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.2) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMsoft (Valence) Structure Functions'
          CALL DFLM4(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMhard (Valence) Structure Functions'
          CALL DFLM5(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMsoft (Gluon) Structure Functions'
          CALL DFLM6(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMhard (Gluon) Structure Functions'
          CALL DFLM7(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMaverage (LO) Structure Functions'
          CALL DFLM8(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLMaverage (DIS) Structure Functions'
          CALL DFLM9(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLM160 (DIS) Structure Functions'
          CALL DFLM1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLM260 (DIS) Structure Functions'
          CALL DFLM2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DFLM360 (DIS) Structure Functions'
          CALL DFLM3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.GT.10) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',   Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.3) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set 1 (MSb) Structure Functions'
          CALL STRMRS1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set 2 (MSb) Structure Functions'
          CALL STRMRS2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set 3 (MSb) Structure Functions'
          CALL STRMRS3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set E Structure Functions'
          CALL STRUCE(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set B Structure Functions'
          CALL STRUCB(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set EP Structure Functions'
          CALL STRCEP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set BP Structure Functions'
          CALL STRCBP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS1 Set E Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL STRH1E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
C       ELSEIF(NSET.EQ.9) THEN
C         IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
C    +             ',   Nset = ',NSET,
C    +             ',  for HMRS1 Set B Structure Functions'
C         IF (FIRST) WRITE(N6,*)
C    +' Warning !!   :  in principle, set has been retracted by the auth
C    +ors !'
C         FIRST = .FALSE.
C         CALL STRH1B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
C         DDSEA = DUSEA
C         DTOP=ZEROD
C         GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for KMRS Set B0 (190) Structure Functions (old)'
          CALL STRKB0O(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS2 Set E Structure Functions'
          CALL STRH2E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS2 Set B Structure Functions'
          CALL STRH2B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.12) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set E+ Structure Functions'
          CALL STRH3EP(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.13) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set E- Structure Functions'
          CALL STRH3EM(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.14) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set E Structure Functions'
          CALL STRH3E(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.15) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set B (190) Structure Functions'
          CALL STRH3B(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.16) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set B (100) Structure Functions'
          CALL STRH3B1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.17) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for HMRS Set B (300) Structure Functions'
          CALL STRH3B3(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.18) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for KMRS Set B- Structure Functions'
          CALL STRKBM(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.19) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for KMRS Set B--R2-SH Structure Functions'
          CALL SFKBMR2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.20) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for KMRS Set B--R5-SH Structure Functions'
          CALL SFKBMR5(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.21) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +            ',  for KMRS Set B0 (190) updated Structure Functions'
          CALL STRKB0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.22) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set B0 (L135) Structure Functions'
          CALL STRKB02(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.23) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set B0 (L160) Structure Functions'
          CALL STRKB03(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.24) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set B0 (L200) Structure Functions'
          CALL STRKB04(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.25) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set B0 (L235) Structure Functions'
          CALL STRKB05(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.26) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set S0 (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.27) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D0 (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.28) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D- (L215) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFMRSDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.29) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set S0p (L230-MSb) Structure Functions'
          CALL SFMSBS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.30) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D0p (L230-MSb) Structure Functions'
          CALL SFMSBD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.31) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D-p (L230-MSb) Structure Functions'
          CALL SFMSBDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.32) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set S0p (L230-DIS) Structure Functions'
          CALL SFDISS0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.33) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D0p (L230-DIS) Structure Functions'
          CALL SFDISD0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.34) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set D-p (L230-DIS) Structure Functions'
          CALL SFDISDM(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.35) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (H) (L230-MSb) Structure Functions'
          CALL SFMSBH(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.36) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (H) (L230-DIS) Structure Functions'
          CALL SFDISH(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.37) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (A) (L230-MSb) Structure Functions'
          CALL SFMSBA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.38) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Fit (A) (L230-MSb) Structure Functions'
          CALL SFMSFA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.39) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L231-MSb) Structure Functions'
          CALL SFMSAP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.40) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Fit (Ap) (L231-MSb) Structure Functions'
          CALL SFMFAP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.41) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (G) (L255-MSb) Structure Functions'
          CALL SFMRSG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.42) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Fit (G) (L255-MSb) Structure Functions'
          CALL SFMSFG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.43) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (A) (L230-MSb) and low Q**2 Structure
     + Functions'
          IF(DQ2.GT.5.D0)
     +    CALL SFMSBA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          IF(DQ2.LE.5.D0)
     +    CALL SFMSQA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.44) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (A) (L230-DIS) and low Q**2 Structure
     + Functions'
          IF(DQ2.GT.5.D0)
     +    CALL SFDISA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          IF(DQ2.LE.5.D0)
     +    CALL SFDSQA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.45) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L150-MSb) Structure Functions'
          CALL SFMS105(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.46) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L201-MSb) Structure Functions'
          CALL SFMS110(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.47) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L266-MSb) Structure Functions'
          CALL SFMS115(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.48) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L344-MSb) Structure Functions'
          CALL SFMS120(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.49) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L435-MSb) Structure Functions'
          CALL SFMS125(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.50) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Ap) (L542-MSb) Structure Functions'
          CALL SFMS130(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.51) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (J) (L344-MSb) Structure Functions'
          CALL SFMRSJ(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.52) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (Jp) (L507-MSb) Structure Functions'
          CALL SFMSJP(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.53) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (R1) (L241-MSb) Structure Functions'
          CALL SFMSR1(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.54) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (R2) (L344-MSb) Structure Functions'
          CALL SFMSR2(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.55) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (R3) (L241-MSb) Structure Functions'
          CALL SFMSR3(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.56) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Set (R4) (L344-MSb) Structure Functions'
          CALL SFMSR4(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.57) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Fit (R1) (L241-MSb) Structure Functions'
          CALL SFMFR1(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.58) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MRS Fit (R2) (L344-MSb) Structure Functions'
          CALL SFMFR2(DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.59) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.4) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set S1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set B1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set B2 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set E1 (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU4(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set 6 (1/2s) (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU5(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set S1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU6(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set B1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU7(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set B2 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU8(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set E1 (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMTU9(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set 6 (1/2s) (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMT10(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for MT Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
        CALL SFMT11(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.12) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 1L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ11(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.13) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 1M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ12(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.14) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 1MS (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ13(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.15) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 1ML (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ14(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.16) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 1D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ15(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.17) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2L (LO) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ21(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.18) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2M (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ22(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.19) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2MS (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ23(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.20) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2MF (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ24(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.21) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2ML (MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ25(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.22) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2D (DIS) Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
       CALL SFCTQ26(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.23) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pL (LO) Structure Functions'
       CALL SRCTQ21(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.24) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pM (MSb) Structure Functions'
       CALL SRCTQ22(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.25) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pMS (MSb) Structure Functions'
       CALL SRCTQ23(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.26) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pMF (MSb) Structure Functions'
       CALL SRCTQ24(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.27) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pML (MSb) Structure Functions'
       CALL SRCTQ25(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.28) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 2pD (DIS) Structure Functions'
       CALL SRCTQ26(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.29) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 3L (LO) Structure Functions'
       CALL SFCTQ31(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.30) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 3M (MSb) Structure Functions'
       CALL SFCTQ32(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.31) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 3D (DIS) Structure Functions'
       CALL SFCTQ33(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.32) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4L (LO) Structure Functions'
       CALL SFCTQ41(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.33) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4D (DIS) Structure Functions'
       CALL SFCTQ42(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.34) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4M (MSb) Structure Functions'
       CALL SFCTQ43(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.35) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4A1 (MSb) Structure Functions'
       CALL SFCTQ44(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.36) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4A2 (MSb) Structure Functions'
       CALL SFCTQ45(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.37) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4A3 (MSb) Structure Functions'
       CALL SFCTQ43(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.38) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4A4 (MSb) Structure Functions'
       CALL SFCTQ46(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.39) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4A5 (MSb) Structure Functions'
       CALL SFCTQ47(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.40) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4HJ (MSb) Structure Functions'
       CALL SFCTQ48(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.41) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for CTEQ Set 4LQ (MSb) Structure Functions'
       CALL SFCTQ49(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.GE.42) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.5) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for old GRV Set HO Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFGRVH(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for old GRV Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*)
     +' Warning !!   :  in principle, set has been retracted by the auth
     +ors !'
          CALL SFGRVL(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV Set HO Structure Functions'
          CALL GRVHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV Set LO Structure Functions'
          CALL GRVLO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-94 Set LO Structure Functions'
          CALL GRVLO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94LO(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-94 Set HO (MSb) Structure Functions'
          CALL GRVHO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94HO(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-94 Set HO (DIS) Structure Functions'
          CALL GRVHO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZSTR,DCHM,DBOT,DTOP,ZGL)
          CALL GRV94DI(DX,DQ,
     +                 DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,ZTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +    ',  for GRSV Set LO polarized Structure Functions (standard)'
          CALL GRSVL1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.9) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +    ',  for GRSV Set LO polarized Structure Functions (valence)'
          CALL GRSVL2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.10) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +    ',  for GRSV Set NLO polarized Structure Functions (standard)'
          CALL GRSVH1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.EQ.11) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +    ',  for GRSV Set NLO polarized Structure Functions (valence)'
          CALL GRSVH2(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 200
        ELSEIF(NSET.GE.12) THEN
          WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.6) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for ABFOW Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for BM Set A Structure Functions'
          CALL SFBMA(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Nucleon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for BM Set B Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for OW-P Set 1 Structure Functions'
          CALL STROWP1(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DGL)
          DDSEA = DUSEA
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for OW-P Set 2 Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for SMRS-P Set 1 Structure Functions'
          CALL SMRSP31(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for SMRS-P Set 2 Structure Functions'
          CALL SMRSP32(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DGL)
          DDSEA = DUSEA
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for SMRS-P Set 3 Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-P Set HO Structure Functions'
          CALL GRVPIHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          DDSEA = DUSEA
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-P Set LO Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for ABFKW-P Set 1 Structure Functions'
          CALL ABFKW1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for ABFKW-P Set 2 Structure Functions'
          CALL ABFKW2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Pion PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for ABFKW-P Set 3 Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DO-G Set LO Structure Functions'
          CALL DOPHO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DO-G Set NLL Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DG-G Set 1 Structure Functions'
          CALL DGPHO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DG-G Set 2 Structure Functions'
          CALL DGPHO2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DG-G Set 3 Structure Functions'
          CALL DGPHO3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for DG-G Set 4 Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for LAC-G Set 1 Structure Functions'
       CALL SFLACG1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for LAC-G Set 2 Structure Functions'
       CALL SFLACG2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for LAC-G Set 3 Structure Functions'
       CALL SFLACG3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
      ELSEIF(NSET.GE.4) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.EQ.4) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GS-G Set HO Structure Functions'
          CALL SFGSHO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GS-G LO Set 1 Structure Functions'
          CALL SFGSLO1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GS-G LO Set 2 Structure Functions'
          CALL SFGSLO2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GS-96-G HO Structure Functions'
          CALL GS96HO(DX,DQ,DUPV,DDNV,DSTR,DCHM,DBOT,DGL)
          DUSEA=DUPV
          DDSEA=DDNV
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GS-96-G LO Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-G Set L-HO Structure Functions'
          CALL GRVGAH0(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-G Set HO Structure Functions'
          CALL GRVGAHO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRV-G Set LO Structure Functions'
          CALL GRVGALO(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) THEN
              WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for GRS-G Set LO Photon Structure Functions'
              WRITE(N6,*) ' !! You MUST call STRUCTP ' 
              WRITE(N6,*) ' !! to get the INPUT parameters right !!!'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for ACFGP-G Set HO Structure Functions'
          CALL ACFGP1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for ACFGP-G Set HO-mc Structure Functions'
          CALL ACFGP2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for AFG-G Set HO Structure Functions'
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
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 1 Structure Functions'
          CALL SFWHI1(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 2 Structure Functions'
          CALL SFWHI2(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 3 Structure Functions'
          CALL SFWHI3(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 4 Structure Functions'
          CALL SFWHI4(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 5 Structure Functions'
          CALL SFWHI5(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
          DBOT=ZEROD
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +          ',  for WHIT-G Set 6 Structure Functions'
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
              WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  for SAS Photon Structure Functions'
              WRITE(N6,*) ' !! You MUST call STRUCTP ' 
              WRITE(N6,*) ' !! to get the INPUT parameters right !!!'
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
  100 IF(FIRST) WRITE(N6,*)
     +' ----------------------------------------------------------------
     +------------------------------------------------'
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
     +------------------------------------------------'
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
