#include "pdf/pilot.h"
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STRUCTP(X,Q2,P2,IP,
     +                             UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C
C  *********************************************************************
C  *                                                                   *
C  *   Main steering routine for all sets of structure functions       *
C  *                                                                   *
C  *                                                                   *
C  *   Input:    X     = x value of parton                             *
C  *             Q2    = QCD scale in GeV**2                           *
C  *             P2    = Virtuality of the Photon in GeV**2            *
C  *             IP    = Flag to evaluate off-shell anomalous photon   *
C  *                     component                                     *
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
        DQ  = SQRT(Q2)
        DQ2 = Q2
        IF (NPTYPE.EQ.3) THEN
            DP2 = P2
            IP2 = IP
        ELSE
            DP2 = ZEROD
            IP2 = 0
        ENDIF
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
C          IF (FIRST) WRITE(N6,*) 
C     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
C       CALL GRVHO(DX,DQ,DUPV,DDNV,DUSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C       DDSEA = DUSEA
C       GOTO 100
          IF (FIRST) WRITE(N6,*) 
     + ' Nucleon PDFs :  MRS Set (G) (L255-MSb) Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFMRSG(DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
          DTOP=ZEROD
          GOTO 100
      ENDIF
C Start with NUCLEON structure functions
      IF(NPTYPE.EQ.1 .OR. NPTYPE.EQ.2) THEN
         CALL STRUCTM(DX,DQ,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
      ENDIF
C Now start PHOTON structure functions
      IF(NPTYPE.EQ.3) THEN
      IF(NGROUP.LE.4) THEN
         CALL STRUCTM(DX,DQ,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
      ELSEIF(NGROUP.EQ.5) THEN
        IF(NSET.LE.3) THEN
         CALL STRUCTM(DX,DQ,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  GRS-G Set LO Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL GRVGALO(DX,DQ,ZUPV,ZDNV,ZUSEA,ZDSEA,ZSTR,DCHM,DBOT,ZGL)
          CALL GRSGALO(DX,DQ2,DP2,
     +                     DUPV,DDNV,DUSEA,DDSEA,DSTR,ZCHM,ZBOT,DGL)
          DTOP=ZEROD
          GOTO 100
        ELSEIF(NSET.GE.5) THEN
          WRITE(N6,*) ' Photon PDFs :  Ngroup = ',NGROUP,
     +             ',   Nset = ',NSET,
     +             ',  Structure Functions not yet exsistent'
          STOP
        ENDIF
      ELSEIF(NGROUP.GE.6 .AND. NGROUP.LE.8) THEN
         CALL STRUCTM(DX,DQ,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
      ELSEIF(NGROUP.EQ.9) THEN
        IF(NSET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 95 Set 1D Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS11(DX,DQ2,DP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.2) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 95 Set 1M Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS12(DX,DQ2,DP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.3) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 95 Set 2D Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS13(DX,DQ2,DP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.4) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 95 Set 2M Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS14(DX,DQ2,DP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.5) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 96 Set 1D Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS21(DX,DQ2,DP2,IP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.6) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 96 Set 1M Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS22(DX,DQ2,DP2,IP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.7) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 96 Set 2D Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS23(DX,DQ2,DP2,IP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NSET.EQ.8) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Photon PDFs :  SaS-G 96 Set 2M Structure Functions'
          IF (FIRST) WRITE(N6,*) 
     + '                Ngroup = ',NGROUP,',   Nset = ',NSET
          CALL SFSAS24(DX,DQ2,DP2,IP2,
     +                DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
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
      END
