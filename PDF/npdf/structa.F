#include "pdf/pilot.h"
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STRUCTA(X,SCALE,ANO,
     +                             UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C
C  *********************************************************************
C  *                                                                   *
C  *   Main steering routine for all sets of structure functions       *
C  *                                                                   *
C  *                                                                   *
C  *   Input:    X     = x value of parton                             *
C  *             SCALE = QCD scale in GeV                              *
C  *             ANO   = atomic mass number of a nucleus               *
C  *                                                                   *
C  *   Output:   UPV   = up valence quark * Ruv-Npdf                   *
C  *             DNV   = down valence quark * Rdv-Npdf                 *
C  *             USEA  = sea (up_bar) * Rub-Npdf                       *
C  *             DSEA  = sea (down_bar) * Rdb-Npdf                     *
C  *             STR   = strange quark * Rs-Npdf                       *
C  *             CHM   = charm quark * Rc-Npdf                         *
C  *             BOT   = bottom quark * Rb-Npdf                        *
C  *             TOP   = top quark * Rt-Npdf                           *
C  *             GL    = gluon * Rg-Npdf                               *
C  *                                                                   *
C  *                                                                   *
C  *   The variables  NPTYPE, NGROUP and NSET  and                     *
C  *   the variables  NATYPE, NAGROUP and NASET  should be,            *
C  *   could be provided by the user via a call to the                 *
C  *   subroutine PDFSET at the initialization phase, where            *
C  *                                                                   *
C  *             NPTYPE = Particle type for proton PDF's               *
C  *                      (number or character string,                 *
C  *                       1,2,3 or 'NU','PI','PH')                    *
C  *                      of desired structure functions set           *
C  *                      (Default: NPTYPE = 1  or  'NU'               *
C  *             NGROUP = author group for proton PDF's                *
C  *                      (number or character string,                 *
C  *                       1 to 7 or i.e. 'DFLM','MRS','MT','GRV',etc) *
C  *                      of desired structure functions set           *
C  *                      (Default: NGROUP = 5  or  'GRV')             *
C  *             NSET   = number of desired structure functions set    *
C  *                      for proton PDF's                             *
C  *                      (Default: NSET = 3)                          *
C  *             NATYPE = 4  for Nuclear Corrections                   *
C  *                      (Default: NATYPE = 4)                        *
C  *             NAGROUP= author group of Nuclear Corrections          *
C  *                      (number or character string, 1 or i.e. 'EKS')*
C  *                      of desired uclear Corrections set            *
C  *                      (Default: NAGROUP = 1)                       *
C  *             NASET  = number of desired Nuclear Corrections set    *
C  *                      (Default: NASET = 1)                         *
C  *                                                                   *
C  *   for each set of structure function.                             *
C  *                                                                   *
C  *   The internal COMMON blocks                                      *
C  *                                                                   *
C  *   COMMON/W50511/  NPTYPE,NGROUP,NSET                              *
C  *   COMMON/W50511a/ NATYPE,NAGROUP,NASET                            *
C  *                                                                   *
C  *   are filled then by the subroutine PDFSET.                       *
C  *                                                                   *
C  *                                                                   *
C  *     Note: STRUCTA returns                                         *
C  *           X * parton distribution function * Nuclear Correction ! *
C  *                                                                   *
C  *                                                                   *
C  *     Author:   H. Plothow-Besch                                    *
C  *               CERN-ETT/TT, CH - 1211 Geneva 23, Switzerland       *
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
#include "pdf/w50511a.inc"
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
         NATYPE  = LATYPE
         NAGROUP = LAGROUP
         NASET   = LNASET
C
         PARM(4) = 'NAtype'
         VAL(4)  = NATYPE
         PARM(5) = 'NAgroup'
         VAL(5)  = NAGROUP
         PARM(6) = 'NAset'
         VAL(6)  = NASET
C
         CALL PDFSET(PARM,VAL)
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
        DX  = X
        DQ  = SCALE
        DQ2 = SCALE*SCALE
        DANO = ANO
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
C      IF    (NPTYPE.EQ.0 .OR. NGROUP.EQ.0 .OR. NSET.EQ.0) THEN
C          IF (FIRST) WRITE(N6,*)
C     + ' Nucleon PDFs :  MRS central-g (L300-MSb) Structure Functions'
C          IF (FIRST) WRITE(N6,*) 
C     + '                 Ngroup = ',NGROUP,',   Nset = ',NSET
C       CALL SFMRS_99A1
C     +               (DX,DQ2,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
C          DTOP=ZEROD
C          GOTO 100
C      ENDIF
C Start with NUCLEAR structure functions
      IF(NATYPE.EQ.4) THEN
      IF(NAGROUP.EQ.1) THEN
        IF(NASET.EQ.1) THEN
          IF (FIRST) WRITE(N6,*) 
     + ' Nuclear PDFs :  EKS98 Structure Functions with Nuclear Correcti
     +ons'
          IF (FIRST) WRITE(N6,*) 
     + '                 NAgroup = ',NAGROUP,',   NAset = ',NASET
       CALL SFEKS98(DX,DQ,DANO,
     +              DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
          GOTO 100
        ELSEIF(NASET.GE.2) THEN
          WRITE(N6,*) ' Nuclear PDFs :  NAgroup = ',NAGROUP,
     +             ',   NAset = ',NASET,
     +             ',   Structure Functions not yet exsistent'
          STOP
        ENDIF
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
C
      END
