#include "isajet/pilot.h"
      SUBROUTINE DIALOG
C
C          Set up user data with interactive commands from standard
C          Fortran input.
C
C          V6.26: Fixed FORCE to allow automatic top decays and to use
C          correct limit for maximum number.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/mbgen.inc"
#include "isajet/force.inc"
#include "isajet/eepar.inc"
#include "isajet/dkytab.inc"
#include "isajet/qcdpar.inc"
#include "isajet/idrun.inc"
#include "isajet/frgpar.inc"
#include "isajet/keys.inc"
#include "isajet/prtout.inc"
#include "isajet/seed.inc"
#include "isajet/types.inc"
#include "isajet/primar.inc"
#include "isajet/jetlim.inc"
#include "isajet/nodcay.inc"
#include "isajet/wcon.inc"
#include "isajet/dylim.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetpar.inc"
#include "isajet/isloop.inc"
#include "isajet/tcpar.inc"
#include "isajet/xmssm.inc"
#include "isajet/sugnu.inc"
#if defined(CERNLIB_PDFLIB)
#include "isajet/w50510.inc"
* Ignoring t=pass
#endif
#if defined(CERNLIB_PDFLIB)
#include "isajet/w50517.inc"
* Ignoring t=pass
#endif
#include "isajet/hcon.inc"
#include "isajet/kkgrav.inc"
#include "isajet/mglims.inc"
C
#include "isajet/const.inc"
C
      LOGICAL SETTYP,DUMY
      CHARACTER*8 TTL(10),WORD,LSTRUC,BLANK,LIST(29)
      CHARACTER*8 WTYP(4),RDID(2)
      CHARACTER*4 YN
      INTEGER IFL,I1,I2,I3,J1,I,IKEY,IJ,J,KK,IDXQK
      INTEGER NSEL,K,KFORCE(5),INDEX,IDG1,IDG2,IDG3,IDG4,IDXQ,IDXLEP
#if defined(CERNLIB_SINGLE)
      REAL    SEED
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION SEED
#endif
      INTEGER NLAP(3,17)
      REAL EMIN,YLIM1,YLIM2,ANGL1,ANGL2,XLIM1,XLIM2,AMW,AMZ
      INTEGER NFTEMP
      CHARACTER*8 HTYPE
#if defined(CERNLIB_PDFLIB)
      CHARACTER*20 PDFPAR(20)
#endif
#if (defined(CERNLIB_PDFLIB))&&(defined(CERNLIB_SINGLE))
      REAL PDFVAL(20)
      REAL DX,DSCALE,DXPDF(-6:6)
#endif
#if (defined(CERNLIB_PDFLIB))&&(defined(CERNLIB_DOUBLE))
      DOUBLE PRECISION PDFVAL(20)
      DOUBLE PRECISION DX,DSCALE,DXPDF(-6:6)
#endif
      REAL AMLIM1,AMLIM2
      INTEGER JLIM1,JLIM2
C
C          Overlapping variable flags.
      DATA NLAP/1,2,3, 1,2,7 ,1,2,8, 1,3,5, 1,3,6, 1,3,7, 1,3,8, 1,5,7,
     X          1,5,8, 1,6,7, 1,6,8, 2,3,7, 2,3,8, 3,5,7, 3,6,7, 3,5,8,
     X          3,6,8/
      DATA BLANK/'        '/
C          Duplicate data to avoid requiring ALDATA
      DATA LIST/'GL','UP','UB','DN','DB','ST','SB','CH','CB','BT','BB',
     $'TP','TB','NUE','ANUE','E-','E+','NUM','ANUM','MU-','MU+',
     $'NUT','ANUT','TAU-','TAU+','ALL','QUARKS','LEPTONS','NUS'/
C
C          Entry
C
      IFL=0
      NFTEMP=0
C
1111  PRINT 11
11    FORMAT(//,10X,' Use CAPITAL LETTERS only for all commands.',//,
     1 10X,' Give title for this run.')
C
C          Read title
C
      READ 1,TTL
    1 FORMAT(10A8)
      WRITE(ITCOM,1) TTL
      IF(TTL(1).EQ.'STOP    ') GOTO 999
C
C          Read energy and no. of events
C
      PRINT*,' total energy, # of events, # to print, # to skip.'
      READ*,ECM,NEVENT,NEVPRT,NJUMP
      WRITE(ITCOM,*)ECM,NEVENT,NEVPRT,NJUMP
C
C          Reset all variables if title is not 'SAME'
C
      IF(TTL(1).NE.'SAME    ') THEN
        DO 20 I=1,10
   20   TITLE(I)=TTL(I)
        CALL RESET
        KEYON=.FALSE.
C
C          Read reaction and set keys and NJET
C
   21   PRINT 22
   22   FORMAT(/,' Choose one of the following reactions:',/,
     $' TWOJET E+E- DRELLYAN MINBIAS SUSY WPAIR HIGGS',
     $' PHOTON TCOLOR WHIGGS')
        READ 3,REAC
    3   FORMAT(A8)
        DO 18 I=1,MXKEYS
18      KEYS(I)=.FALSE.
        KEYON=.FALSE.
C
        IF(REAC.EQ.'TWOJET  ') THEN
          KEYS(1)=.TRUE.
          IKEY=1
          NJET=2
        ELSEIF(REAC.EQ.'E+E-    ') THEN
          KEYS(2)=.TRUE.
          IKEY=2
          NJET=2
        ELSEIF(REAC.EQ.'DRELLYAN') THEN
          KEYS(3)=.TRUE.
          IKEY=3
          NJET=3
        ELSEIF(REAC.EQ.'MINBIAS ') THEN
          KEYS(4)=.TRUE.
          IKEY=4
          NJET=0
        ELSEIF(REAC.EQ.'SUPERSYM'.OR.REAC.EQ.'SUSY    ') THEN
          KEYS(5)=.TRUE.
          IKEY=5
          NJET=2
        ELSEIF(REAC.EQ.'WPAIR   ') THEN
          KEYS(6)=.TRUE.
          IKEY=6
          NJET=2
        ELSEIF(REAC.EQ.'HIGGS   ') THEN
          KEYS(7)=.TRUE.
          IKEY=7
          NJET=2
       ELSEIF(REAC.EQ.'PHOTON  ') THEN
         KEYS(8)=.TRUE.
         IKEY=8
         NJET=2
        ELSEIF(REAC.EQ.'TCOLOR  ') THEN
          KEYS(9)=.TRUE.
          IKEYS=9
          NJET=2
        ELSEIF(REAC.EQ.'WHIGGS  ') THEN
          KEYS(10)=.TRUE.
          IKEY=10
          NJET=2
        ELSE
          PRINT* ,' You must choose a valid reaction, try again.'
          GOTO 21
        ENDIF
C
        WRITE(ITCOM,3) REAC
      ENDIF
C
      EMIN=1.0
      SCM=ECM**2
      HALFE=ECM/2
      YLIM2=ALOG(ECM)
      YLIM1=-YLIM2
      NSEL=0
      ANGL1=0
      ANGL2=2*PI
      XLIM1=-1.0
      XLIM2=1.0
C
C          Read keyword. For each recognized keyword read corresponding
C          variables and set LOC flag.
C
      NSEL=0
1000  PRINT 99
99    FORMAT(//,
     $'      Choose among the following:',/,
     $'      Jet limits:',/,
     $' JETTYPE1 JETTYPE2 JETTYPE3 P PHI PT TH X Y WMODE1 WMODE2',/,
     $'      Drell Yan/Higgs limits:',/,
     $' HTYPE PHIW QMH QMW QTW THW WTYPE XW YW',/,
     $'      Decays:',/,
     $' FORCE FORCE1 NODECAY NOETA NOEVOLVE NOFRGMNT NOPI0',/,
     $'      Physics parameters:',/,
     $' CUTJET CUTOFF FRAGMENT GAUGINO HMASS HMASSES LAMBDA',/,
     $' MSSMA MSSMB MSSMC MSSMD MSSME SIGQT SIN2W SLEPTON SQUARK',/,
     $' SUGRA TCMASS TMASS WMASS XGEN',/,
     $'      Other:',/,
     $' BEAMS EPOL NPOMERON NSIGMA NTRIES PDFLIB SEED STRUC WFUDGE',/,
     $' WMMODE WPMODE Z0MODE',/,
     $'      Terminate with END. Use HELP for help.',/)
C
100   PRINT* ,' Give a variable name.'
      READ 3,WORD
      IF(WORD.EQ.'HELP    ') GO TO 1000
      NSEL=NSEL+1
C
C          Keyword END
      IF(WORD.EQ.'END     ') THEN
C
C              End of run
C
          WRITE(ITCOM,3) WORD
          PRINT 889
  889 FORMAT(//,10X,' NEW RUN',/,' If you are finished answer STOP.',/,
     1  ' If you answer SAME the parameters from previous run will be',
     2  ' used,',/,' unless you request otherwise.')
        GO TO 101
      ENDIF
C
C          Keyword P
      IF(WORD.EQ.'P       ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,EMIN,HALFE
        READ*, (PMIN(K),PMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(PMIN(K),PMAX(K),K=1,NJET)
        LOC(1)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword Y
      IF(WORD.EQ.'Y       ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,YLIM1,YLIM2
        READ*, (YJMIN(K),YJMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(YJMIN(K),YJMAX(K),K=1,NJET)
        LOC(2)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword X
      IF(WORD.EQ.'X       ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,XLIM1,XLIM2
        READ*, (XJMIN(K),XJMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(XJMIN(K),XJMAX(K),K=1,NJET)
        LOC(3)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword PHI
      IF(WORD.EQ.'PHI     ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,ANGL1,ANGL2
        READ*, (PHIMIN(K),PHIMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(PHIMIN(K),PHIMAX(K),K=1,NJET)
        LOC(4)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword TH
        IF(WORD.EQ.'TH      '.OR.WORD.EQ.'THETA   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,ANGL1,PI
        READ*, (THMIN(K),THMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(THMIN(K),THMAX(K),K=1,NJET)
        LOC(5)=NSEL
        LOC(6)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword PT
      IF(WORD.EQ.'PT      '.OR.WORD.EQ.'PPERP   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 562
        PRINT 563,EMIN,HALFE
        READ*, (PTMIN(K),PTMAX(K),K=1,NJET)
        WRITE(ITCOM,*)(PTMIN(K),PTMAX(K),K=1,NJET)
        LOC(7)=NSEL
        LOC(8)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NODECAY
      IF(WORD.EQ.'NODECAY ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571,NODCAY
571     FORMAT(L1)
        WRITE(ITCOM,571) NODCAY
        LOC(9)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NOETA
      IF(WORD.EQ.'NOETA   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571, NOETA
        WRITE(ITCOM,571) NOETA
        LOC(10)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NOPI0
      IF(WORD.EQ.'NOPI0   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571, NOPI0
        WRITE(ITCOM,571)  NOPI0
        LOC(11)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword BEAMS
      IF(WORD.EQ.'BEAMS   ') THEN
        WRITE(ITCOM,3) WORD
   58   PRINT*,' Select each beam enclosed in single quotes.'
        PRINT*,' Allowed names are: P  N  AP  AN'
        READ *,RDID(1),RDID(2)
        IDIN(1)=0
        IDIN(2)=0
        DO 583 K=1,2
        IF(RDID(K).EQ.'P       ') IDIN(K)=+1120
        IF(RDID(K).EQ.'AP      ') IDIN(K)=-1120
        IF(RDID(K).EQ.'N       ') IDIN(K)=+1220
        IF(RDID(K).EQ.'AN      ') IDIN(K)=-1220
583     CONTINUE
        IF(IDIN(1)*IDIN(2).EQ.0) GOTO 58
        LOC(12)=NSEL
        WRITE(ITCOM,4) RDID(1),RDID(2)
        GO TO 1000
      ENDIF
C
C          Keyword FRAGMENT
      IF(WORD.EQ.'FRAGMENT') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, FRPAR
        WRITE(ITCOM,*) FRPAR
        LOC(13)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SEED
      IF(WORD.EQ.'SEED    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, SEED
        WRITE(ITCOM,*) SEED
        CALL RANFST(SEED)
        LOC(14)=NSEL
        GO TO 1000
      ENDIF
C
C          Keywords JETTYPE1, JETTYPE2, JETTYPE3
      IF(WORD.EQ.'JETTYPE1'.OR.WORD.EQ.'JETTYPE2'.OR.WORD.EQ.
     $'JETTYPE3') THEN
        WRITE(ITCOM,3) WORD
C          Print allowed types
60      IF(KEYS(6)) THEN
          PRINT 6001
          PRINT 6002
6001      FORMAT(' Give a list of jet types enclosed in single quotes ',
     $    'and separated by commas.',/,' Terminate list with a /')
6002      FORMAT(' The recognized W types are: W+, W-, Z0, GM and ALL')
        ELSE
          PRINT 6001
          PRINT 6003, LIST
6003      FORMAT(' The recognized jet types are:',/,4(1X,8A8,/),1X,5A8)
        ENDIF
        IF(KEYS(1)) THEN
          PRINT 6004
6004      FORMAT(' A fourth generation is also allowed:',/,
     $    ' X, XB for up, Y, YB for down.',/,
     $    ' The fourth generation is normally turned off.')
        ENDIF
        IF(WORD.EQ.'JETTYPE1') IJ=1
        IF(WORD.EQ.'JETTYPE2') IJ=2
        IF(WORD.EQ.'JETTYPE3') IJ=3
        READ*,(JETYP(K,IJ),K=1,25)
        DO 61 K=1,25
61      IF(JETYP(K,IJ).NE.BLANK) NJTTYP(IJ)=NJTTYP(IJ)+1
C          Check that only legal jet types are in the list
        IF(SETTYP(0)) GOTO 60
        WRITE(ITCOM,4)(JETYP(K,IJ),K=1,NJTTYP(IJ))
    4   FORMAT(1H',A,2H',)
        IF(NJTTYP(IJ).LT.25) WRITE(ITCOM,*) '/'
        LOC(15)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SIN2W
      IF(WORD.EQ.'SIN2W   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, SIN2W
        WRITE(ITCOM,*) SIN2W
        LOC(17)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword TMASS
      IF(WORD.EQ.'TMASS   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, AMLEP(6),AMLEP(7),AMLEP(8)
        WRITE(ITCOM,*) AMLEP(6),AMLEP(7),AMLEP(8)
        LOC(18)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword QMH (note that it uses same variable as QMW)
      IF(WORD.EQ.'QMH     ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum MASS for Higgs.'
        READ*,QMIN,QMAX
        WRITE(ITCOM,*)QMIN,QMAX
        LOC(19)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword QMW
      IF(WORD.EQ.'QMW     ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum MASS for W(Z0).'
        READ*,QMIN,QMAX
        WRITE(ITCOM,*)QMIN,QMAX
        LOC(19)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword QTW
      IF(WORD.EQ.'QTW     ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum PT for W(Z0).'
        PRINT 563,EMIN,HALFE
        PRINT*,' Fix QTW to zero for standard Drell-Yan.'
        READ*, QTMIN,QTMAX
        WRITE(ITCOM,*) QTMIN,QTMAX
        LOC(20)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword YW
      IF(WORD.EQ.'YW      ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum Y for W(Z0).'
        PRINT 563,YLIM1,YLIM2
        READ*, YWMIN,YWMAX
        WRITE(ITCOM,*)YWMIN,YWMAX
        LOC(21)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword XW
      IF(WORD.EQ.'XW      ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum X for W(Z0).'
        PRINT 563,XLIM1,XLIM2
        READ*, XWMIN,XWMAX
        WRITE(ITCOM,*)XWMIN,XWMAX
        LOC(22)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword THW
      IF(WORD.EQ.'THW     ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum THETA for W(Z0).'
        PRINT 563,ANGL1,PI
        READ*, THWMIN,THWMAX
        WRITE(ITCOM,*)THWMIN,THWMAX
        LOC(23)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword PHIW
      IF(WORD.EQ.'PHIW    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give minimum and maximum PHI for W(Z0).'
        PRINT 563,ANGL1,ANGL2
        READ*, PHWMIN,PHWMAX
        WRITE(ITCOM,*)PHWMIN,PHWMAX
        LOC(24)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NONUNU
      IF(WORD.EQ.'NONUNU  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571,NONUNU
        WRITE(ITCOM,571) NONUNU
        LOC(25)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword WTYPE
      IF(WORD.EQ.'WTYPE   ') THEN
        WRITE(ITCOM,3) WORD
   92   PRINT*,' Choose from one set of bosons: W+,W- or Z0,GM'
        PRINT*,' Enclose each in single quotes and terminate with / '
        DO 191 J=1,4
        WTYP(J)=BLANK
  191   GODY(J)=.FALSE.
        READ*,WTYP
        DO 192 K=1,4
        IF(WTYP(K).EQ.'GM      ') GODY(1)=.TRUE.
        IF(WTYP(K).EQ.'W+      ') GODY(2)=.TRUE.
        IF(WTYP(K).EQ.'W-      ') GODY(3)=.TRUE.
        IF(WTYP(K).EQ.'Z0      ') GODY(4)=.TRUE.
  192   CONTINUE
        IF(GODY(1)) JWTYP=1
        IF(GODY(2).OR.GODY(3)) JWTYP=3
        IF(GODY(4)) JWTYP=4
        IF((GODY(3).OR.GODY(2)).AND.(GODY(1).OR.GODY(4)).OR.JWTYP
     1  .EQ.0) THEN
          PRINT 2003
          GOTO 92
        ELSE
          WRITE(ITCOM,4) WTYP
          LOC(26)=NSEL
        ENDIF
        GO TO 1000
      ENDIF
C
C          Keyword LAMBDA
      IF(WORD.EQ.'LAMBDA  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, ALAM
        WRITE(ITCOM,*) ALAM
        ALAM2=ALAM**2
        LOC(27)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NTRIES
      IF(WORD.EQ.'NTRIES  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,NTRIES
        WRITE(ITCOM,*)NTRIES
        LOC(28)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword CUTOFF
      IF(WORD.EQ.'CUTOFF  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,CUTOFF,CUTPOW
        WRITE(ITCOM,*)CUTOFF,CUTPOW
        LOC(29)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword XGEN
      IF(WORD.EQ.'XGEN    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,XGEN
        WRITE(ITCOM,*)XGEN
        LOC(30)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SIGQT
      IF(WORD.EQ.'SIGQT   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,SIGQT
        WRITE(ITCOM,*)SIGQT
        LOC(31)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword CUTJET
      IF(WORD.EQ.'CUTJET  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, CUTJET
        WRITE(ITCOM,*) CUTJET
        LOC(32)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword WFUDGE
      IF(WORD.EQ.'WFUDGE  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, WFUDGE
        WRITE(ITCOM,*) WFUDGE
        LOC(50)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword STRUC
      IF(WORD.EQ.'STRUC   ') THEN
        WRITE(ITCOM,3) WORD
        ISTRUC=0
   69   PRINT*,' Choose a structure function: OWENS, BAIER or EICHTEN'
        PRINT*,' Enclose in single quotes'
        READ*,LSTRUC
        IF(LSTRUC.EQ.'OWENS   ') ISTRUC=1
        IF(LSTRUC.EQ.'BAIER   ') ISTRUC=2
        IF(LSTRUC.EQ.'EICHTEN ') ISTRUC=3
        IF(LSTRUC.EQ.'DUKE    '.OR.LSTRUC.EQ.'DO      ') ISTRUC=4
        IF(LSTRUC.EQ.'CTEQ2L  ') ISTRUC=5
        IF(LSTRUC.EQ.'CTEQ    '.OR.LSTRUC.EQ.'CTEQ3L  ') ISTRUC=6
        IF(ISTRUC.EQ.0) GOTO 69
        WRITE(ITCOM,*) LSTRUC
        LOC(33)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NPOMERON
      IF(WORD.EQ.'NPOMERON') THEN
        WRITE(ITCOM,3) WORD
   93   PRINT*,' Give min and max; allowed range is 1 to 20.'
        READ*,MNPOM,MXPOM
        IF(MNPOM.LT.1.OR.MNPOM.GT.MXPOM.OR.MXPOM.GT.LIMPOM)
     1  GO TO 93
        WRITE(ITCOM,*)MNPOM,MXPOM
        LOC(34)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword FORCE
      IF(WORD.EQ.'FORCE   ') THEN
        PRINT*,' Force decay of a particle and its anti-particle.'
        PRINT*,' Use FORCE1 if you want to decay only one of them.'
        PRINT*,' Give a particle IDENT and IDENTs for decay mode.'
        PRINT*,' Max. number of decay products is 5.'
        NFTEMP=NFTEMP+2
        IF(NFTEMP.LE.MXFORC) THEN
          WRITE(ITCOM,3) WORD
          DO 72 K=1,5
  72      MFORCE(K,NFTEMP)=0
          READ*,IFORCE(NFTEMP),(MFORCE(K,NFTEMP),K=1,5)
          WRITE(ITCOM,*)IFORCE(NFTEMP),(MFORCE(K,NFTEMP),K=1,5)
          LOC(35)=NSEL
        ELSE
          PRINT*,' You exceeded the maximum allowed forced decays.'
        ENDIF
        GO TO 1000
      ENDIF
C
C          Keyword FORCE1
      IF(WORD.EQ.'FORCE1  ') THEN
        PRINT*,' Force decay of a particle.'
        PRINT*,' Give a particle IDENT and IDENTs for decay mode.'
        PRINT*,' Max. number of decay products is 5.'
        NFTEMP=NFTEMP+1
        IF(NFTEMP.LE.MXFORC) THEN
          WRITE(ITCOM,3) WORD
          DO 73 K=1,5
  73      MFORCE(K,NFTEMP)=0
          READ*,IFORCE(NFTEMP),(MFORCE(K,NFTEMP),K=1,5)
          WRITE(ITCOM,*)IFORCE(NFTEMP),(MFORCE(K,NFTEMP),K=1,5)
          LOC(35)=NSEL
        ELSE
          PRINT*,' You exceeded the maximun allowed forced decays.'
        ENDIF
        GO TO 1000
      ENDIF
C
C          Keyword HMASSES - also see HMASS
      IF(WORD.EQ.'HMASSES ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,
     $  ' Give Higgs masses (HIGGS,H10,H20,H30,H40,H1+,H2+,H1++,H2++)'
        CALL FLAVOR(80,I1,I2,I3,J1,INDEX)
        READ*,(AMLEP(INDEX+K),K=1,9)
        WRITE(ITCOM,*)(AMLEP(INDEX+K),K=1,9)
        LOC(36)=NSEL
        GO TO 1000
      ENDIF
C
C          Keywords WMODE1,WMODE2
      IF(WORD.EQ.'WMODE1  '.OR.WORD.EQ.'WMODE2  ') THEN
        WRITE(ITCOM,3) WORD
95      PRINT 6001
        PRINT 6003, LIST
        IF(WORD.EQ.'WMODE1  ') IJ=1
        IF(WORD.EQ.'WMODE2  ') IJ=2
        READ*,(WWTYP(K,IJ),K=1,25)
        DO 372 K=1,25
372     IF(WWTYP(K,IJ).NE.BLANK) NWWTYP(IJ)=NWWTYP(IJ)+1
        IF(SETTYP(0)) GOTO 95
        WRITE(ITCOM,4)(WWTYP(K,IJ),K=1,NWWTYP(IJ))
        IF(NWWTYP(IJ).LT.25) WRITE(ITCOM,*) '/'
        LOC(37)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NOEVOLVE
      IF(WORD.EQ.'NOEVOLVE') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571,NOEVOL
        WRITE(ITCOM,571) NOEVOL
        LOC(38)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NOHADRON
      IF(WORD.EQ.'NOHADRON') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571,NOEVOL
        WRITE(ITCOM,571) NOEVOL
        LOC(39)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword GAUGINO
      IF(WORD.EQ.'GAUGINO ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give masses (GLSS,GMSS,W+SS,W-SS)'
        CALL FLAVOR(29,I1,I2,I3,J1,IDG1)
        CALL FLAVOR(30,I1,I2,I3,J1,IDG2)
        CALL FLAVOR(39,I1,I2,I3,J1,IDG3)
        CALL FLAVOR(40,I1,I2,I3,J1,IDG4)
        READ*, AMLEP(IDG1),AMLEP(IDG2),AMLEP(IDG3),AMLEP(IDG4)
        WRITE(ITCOM,*) AMLEP(IDG1),AMLEP(IDG2),AMLEP(IDG3),AMLEP(IDG4)
        LOC(40)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SQUARK
      IF(WORD.EQ.'SQUARK  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give masses (UPSS,DNSS,STSS,CHSS,BTSS,TPSS).'
        CALL FLAVOR(21,I1,I2,I3,J1,IDXQK)
        READ*, (AMLEP(IDXQK+K-1),K=1,6)
        WRITE(ITCOM,*) (AMLEP(IDXQK+K-1),K=1,6)
        LOC(41)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SLEPTON
      IF(WORD.EQ.'SLEPTON ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,' Give masses (NUESS,E-SS,NUMSS,MU-SS,NUTSS,T-SS).'
        CALL FLAVOR(31,I1,I2,I3,J1,IDXLEP)
        READ*, (AMLEP(IDXLEP+K-1),K=1,6)
        WRITE(ITCOM,*) (AMLEP(IDXLEP+K-1),K=1,6)
        LOC(42)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NSIGMA
      IF(WORD.EQ.'NSIGMA  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,NSIGMA
        WRITE(ITCOM,*) NSIGMA
        LOC(43)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword XGENSS
      IF(WORD.EQ.'XGENSS  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*, XGENSS(9),(XGENSS(KK),KK=1,8)
        WRITE(ITCOM,*) XGENSS(9),(XGENSS(KK),KK=1,8)
        LOC(44)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword HMASS - just standard Higgs
      IF(WORD.EQ.'HMASS   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        CALL FLAVOR(81,I1,I2,I3,J1,INDEX)
        READ *, AMLEP(INDEX)
        WRITE(ITCOM,*) AMLEP(INDEX)
        LOC(45)=NSEL
        GO TO 1000
      ENDIF
C
C          Keywords WPMODE, WMMODE, Z0MODE
      IF(WORD.EQ.'WPMODE  '.OR.WORD.EQ.'WMMODE  '
     $.OR.WORD.EQ.'Z0MODE  ') THEN
        IF(WORD.EQ.'WPMODE  ') IJ=1
        IF(WORD.EQ.'WMMODE  ') IJ=2
        IF(WORD.EQ.'Z0MODE  ') IJ=3
        WRITE(ITCOM,3) WORD
        PRINT 6001
        PRINT 6003, LIST
        READ *, (WMODES(K,IJ),K=1,25)
        DO 463 K=1,25
463     IF(WMODES(K,IJ).NE.BLANK) NWMODE(IJ)=NWMODE(IJ)+1
        WRITE(ITCOM,*) (WMODES(K,IJ),K=1,NWMODE(IJ))
        LOC(46)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword WMASS
      IF(WORD.EQ.'WMASS   ') THEN
        PRINT 584
        READ*, AMW,AMZ
        WRITE(ITCOM,*) AMW,AMZ
        WMASS(1)=0.
        WMASS(2)=AMW
        WMASS(3)=AMW
        WMASS(4)=AMZ
        CALL FLAVOR(80,I1,I2,I3,J,INDEX)
        AMLEP(INDEX)=AMW
        CALL FLAVOR(90,I1,I2,I3,J,INDEX)
        AMLEP(INDEX)=AMZ
        LOC(47)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NEVOLVE
      IF(WORD.EQ.'NEVOLVE ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,NEVOLV
        WRITE(ITCOM,*) NEVOLV
        LOC(48)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword NHADRON
      IF(WORD.EQ.'NHADRON ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,NFRGMN
        WRITE(ITCOM,*)NFRGMN
        LOC(49)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword TCMASS
      IF(WORD.EQ.'TCMASS  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 584
        READ*,TCMRHO,TCGRHO
        WRITE(ITCOM,*) TCMRHO,TCGRHO
        LOC(50)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MSSMA
      IF(WORD.EQ.'MSSMA   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 585
        READ *, XGLSS,XMUSS,XHASS,XTBSS
        WRITE(ITCOM,*) XGLSS,XMUSS,XHASS,XTBSS
        LOC(51)=NSEL
        GOMSSM=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword MSSMB
      IF(WORD.EQ.'MSSMB   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 586
        READ*, XQ1SS,XDRSS,XURSS,XL1SS,XERSS
        WRITE(ITCOM,*) XQ1SS,XDRSS,XURSS,XL1SS,XERSS
        LOC(52)=NSEL
        GOMSSM=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword MSSMC
      IF(WORD.EQ.'MSSMC   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 587
        READ*, XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
        WRITE(ITCOM,*)XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
        LOC(53)=NSEL
        GOMSSM=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword PDFLIB: parameters for PDFLIB
#if defined(CERNLIB_PDFLIB)
      IF(WORD.EQ.'PDFLIB  ') THEN
        WRITE(ITCOM,3) WORD
        DO 541 I=1,20
          PDFPAR(I)='                    '
          PDFVAL(I)=0
541     CONTINUE
        PRINT 588
        READ*, (PDFPAR(I),PDFVAL(I),I=1,20)
        DO 542 I=1,20
          IF(PDFPAR(I).NE.'                    ') THEN
            WRITE(ITCOM,5401) PDFPAR(I)
5401        FORMAT(1H',A,1H')
            WRITE(ITCOM,*) PDFVAL(I)
          ENDIF
542     CONTINUE
        WRITE(ITCOM,*) '/'
        ISTRUC=-999
        LOC(54)=NSEL
        GO TO 1000
      ENDIF
#endif
C
C          Keyword SUGRA
      IF(WORD.EQ.'SUGRA   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 589
        READ*, XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
        WRITE(ITCOM,*) XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
        LOC(55)=NSEL
        GOMSSM=.TRUE.
        GOSUG=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword HTYPE
      IF(WORD.EQ.'HTYPE   ') THEN
        PRINT*,' Enter Higgs type (HL0, HH0, HA0) in single quotes'
        READ*, HTYPE
        IHTYPE=0
        IF(HTYPE.EQ.'HL0     ') IHTYPE=82
        IF(HTYPE.EQ.'HH0     ') IHTYPE=83
        IF(HTYPE.EQ.'HA0     ') IHTYPE=84
        IF(IHTYPE.EQ.0) THEN
          PRINT*, 'Invalid Higgs type - try again'
          GO TO 1000
        ENDIF
        WRITE(ITCOM,3) WORD
        WRITE(ITCOM,*) HTYPE
        LOC(56)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword EPOL
      IF(WORD.EQ.'EPOL    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 590
        READ*,PLEM,PLEP
        WRITE(ITCOM,*) PLEM,PLEP
        LOC(57)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MSSMD
      IF(WORD.EQ.'MSSMD   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 591
        READ*,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
        WRITE(ITCOM,*) XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
        LOC(58)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MSSME
      IF(WORD.EQ.'MSSME   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 592
        READ*, XM1SS,XM2SS
        WRITE(ITCOM,*) XM1SS,XM2SS
        LOC(59)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword GMSB
      IF(WORD.EQ.'GMSB    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 593
        READ*, XLAMGM,XMESGM,XN5GM,XTGBSU,XSMUSU,XCMGV
        WRITE(ITCOM,*) XLAMGM,XMESGM,XN5GM,XTGBSU,XSMUSU,XCMGV
        LOC(60)=NSEL
        GOMSSM=.TRUE.
        GOGMSB=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword NUSUG1: optional GUT scale gaugino masses
      IF(WORD.EQ.'NUSUG1   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 594
        READ*, XNUSUG(1),XNUSUG(2),XNUSUG(3)
        WRITE(ITCOM,*) XNUSUG(1),XNUSUG(2),XNUSUG(3)
        LOC(61)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG2: optional GUT scale A terms
      IF(WORD.EQ.'NUSUG2   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 595
        READ*, XNUSUG(6),XNUSUG(5),XNUSUG(4)
        WRITE(ITCOM,*) XNUSUG(6),XNUSUG(5),XNUSUG(4)
        LOC(62)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG3: optional GUT scale Higgs masses
      IF(WORD.EQ.'NUSUG3   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 596
        READ*, XNUSUG(7),XNUSUG(8)
        WRITE(ITCOM,*) XNUSUG(7),XNUSUG(8)
        LOC(63)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG4: optional GUT scale 1st/2nd gen. masses
      IF(WORD.EQ.'NUSUG4   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 597
        READ*, XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10)
     $,XNUSUG(9)
        WRITE(ITCOM,*) XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10)
     $,XNUSUG(9)
        LOC(64)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG5: optional GUT scale 3rd gen. masses
      IF(WORD.EQ.'NUSUG5   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 598
        READ*, XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15)
     $,XNUSUG(14)
        WRITE(ITCOM,*) XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15)
     $,XNUSUG(14)
        LOC(65)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOGRAV: No gravitino decays
      IF(WORD.EQ.'NOGRAV  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 572
        READ 571,NOGRAV
        WRITE(ITCOM,571) NOGRAV
        LOC(66)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MGVTNO: Sets the gravitino mass
      IF(WORD.EQ.'MGVTNO  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 599
        READ*, XMGVTO
        WRITE(ITCOM,*) XMGVTO
        LOC(67)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword AL3UNI: Impose alpha_s unification at M_GUT
      IF(WORD.EQ.'AL3UNI  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 600
        READ*, AL3UNI
        WRITE(ITCOM,*) AL3UNI
        LOC(68)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword GMSB2: additional GMSB parameters
      IF(WORD.EQ.'GMSB2   ') THEN
       WRITE(ITCOM,3) WORD
       PRINT 601
       READ*, XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
       WRITE(ITCOM,*) XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
       LOC(69)=NSEL
       GO TO 1000
      ENDIF
C
C          Keyword EEBREM: do bremsstralung
      IF(WORD.EQ.'EEBREM  ') THEN
       WRITE(ITCOM,3) WORD
       PRINT 602
       READ*, RSHMIN,RSHMAX
       WRITE(ITCOM,*) RSHMIN,RSHMAX
       IBREM=.TRUE.
       LOC(70)=NSEL
       GO TO 1000
      ENDIF
C
C          Keyword EEBEAM: do beamstralung
      IF(WORD.EQ.'EEBEAM   ') THEN
       WRITE(ITCOM,3) WORD
       PRINT 603
       READ*, RSHMIN,RSHMAX,UPSLON,SIGZ
       WRITE(ITCOM,*) RSHMIN,RSHMAX,UPSLON,SIGZ
       IBREM=.TRUE.
       IBEAM=.TRUE.
       LOC(71)=NSEL
       GO TO 1000
      ENDIF
C
C          Keyword QMKKG (QMW for EXTRADIM)
      IF(WORD.EQ.'QMKKG   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,'KK graviton mass limits?'
        READ*, QMIN,QMAX
        WRITE(ITCOM,*) QMIN,QMAX
        LOC(19)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword QTKKG (QTW for EXTRADIM)
      IF(WORD.EQ.'QTKKG   ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,'KK graviton qt limits?'
        READ*, QTMIN,QTMAX
        WRITE(ITCOM,*) QTMIN,QTMAX
        LOC(20)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword EXTRAD for EXTRADIM
      IF(WORD.EQ.'EXTRAD  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,'KK delta n, mass, UVCUT flag?'
        READ*, NEXTRAD,MASSD,UVCUT
        WRITE(ITCOM,*) NEXTRAD,MASSD,UVCUT
        LOC(72)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MIJLIM
      IF(WORD.EQ.'MIJLIM  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,'Dijet pair mass limits i,j,min,max?'
        READ*, JLIM1,JLIM2,AMLIM1,AMLIM2
        WRITE(ITCOM,*) JLIM1,JLIM2,AMLIM1,AMLIM2
        LOC(73)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword MTOT
      IF(WORD.EQ.'MTOT    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT*,'Total mass limits?'
        READ*, EHMGMN,EHMGMX
        WRITE(ITCOM,*) EHMGMN,EHMGMX
        LOC(74)=NSEL
        GO TO 1000
      ENDIF
C
C          Keyword SUGRHN: enter right-neutrino parameters
      IF(WORD.EQ.'SUGRHN   ') THEN
       WRITE(ITCOM,3) WORD
       PRINT 604
       READ*, XMN3NR,XMAJNR,XANSS,XNRSS
       WRITE(ITCOM,*) XMN3NR,XMAJNR,XANSS,XNRSS
       LOC(75)=NSEL
       GO TO 1000
      ENDIF
C
C          Keyword AMSB
      IF(WORD.EQ.'AMSB    ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 605
        READ*, XM0SU,XMHSU,XTGBSU,XSMUSU
        WRITE(ITCOM,*) XM0SU,XMHSU,XTGBSU,XSMUSU
        LOC(76)=NSEL
        GOMSSM=.TRUE.
        GOSUG=.TRUE.
        GOAMSB=.TRUE.
        GO TO 1000
      ENDIF
C
C          Keyword SSBCSC
      IF(WORD.EQ.'SSBCSC  ') THEN
        WRITE(ITCOM,3) WORD
        PRINT 606
        READ*, XSBCS
        WRITE(ITCOM,*) XSBCS
        LOC(77)=NSEL
        GO TO 1000
      ENDIF
C
C          None of the above
C
      PRINT 99
      GO TO 1000
C
  562 FORMAT(' Give jet limits: min1,max1,min2,max2, etc.',
     1' Terminate list with a /.',/,
     2' To fix a variable give only minimum value.',/)
  563 FORMAT(' Allowed range is ',F9.2,' - ',F9.2)
  572 FORMAT('  T or F?')
  584 FORMAT(' Value?')
  585 FORMAT(' Give M(gluino), mu, M(ha), tan(beta)')
  586 FORMAT(' Give M(ul), M(dr), M(ur), M(el), M(er), ')
  587 FORMAT(' Give M(tl), M(br), M(tr), M(taul), M(taur), A(t),',
     $' A(b), A(tau)')
  588 FORMAT(' Give PDFLIB keyword, value, ...')
  589 FORMAT(' Give M_0, M_(1/2), A_0, tan(beta), sgn(mu)')
  590 FORMAT(' Give polarizations P_L(e-), P_L(e+)')
  591 FORMAT(' Give optional M(cl), M(sr), M(cr), M(mul), M(mur), ')
  592 FORMAT(' Give optional M_1, M_2')
  593 FORMAT(' Give LAM, M_MES, XN5, tan(beta), sgn(mu), C_mgv')
  594 FORMAT(' Give optional GUT scale M_1, M_2, M_3')
  595 FORMAT(' Give optional GUT scale A_t, A_b, A_tau')
  596 FORMAT(' Give optional GUT scale m_Hd, m_Hu')
  597 FORMAT(' Give optional M(ul), M(dr), M(ur), M(el), M(er), ')
  598 FORMAT(' Give optional M(tl), M(br), M(tr), M(Ll), M(Lr), ')
  599 FORMAT(' Enter optional m(gravitino), ')
  600 FORMAT(' Enter .TRUE. if alpha_s unification is wanted:, ')
  601 FORMAT(' Enter Rsl,dmH_d,dmH_u,RDY,n5_1,n5_2,n5_3: ')
  602 FORMAT(' Enter RSHMIN,RSHMAX for bremsstrahlung: ')
  603 FORMAT(' Enter RSHMIN,RSHMAX,UPSLON,SIGZ for beamstrahlung: ')
  604 FORMAT(' Enter M_N3, M_MAJ, A_N, M_NRSS for RHN model: ')
  605 FORMAT(' Give M_0, M_(3/2), tan(beta), sgn(mu)')
  606 FORMAT(' Enter Q_max= max scale choice for SUSY BCs')
C
C          CHECK FOR ERRORS AND OVERLAPPING VARIABLES.
C          PRINT ERROR MESSAGE OR WARNINGS
  101 CONTINUE
      IF(LOC(2)*LOC(5).NE.0) PRINT 2001
      IF(LOC(2)*LOC(6).NE.0) PRINT 2001
      IF(LOC(15).NE.0.OR.LOC(37).NE.0) DUMY=SETTYP(0)
      DO 120 I=1,17
      I1=NLAP(1,I)
      I2=NLAP(2,I)
      I3=NLAP(3,I)
      IF(LOC(I1)*LOC(I2)*LOC(I3).NE.0) PRINT 1001
  120 CONTINUE
 1001 FORMAT(//,2X,'YOU HAVE GIVEN LIMITS FOR AN OVERLAPPING SET OF',
     C  ' VARIABLES. SET MINIMIZING PPERP INTERVAL WILL BE USED')
 2001 FORMAT(//,2X,' YOU CANNOT GIVE LIMITS FOR BOTH THETA AND Y. MAKE',
     C  ' UP YOUR MIND. JOB TERMINATED')
 2003 FORMAT(/,'  YOU CANNOT RUN WS AND Z0 OR GAMMAS AT THE SAME TIME')
 2004 FORMAT(//,'  PARAMETER OUT OF RANGE. ')
      GOTO 1111
  999 REWIND ITCOM
      RETURN
      END
