#include "isajet/pilot.h"
      SUBROUTINE READIN(IFL)
C
C       Read in user data and execute SETTYP if appropriate values 
C       are set. IFL return values:
C       IFL = 0       Good parameter set
C       IFL = 1001    Stop
C       IFL > 0       Error. Program will continue reading data but
C                     will exit when END or unrecognizable keyword
C                     is found.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/mbgen.inc"
#include "isajet/force.inc"
#include "isajet/dkytab.inc"
#include "isajet/qcdpar.inc"
#include "isajet/eepar.inc"
#include "isajet/idrun.inc"
#include "isajet/frgpar.inc"
#include "isajet/keys.inc"
#include "isajet/kkgrav.inc"
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
#include "isajet/mglims.inc"
C
      LOGICAL SETTYP,DUMY
      CHARACTER*8 TTL(10),WORD,LSTRUC,BLANK
      CHARACTER*8 WTYP(4),RDID(2)
      CHARACTER*40 V,VISAJE
      INTEGER NLAP(3,17)
      INTEGER IDANTI,ID,IDB
      INTEGER IFL,I1,I2,I3,J1,I,IKEY,IJ,J,KK,IDABS
      INTEGER IDXQKL,IDXQKR
      INTEGER NSEL,K,KFORCE(5),INDEX,IDG1,IDG2,IDG3,IDG4,IDXLEP
      REAL AMW,AMZ
      CHARACTER*8 HTYPE
      INTEGER JLIM1,JLIM2
      REAL AMLIM1,AMLIM2
#if defined(CERNLIB_SINGLE)
      REAL    SEED
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION SEED
#endif
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
C
C          Overlapping variable flags.
      DATA NLAP/1,2,3, 1,2,7 ,1,2,8, 1,3,5, 1,3,6, 1,3,7, 1,3,8, 1,5,7,
     X          1,5,8, 1,6,7, 1,6,8, 2,3,7, 2,3,8, 3,5,7, 3,6,7, 3,5,8,
     X          3,6,8/
      DATA BLANK/'        '/
C
C          Entry
      IFL=0
      V=VISAJE()
      WRITE(ITLIS,10) V
10    FORMAT('1',//5X,'*****  ',A40,'  *****')
      WRITE(ITLIS,11)
   11 FORMAT(////30X,' COMMANDS READ BY READIN')
C
C          Read title
C
      READ(ITCOM,1) TTL
    1 FORMAT(10A8)
      WRITE(ITLIS,2) TTL
    2 FORMAT(' ',10A8)
      IF(TTL(1).EQ.'STOP    ') THEN
        IFL=1001
        RETURN
      ENDIF
C
C          Read energy and no. of events
C
      READ(ITCOM,*) ECM,NEVENT,NEVPRT,NJUMP
      WRITE(ITLIS,*) ECM,NEVENT,NEVPRT,NJUMP
C
C          Reset all variables and set process if title is not 'SAME'
C
      IF(TTL(1).NE.'SAME    ') THEN
        DO 20 I=1,10
   20   TITLE(I)=TTL(I)
        CALL RESET
        KEYON=.FALSE.
C          Read reaction
        READ(ITCOM,3) REAC
    3   FORMAT(A8)
        WRITE(ITLIS,4) REAC
   4    FORMAT(1X,A8)
        DO 18 I=1,MXKEYS
18      KEYS(I)=.FALSE.
        KEYON=.FALSE.
C          Set KEYS and NJET
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
        ELSEIF(REAC.EQ.'EXTRADIM') THEN
          KEYS(11)=.TRUE.
          IKEY=11
          NJET=3
        ELSEIF(REAC.EQ.'ZJJ     ') THEN
          KEYS(12)=.TRUE.
          IKEY=12
          NJET=3
        ELSE
          KEYON=.FALSE.
  890     WRITE(ITLIS,1999)
          IFL=9
          RETURN
        ENDIF
      ENDIF
C
      SCM=ECM**2
      HALFE=ECM/2
      NSEL=0
C
C          Read keyword. For each recognized keyword read corresponding
C          variables and set LOC flag.
C
      NSEL=0
100   CONTINUE
      READ(ITCOM,3) WORD
      WRITE(ITLIS,4) WORD
      NSEL=NSEL+1
C
C          Keyword END
      IF(WORD.EQ.'END     ') THEN
C          Check for previous error
        IF(IFL.NE.0) RETURN
C          Check inconsistent limits
        IF(LOC(2)*LOC(5).NE.0.OR.LOC(2)*LOC(6).NE.0) THEN
          WRITE(ITLIS,2001)
          IFL=11
        ENDIF
C          Set and check jet types
        IF(LOC(15).NE.0.OR.LOC(37).NE.0.OR.LOC(46).NE.0) THEN
          IF(SETTYP(0)) THEN
            WRITE(ITLIS,2006)
            IFL=12
          ENDIF
        ENDIF
C          Check MSSM/SUGRA conflict
        IF((LOC(51).NE.0.OR.LOC(52).NE.0.OR.LOC(53).NE.0).AND.
     $  LOC(55).NE.0) THEN
          WRITE(ITLIS,2007)
          IFL=29
        ENDIF
C          Check overlapping limits
        DO 120 I=1,17
          I1=NLAP(1,I)
          I2=NLAP(2,I)
          I3=NLAP(3,I)
          IF(LOC(I1)*LOC(I2)*LOC(I3).NE.0) WRITE(ITLIS,1001)
  120   CONTINUE
C          Setup PDFLIB
#if defined(CERNLIB_PDFLIB)
        IF(ISTRUC.EQ.-999) THEN
          WRITE(ITLIS,1200)
1200      FORMAT(//
     $    '1********************************'/
     $    ' *                              *'/
     $    ' * INITIALIZE PDFLIB FOR ISAJET *'/
     $    ' *                              *'/
     $    ' ********************************'/)
          N6=ITLIS
          IFLPRT=2
          CALL PDFSET(PDFPAR,PDFVAL)
          CALL PFTOPDG(0.5D0,1.0D2,DXPDF)
          IFLPRT=0
        ENDIF
#endif
C          Check EXTRADIM parameters are set
        IF(KEYS(11).AND.LOC(72).EQ.0) THEN
          WRITE(ITLIS,*) 'YOU FORGOT TO SET EXTRAD PARAMETERS'
          IFL=72
        ENDIF  
C
        RETURN
      ENDIF
C
C          Keyword P
      IF(WORD.EQ.'P       ') THEN
        READ(ITCOM,*)  (PMIN(K),PMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (PMIN(K),PMAX(K),K=1,NJET)
        LOC(1)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword Y
      IF(WORD.EQ.'Y       ') THEN
        READ(ITCOM,*)  (YJMIN(K),YJMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (YJMIN(K),YJMAX(K),K=1,NJET)
        LOC(2)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword X
      IF(WORD.EQ.'X       ') THEN
        READ(ITCOM,*)  (XJMIN(K),XJMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (XJMIN(K),XJMAX(K),K=1,NJET)
        LOC(3)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword PHI
      IF(WORD.EQ.'PHI     ') THEN
        READ(ITCOM,*)  (PHIMIN(K),PHIMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (PHIMIN(K),PHIMAX(K),K=1,NJET)
        LOC(4)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword TH
      IF(WORD.EQ.'TH      '.OR.WORD.EQ.'THETA   ') THEN
        READ(ITCOM,*)  (THMIN(K),THMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (THMIN(K),THMAX(K),K=1,NJET)
        LOC(5)=NSEL
        LOC(6)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword PT
      IF(WORD.EQ.'PT      '.OR.WORD.EQ.'PPERP   ') THEN
        READ(ITCOM,*)  (PTMIN(K),PTMAX(K),K=1,NJET)
        WRITE(ITLIS,*) (PTMIN(K),PTMAX(K),K=1,NJET)
        LOC(7)=NSEL
        LOC(8)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NODECAY
      IF(WORD.EQ.'NODECAY ') THEN
        READ(ITCOM,571) NODCAY
571     FORMAT(L1)
        WRITE(ITLIS,572) NODCAY
572     FORMAT(' ',L1)
        LOC(9)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOETA
      IF(WORD.EQ.'NOETA   ') THEN
        READ(ITCOM,571) NOETA
        WRITE(ITLIS,572) NOETA
        LOC(10)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOPI0
      IF(WORD.EQ.'NOPI0   ') THEN
        READ(ITCOM,571) NOPI0
        WRITE(ITLIS,572) NOPI0
        LOC(11)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword BEAMS
      IF(WORD.EQ.'BEAMS   ') THEN
        READ(ITCOM,*) RDID(1),RDID(2)
        WRITE(ITLIS,*) RDID(1),RDID(2)
        IDIN(1)=0
        IDIN(2)=0
        DO 123 K=1,2
          IF(RDID(K).EQ.'P       ') IDIN(K)=+1120
          IF(RDID(K).EQ.'AP      ') IDIN(K)=-1120
          IF(RDID(K).EQ.'N       ') IDIN(K)=+1220
          IF(RDID(K).EQ.'AN      ') IDIN(K)=-1220
123     CONTINUE
        IF(IDIN(1)*IDIN(2).EQ.0) THEN
          WRITE(ITLIS,2002)
          IFL=13
        ENDIF
        LOC(12)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword FRAGMENT
      IF(WORD.EQ.'FRAGMENT') THEN
        READ(ITCOM,*)  FRPAR
        WRITE(ITLIS,*)  FRPAR
        LOC(13)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SEED
      IF(WORD.EQ.'SEED    ') THEN
        READ(ITCOM,*) SEED
        WRITE(ITLIS,*) SEED
        CALL RANFST(SEED)
        WRITE(XSEED,'(E24.15)') SEED
        LOC(14)=NSEL
        GO TO 100
      ENDIF
C
C          Keywords JETTYPE1, JETTYPE2, JETTYPE3, ...
C          (Yes, this is ugly)
      IF(WORD.EQ.'JETTYPE1'.OR.WORD.EQ.'JETTYPE2'.OR.
     $WORD.EQ.'JETTYPE3'.OR.WORD.EQ.'JETTYPE4'.OR.
     $WORD.EQ.'JETTYPE5'.OR.WORD.EQ.'JETTYPE6'.OR.
     $WORD.EQ.'JETTYPE7'.OR.WORD.EQ.'JETTYPE8') THEN
        IF(WORD.EQ.'JETTYPE1') IJ=1
        IF(WORD.EQ.'JETTYPE2') IJ=2
        IF(WORD.EQ.'JETTYPE3') IJ=3
        IF(WORD.EQ.'JETTYPE4') IJ=4
        IF(WORD.EQ.'JETTYPE5') IJ=5
        IF(WORD.EQ.'JETTYPE6') IJ=6
        IF(WORD.EQ.'JETTYPE7') IJ=7
        IF(WORD.EQ.'JETTYPE8') IJ=8
        DO 151 K=1,30
          JETYP(K,IJ)=BLANK
151     CONTINUE
        READ(ITCOM,*) (JETYP(K,IJ),K=1,30)
        DO 152 K=1,25
152     IF(JETYP(K,IJ).NE.BLANK) NJTTYP(IJ)=NJTTYP(IJ)+1
        WRITE(ITLIS,*) (JETYP(K,IJ),K=1,NJTTYP(IJ))
        LOC(15)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SIN2W
      IF(WORD.EQ.'SIN2W   ') THEN
        READ(ITCOM,*)  SIN2W
        WRITE(ITLIS,*) SIN2W
        LOC(17)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword TMASS
      IF(WORD.EQ.'TMASS   ') THEN
        READ(ITCOM,*)  AMLEP(6),AMLEP(7),AMLEP(8)
        WRITE(ITLIS,*) AMLEP(6),AMLEP(7),AMLEP(8)
        LOC(18)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword QMH
      IF(WORD.EQ.'QMH     ') THEN
        READ(ITCOM,*) QMIN,QMAX
        WRITE(ITLIS,*) QMIN,QMAX
        LOC(19)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword QMW
      IF(WORD.EQ.'QMW     ') THEN
        READ(ITCOM,*) QMIN,QMAX
        WRITE(ITLIS,*) QMIN,QMAX
        LOC(19)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword QTW
      IF(WORD.EQ.'QTW     ') THEN
        READ(ITCOM,*)  QTMIN,QTMAX
        WRITE(ITLIS,*) QTMIN,QTMAX
        LOC(20)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword YW
      IF(WORD.EQ.'YW      ') THEN
        READ(ITCOM,*)  YWMIN,YWMAX
        WRITE(ITLIS,*) YWMIN,YWMAX
        LOC(21)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword XW
      IF(WORD.EQ.'XW      ') THEN
        READ(ITCOM,*)  XWMIN,XWMAX
        WRITE(ITLIS,*) XWMIN,XWMAX
        LOC(22)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword THW
      IF(WORD.EQ.'THW     ') THEN
        READ(ITCOM,*)  THWMIN,THWMAX
        WRITE(ITLIS,*) THWMIN,THWMAX
        LOC(23)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword PHIW
      IF(WORD.EQ.'PHIW    ') THEN
        READ(ITCOM,*)  PHWMIN,PHWMAX
        WRITE(ITLIS,*) PHWMIN,PHWMAX
        LOC(24)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NONUNU
      IF(WORD.EQ.'NONUNU  ') THEN
        READ(ITCOM,571) NONUNU
        WRITE(ITLIS,572) NONUNU
        LOC(25)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword WTYPE
      IF(WORD.EQ.'WTYPE   ') THEN
        DO 261 J=1,4
          WTYP(J)=BLANK
          GODY(J)=.FALSE.
261     CONTINUE
        READ(ITCOM,*) WTYP
        WRITE(ITLIS,*) WTYP
        DO 262 K=1,4
          IF(WTYP(K).EQ.'GM      ') GODY(1)=.TRUE.
          IF(WTYP(K).EQ.'W+      ') GODY(2)=.TRUE.
          IF(WTYP(K).EQ.'W-      ') GODY(3)=.TRUE.
          IF(WTYP(K).EQ.'Z0      ') GODY(4)=.TRUE.
  262   CONTINUE
        IF(GODY(1)) JWTYP=1
        IF(GODY(2).OR.GODY(3)) JWTYP=3
        IF(GODY(4)) JWTYP=4
        IF((GODY(2).OR.GODY(3)).AND.(GODY(1).OR.GODY(4))) THEN
          WRITE(ITLIS,2003)
          IFL=13
        ENDIF
        LOC(26)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword LAMBDA
      IF(WORD.EQ.'LAMBDA  ') THEN
        READ(ITCOM,*)  ALAM
        WRITE(ITLIS,*)  ALAM
        ALAM2=ALAM**2
        LOC(27)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NTRIES
      IF(WORD.EQ.'NTRIES  ') THEN
        READ(ITCOM,*) NTRIES
        WRITE(ITLIS,*) NTRIES
        LOC(28)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword CUTOFF
      IF(WORD.EQ.'CUTOFF  ') THEN
        READ(ITCOM,*) CUTOFF,CUTPOW
        WRITE(ITLIS,*) CUTOFF,CUTPOW
        LOC(29)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword XGEN
      IF(WORD.EQ.'XGEN    ') THEN
        READ(ITCOM,*) XGEN
        WRITE(ITLIS,*) XGEN
        LOC(30)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SIGQT
      IF(WORD.EQ.'SIGQT   ') THEN
        READ(ITCOM,*) SIGQT
        WRITE(ITLIS,*) SIGQT
        LOC(31)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword CUTJET
      IF(WORD.EQ.'CUTJET  ') THEN
        READ(ITCOM,*)  CUTJET
        WRITE(ITLIS,*)  CUTJET
        LOC(32)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword WFUDGE
      IF(WORD.EQ.'WFUDGE  ') THEN
        READ(ITCOM,*)  WFUDGE
        WRITE(ITLIS,*)  WFUDGE
        LOC(50)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword STRUC
      IF(WORD.EQ.'STRUC   ') THEN
        ISTRUC=0
        READ(ITCOM,*) LSTRUC
        WRITE(ITLIS,4) LSTRUC
        IF(LSTRUC.EQ.'OWENS   ') ISTRUC=1
        IF(LSTRUC.EQ.'BAIER   ') ISTRUC=2
        IF(LSTRUC.EQ.'EICHTEN '.OR.LSTRUC.EQ.'EHLQ    ') ISTRUC=3
        IF(LSTRUC.EQ.'DUKE    '.OR.LSTRUC.EQ.'DO      ') ISTRUC=4
        IF(LSTRUC.EQ.'CTEQ2L  ') ISTRUC=5
        IF(LSTRUC.EQ.'CTEQ    '.OR.LSTRUC.EQ.'CTEQ3L  ') ISTRUC=6
        IF(ISTRUC.EQ.0) THEN
          WRITE(ITLIS,2002)
          IFL=13
        ENDIF
        LOC(33)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NPOMERON
      IF(WORD.EQ.'NPOMERON') THEN
        READ(ITCOM,*) MNPOM,MXPOM
        WRITE(ITLIS,*) MNPOM,MXPOM
        IF(MNPOM.LT.1.OR.MNPOM.GT.MXPOM.OR.MXPOM.GT.LIMPOM) THEN
          WRITE(ITLIS,2004)
          IFL=14
        ENDIF
        LOC(34)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword FORCE
      IF(WORD.EQ.'FORCE   ') THEN
        NFORCE=NFORCE+1
        IF(NFORCE.GT.MXFORC-1) THEN
          WRITE(ITLIS,2004)
          IFL=14
        ENDIF
        DO 351 K=1,5
351     KFORCE(K)=0
        READ(ITCOM,*) IFORCE(NFORCE),(KFORCE(K),K=1,5)
        CALL ORDER(IFORCE(NFORCE),KFORCE,MFORCE(1,NFORCE),
     $  MEFORC(NFORCE))
        WRITE(ITLIS,*) IFORCE(NFORCE),(MFORCE(K,NFORCE),K=1,5)
        ID=IFORCE(NFORCE)
        IDABS=IABS(ID)
        IF(IDABS.LT.6) THEN
          WRITE(ITLIS,2005)
          IFL=15
        ENDIF
        IDB=IDANTI(ID)
        IF(IDB.NE.ID) THEN
          IFORCE(NFORCE+1)=IDB
          DO 352 K=1,5
  352     MFORCE(K,NFORCE+1)=IDANTI(MFORCE(K,NFORCE))
          NFORCE=NFORCE+1
        ENDIF
        LOC(35)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword FORCE1
      IF(WORD.EQ.'FORCE1  ') THEN
        NFORCE=NFORCE+1
        IF(NFORCE.GT.MXFORC) THEN
          WRITE(ITLIS,2004)
          IFL=14
        ENDIF
        DO 353 K=1,5
353     KFORCE(K)=0
        READ(ITCOM,*) IFORCE(NFORCE),(KFORCE(K),K=1,5)
        CALL ORDER(IFORCE(NFORCE),KFORCE,MFORCE(1,NFORCE),
     $  MEFORC(NFORCE))
        WRITE(ITLIS,*) IFORCE(NFORCE),(MFORCE(K,NFORCE),K=1,5)
        IF(IABS(IFORCE(NFORCE)).LT.6) THEN
          WRITE(ITLIS,2005)
          IFL=15
        ENDIF
        LOC(35)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword HMASSES - also see HMASS
      IF(WORD.EQ.'HMASSES ') THEN
        CALL FLAVOR(80,I1,I2,I3,J1,INDEX)
        READ(ITCOM,*) (AMLEP(INDEX+K),K=1,9)
        WRITE(ITLIS,*) (AMLEP(INDEX+K),K=1,9)
        LOC(36)=NSEL
        GO TO 100
      ENDIF
C
C          Keywords WMODE1,WMODE2
      IF(WORD.EQ.'WMODE1  '.OR.WORD.EQ.'WMODE2  ') THEN
        IF(WORD.EQ.'WMODE1  ') IJ=1
        IF(WORD.EQ.'WMODE2  ') IJ=2
        READ(ITCOM,*) (WWTYP(K,IJ),K=1,25)
        DO 372 K=1,25
372     IF(WWTYP(K,IJ).NE.BLANK) NWWTYP(IJ)=NWWTYP(IJ)+1
        WRITE(ITLIS,*) (WWTYP(K,IJ),K=1,NWWTYP(IJ))
        LOC(37)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOEVOLVE
      IF(WORD.EQ.'NOEVOLVE') THEN
        READ (ITCOM,571) NOEVOL
        WRITE(ITLIS,572) NOEVOL
        LOC(38)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOHADRON
      IF(WORD.EQ.'NOHADRON') THEN
        READ (ITCOM,571) NOHADR
        WRITE(ITLIS,572) NOHADR
        LOC(39)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword GAUGINO
      IF(WORD.EQ.'GAUGINO ') THEN
        CALL FLAVOR(29,I1,I2,I3,J1,IDG1)
        CALL FLAVOR(30,I1,I2,I3,J1,IDG2)
        CALL FLAVOR(39,I1,I2,I3,J1,IDG3)
        CALL FLAVOR(40,I1,I2,I3,J1,IDG4)
        READ(ITCOM,*) AMLEP(IDG1),AMLEP(IDG2),AMLEP(IDG3),AMLEP(IDG4)
        WRITE(ITLIS,*) AMLEP(IDG1),AMLEP(IDG2),AMLEP(IDG3),AMLEP(IDG4)
        LOC(40)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SQUARK
      IF(WORD.EQ.'SQUARK  ') THEN
        CALL FLAVOR(21,I1,I2,I3,J1,IDXQKL)
        READ(ITCOM,*) (AMLEP(IDXQKL+K-1),K=1,6)
        WRITE(ITLIS,*) (AMLEP(IDXQKL+K-1),K=1,6)
        CALL FLAVOR(41,I1,I2,I3,J1,IDXQKR)
        DO 411 K=1,6
          AMLEP(IDXQKR+K-1)=AMLEP(IDXQKL+K-1)
411     CONTINUE
        LOC(41)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SLEPTON
      IF(WORD.EQ.'SLEPTON ') THEN
        CALL FLAVOR(31,I1,I2,I3,J1,IDXLEP)
        READ(ITCOM,*) (AMLEP(IDXLEP+K-1),K=1,6)
        WRITE(ITLIS,*) (AMLEP(IDXLEP+K-1),K=1,6)
        LOC(42)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NSIGMA
      IF(WORD.EQ.'NSIGMA  ') THEN
        READ(ITCOM,*) NSIGMA
        WRITE(ITLIS,*) NSIGMA
        LOC(43)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword XGENSS
      IF(WORD.EQ.'XGENSS  ') THEN
        READ(ITCOM,*) XGENSS(9),(XGENSS(KK),KK=1,8)
        WRITE(ITLIS,*) XGENSS(9),(XGENSS(KK),KK=1,8)
        LOC(44)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword HMASS - just standard Higgs
      IF(WORD.EQ.'HMASS   ') THEN
        CALL FLAVOR(81,I1,I2,I3,J1,INDEX)
        READ(ITCOM,*) AMLEP(INDEX)
        WRITE(ITLIS,*) AMLEP(INDEX)
        LOC(45)=NSEL
        GO TO 100
      ENDIF
C
C          Keywords WPMODE, WMMODE, Z0MODE
      IF(WORD.EQ.'WPMODE  '.OR.WORD.EQ.'WMMODE  '
     $.OR.WORD.EQ.'Z0MODE  ') THEN
        IF(WORD.EQ.'WPMODE  ') IJ=1
        IF(WORD.EQ.'WMMODE  ') IJ=2
        IF(WORD.EQ.'Z0MODE  ') IJ=3
        READ(ITCOM,*) (WMODES(K,IJ),K=1,25)
        DO 463 K=1,25
463     IF(WMODES(K,IJ).NE.BLANK) NWMODE(IJ)=NWMODE(IJ)+1
        WRITE(ITLIS,*) (WMODES(K,IJ),K=1,NWMODE(IJ))
        LOC(46)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword WMASS
      IF(WORD.EQ.'WMASS   ') THEN
        READ(ITCOM,*) AMW,AMZ
        WRITE(ITLIS,*) AMW,AMZ
        WMASS(1)=0.
        WMASS(2)=AMW
        WMASS(3)=AMW
        WMASS(4)=AMZ
        CALL FLAVOR(80,I1,I2,I3,J,INDEX)
        AMLEP(INDEX)=AMW
        CALL FLAVOR(90,I1,I2,I3,J,INDEX)
        AMLEP(INDEX)=AMZ
        LOC(47)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NEVOLVE
      IF(WORD.EQ.'NEVOLVE ') THEN
        READ(ITCOM,*) NEVOLV
        WRITE(ITLIS,*) NEVOLV
        LOC(48)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NHADRON
      IF(WORD.EQ.'NHADRON ') THEN
        READ(ITCOM,*) NFRGMN
        WRITE(ITLIS,*) NFRGMN
        LOC(49)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword TCMASS
      IF(WORD.EQ.'TCMASS  ') THEN
        READ(ITCOM,*) TCMRHO,TCGRHO
        WRITE(ITLIS,*) TCMRHO,TCGRHO
        LOC(50)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MSSMA: gluino, mu, mha, tanb
      IF(WORD.EQ.'MSSMA   ') THEN
        READ(ITCOM,*) XGLSS,XMUSS,XHASS,XTBSS
        WRITE(ITLIS,*) XGLSS,XMUSS,XHASS,XTBSS
        GOMSSM=.TRUE.
        LOC(51)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MSSMB: 1st generation soft terms
      IF(WORD.EQ.'MSSMB   ') THEN
        READ(ITCOM,*) XQ1SS,XDRSS,XURSS,XL1SS,XERSS
        WRITE(ITLIS,*) XQ1SS,XDRSS,XURSS,XL1SS,XERSS
        LOC(52)=NSEL
        GOMSSM=.TRUE.
        GO TO 100
      ENDIF
C
C          Keyword MSSMC: 3rd generation soft terms
      IF(WORD.EQ.'MSSMC   ') THEN
        READ(ITCOM,*)XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
       WRITE(ITLIS,*)XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
        LOC(53)=NSEL
        GOMSSM=.TRUE.
        GO TO 100
      ENDIF
C
C          Keyword PDFLIB: parameters for PDFLIB
#if defined(CERNLIB_PDFLIB)
      IF(WORD.EQ.'PDFLIB  ') THEN
        DO 541 I=1,20
          PDFPAR(I)='                    '
          PDFVAL(I)=0
541     CONTINUE
        READ(ITCOM,*) (PDFPAR(I),PDFVAL(I),I=1,20)
        DO 542 I=1,20
          IF(PDFPAR(I).NE.'                    ') THEN
            WRITE(ITLIS,*) PDFPAR(I),PDFVAL(I)
          ENDIF
542     CONTINUE
        ISTRUC=-999
        LOC(54)=NSEL
        GO TO 100
      ENDIF
#endif
C
C          Keyword SUGRA
      IF(WORD.EQ.'SUGRA   ') THEN
        READ(ITCOM,*) XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
        WRITE(ITLIS,*) XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
        LOC(55)=NSEL
        GOMSSM=.TRUE.
        GOSUG=.TRUE.
        GO TO 100
      ENDIF
C
C          Keyword HTYPE
      IF(WORD.EQ.'HTYPE   ') THEN
        READ(ITCOM,*) HTYPE
        WRITE(ITLIS,*) HTYPE
        LOC(56)=NSEL
        IHTYPE=0
        IF(HTYPE.EQ.'HL0     ') IHTYPE=82
        IF(HTYPE.EQ.'HH0     ') IHTYPE=83
        IF(HTYPE.EQ.'HA0     ') IHTYPE=84
        IF(IHTYPE.EQ.0) THEN
          WRITE(ITLIS,2000) HTYPE
          IFL=16
        ENDIF
        GO TO 100
      ENDIF
C
C          Keyword EPOL
      IF(WORD.EQ.'EPOL  ') THEN
        READ(ITCOM,*) PLEM,PLEP
        WRITE(ITLIS,*) PLEM,PLEP
        LOC(57)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MSSMD: optional 2nd geenration soft terms
      IF(WORD.EQ.'MSSMD   ') THEN
        READ(ITCOM,*) XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
        WRITE(ITLIS,*) XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
        LOC(58)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MSSME: optional U(1) and SU(2) gaugino masses
      IF(WORD.EQ.'MSSME   ') THEN
        READ(ITCOM,*) XM1SS,XM2SS
        WRITE(ITLIS,*) XM1SS,XM2SS
        LOC(59)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword GMSB: gauge-mediated SUSY breaking model
      IF(WORD.EQ.'GMSB    ') THEN
        READ(ITCOM,*) XLAMGM,XMESGM,XN5GM,XTGBSU,XSMUSU,XCMGV
        WRITE(ITLIS,*) XLAMGM,XMESGM,XN5GM,XTGBSU,XSMUSU,XCMGV
        GOMSSM=.TRUE.
        GOGMSB=.TRUE.
        LOC(60)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG1: optional GUT scale gaugino masses
      IF(WORD.EQ.'NUSUG1   ') THEN
        READ(ITCOM,*) XNUSUG(1),XNUSUG(2),XNUSUG(3)
        WRITE(ITLIS,*) XNUSUG(1),XNUSUG(2),XNUSUG(3)
        LOC(61)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG2: optional GUT scale A terms
      IF(WORD.EQ.'NUSUG2   ') THEN
        READ(ITCOM,*) XNUSUG(6),XNUSUG(5),XNUSUG(4)
        WRITE(ITLIS,*) XNUSUG(6),XNUSUG(5),XNUSUG(4)
        LOC(62)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG3: optional GUT scale Higgs masses
      IF(WORD.EQ.'NUSUG3   ') THEN
        READ(ITCOM,*) XNUSUG(7),XNUSUG(8)
        WRITE(ITLIS,*) XNUSUG(7),XNUSUG(8)
        LOC(63)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG4: optional GUT scale 1st/2nd gen. masses
      IF(WORD.EQ.'NUSUG4   ') THEN
        READ(ITCOM,*) XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10)
     $,XNUSUG(9)
        WRITE(ITLIS,*) XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10)
     $,XNUSUG(9)
        LOC(64)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NUSUG5: optional GUT scale 3rd gen. masses
      IF(WORD.EQ.'NUSUG5   ') THEN
        READ(ITCOM,*) XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15)
     $,XNUSUG(14)
        WRITE(ITLIS,*) XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15)
     $,XNUSUG(14)
        LOC(65)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword NOGRAV: No gravitino decays
      IF(WORD.EQ.'NOGRAV  ') THEN
        READ(ITCOM,571) NOGRAV
        WRITE(ITLIS,572) NOGRAV
        LOC(66)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MGVTNO: Sets the gravitino mass
      IF(WORD.EQ.'MGVTNO  ') THEN
        READ(ITCOM,*) XMGVTO
        WRITE(ITLIS,*) XMGVTO
        LOC(67)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword AL3UNI: Impose alpha_s unification at M_GUT
      IF(WORD.EQ.'AL3UNI  ') THEN
        READ(ITCOM,571) AL3UNI
        WRITE(ITLIS,572) AL3UNI
        LOC(68)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword GMSB2: additional GMSB parameters
      IF(WORD.EQ.'GMSB2   ') THEN
        READ(ITCOM,*) XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
        WRITE(ITLIS,*) XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
        LOC(69)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword EEBREM: invoke bremsstrahlung in e+e- reactions
      IF(WORD.EQ.'EEBREM  ') THEN
        READ(ITCOM,*) RSHMIN,RSHMAX
        WRITE(ITLIS,*) RSHMIN,RSHMAX
        IBREM=.TRUE.
        LOC(70)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword EEBEAM: invoke beamstrahlung in e+e- reactions
      IF(WORD.EQ.'EEBEAM  ') THEN
        READ(ITCOM,*) RSHMIN,RSHMAX,UPSLON,SIGZ
        WRITE(ITLIS,*) RSHMIN,RSHMAX,UPSLON,SIGZ
        IBREM=.TRUE.
        IBEAM=.TRUE.
        LOC(71)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword QMKKG (QMW for EXTRADIM)
      IF(WORD.EQ.'QMKKG   ') THEN
        READ(ITCOM,*) QMIN,QMAX
        WRITE(ITLIS,*) QMIN,QMAX
        LOC(19)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword QTKKG (QTW for EXTRADIM)
      IF(WORD.EQ.'QTKKG   ') THEN
        READ(ITCOM,*)  QTMIN,QTMAX
        WRITE(ITLIS,*) QTMIN,QTMAX
        LOC(20)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword EXTRAD for EXTRADIM
      IF(WORD.EQ.'EXTRAD  ') THEN
        READ(ITCOM,*) NEXTRAD,MASSD,UVCUT
        WRITE(ITLIS,*) NEXTRAD,MASSD,UVCUT
        LOC(72)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MIJLIM
      IF(WORD.EQ.'MIJLIM  ') THEN
        READ(ITCOM,*) JLIM1,JLIM2,AMLIM1,AMLIM2
        WRITE(ITLIS,*) JLIM1,JLIM2,AMLIM1,AMLIM2
        IF(JLIM1.EQ.0.AND.JLIM2.EQ.0) THEN
          DO 720 I=1,NJET
            DO 721 J=1,NJET
              AMIJMN(I,J)=AMLIM1
              AMIJMX(I,J)=AMLIM2
721         CONTINUE
720       CONTINUE
        ELSEIF(JLIM1.GT.0.AND.JLIM1.LE.NJET.AND.JLIM2.GT.0.AND.
     $  JLIM2.LE.NJET) THEN
          AMIJMN(JLIM1,JLIM2)=AMLIM1
          AMIJMN(JLIM2,JLIM1)=AMLIM1
          AMIJMX(JLIM1,JLIM2)=AMLIM2
          AMIJMX(JLIM2,JLIM1)=AMLIM2
        ELSE
          WRITE(ITLIS,2008)
          IFL=73
        ENDIF
        LOC(73)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword MTOT
      IF(WORD.EQ.'MTOT    ') THEN
        READ(ITCOM,*) EHMGMN,EHMGMX
        WRITE(ITLIS,*) EHMGMN,EHMGMX
        LOC(74)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword SUGRHN
      IF(WORD.EQ.'SUGRHN  ') THEN
        READ(ITCOM,*) XMN3NR,XMAJNR,XANSS,XNRSS
        WRITE(ITLIS,*) XMN3NR,XMAJNR,XANSS,XNRSS
        LOC(75)=NSEL
        GO TO 100
      ENDIF
C
C          Keyword AMSB
      IF(WORD.EQ.'AMSB    ') THEN
        READ(ITCOM,*) XM0SU,XMHSU,XTGBSU,XSMUSU
        WRITE(ITLIS,*) XM0SU,XMHSU,XTGBSU,XSMUSU
        LOC(76)=NSEL
        GOMSSM=.TRUE.
        GOSUG=.TRUE.
        GOAMSB=.TRUE.
        GO TO 100
      ENDIF
C
C          Keyword SSBCSC
      IF(WORD.EQ.'SSBCSC  ') THEN
        READ(ITCOM,*) XSBCS
        WRITE(ITLIS,*) XSBCS
        LOC(77)=NSEL
        GO TO 100
      ENDIF
C
C          None of the above
C
      WRITE(ITLIS,2000) WORD
      IFL=10
      RETURN
C
C          Error message or warnings
C
 1001     FORMAT(//2X,'YOU HAVE GIVEN LIMITS FOR AN OVERLAPPING SET', 
     $  ' OF VARIABLES. SET MINIMIZING PPERP INTERVAL WILL BE USED.')
 1999     FORMAT(//' YOU FORGOT TO SELECT A PROCESS FOR GENERATION.'
     $    /' AVAILABLE AT PRESENT ARE ',
     $    /' TWOJET  E+E-  DRELLYAN  MINBIAS  WPAIR  SUPERSYM,' 
     $    /' HIGGS  PHOTON  TCOLOR')
 2000 FORMAT(//2X,A8,' IS NOT A RECOGNIZABLE PARAMETER. JOB TERMINATED')
 2001     FORMAT(//2X,' YOU CANNOT GIVE LIMITS FOR BOTH THETA AND Y.',
     $    ' MAKE UP YOUR MIND. JOB TERMINATED.')
 2002     FORMAT(/'  WHAT IS THAT SUPPOSED TO BE')
 2003     FORMAT(/'  YOU CANNOT RUN WS AND Z0 OR GAMMAS AT THE',
     $    ' SAME TIME. JOB TERMINATED')
 2004 FORMAT(//'  PARAMETER OUT OF RANGE. JOB TERMINATED.')
 2005 FORMAT(//'  YOU CANNOT FORCE DECAY OF A QUARK, YOU MUST CHOSE '
     $,' A PARTICLE')
 2006     FORMAT(//2X,' INVALID JETTYPE VALUES. JOB TERMINATED.')
 2007 FORMAT(//2X,'YOU CANNOT USE MSSM AND SUGRA SIMULTANEOUSLY')
 2008 FORMAT(//2X,'INVALID JET-JET MASS LIMITS. JOB TERMINATED.')
C
      END
