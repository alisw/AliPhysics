#include "isajet/pilot.h"
      SUBROUTINE PRTLIM
C
C          Print initial conditions and limits for generating jets
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/force.inc"
#include "isajet/mbgen.inc"
#include "isajet/qcdpar.inc"
#include "isajet/qlmass.inc"
#include "isajet/wcon.inc"
#include "isajet/jetpar.inc"
#include "isajet/keys.inc"
#include "isajet/kkgrav.inc"
#include "isajet/frgpar.inc"
#include "isajet/nodcay.inc"
#include "isajet/prtout.inc"
#include "isajet/seed.inc"
#include "isajet/types.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/ptpar.inc"
#include "isajet/idrun.inc"
#include "isajet/dylim.inc"
#include "isajet/hcon.inc"
#include "isajet/isloop.inc"
#include "isajet/xmssm.inc"
C
      INTEGER I,II,K,NPRT,I1,I2,I3,J1,INDEX,IQ,KK,KKK,NN,N0J
      REAL AMASS
      CHARACTER*8 LSTRUC(6),LMODE(5),STUF(6),IDFMT(2)
      CHARACTER*8 WTITL(4)
      CHARACTER*8 LABEL,L0
      CHARACTER*8 BLANK
      CHARACTER*40 V,VISAJE
      REAL AM(6),AML(6)
      INTEGER NPRSS
      PARAMETER (NPRSS=32)
      INTEGER IDPRSS(NPRSS)
      REAL AMPRSS(NPRSS)
      CHARACTER*8 LPRSS(NPRSS)
      DATA LSTRUC/'OWENS','BAIER','EICHTEN','DUKE','CTEQ2L','CTEQ3L'/
      DATA WTITL/'GM','W+','W-','Z0'/
      DATA BLANK/' '/
      DATA IDPRSS/21,22,23,24,25,26,41,42,43,44,45,46,
     $31,32,33,34,35,36,52,54,56,
     $29,30,40,50,60,39,49,82,83,84,86/
C
C          Print version
      V=VISAJE()
      WRITE(ITLIS,1000) V
1000  FORMAT('1',44('*')/' *',42X,'*'/
     C  ' * ',A40,' *'/
     C  ' *',42X,'*'/' ',44('*')/)
C
C          Print title, reaction, energy, number, run id
      WRITE(ITLIS,1010) TITLE
1010  FORMAT(/11X,10A8)
      WRITE(ITLIS,1020) NJET
1020  FORMAT(/2X,'NUMBER OF JETS TO BE GENERATED PER EVENT',I3)
      DO 100 I=1,2
        IDFMT(I)=LABEL(IDIN(I))
100   CONTINUE
      WRITE(ITLIS,1030) IDFMT(1),IDFMT(2),ECM
1030  FORMAT(/2X,A8,' ON ',A8,' AT COM ENERGY',E15.4)
      WRITE(ITLIS,1040) REAC,NEVENT
1040  FORMAT(/2X,'NUMBER OF ',A8,' EVENTS TO BE GENERATED',I10)
      IF(NEVPRT.GT.0) WRITE(ITLIS,1050) NEVPRT,NJUMP
1050  FORMAT(/2X,'PRINT A MAXIMUM OF ',I6,
     C  ' EVENTS SKIPPING ',I6,' EVENTS AT A TIME')
      WRITE(ITLIS,1060) IDG
1060  FORMAT(/2X,'RUN ID   ',2I10.6)
C
C          Print W/Higgs parameters
C
      IF(KEYS(3).OR.KEYS(7).OR.KEYS(11)) THEN
        IF(KEYS(3)) THEN
          II=0
          DO 200 I=1,4
            IF(.NOT.GODY(I)) GOTO 200
            II=II+1
            STUF(II)=WTITL(I)
200       CONTINUE
        ELSEIF(KEYS(11)) THEN
	  II=1
	  STUF(II)='GRAV'
        ELSE
          II=1
          STUF(II)='HIGGS'
          IF(IHTYPE.EQ.82) STUF(II)='HL0'
          IF(IHTYPE.EQ.83) STUF(II)='HH0'
          IF(IHTYPE.EQ.84) STUF(II)='HA0'
        ENDIF
        WRITE(ITLIS,2000) (STUF(K),K=1,II)
2000    FORMAT(//10X,'PARAMETERS FOR',4(2X,A8))
        WRITE(ITLIS,2010) QMIN,QMAX
2010    FORMAT(' MASS LIMITS',15X,2E15.4)
        WRITE(ITLIS,3020) QTMIN,QTMAX
        WRITE(ITLIS,3030) THWMIN,THWMAX
        WRITE(ITLIS,3040) PHWMIN,PHWMAX
        WRITE(ITLIS,3050) YWMIN,YWMAX
        WRITE(ITLIS,3060) XWMIN,XWMAX
        NPRT=0
        IF(FIXQT) THEN
          NPRT=NPRT+1
          STUF(NPRT)='QTW     '
        ENDIF
        IF(FIXQM) THEN
          NPRT=NPRT+1
          STUF(NPRT)='QMW     '
        ENDIF
        IF(FIXYW) THEN
          NPRT=NPRT+1
          STUF(NPRT)='YW      '
        ENDIF
        IF(FIXXW) THEN
          NPRT=NPRT+1
          STUF(NPRT)='XW      '
        ENDIF
        IF(FIXPHW) THEN
          NPRT=NPRT+1
          STUF(NPRT)='PHW     '
        ENDIF
        IF(NPRT.NE.0) WRITE(ITLIS,2070) (STUF(K),K=1,NPRT)
      ENDIF
C
C          Print jet parameters
C
      DO 300 I=1,6
300   STUF(I)=BLANK
      IF(KEYS(11)) THEN
        N0J=3
      ELSE
        N0J=1
      ENDIF
      DO 310 I=N0J,NJET
        NPRT=0
        WRITE(ITLIS,3000) I
3000    FORMAT(//10X,'JET NO.',I3,/)
        WRITE(ITLIS,3010) PMIN(I),PMAX(I)
3010    FORMAT(' MOMENTUM LIMITS',11X,2E15.4)
        WRITE(ITLIS,3020) PTMIN(I),PTMAX(I)
3020    FORMAT(' TRANSVERSE MOMENTUM LIMITS',2E15.4)
        WRITE(ITLIS,3030) THMIN(I),THMAX(I)
3030    FORMAT(' THETA LIMITS',14X,2E15.4)
        WRITE(ITLIS,3040) PHIMIN(I),PHIMAX(I)
3040    FORMAT(' PHI LIMITS',16X,2E15.4)
        WRITE(ITLIS,3050) YJMIN(I),YJMAX(I)
3050    FORMAT(' RAPIDITY (Y) LIMITS',7X,2E15.4)
        WRITE(ITLIS,3060) XJMIN(I),XJMAX(I)
3060    FORMAT(' FEYNMAN X LIMITS',10X,2E15.4)
        IF(.NOT.GOALL(I)) WRITE(ITLIS,3070) (JETYP(K,I),K=1,NJTTYP(I))
3070    FORMAT(' JET TYPE',23X,A8,1X,A8,1X,A8,1X,A8,1X,A8)
        IF((KEYS(6).OR.KEYS(7)).AND..NOT.ALLWW(I))
     $  WRITE(ITLIS,3080) (WWTYP(K,I),K=1,NWWTYP(I))
3080    FORMAT(' DECAY MODES',20X,A8,1X,A8,1X,A8,1X,A8,1X,A8)
        IF(FIXP(I)) THEN
          NPRT=NPRT+1
          STUF(NPRT)='P       '
        ENDIF
        IF(FIXPT(I)) THEN
          NPRT=NPRT+1
          STUF(NPRT)='PT      '
        ENDIF
        IF(FIXYJ(I)) THEN
          NPRT=NPRT+1
          STUF(NPRT)='TH      '
          NPRT=NPRT+1
          STUF(NPRT)='Y      '
        ENDIF
        IF(FIXPHI(I)) THEN
          NPRT=NPRT+1
          STUF(NPRT)='PHI     '
        ENDIF
        IF(FIXXJ(I)) THEN
          NPRT=NPRT+1
          STUF(NPRT)='X       '
        ENDIF
        IF(NPRT.EQ.0) GOTO 310
        WRITE(ITLIS,2070) (STUF(K),K=1,NPRT)
2070    FORMAT(/5X,'FOLLOWING PARAMETERS HAVE BEEN FIXED ',6A8)
310   CONTINUE
C
C          Print structure functions, QCD parameters, W parameters,
C          and other parameters changed from their default values.
C
      IF(KEYS(1).OR.KEYS(3).OR.KEYS(5).OR.KEYS(6).OR.KEYS(7)
     $.OR.KEYS(8).OR.KEYS(9).OR.KEYS(11).OR.KEYS(12)) THEN
        IF(ISTRUC.GT.0) THEN
          WRITE(ITLIS,4000) LSTRUC(ISTRUC)
4000      FORMAT(//1X,A8,' STRUCTURE FUNCTIONS')
        ELSEIF(ISTRUC.EQ.-999) THEN
          WRITE(ITLIS,4001)
4001      FORMAT(//1X,' PDFLIB STRUCTURE FUNCTIONS')
        ENDIF
      ENDIF
      WRITE(ITLIS,4010) ALAM,CUTJET
4010  FORMAT(//' QCD LAMBDA =',F10.4,10X,'JET CUTOFF MASS =',F10.3)
      WRITE(ITLIS,4020) AMLEP(6),AMLEP(7),AMLEP(8)
4020  FORMAT(/' HEAVY QUARK MASSES =',3F8.2)
      IF(LOC(36).NE.0) THEN
        CALL FLAVOR(80,I1,I2,I3,J1,INDEX)
        WRITE(ITLIS,4030) (AMLEP(INDEX+K),K=1,9)
4030    FORMAT(/' HIGGS MASSES =',6F8.2/15X,3F8.2)
      ENDIF
C          Supersymmetry
      IF(KEYS(5).AND..NOT.GOMSSM) THEN
        DO 410 IQ=1,6
          AM(IQ)=AMASS(20+IQ)
          AML(IQ)=AMASS(30+IQ)
410     CONTINUE
        WRITE(ITLIS,4040) (AM(KK),KK=1,6)
4040    FORMAT(//' SQUARK  MASSES ',7F9.2)
        WRITE(ITLIS,4050) (AML(KK),KK=1,6)
4050    FORMAT(' SLEPTON MASSES ',7F9.2)
        AM(1)=AMASS(29)
        AM(2)=AMASS(30)
        AM(3)=AMASS(39)
        AM(4)=AMASS(40)
        WRITE(ITLIS,4060) (AM(KK),KK=1,4)
4060    FORMAT(' GAUGINO MASSES (WITH SIGNS) ',7F9.2)
        IF(LOC(44).NE.0) WRITE(ITLIS,4070) XGENSS(9),
     $  (XGENSS(KK),KK=1,8) 
4070    FORMAT(/' FRAG. PARAM. XGENSS = ',9F8.3)
      ENDIF
C          MSSM model
      IF(GOMSSM) THEN
        DO 420 I=1,NPRSS
          LPRSS(I)=LABEL(IDPRSS(I))
          AMPRSS(I)=AMASS(IDPRSS(I))
420     CONTINUE
        WRITE(ITLIS,4100) (LPRSS(I),AMPRSS(I),I=1,NPRSS)
4100    FORMAT(/' MSSM MASSES (WITHOUT SIGNS):'/
     $  10(' M(',A5,') = ',F10.3,5X,'M(',A5,') = ',F10.3,
     $  5X,'M(',A5,') = ',F10.3/),
     $  ' M(',A5,') = ',F10.3,5X,'M(',A5,') = ',F10.3)
        IF(.NOT.GOSUG) THEN
          WRITE(ITLIS,4110) XTBSS,XMUSS,XATSS
4110      FORMAT(/' OTHER MSSM PARAMETERS:'/
     $    ' TAN(BETA) = ',F10.3,5X,'MU = ',F10.3,5X,'A_t = ',F10.3)
        ENDIF
        IF(GOSUG) THEN
          WRITE(ITLIS,4120) XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
4120      FORMAT(/' DETERMINED FROM SUGRA INPUT:'/
     $    ' M_0       =',F10.3,'   M_(1/2) =',F10.3,'   A_0 =',F10.3/
     $    ' TAN(BETA) =',F10.4,'   SGN(MU) =',F6.1)
        ENDIF
      ENDIF
C
C          Other parameters
C
      IF(LOC(30).NE.0) WRITE(ITLIS,4200) XGEN
4200  FORMAT(/' FRAGMENTATION PARAMETER XGEN =',6F8.3)
      IF(LOC(31).NE.0) WRITE(ITLIS,4210) SIGQT
4210  FORMAT(/' FRAGMENTATION PARAMETER SIGQT =',F8.3)
      IF(KEYS(2).OR.KEYS(3).OR.KEYS(6).OR.KEYS(7)) THEN
        WRITE(ITLIS,4220) SIN2W,WMASS(2),WMASS(4),WGAM(2),WGAM(4)
4220    FORMAT(//' WEINBERG MODEL',5X,'SIN**2(THETA-W)=',F8.4/
     $  ' MASSES = ',F8.2,',',F8.2,'   WIDTHS = ',F8.3,',',F8.3)
      ENDIF
      IF(KEYS(3)) THEN
        IF(.NOT.STDDY) WRITE(ITLIS,4230) CUTOFF,CUTPOW
4230    FORMAT(/' CUTOFF FUNCTION IS QT**2=',E11.4,'*Q**',E11.4)
        IF(LOC(50).NE.0) WRITE(ITLIS,4240) WFUDGE
4240    FORMAT(/' W fudge factor (WFUDGE) = ',F8.3)
      ENDIF
      IF(KEYS(7)) THEN
        WRITE(ITLIS,4250) HMASS,HGAM
4250    FORMAT(/' HIGGS MASS = ',F8.2,'  WIDTH = ',F8.3)
      ENDIF
      WRITE(ITLIS,4260) XSEED
4260  FORMAT(/' SEED FOR RANDOM NUMBER GENERATOR = ',A24)
      IF(LOC(13).NE.0) WRITE(ITLIS,4270) FRPAR
4270  FORMAT(//'  FRAGMENTATION PARAMETERS ',8E11.3,2(/,26X,8E11.3))
      IF(LOC(34).NE.0) WRITE(ITLIS,4280) MNPOM,MXPOM
4280  FORMAT(//'  NUMBER OF POMERONS =',I4,' TO',I4)
C
C          Print KKG parameters
      IF(KEYS(11)) THEN
        WRITE(ITLIS,4291) NEXTRAD
        WRITE(ITLIS,4292) MASSD
        WRITE(ITLIS,4296) SURFD
        WRITE(ITLIS,4295) KKGSD
        WRITE(ITLIS,4297) UVCUT
4291    FORMAT(//' NB EXTRA-DIMENSIONS',7X,I4)
4292    FORMAT(' SCALE M_D',15X,E15.4)
4295    FORMAT(' KKGSD FACTOR',12X,E15.4)
4296    FORMAT(' SD SURFACE',14X,E15.4)
4297    FORMAT(' UV CUTOFF',17X,I4)
      ENDIF
C
C          Decay switches
C
      IF(LOC(46).NE.0) THEN
        WRITE(ITLIS,4300)
4300    FORMAT(//'  SECONDARY W DECAY MODES:')
        NN=MAX0(NWMODE(1),1)
        WRITE(ITLIS,4310) (WMODES(KKK,1),KKK=1,NN)
4310    FORMAT('  W+ --> '
     $  ,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8)
        NN=MAX0(NWMODE(2),1)
        WRITE(ITLIS,4320) (WMODES(KKK,2),KKK=1,NN)
4320    FORMAT('  W- --> '
     $  ,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8)
        WRITE(ITLIS,4330) (WMODES(KKK,3),KKK=1,NN)
4330    FORMAT('  Z0 --> '
     $  ,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8,1X,A8)
      ENDIF
C
      IF(NODCAY) WRITE(ITLIS,4400)
4400  FORMAT(//'  NO DECAYS OF HADRONS WILL BE GENERATED')
      IF(NOETA) WRITE(ITLIS,4410)
4410  FORMAT(//'  NO DECAYS OF ETAS WILL BE GENERATED')
      IF(NOPI0) WRITE(ITLIS,4420)
4420  FORMAT(//'  NO DECAYS OF PI0S WILL BE GENERATED')
      IF(NONUNU) WRITE(ITLIS,4430)
4430  FORMAT(//'  NO DECAYS Z0---->NU+NU WILL BE GENERATED')
      IF(NOEVOL) WRITE(ITLIS,4440)
4440  FORMAT(//'  NO QCD JET EVOLUTION WILL BE DONE')
      IF(NOHADR) WRITE(ITLIS,4450)
4450  FORMAT(//'  NO JET HADRONIZATION WILL BE DONE')
      IF(GOGMSB.AND.NOGRAV) WRITE(ITLIS,4460)
4460  FORMAT(//'  NO GRAVITINO DECAYS WILL BE GENERATED')
C
C          Print forced decay modes and M.E. flag
C
      IF(NFORCE.NE.0) THEN
        WRITE(ITLIS,4500)
4500    FORMAT(//7X,'FORCED DECAY MODES (MATRIX ELEMENT FLAGS)'/
     $  2X,'PART',6X,'DECAY MODE')
        DO 450 I=1,NFORCE
          IF(IFORCE(I).EQ.0) GOTO 450
          L0=LABEL(IFORCE(I))
          DO 451 K=1,5
            LMODE(K)=BLANK
            IF(MFORCE(K,I).EQ.0) GO TO 451
            LMODE(K)=LABEL(MFORCE(K,I))
451       CONTINUE
          WRITE(ITLIS,4510) L0,(LMODE(K),K=1,5),MEFORC(I)
4510      FORMAT(2X,6A10,'(M.E. =',I5,')')
450     CONTINUE
      ENDIF
C
C          Print multiple evolution/fragmentation information
C
      IF(NEVOLV.NE.1.OR.NFRGMN.NE.1) THEN
        WRITE(ITLIS,4600) NEVOLV,NFRGMN
4600    FORMAT(//,
     $  '  MULTIPLE EVOLUTION AND FRAGMENTATION VERSION'/
     $  '  EVENTS WILL BE EVOLVED',I6,' TIMES'/
     $  '  AND FRAGMENTED        ',I6,' TIMES'/)
      ENDIF
C
      RETURN
      END
