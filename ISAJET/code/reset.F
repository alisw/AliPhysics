#include "isajet/pilot.h"
      SUBROUTINE RESET
C          RESET ALL USER DEFINED VARIABLES
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/dylim.inc"
#include "isajet/frgpar.inc"
#include "isajet/hcon.inc"
#include "isajet/jetlim.inc"
#include "isajet/jetpar.inc"
#include "isajet/nodcay.inc"
#include "isajet/primar.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/seed.inc"
#include "isajet/sspar.inc"
#include "isajet/tcpar.inc"
#include "isajet/types.inc"
#include "isajet/wcon.inc"
#include "isajet/force.inc"
#include "isajet/mbgen.inc"
#include "isajet/isloop.inc"
#include "isajet/limevl.inc"
#include "isajet/xmssm.inc"
#include "isajet/eepar.inc"
#include "isajet/mglims.inc"
C
      INTEGER I,I1,I2,I3,J1,INDEX,J,K
      REAL UNDEF,AMASS
      CHARACTER*8 BLANK
      DATA BLANK/'        '/
      DATA UNDEF/-1.E9/
C          RESET DYLIM
      DO 110 I=1,12
      BLIM1(I)=UNDEF
      SETLMQ(I)=.TRUE.
110   CONTINUE
C          RESET FRGPAR
      PUD=.43
      PBARY=.10
      SIGQT=.35
      PEND=.14
      XGEN(1)=.96
      XGEN(2)=3.
      XGEN(3)=0.
      XGEN(4)=.8
      XGEN(5)=.5
      XGEN(6)=.5
      XGEN(7)=.5
      XGEN(8)=.5
      DO 111 K=1,9
111   XGENSS(K)=.5
      PSPIN1(1)=.5
      PSPIN1(2)=.5
      PSPIN1(3)=.5
      PSPIN1(4)=.75
      PSPIN1(5)=.75
      PSPIN1(6)=.75
      PSPIN1(7)=.75
      PSPIN1(8)=.75
      PMIXX1(1)=.25
      PMIXX1(2)=.25
      PMIXX1(3)=.5
      PMIXX1(4)=0.
      PMIXX1(5)=.5
      PMIXX1(6)=1.
      PMIXX2(1)=.5
      PMIXX2(2)=.5
      PMIXX2(3)=1.
      PMIXX2(4)=0.
      PMIXX2(5)=0.
      PMIXX2(6)=1.
C          RESET ISLOOP
      NEVOLV=1
      NFRGMN=1
C          RESET JETLIM
      DO 120 I=1,12*MXLIM
        BLIMS(I)=UNDEF
        SETLMJ(I)=.TRUE.
120   CONTINUE
C          RESET NODCAY
      NODCAY=.FALSE.
      NOETA=.FALSE.
      NOPI0=.FALSE.
      NONUNU=.FALSE.
      NOEVOL=.FALSE.
      NOHADR=.FALSE.
      NOGRAV=.FALSE.
C          RESET PRIMAR
      IDIN(1)=1120
      IDIN(2)=1120
      NTRIES=1000
      NSIGMA=20
C          RESET QCDPAR
      ALAM=.2
      ALAM2=ALAM**2
      CUTJET=6.
      ISTRUC=6
C          RESET QLMASS
      AMLEP(6)=175.
      AMLEP(7)=-1.
      AMLEP(8)=-1.
      DO 125 I=1,9
      CALL FLAVOR(80+I,I1,I2,I3,J1,INDEX)
125   AMLEP(INDEX)=0.
      CALL FLAVOR(29,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=100.
      CALL FLAVOR(30,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=0.
      CALL FLAVOR(39,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=100.
      CALL FLAVOR(40,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=100.
      DO 126 I=1,6
      CALL FLAVOR(20+I,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=100.+AMASS(I)
      CALL FLAVOR(30+I,I1,I2,I3,J1,INDEX)
      AMLEP(INDEX)=100.+AMASS(I+10)
126   CONTINUE
C          RESET Q1Q2
      DO 130 I=1,MXGOQ
      DO 130 J=1,MXGOJ
130   GOQ(I,J)=.TRUE.
      DO 131 I=1,MXGOJ
131   GOALL(I)=.TRUE.
      GODY(1)=.TRUE.
      GODY(2)=.FALSE.
      GODY(3)=.FALSE.
      GODY(4)=.TRUE.
      DO 132 I=1,2
      ALLWW(I)=.TRUE.
      DO 132 J=1,25
132   GOWW(J,I)=.TRUE.
      DO 133 I=1,3
      DO 133 J=1,25
133   GOWMOD(J,I)=.TRUE.
C          RESET TCPAR
      TCMRHO=1000.
      TCGRHO=100.
C          RESET TYPES
      DO 140 I=1,NTYP
140   LOC(I)=0
      DO 141 I=1,MXTYPE
      NJTTYP(I)=0
      JETYP(1,I)='ALL     '
      DO 141 K=2,30
141   JETYP(K,I)=BLANK
      JWTYP=4
      DO 142 I=1,2
      NWWTYP(I)=0
      WWTYP(1,I)='ALL     '
      DO 142 K=2,4
142   WWTYP(K,I)=BLANK
      DO 143 I=1,3
      NWMODE(I)=0
      WMODES(1,I)='ALL     '
      DO 143 K=2,30
143   WMODES(K,I)=BLANK
C          RESET WCON
      SIN2W=.232
      WMASS(2)=80.2
      WMASS(3)=WMASS(2)
      WMASS(4)=91.19
      CALL FLAVOR(80,I1,I2,I3,J,INDEX)
      AMLEP(INDEX)=WMASS(2)
      CALL FLAVOR(90,I1,I2,I3,J,INDEX)
      AMLEP(INDEX)=WMASS(4)
      CUTOFF=.200
      CUTPOW=1.0
      WFUDGE=1.75
C          RESET MBGEN
      MNPOM=1
      MXPOM=LIMPOM
C          RESET FORCE
      NFORCE=0
C
C          RESET QCD EVOLUTION CUTS
      USELIM=.FALSE.
      CONCUT=1.0
C
C          RESET SSPAR
      AMGVSS=1.E20
C
C          RESET XMSSM
      GOMSSM=.FALSE.
      GOSUG=.FALSE.
      GOAMSB=.FALSE.
      XM1SS=1.E20
      XM2SS=1.E20
      XMAJNR=1.E20
      XMGVTO=1.E20
C
C          RESET HCON
      IHTYPE=0
C
C          RESET EEPAR
      PLEP=0.
      PLEM=0.
C
C          RESET MGLIMS
      EHMGMN=-1.E9
      EHMGMX=-1.E9
      YHMGMN=-1.E9
      YHMGMX=-1.E9
      DO 150 I=1,MXLIM
        DO 151 J=1,MXLIM
          AMIJMN(I,J)=-1.E9
          AMIJMX(I,J)=-1.E9
151     CONTINUE
150   CONTINUE
C
      RETURN
      END
