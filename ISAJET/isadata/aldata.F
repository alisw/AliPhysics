#include "isajet/pilot.h"
      BLOCK DATA ALDATA
C          INITIALIZE ALL COMMON BLOCKS
C.......................................................................
C          WARNINGS: MANY VARIABLES SET IN ALDATA ARE ALSO SET BY      .
C          SUBROUTINE RESET.                                           .
C                                                                      .
C          ALDATA SHOULD ALWAYS BE LOADED WHEN USING ISAJET OR WHEN    .
C          READING AN ISAJET TAPE.                                     .
C.......................................................................
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/lstprt.inc"
#include "isajet/dkytab.inc"
#include "isajet/dylim.inc"
#include "isajet/eepar.inc"
#include "isajet/frgpar.inc"
#include "isajet/idrun.inc"
#include "isajet/jetlim.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetset.inc"
#include "isajet/jetsig.inc"
#include "isajet/limevl.inc"
#include "isajet/mbpar.inc"
#include "isajet/nodcay.inc"
#include "isajet/partcl.inc"
#include "isajet/primar.inc"
#include "isajet/prtout.inc"
#include "isajet/qcdpar.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/seed.inc"
#include "isajet/sspar.inc"
#include "isajet/tcpar.inc"
#include "isajet/totals.inc"
#include "isajet/types.inc"
#include "isajet/wcon.inc"
#include "isajet/mbgen.inc"
#include "isajet/force.inc"
#include "isajet/zevel.inc"
#include "isajet/final.inc"
#include "isajet/keys.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
#include "isajet/sugnu.inc"
#include "isajet/isapw.inc"
#include "isajet/sstype.inc"
#include "isajet/listss.inc"
C
      INTEGER III,JJJ
      INTEGER MXGOQJ
      PARAMETER (MXGOQJ=MXGOJ*MXGOQ)
      INTEGER MXGOWJ
      PARAMETER (MXGOWJ=25*MXGOJ)
      INTEGER MXT29
      PARAMETER (MXT29=29*MXTYPE)
C          SUSY IDENT codes from /SSTYPE/
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
C
C          DATA FOR IDRUN
C          IDVER=100*VERSION+CYCLE
C     DATA IDVER/600/
C
C          DATA FOR ITAPES
      DATA ITDKY,ITEVT,ITCOM,ITLIS/1,2,5,6/
C
C          DATA FOR QLMASS
C          AMLEP LABELED BY INDEX...SEE FLAVOR
C          SETW RESETS W+- AND Z0 MASSES
      DATA AMLEP/.3,.3,.5,1.6,5.0,175.,-1.,-1.,0.,0.,
     $0.,.511003E-3,0.,.105661,0.,1.777,3*-1.,.49767,.49767,
     $79*0./
      DATA NQLEP,NMES,NBARY/61,2,2/
C
C          DATA FOR PJETS
      DATA IDJETS/MXJETS*0/,IDENTW/0/
C
C          DATA FOR PINITS
      DATA IDINIT/2*0/
C
C          DATA FOR LSTPRT
      DATA LSTPRT/0/
C
C          DATA FOR MBPAR
      DATA PUD0/.45/,PJSPN,PISPN/2*.5/,SIGQT0/.35/,XGEN0/.9,1./,PMIX01/
     $.25,.25,.5,0.,.5,1./,PMIX02/.5,.5,1.,0.,0.,1./
      DATA PBARY0/.075/
C
C          DATA FOR MBGEN
      DATA MNPOM,MXPOM/1,LIMPOM/
C
C          DATA FOR SEED
      DATA XSEED/'0'/
C
C          DATA FOR TCPAR
      DATA TCMRHO,TCGRHO/1000.,100./
C
C          DATA FOR FRGPAR
C          F(X)=1-XGEN(1)+XGEN(1)*(XGEN(2)+1)*(1-X)**XGEN(2) FOR U,D,S
C          PETERSON FRAGMENTATION, EPSILON=XGEN(I)*M(I)**2 FOR C,B,T
      DATA PUD,PBARY/.43,.10/
      DATA SIGQT,PEND/.35,.14/
      DATA XGEN/.96,3.,0.,.8,.5,.5,.5,.5/
      DATA PSPIN1/.5,.5,.5,.75,.75,.75,.75,.75/
      DATA PMIX1/.25,.25,.5,0.,.5,1./,PMIX2/.5,.5,1.,0.,0.,1./
      DATA XGENSS/9*.5/
C
C          DATA FOR JETLIM
      DATA BLIMS/MXLX12*-1.E9/
C
C          DATA FOR NODCAY
      DATA NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR/6*.FALSE./
      DATA NOGRAV/.FALSE./
C
C          DATA FOR TYPES
      DATA LOC/100*0/,NTYP/100/
      DATA NJTTYP/MXTYPE*0/
      DATA (JETYP(1,JJJ),JJJ=1,MXTYPE)/MXTYPE*'ALL     '/,
     $((JETYP(III,JJJ),III=2,30),JJJ=1,MXTYPE)/MXT29*'        '/
      DATA NWWTYP/2*0/
      DATA (WWTYP(1,JJJ),JJJ=1,2)/2*'ALL     '/,
     $((WWTYP(III,JJJ),III=2,30),JJJ=1,2)/58*'        '/
      DATA JWTYP/4/
C
C          DATA FOR PRIMAR
      DATA IDIN/1120,1120/
      DATA NTRIES/1000/
      DATA NSIGMA/20/
C
C          DATA FOR DKYTAB
      DATA LOOK/MXLOOK*0/
      DATA CBR/MXDKY*0./
      DATA MODE/MXDKY*0,MXDKY*0,MXDKY*0,MXDKY*0,MXDKY*0/
C
C          DATA FOR Q1Q2
      DATA GOQ/MXGOQJ*.TRUE./
      DATA GOALL/MXGOJ*.TRUE./
      DATA GODY/.TRUE.,.FALSE.,.FALSE.,.TRUE./
      DATA GOWW/50*.TRUE./,ALLWW/2*.TRUE./
      DATA GOWMOD/MXGOWJ*.TRUE./
C
C          DATA FOR WCON
      DATA MATCH/
     $0,3,2,5,4,7,6,9,8,11,10,13,12,0,0,17,16,0,0,21,20,0,0,25,24,
     $0,5,0,0,2,0,8,7,0,0,12,11,0,17,0,0,14,21,0,0,18,25,0,0,22,
     $0,0,4,3,0,9,0,0,6,13,0,0,10,0,16,15,0,0,20,19,0,0,24,23,0,
     $0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24/
      DATA CUTOFF,CUTPOW/.200,1.0/
      DATA WMASS/0.,80.2,80.2,91.19/
      DATA WFUDGE/1.85/
C
C          DATA FOR TOTALS
      DATA NKINPT,NWGEN,NKEEP/3*0/,SUMWT/0./
C
C          DATA FOR DYLIM
      DATA BLIM1/12*-1.E9/
C
C          DATA FOR EEPAR
      DATA PLEP/0./,PLEM/0./,IBREM/.FALSE./,IBEAM/.FALSE./
C
C          DATA FOR PARTCL
      DATA NPTCL/0/
C
C            DATA FOR PRTOUT
      DATA NEVPRT,NJUMP/1,1/
C
C          DATA FOR JETSET
      DATA NJSET/0/
C
C          DATA FOR QCDPAR
      DATA ALAM,ALAM2/.2,.04/,CUTJET/6./,ISTRUC/6/
C
C          DATA FOR FORCE
      DATA NFORCE/0/
C
C          DATA FOR NRECS
      DATA NRECS/0/
C
C          DATA FOR KEYS
      DATA KEYS/MXKEYS*.FALSE./
C
C          DATA FOR MATCHH
      DATA MATCHH/
     $1,3,2,5,4,7,6,9,8,11,10,13,12,
     $15,14,17,16,19,18,21,20,23,22,25,24,
     $26,28,27,29/
      DATA USELIM/.FALSE./
      DATA CONCUT/1.0/
C
C          DATA FOR SSPAR
      DATA AMGVSS/1.E20/
C
C          DATA FOR XMSSM
      DATA GOMSSM/.FALSE./,GOSUG/.FALSE./,GOGMSB/.FALSE./
      DATA GOAMSB/.FALSE./
      DATA AL3UNI/.FALSE./
      DATA XM1SS,XM2SS/1.E20,1.E20/
      DATA XMGVTO/1.E20/
      DATA XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS/1.E20,1.E20,1.E20,1.E20,1.E20/
      DATA XRSLGM,XDHDGM,XDHUGM,XDYGM/1.,0.,0.,0./
      DATA XN51GM,XN52GM,XN53GM/0.,0.,0./
      DATA XMN3NR/0./,XMAJNR/1.E20/,XANSS/0./,XNRSS/0./,XSBCS/0./
C          DATA FOR SUGNU
      DATA XNUSUG/18*1.E20/
C
C          DATA FOR ISAPW
      DATA ISAPW1/'ALDATA REQUIRED BY FORTRAN G,H'/
C
C          DATA FOR LISTSS
      DATA LISTSS/ISGL,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2,
     $ISW1,MSW1,ISW2,MSW2,ISZ1,ISZ2,ISZ3,ISZ4,
     $ISNEL,MSNEL,ISEL,MSEL,ISNML,MSNML,ISMUL,MSMUL,ISNTL,MSNTL,
     $ISTAU1,MSTAU1,ISER,MSER,ISMUR,MSMUR,ISTAU2,MSTAU2,
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,82,83,84,86,-86/
C
      END
