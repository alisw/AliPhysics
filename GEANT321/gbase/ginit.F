*
* $Id$
*
* $Log$
* Revision 1.4  1996/09/30 14:25:31  ravndal
* Windows NT related modifications
*
* Revision 1.3  1996/04/15 14:14:55  ravndal
* PATCHY does not set date and time any more
*
* Revision 1.2  1996/03/13 17:20:42  ravndal
* Mod's for the parallel version
*
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 22/02/95  16.02.44  by  S.Giani
*-- Author :
      SUBROUTINE GINIT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       GEANT initialisation routine                             *
C.    *                                                                *
C.    *         IFINIT(1)=1 free                                       *
C.    *               (2)=1 if GZINIT   "  "                           *
C.    *               (3)=1 if GLUND or GLUNDI have been called        *
C.    *               (4)=1 if GHEINI or GPGHEI have been called       *
C.    *               (5)=1 if GHCASC has been called                  *
C.    *               (6)=1 if GLUDKY has been called                  *
C.    *               (7)=1 if GTAU  "   "                             *
C.    *               (8)=1 if GPRELA                                  *
C.    *               (9)=1 if GPCXYZ                                  *
C.    *              (10)=1 if GDRAW                                   *
C.    *              (11)=1 if INIT_GMR                                *
C.    *              (12)=1 if GET_GEANT_STRUCTURE                     *
C.    *              (13)=1 if GPIONS                                  *
C.    *              (14)=1 if GDINIT                                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gclist.inc"
#include "geant321/gcsets.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcpush.inc"
#include "geant321/gctime.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcparm.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcopti.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcrz.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
*
*     COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
*     COMMON/GCONSX/EMMU,PMASS,AVO
      COMMON/GCONST/CONS1(8)
      COMMON/GCONSX/CONS2(3)
#if defined(CERNLIB_USRJMP)
      EXTERNAL      GUDCAY, GUDIGI, GUDTIM, GUFLD , GUHADR, GUIGET,
     +              GUINME, GUINTI, GUKINE, GUNEAR, GUOUT , GUPHAD,
     +              GUSKIP, GUSTEP, GUSWIM, GUTRAK, GUTREV, GUVIEW,
     +              GUPARA
#endif
C
C             ZEBRA system common blocks
C
      COMMON /ZHEADP/IQHEAD(20),IQDATE,IQTIME,IQPAGE,NQPAGE(4)
      COMMON /ZMACH/ NQBITW,NQBITC,NQCHAW,NQLNOR,NQLMAX,NQLPTH,NQRMAX
     +,              IQLPCT,IQNIL
      COMMON /ZSTATE/QVERSN,NQPHAS,IQDBUG,NQDCUT,NQWCUT,NQERR
     +,              NQLOGD,NQLOGM,NQLOCK,NQDEVZ,NQAUGM(6)
      COMMON /ZUNIT/ IQREAD,IQPRNT,IQPR2,IQLOG,IQPNCH,IQTTIN,IQTYPE
C
      CHARACTER*4 MEC(MAXMEC),MEC1(MAXME1),DFLT(2)
      CHARACTER*8 CHVERS,CHDATE

      CHARACTER*4 CDUMMY
#if defined(CERNLIB_MONITOR)
      CHARACTER*32 CHINFO
#endif
      PARAMETER (LEFTM1=MAXME1-9)
      SAVE LOAD
      DATA LOAD/0/
      DATA  DFLT /'    ','XXXX'/
      DATA MEC/'NEXT','MULS','LOSS','FIEL','DCAY','PAIR','COMP','PHOT'
     +        ,'BREM','DRAY','ANNI','HADR','ECOH','EVAP','FISS','ABSO'
     +        ,'ANNH','CAPT','EINC','INHE','MUNU','TOFM','PFIS','SCUT'
     +        ,'RAYL','PARA','PRED','LOOP','NULL','STOP'/
      DATA MEC1/'LABS','LREF','SMAX','SCOR','CKOV','REFL','REFR',
     +          'SYNC','STRA',LEFTM1*'    '/
C.
C.    ------------------------------------------------------------------
C.
C
#if defined (CERNLIB_PARA)
      CALL GPINIT()
#endif
      IDATQQ = 0
      ITIMQQ = 0
      CALL GETVER(CHVERS,CHDATE,GVERSN)
      ZVERSN = QVERSN
      IGDATE = IQDATE
      IGTIME = IQTIME
      LIN    = IQREAD
      IF(IQTTIN.NE.0) LIN=IQTTIN
      LOUT   = IQPRNT
C
      WRITE (CHMAIL,10100) CHVERS,CHDATE
      CALL GMAIL(0,0)
*
      GVERSC = 0.0
#include "geant321/gversc.inc"
* Ignoring t=pass
      IF (GVERSC.NE.0.0) THEN
         WRITE (CHMAIL,10200) GVERSC
         CALL GMAIL(0,0)
      ENDIF
#if defined(CERNLIB_MONITOR)
*
      WRITE(CHINFO,10000) GVERSN, GVERSC
10000 FORMAT(' Version/Cradle: ',F7.4,'/',F7.4)
      CALL GEAMON(1,CHINFO)
#endif
C      Since conversion to CVS library handling, this
C      is ommitted because PATCHY does not anymore 
C      set IDATQQ, ITIMQQ
C	
C      WRITE(CHMAIL,10300)IDATQQ,ITIMQQ
C      CALL GMAIL(0,1)
C
#if defined(CERNLIB_USRJMP)
      JUDCAY = JUMPAD(GUDCAY)
      JUDIGI = JUMPAD(GUDIGI)
* GUDTIM is a function
      JUFLD  = JUMPAD(GUFLD)
      JUHADR = JUMPAD(GUHADR)
      JUIGET = JUMPAD(GUIGET)
      JUINME = JUMPAD(GUINME)
      JUINTI = JUMPAD(GUINTI)
      JUKINE = JUMPAD(GUKINE)
      JUNEAR = JUMPAD(GUNEAR)
      JUOUT  = JUMPAD(GUOUT)
      JUPHAD = JUMPAD(GUPHAD)
      JUSKIP = JUMPAD(GUSKIP)
      JUSTEP = JUMPAD(GUSTEP)
      JUSWIM = JUMPAD(GUSWIM)
      JUTRAK = JUMPAD(GUTRAK)
      JUTREV = JUMPAD(GUTREV)
      JUVIEW = JUMPAD(GUVIEW)
      JUPARA = JUMPAD(GUPARA)
#endif
C
      CALL FFINIT(0)
      NBIT   = NQBITW
      CALL UCTOH(DFLT,IDFLT,4,4)
      CONS1( 1) = PI
      CONS1( 2) = TWOPI
      CONS1( 3) = PIBY2
      CONS1( 4) = DEGRAD
      CONS1( 5) = RADDEG
      CONS1( 6) = CLIGHT
      CONS1( 7) = BIG
      CONS1( 8) = EMASS
      CONS2( 1) = EMMU
      CONS2( 2) = PMASS
      CONS2( 3) = AVO
      DO 10 J=1,MXGKIN
         TOFD(J)  = 0.
         IFLGK(J) = 0
   10 CONTINUE
C
      IGAUTO= 1
      IPAIR = 1
      ICOMP = 1
      IPHOT = 1
      IRAYL = 0
      IBREM = 1
      IHADR = 1
      IANNI = 1
      IDRAY = 1
      IMUNU = 1
      IPFIS = 0
      IDCAY = 1
      ILOSS = 2
      IMULS = 1
      ILABS = -1
      ITCKOV = 0
      ISYNC = 0
      ISTRA = 0
C
      IABAN = 1
      DPHYS1 = IABAN
C
      CUTGAM = 0.001
      CUTELE = 0.001
      CUTHAD = 0.01
      CUTNEU = 0.01
      CUTMUO = 0.01
      TOFMAX = BIG
      DO 20 J=1,5
         GCUTS(J) = 0.
   20 CONTINUE
C
C               The following cuts can be changed by data card CUTS
C               If they are now changed, then the routine GPHYSI
C               will change them respectively to
C               BCUTE=CUTGAM,BCUTM=CUTGAM, DCUTE=CUTELE, DCUTM=CUTELE
C               and PPCUTM=4.*EMASS
C
      DCUTE = BIG
      DCUTM = BIG
      BCUTE = BIG
      BCUTM = BIG
      PPCUTM= BIG
      ISTPAR= 1
      IOPTIM= 0
C
      NCVERT =  5
      NCKINE =  50
      NCJXYZ = 100
      NPVERT = 5
      NPKINE = 10
      NPJXYZ = 200
C
      IKINE = 0
      DO 30 J=1,10
         PKINE(J) = BIG
   30 CONTINUE
      CALL VZERO (IHSET,26)
      CALL VZERO (NHSTA,9)
      CALL VFILL (LHSTA, 180, IDFLT)
      CALL VFILL (LRGET,  40, IDFLT)
      CALL VZERO (NUNITS,6)
      CALL VZERO (IDEBUG,42)
      CALL VZERO (NMATE,9)
      CALL VZERO (NLEVEL,306)
      NALIVE = 0
      NTMSTO = 0
      NJTMAX = 0
      NJTMIN = 0
      NTSTKP = 500
      NTSTKS = 100
*
*-------- Scan parameters defaults
      SCANFL = .FALSE.
      NPHI   = 90
      PHIMIN = 0.
      PHIMAX = 360.
      IPHI1  = 1
      IPHIL  = NPHI
      NTETA  = 90
      TETMID(1) = -10.
      TETMID(2) = 0.
      TETMID(3) = -1.
      TETMAD(1) =  10.
      TETMAD(2) = 180.
      TETMAD(3) = 1.
      MODTET = 1
      CALL VFILL (ISLIST,MSLIST,IDFLT)
      CALL UCTOH(DFLT(2),ISLIST(1),4,4)
      NSLIST = 1
      VSCAN(1) = 0.
      VSCAN(2) = 0.
      VSCAN(3) = 0.
      FACTX0 = 100.
      FACTL  = 10.
      FACTSF = 100.
      FACTR  = 100.
*---      Parametrization cut=0 means no parametrization
      IPARAM  = 0
      DO 40 J=1,5
         PACUTS(J) = 0.
   40 CONTINUE
*---      Size for the primary parametrization stak
      MPSTAK = 1000
*---      Number of particles generated for every shower
      NPGENE = 20
*-------- Scan parameters defaults
C
      RZTAGS(1)='STRUCTUR'
      RZTAGS(2)='TRIG-NR '
      RZTAGS(3)='RUNG-NRT '
      RZTAGS(4)='USER-ID '
      NRGET  = 0
      NRSAVE = 0
      NRECRZ = 1000
C
      IPAOLD =-1
      NUMOLD = 0
C
      IEVENT = 0
      IDRUN  = 1
      NHEAD  = 10
      NTMED  = 100
      NMATE  = 100
      NROTM  = 100
      NPART  = 100
      NEVENT = 10000000
C
      TIMINT = 0.
      TIMEND = 1.
      ITIME  = 1
C
      CALL UCTOH(MEC,NAMEC,4,MAXMEC*4)
      CALL UCTOH(MEC1,NAMEC1,4,MAXME1*4)
      MAXNST=10000
C
C             Constants for energy loss and physics processes
C
      UPWGHT=1.
      NEKBIN=90
      NEK1=NEKBIN+1
      EKMIN=1.E-5
      EKMAX=1.E+4
C
C             Initialize Random number generator
C
      NRNDM(1) = 0
      NRNDM(2) = 0
      CALL GRNDMQ(0,0,1,' ')
#if defined (CERNLIB_PARA)
      CALL GPDEFRNG(1)
#endif
C
C             Constants for multiple scattering (GMUL)
C
      DXM=TWOPI/100.
      XM=-0.5*DXM
      SQ=-0.0099999
      DO 50 I=1,101
         SQ=SQ+0.01
         IF(I.LT.101)SQRMUL(I)=SQRT(-2.*LOG(SQ))
         XM=XM+DXM
         SINMUL(I)=SIN(XM)
   50 COSMUL(I)=COS(XM)
      SQRMUL(101)=0.01
C
C               This piece of code to force loading of default
C               routines from the GEANG file on some machines
C               like VAX.
C
      IF(LOAD.NE.0)THEN
         CALL GWORK (IP1)
#if !defined(CERNLIB_USRJMP)
         CALL GUDCAY
         CALL GUDIGI
         CALL GUDTIM( P2, P3,IP4, P4)
         CALL GUFLD ( P1, P2)
         CALL GUHADR
         CALL GUIGET(IP1,IP2,IP3)
         CALL GUINME( P1, P2, P3,IP4)
         CALL GUINTI
         CALL GUKINE
         CALL GUNEAR(IP1,IP2, P3,IP4)
         CALL GUOUT
         CALL GUPHAD
         CALL GUSKIP(IP1)
         CALL GUSTEP
         CALL GUSWIM( P1, P2, P3, P4)
         CALL GUTRAK
         CALL GUTREV
         CALL GUVIEW(IP1,IP2,CDUMMY,IP4)
         CALL GUPARA
#endif
      ENDIF
C
10100 FORMAT('1*****  GEANT Version ',A8,' Released on ',A8)
10200 FORMAT('0*****  Correction Cradle Version ',F7.4)
10300 FORMAT(' *****  Library compiled on ',I6,' at ',I4,' *****')
      END
