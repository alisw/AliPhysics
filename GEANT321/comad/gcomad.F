*CMZ :          22/01/99  15.54.06  by  Rene Brun
*-- Author :    Federico Carminati   28/11/98
      SUBROUTINE GCOMAD(CHCOMM,IADD)
*#include "geant321/gcasho.inc"
      COMMON/GCASHO/ZMED,AMED,DMED,E0MED,ZSMED(50),ESMED(50),ALFA,
     *             STEPALT,PLIN,PLOG,BE2,PLASM,TRNSMA,
*     *             STEP,PLIN,PLOG,BE2,PLASM,TRNSMA,
     *             BOSC(50),AOSC(50),EOSC(50),ZOSC(50),EMEAN,
     *             CMGO(2000),EMGO,EMGOMI,
     *             NSMED,IOSC(50),NOSC,NMGO,NMGOMA
C
      REAL         ZMED,AMED,DMED,E0MED,ZSMED,ESMED,ALFA,
     *             STEPALT,PLIN,PLOG,BE2,PLASM,TRNSMA,
*     *             STEP,PLIN,PLOG,BE2,PLASM,TRNSMA,
     *             BOSC,AOSC,EOSC,ZOSC,EMEAN,CMGO,EMGO,EMGOMI
      INTEGER      NSMED,IOSC,NOSC,NMGO,NMGOMA
C
#include "geant321/gcbank.inc"
#include "geant321/gccurs.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcfdim.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gchil2.inc"
#include "geant321/gchiln.inc"
#include "geant321/gchvir.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcjump.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gclist.inc"
#include "geant321/gclund.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcmutr.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcnum.inc"
*#include "geant321/gcomis.inc"
*      COMMON/GCOMIS/JUINIT,JUGEOM,JUKINE,JUSTEP,JUOUT,JULAST
      COMMON/GCOMIS/KUINIT,KUGEOM,KUKINE,KUSTEP,KUOUT,KULAST
      DIMENSION JPCOMS(6)
      EQUIVALENCE (JPCOMS,JUINIT)
*
#include "geant321/gconst.inc"
#include "geant321/gcopti.inc"
#include "geant321/gcpara.inc"
#include "geant321/gcparm.inc"
#include "geant321/gcpmxz.inc"
#include "geant321/gcphnr.inc"
#include "geant321/gcphpr.inc"
#include "geant321/gcphrd.inc"
#include "geant321/gcphxs.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcpixe.inc"
#include "geant321/gcpoly.inc"
#include "geant321/gcpush.inc"
#include "geant321/gcrayt.inc"
#include "geant321/gcrz.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcsets.inc"
#include "geant321/gcshpt.inc"
*#include "geant321/gcspee.inc"
*
*      REAL S1,S2,S3,SS1,SS2,SS3,SRAGMX,SRAGMN,
      REAL S1,S2ALT,S3,SS1,SS2,SS3,SRAGMX,SRAGMN,
     +     RAINT1,RAINT2,RMIN1,RMIN2,RMAX1,RMAX2
      INTEGER ISCOP,NTIM,NTFLAG,IOLDCU,ITSTCU,ISUBLI,IPORLI
      INTEGER LPASS,JPORJJ,LEP,JSC
*
*      COMMON/GCSPEE/S1,S2,S3,SS1,SS2,SS3,LEP,IPORLI,ISUBLI,
      COMMON/GCSPEE/S1,S2ALT,S3,SS1,SS2,SS3,LEP,IPORLI,ISUBLI,
     +              SRAGMX,SRAGMN,RAINT1,RAINT2,RMIN1,RMIN2,
     +              RMAX1,RMAX2,JPORJJ,ITSTCU,IOLDCU,ISCOP,
     +              NTIM,NTFLAG,LPASS,JSC
*

#include "geant321/gcstak.inc"
#include "geant321/gcstra.inc"
#include "geant321/gctime.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcursb.inc"
#include "geant321/gcvdma.inc"
#include "geant321/gcvol1.inc"
#include "geant321/gcvol2.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcxlun.inc"
#include "geant321/gfkdis.inc"
#include "geant321/gsecti.inc"
      COMMON / QUEST / IQUEST(100)
*
      CHARACTER *(*) CHCOMM
#if defined(CERNLIB_DECS)
      INTEGER*8 GCADDI, GCADDF, GCADDL, GCADDC, IADD
#else          
      INTEGER GCADDI, GCADDF, GCADDL, GCADDC
#endif
*
      IF(CHCOMM.EQ.'GCASHO') THEN
         IADD=GCADDF(ZMED)
      ELSE IF(CHCOMM.EQ.'IQ') THEN
         IADD=GCADDI(LQ(8))
      ELSE IF(CHCOMM.EQ.'LQ') THEN
         IADD=GCADDF(FENDQ(16))
      ELSE IF(CHCOMM.EQ.'GCCURS') THEN
         IADD=GCADDI(INTFLA)
      ELSE IF(CHCOMM.EQ.'GCCUTS') THEN
         IADD=GCADDF(CUTGAM)
      ELSE IF(CHCOMM.EQ.'GCDRAW') THEN
         IADD=GCADDI(NUMNOD)
      ELSE IF(CHCOMM.EQ.'GCFDIM') THEN
         IADD=GCADDI(NPROC)
      ELSE IF(CHCOMM.EQ.'GCFLAG') THEN
         IADD=GCADDI(IDEBUG)
      ELSE IF(CHCOMM.EQ.'GCGOBJ') THEN
         IADD=GCADDI(IST)
      ELSE IF(CHCOMM.EQ.'GCBLIM') THEN
         IADD=GCADDI(IHOLE)
      ELSE IF(CHCOMM.EQ.'GCHIL2') THEN
         IADD=GCADDI(LARETT)
      ELSE IF(CHCOMM.EQ.'GCHILN') THEN
         IADD=GCADDI(LARECG)
      ELSE IF(CHCOMM.EQ.'GCHVIR') THEN
         IADD=GCADDI(JVIRT)
      ELSE IF(CHCOMM.EQ.'GCJLOC') THEN
         IADD=GCADDI(NJLOC)
      ELSE IF(CHCOMM.EQ.'GCJLCK') THEN
         IADD=GCADDI(NJLCK)
      ELSE IF(CHCOMM.EQ.'GCJUMP') THEN
         IADD=GCADDI(JUDCAY)
      ELSE IF(CHCOMM.EQ.'GCKINE') THEN
         IADD=GCADDI(IKINE)
      ELSE IF(CHCOMM.EQ.'GCKING') THEN
         IADD=GCADDI(KCASE)
      ELSE IF(CHCOMM.EQ.'GCKIN2') THEN
         IADD=GCADDI(NGPHOT)
      ELSE IF(CHCOMM.EQ.'GCKIN3') THEN
         IADD=GCADDF(GPOS)
      ELSE IF(CHCOMM.EQ.'GCLINK') THEN
         IADD=GCADDI(JDIGI)
      ELSE IF(CHCOMM.EQ.'GCLIST') THEN
         IADD=GCADDI(NHSTA)
      ELSE IF(CHCOMM.EQ.'GCLUND') THEN
         IADD=GCADDI(IFLUND)
      ELSE IF(CHCOMM.EQ.'GCMATE') THEN
         IADD=GCADDI(NMAT)
      ELSE IF(CHCOMM.EQ.'GCMULO') THEN
         IADD=GCADDF(SINMUL)
      ELSE IF(CHCOMM.EQ.'GCMUTR') THEN
         IADD=GCADDI(NCVOLS)
      ELSE IF(CHCOMM.EQ.'GCMUTC') THEN
         IADD=GCADDC(GNASH)
      ELSE IF(CHCOMM.EQ.'GCMZFO') THEN
         IADD=GCADDI(IOMATE)
      ELSE IF(CHCOMM.EQ.'GCNUM') THEN
         IADD=GCADDI(NMATE)
      ELSE IF(CHCOMM.EQ.'GCOMIS') THEN
         IADD=GCADDI(KUINIT)
      ELSE IF(CHCOMM.EQ.'GCONST') THEN
         IADD=GCADDF(PI)
      ELSE IF(CHCOMM.EQ.'GCOPTI') THEN
         IADD=GCADDI(IOPTIM)
      ELSE IF(CHCOMM.EQ.'GCPARA') THEN
         IADD=GCADDF(EPSIX0)
      ELSE IF(CHCOMM.EQ.'GCPARM') THEN
         IADD=GCADDI(IPARAM)
      ELSE IF(CHCOMM.EQ.'GCPHNR') THEN
         IADD=GCADDI(IGNRFN)
      ELSE IF(CHCOMM.EQ.'GCPHPR') THEN
         IADD=GCADDF(GFLUPR)
      ELSE IF(CHCOMM.EQ.'GCPHRD') THEN
         IADD=GCADDF(GPHRAT)
      ELSE IF(CHCOMM.EQ.'GCPXRN') THEN
         IADD=GCADDC(CRNGUP)
      ELSE IF(CHCOMM.EQ.'GCPXCF') THEN
         IADD=GCADDF(COFS)
      ELSE IF(CHCOMM.EQ.'GCPHYS') THEN
         IADD=GCADDI(IPAIR)
      ELSE IF(CHCOMM.EQ.'GCPHLT') THEN
         IADD=GCADDI(ILABS)
      ELSE IF(CHCOMM.EQ.'GCPIXE') THEN
         IADD=GCADDI(LIMPRE)
      ELSE IF(CHCOMM.EQ.'GCPOLY') THEN
         IADD=GCADDI(IZSEC)
      ELSE IF(CHCOMM.EQ.'GCPUSH') THEN
         IADD=GCADDI(NCVERT)
      ELSE IF(CHCOMM.EQ.'GCRAYT') THEN
         IADD=GCADDI(INTEN)
      ELSE IF(CHCOMM.EQ.'GCRZ1') THEN
         IADD=GCADDI(NRECRZ)
      ELSE IF(CHCOMM.EQ.'GCRZ2') THEN
         IADD=GCADDC(RZTAGS)
      ELSE IF(CHCOMM.EQ.'GCSCAL') THEN
         IADD=GCADDI(ISLINK)
      ELSE IF(CHCOMM.EQ.'GCSCAN') THEN
         IADD=GCADDL(SCANFL)
      ELSE IF(CHCOMM.EQ.'GCSCAC') THEN
         IADD=GCADDC(SFIN)
      ELSE IF(CHCOMM.EQ.'GCSETS') THEN
         IADD=GCADDI(IHSET)
      ELSE IF(CHCOMM.EQ.'GCSHPT') THEN
         IADD=GCADDI(NSHLST)
      ELSE IF(CHCOMM.EQ.'GCSPEE') THEN
         IADD=GCADDF(S1)
      ELSE IF(CHCOMM.EQ.'GCSTAK') THEN
         IADD=GCADDI(NJTMAX)
      ELSE IF(CHCOMM.EQ.'GCSTRA') THEN
         IADD=GCADDF(EMAX)
      ELSE IF(CHCOMM.EQ.'GCTIME') THEN
         IADD=GCADDF(TIMINT)
      ELSE IF(CHCOMM.EQ.'GCTMED') THEN
         IADD=GCADDI(NUMED)
      ELSE IF(CHCOMM.EQ.'GCTRAK') THEN
         IADD=GCADDF(VECT)
      ELSE IF(CHCOMM.EQ.'GCTPOL') THEN
         IADD=GCADDF(POLAR)
      ELSE IF(CHCOMM.EQ.'GCUNIT') THEN
         IADD=GCADDI(LIN)
      ELSE IF(CHCOMM.EQ.'GCURSB') THEN
         IADD=GCADDI(NUMNDS)
      ELSE IF(CHCOMM.EQ.'GCVDMA') THEN
         IADD=GCADDI(NVMANY)
      ELSE IF(CHCOMM.EQ.'GCVOL1') THEN
         IADD=GCADDI(NLEVL1)
      ELSE IF(CHCOMM.EQ.'GCVOL2') THEN
         IADD=GCADDI(NLEVE2)
      ELSE IF(CHCOMM.EQ.'GCVOLU') THEN
         IADD=GCADDI(NLEVEL)
      ELSE IF(CHCOMM.EQ.'GCXLUN') THEN
         IADD=GCADDI(LUNIT)
      ELSE IF(CHCOMM.EQ.'GFKDIS') THEN
         IADD=GCADDF(ZINE)
      ELSE IF(CHCOMM.EQ.'GSECTI') THEN
         IADD=GCADDF(AIEL)
      ELSE IF(CHCOMM.EQ.'GCTLIT') THEN
         IADD=GCADDF(THRIND)
      ELSE IF(CHCOMM.EQ.'QUEST') THEN
         IADD=GCADDI(IQUEST)
      ELSE IF(CHCOMM.EQ.'GCBANK') THEN
         IADD=GCADDI(NZEBRA)
      ELSE
         WRITE(6,123) CHCOMM
         IADD=0
      ENDIF
*
 123  FORMAT(' *** GCOMAD: cannot load common ',A)
      END
