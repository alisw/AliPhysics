*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/30 11:25:57  ravndal
* Implicit Fortran data type convention overlooked
*
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 17/01/95  18.01.49  by  S.Giani
*-- Author :
      SUBROUTINE GXGEOM
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Geometry commands                                         *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                  P.Zanarini  **********                        *
C.    *                  N.Hoimyr 1992    **********                   *
C.    *                  S.Giani  1992    **********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/pawc.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#include "geant321/gclist.inc"
#include "geant321/gconsp.inc"
      CHARACTER*4 CHNAME,CHISH,CHPAR,IONLY,CYESNO,CHAX
      CHARACTER*32 CHPATL,VNAME
      CHARACTER*24 FNAME
      CHARACTER*4 ANAME
      CHARACTER*24 INST,SITE,DEPT,RESP
      DIMENSION ARRAY(50),UBUF(1)
      DATA UBUF/1./
C.
C.    ------------------------------------------------------------------
C.
      IWKSTY = IGIWTY(1)
      CALL KUPATL(CHPATL,NPAR)
*
      IF(CHPATL.EQ.'OPTI') THEN
         CALL KUGETI(IOPT)
         CALL GOPTIM(IOPT)
      ELSEIF (CHPATL.EQ.'SVOLU') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETC(CHISH,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETI(NP)
         CALL KUGETV(VNAME,LPAR,LLL)
         CALL GSVOLU(CHNAME,CHISH,NUMED,QQ(LPAR),NP,IVOLU)
*
      ELSEIF (CHPATL.EQ.'SPOS') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(N)
         CALL KUGETC(CHPAR,NCH)
         CALL KUGETR(X0)
         CALL KUGETR(Y0)
         CALL KUGETR(Z0)
         CALL KUGETI(IROT)
         CALL KUGETC(IONLY,NCH)
         CALL GSPOS(CHNAME,N,CHPAR,X0,Y0,Z0,IROT,IONLY)
*
      ELSEIF (CHPATL.EQ.'SDVN') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETC(CHPAR,NCH)
         CALL KUGETI(NDIV)
         CALL KUGETC(CHAX,NCH)
         IF (CHAX.EQ.'X'.OR.CHAX.EQ.'1') THEN
            IAX=1
         ELSEIF (CHAX.EQ.'Y'.OR.CHAX.EQ.'2') THEN
            IAX=2
         ELSEIF (CHAX.EQ.'Z'.OR.CHAX.EQ.'3') THEN
            IAX=3
         ENDIF
         IF(IAX.LE.3.AND.IAX.GE.1) THEN
            CALL GSDVN(CHNAME,CHPAR,NDIV,IAX)
         ELSE
            WRITE(CHMAIL,10000)
10000       FORMAT(' *** GXGEOM *** Wrong value of IAX')
            CALL GMAIL(0,0)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'PVOLU') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPVOLX(NUMB)
         ELSE
            CALL GPVOLU(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'SROTM') THEN
         CALL KUGETI(N)
         CALL KUGETR(THETA1)
         CALL KUGETR(PHI1)
         CALL KUGETR(THETA2)
         CALL KUGETR(PHI2)
         CALL KUGETR(THETA3)
         CALL KUGETR(PHI3)
         CALL GSROTM(N,THETA1,PHI1,THETA2,PHI2,THETA3,PHI3)
*
      ELSEIF (CHPATL.EQ.'PROTM') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPROTX(NUMB)
         ELSE
            CALL GPROTM(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'STMED') THEN
         NMED=1
         NMAT=1
         IFIELD=0
         EPSIL=0.01
         ISVOL=0
         FIELDM=0.
         TMAXFD=0.01
         STEMAX=BIG
         DEEMAX=0.01
         STMIN=0.1
         CALL KUGETI(NMED)
         CALL KUGETC(VNAME,NCH)
         CALL KUGETI(NMAT)
         CALL KUGETI(ISVOL)
         CALL KUGETI(IFIELD)
         CALL KUGETR(FIELDM)
         CALL KUGETR(TMAXFD)
         CALL KUGETR(STEMAX)
         CALL KUGETR(DEEMAX)
         CALL KUGETR(EPSIL)
         CALL KUGETR(STMIN)
         CALL GSTMED(NMED,VNAME,NMAT,ISVOL,IFIELD,FIELDM,TMAXFD,
     +   STEMAX,DEEMAX,EPSIL,STMIN,UBUF,0)
*
      ELSEIF (CHPATL.EQ.'PTMED') THEN
         CALL KUGETI(NUMB)
         IF(IWKSTY.GE.1.AND.IWKSTY.LE.10) THEN
            CALL GPTMEX(NUMB)
         ELSE
            CALL GPTMED(NUMB)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'EDITV') THEN
         CALL KUGETI(NUM)
         IF(NUM.LE.0)THEN
            CALL GGCLOS
         ELSE
            CALL GEDITV(NUM)
         ENDIF
*
      ELSEIF(CHPATL.EQ.'SBOX') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(HALFX)
         CALL KUGETR(HALFY)
         CALL KUGETR(HALFZ)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=HALFX
         ARRAY(2)=HALFY
         ARRAY(3)=HALFZ
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=3
         ENDIF
         CALL GSVOLU(CHNAME,'BOX ',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'STRD1') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(HLFDWX)
         CALL KUGETR(HLFUPX)
         CALL KUGETR(HALFY)
         CALL KUGETR(HALFZ)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=HLFDWX
         ARRAY(2)=HLFUPX
         ARRAY(3)=HALFY
         ARRAY(4)=HALFZ
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=4
         ENDIF
         CALL GSVOLU(CHNAME,'TRD1',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'STRD2') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(HLFDWX)
         CALL KUGETR(HLFUPX)
         CALL KUGETR(HLFDWY)
         CALL KUGETR(HLFUPY)
         CALL KUGETR(HALFZ)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=HLFDWX
         ARRAY(2)=HLFUPX
         ARRAY(3)=HLFDWY
         ARRAY(4)=HLFUPY
         ARRAY(5)=HALFZ
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=5
         ENDIF
         CALL GSVOLU(CHNAME,'TRD2',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'STUBE') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(XINRAD)
         CALL KUGETR(OUTRAD)
         CALL KUGETR(HALFZ)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=XINRAD
         ARRAY(2)=OUTRAD
         ARRAY(3)=HALFZ
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=3
         ENDIF
         CALL GSVOLU(CHNAME,'TUBE',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'STUBS') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(XINRAD)
         CALL KUGETR(OUTRAD)
         CALL KUGETR(HALFZ)
         CALL KUGETR(SPHI)
         CALL KUGETR(EPHI)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=XINRAD
         ARRAY(2)=OUTRAD
         ARRAY(3)=HALFZ
         ARRAY(4)=SPHI
         ARRAY(5)=EPHI
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=5
         ENDIF
         CALL GSVOLU(CHNAME,'TUBS',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'SCONE') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(XINRDW)
         CALL KUGETR(OUTRDW)
         CALL KUGETR(XINRUP)
         CALL KUGETR(OUTRUP)
         CALL KUGETR(HALFZ)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=XINRDW
         ARRAY(2)=OUTRDW
         ARRAY(3)=XINRUP
         ARRAY(4)=OUTRUP
         ARRAY(5)=HALFZ
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=5
         ENDIF
         CALL GSVOLU(CHNAME,'CONE',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'SCONS') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(XINRDW)
         CALL KUGETR(OUTRDW)
         CALL KUGETR(XINRUP)
         CALL KUGETR(OUTRUP)
         CALL KUGETR(HALFZ)
         CALL KUGETR(SPHI)
         CALL KUGETR(EPHI)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=XINRDW
         ARRAY(2)=OUTRDW
         ARRAY(3)=XINRUP
         ARRAY(4)=OUTRUP
         ARRAY(5)=HALFZ
         ARRAY(6)=SPHI
         ARRAY(7)=EPHI
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=7
         ENDIF
         CALL GSVOLU(CHNAME,'CONS',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'SSPHE') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(XINRAD)
         CALL KUGETR(OUTRAD)
         CALL KUGETR(SPHI)
         CALL KUGETR(EPHI)
         CALL KUGETR(STHETA)
         CALL KUGETR(ETHETA)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=XINRAD
         ARRAY(2)=OUTRAD
         ARRAY(3)=SPHI
         ARRAY(4)=EPHI
         ARRAY(5)=STHETA
         ARRAY(6)=ETHETA
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=6
         ENDIF
         CALL GSVOLU(CHNAME,'SPHE',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF(CHPATL.EQ.'SPARA') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETI(NUMED)
         CALL KUGETR(HALFX)
         CALL KUGETR(HALFY)
         CALL KUGETR(HALFZ)
         CALL KUGETR(AXIS)
         CALL KUGETR(PHI)
         CALL KUGETR(THETA)
         CALL KUGETC(CYESNO,NCHAR)
         ARRAY(1)=HALFX
         ARRAY(2)=HALFY
         ARRAY(3)=HALFZ
         ARRAY(4)=AXIS
         ARRAY(5)=PHI
         ARRAY(6)=THETA
         IF(CYESNO.EQ.'YES')THEN
            NUMP=0
         ELSE
            NUMP=6
         ENDIF
         CALL GSVOLU(CHNAME,'PARA',NUMED,ARRAY,NUMP,IVOLU)
*
      ELSEIF (CHPATL.EQ.'CADINT') THEN
         CALL KUGETS(FNAME,NCH)
         CALL KUGETC(ANAME,NCH)
         CALL KUGETI(NBINS)
         CALL KUGETI(LUNIT)
         CALL KUGETI(LUNIT2)
         CALL KUGETS(INST,NCH)
         CALL KUGETS(SITE,NCH)
         CALL KUGETS(DEPT,NCH)
         CALL KUGETS(RESP,NCH)
         CALL GTXSET(FNAME,ANAME,NBINS,LUNIT,LUNIT2,
     +INST,SITE,DEPT,RESP)
*
      ELSEIF (CHPATL.EQ.'WEUCLID') THEN
         CALL KUGETI(LUN)
         CALL KUGETS(FNAME,NCH)
         CALL KUGETC(CHNAME,NCH1)
         CALL KUGETI(NUMBER)
         CALL KUGETI(NLEVEL)
         CALL GWEUCL (LUN,FNAME(1:NCH),CHNAME,NUMBER,NLEVEL)
*
      ELSEIF (CHPATL.EQ.'REUCLID') THEN
         CALL KUGETI(LUN)
         CALL KUGETS(FNAME,NCH)
         CALL GREUCL (LUN,FNAME(1:NCH))
*
      ENDIF
*
      END
