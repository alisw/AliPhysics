*
* $Id$
*
* $Log$
* Revision 1.2  1996/09/30 13:40:37  ravndal
* IABAN initialization
*
* Revision 1.1.1.1  1995/10/24 10:20:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 22/02/95  08.50.24  by  S.Ravndal
*-- Author :
      SUBROUTINE GFFGO
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to define and read GEANT/FFREAD data cards       *
C.    *      If user data cards have been defined via FFKEY            *
C.    *      they will be read as well                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gclist.inc"
#include "geant321/gcparm.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctime.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcopti.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcrz.inc"
C.
C.    ------------------------------------------------------------------
C.
      IABAN = NINT(DPHYS1)
      CALL FFKEY ('ANNI',IANNI , 1,'INTEGER')
      CALL FFKEY ('BREM',IBREM , 1,'INTEGER')
      CALL FFKEY ('COMP',ICOMP , 1,'INTEGER')
      CALL FFKEY ('CUTS',CUTGAM,16,'REAL')
      CALL FFKEY ('DEBU',IDEMIN, 3,'INTEGER')
      CALL FFKEY ('DCAY',IDCAY , 1,'INTEGER')
      CALL FFKEY ('DRAY',IDRAY , 1,'INTEGER')
      CALL FFKEY ('LABS',ILABS , 1,'INTEGER')
      CALL FFKEY ('GEOM',LGEOM ,20,'MIXED')
      CALL FFKEY ('GET ',LGET  ,20,'MIXED')
      CALL FFKEY ('HADR',IHADR , 1,'INTEGER')
      CALL FFKEY ('HSTA',LHSTA ,20,'MIXED')
      CALL FFKEY ('KINE',IKINE ,11,'MIXED')
      CALL FFKEY ('LOSS',ILOSS , 1,'INTEGER')
      CALL FFKEY ('MULS',IMULS , 1,'INTEGER')
      CALL FFKEY ('MUNU',IMUNU , 1,'INTEGER')
      CALL FFKEY ('PAIR',IPAIR , 1,'INTEGER')
*     CALL FFKEY ('PATR',NJTMAX, 4,'INTEGER')
      CALL FFKEY ('SORD',ISTORD, 1,'INTEGER')
      CALL FFKEY ('PFIS',IPFIS , 1,'INTEGER')
      CALL FFKEY ('PHOT',IPHOT , 1,'INTEGER')
      CALL FFKEY ('PLOT',LPLOT ,20,'MIXED')
      CALL FFKEY ('PRIN',LPRIN ,20,'MIXED')
      CALL FFKEY ('RAYL',IRAYL , 1,'INTEGER')
      CALL FFKEY ('RGET',LRGET ,20,'MIXED')
      CALL FFKEY ('RSAV',LRSAVE,20,'MIXED')
      CALL FFKEY ('RNDM',NRNDM , 2,'INTEGER')
      CALL FFKEY ('RUNG',IDRUN , 2,'INTEGER')
      CALL FFKEY ('SAVE',LSAVE ,20,'MIXED')
      CALL FFKEY ('SETS',LSETS ,20,'MIXED')
      CALL FFKEY ('STAT',LSTAT ,20,'MIXED')
      CALL FFKEY ('SWIT',ISWIT ,10,'INTEGER')
      CALL FFKEY ('TIME',TIMINT, 3,'MIXED')
      CALL FFKEY ('TRIG',NEVENT, 1,'INTEGER')
      CALL FFKEY ('VIEW',LVIEW ,20,'MIXED')
      CALL FFKEY ('ERAN',EKMIN , 3,'MIXED')
      CALL FFKEY ('AUTO',IGAUTO, 1,'INTEGER')
      CALL FFKEY ('OPTI',IOPTIM, 1,'INTEGER')
      CALL FFKEY ('CKOV',ITCKOV, 1,'INTEGER')
      CALL FFKEY ('SYNC',ISYNC,  1,'INTEGER')
      CALL FFKEY ('STRA',ISTRA,  1,'INTEGER')
      CALL FFKEY ('ABAN',IABAN,  1,'INTEGER')
*--------------- SCAN/Parametrize cards
      CALL FFKEY ('SCAN',SCANFL, 8,'MIXED')
      CALL FFKEY ('SCAL',ISLIST,MSLIST,'MIXED')
      CALL FFKEY ('SCAP',VSCAN , 6,'REAL')
      CALL FFKEY ('PCUT',IPARAM,  6,'MIXED')
      CALL FFKEY ('PNUM',MPSTAK,  2,'MIXED')
C
C             Now read data cards
C
      CALL FFGO
C
C             Get some parameters from the data cards
C
      CALL GETNUM (LHSTA ,NHSTA)
      CALL GETNUM (LGET  ,NGET)
      CALL GETNUM (LSAVE ,NSAVE)
      CALL GETNUM (LRGET ,NRGET)
      CALL GETNUM (LRSAVE,NRSAVE)
      CALL GETNUM (LSETS ,NSETS)
      CALL GETNUM (LPRIN ,NPRIN)
      CALL GETNUM (LGEOM ,NGEOM)
      CALL GETNUM (LVIEW ,NVIEW)
      CALL GETNUM (LPLOT ,NPLOT)
      CALL GETNUM (LSTAT ,NSTAT)
      CALL GETNUM (ISLIST,NSLIST)
*
      IF(ILABS.NE.0) THEN
         IF(ITCKOV.NE.0) ILABS=1
      ENDIF
*
      IF (NINT(DPHYS1).NE.IABAN) THEN
          DPHYS1 = IABAN
      END IF
*
      IF(MODTET.LT.1.OR.MODTET.GT.MAXMDT) THEN
         WRITE(LOUT,10000) MODTET
10000    FORMAT(' MODTET = ',I2,' out of range - 1 assumed')
         MODTET=1
      ENDIF
      TETMIN = TETMID(MODTET)
      TETMAX = TETMAD(MODTET)
      IF(IDEVT.GT.0)IDEVT=IDEVT-1
C
C             Set IDEBUG flag for initialisation phase
C
      IF(IDEMIN.LT.0.AND.IDEMAX.GE.0)THEN
         IDEMIN=-IDEMIN
         IDEBUG=1
      ENDIF
C
C             Initialise the random number generator
C
      IF(NRNDM(2).NE.0)THEN
         CALL GRNDMQ(NRNDM(1),NRNDM(2),1,'S')
      ELSEIF(NRNDM(1).GT.0)THEN
         ISEQ=NRNDM(1)
         CALL GRNDMQ(NRNDM(1),NRNDM(2),ISEQ,'Q')
         CALL GRNDMQ(NRNDM(1),NRNDM(2),ISEQ,'S')
      ENDIF
C
  999 END
