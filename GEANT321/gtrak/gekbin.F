*
* $Id$
*
* $Log$
* Revision 1.2  1998/01/21 16:44:37  gunter
* Add protection of iekbin getting too big due to numerical inacurecies
*
* Revision 1.1.1.1  1995/10/24 10:21:40  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GEKBIN
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to find bin number in kinetic energy table       *
C.    *       stored in ELOW(NEKBIN)                                   *
C.    *                                                                *
C.    *    ==>Called by : GBREME,GBREMM,GCOMP,GDRAY,GHEISH,GPAIRM,     *
C.    *                   GLTRAC,GTELEC,GTHADR,GTMUON                  *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctrak.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
C.
C.    ------------------------------------------------------------------
C.
      EK=GEKIN
      IF (ITRTYP.EQ.4.OR.ITRTYP.EQ.8) EK=GEKIN*PMASS/AMASS
      IF (EK.LT.ELOW(1)) THEN
         IEKBIN=1
         GEKRAT=0.
      ELSEIF (EK.GE.ELOW(NEK1)) THEN
         IEKBIN=NEKBIN
         GEKRAT=1.
      ELSE
         IEKBIN=GEKA*LOG10(EK) + GEKB
         IF (IEKBIN.GT.NEKBIN ) IEKBIN=NEKBIN
         GEKRAT=(EK-ELOW(IEKBIN))/(ELOW(IEKBIN+1)-ELOW(IEKBIN))
         GEKRAT=MIN(1.,MAX(GEKRAT,0.))
      ENDIF
C
      END
