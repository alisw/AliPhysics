*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GNCKOV
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   This routine calculates the number of Cerenkov photons       *
C.    *   produced per centimeter in the current medium.               *
C.    *                                                                *
C.    *   ==>Called by : GRANGI                                        *
C.    *      Authors     R.Jones, F.Carminati ********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
*
      PARAMETER (RFACT=369.81E9)
*
*   ------------------------------------------------------------------
*
* *** See whether we generate at least one photon
*
      THRIND = GETOT/VECT(7)
      IF(Q(JINDEX+NPCKOV).LT.THRIND) THEN
         DP   = 0.
         GE   = 0.
      ELSEIF(Q(JINDEX+1).GE.THRIND) THEN
         PMIN = Q(JTCKOV+2)
         DP   = Q(JTCKOV+NPCKOV+1)-PMIN
         GE   = Q(JCURIN+NPCKOV)
         JMIN = 1
      ELSE
         JMIN = 1
         JMAX = NPCKOV
   10    JMED = (JMIN+JMAX)/2
         IF(Q(JINDEX+JMED).LT.THRIND) THEN
            JMIN = JMED
         ELSE
            JMAX = JMED
         ENDIF
         IF(JMAX-JMIN.GT.1) GO TO 10
         RATIO =
     +   (THRIND-Q(JINDEX+JMIN))/(Q(JINDEX+JMIN+1)-Q(JINDEX+JMIN))
         RATI1 = 1.-RATIO
         PMIN  = Q(JTCKOV+JMIN+1)*RATI1+Q(JTCKOV+JMIN+2)*RATIO
         DP    = Q(JTCKOV+NPCKOV+1)-PMIN
         GEMIN = Q(JCURIN+JMIN)*RATI1+Q(JCURIN+JMIN+1)*RATIO
         GE    = Q(JCURIN+NPCKOV)-GEMIN
      ENDIF
      DNDL = RFACT*(CHARGE**2)*(DP-GE*THRIND**2)
 999  END
