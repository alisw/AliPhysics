*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE GMUNUI
C
C *** CALCULATION OF MUON-NUCLEUS CROSS SECTIONS IN CURRENT MATERIAL ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C CALLED BY : GPHYSI
C ORIGIN : F.CARMINATI, H.FESEFELDT
C
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
C
C
C --- GET MUON MASS FROM GEANT PARTICLE DATA BANK ---
      LNVE=LQ(JPART-5)
      HMAMU=Q(LNVE+7)
C
C --- Infinite mean free patch for EKIN < 1 GeV (HJT/NVE 27-sep-1990) ---
      SIG=0.
      IF(JMUNU.EQ.0)GO TO 9999
      IF (ELOW(IEKBIN) .LT. 1.) GO TO 90
C
      IF(Z.LT.1.)GO TO 90
      EN     = ELOW(IEKBIN)+HMAMU
      SIG    = 0.0003
      IF(EN.GT.30.) SIG = SIG*(EN/30.)**0.25
      SIG    = SIG*DENS*AVO*1E-3
C
  90  IF(SIG.GT.0.)THEN
         Q(JMUNU+IEKBIN)=1./SIG
      ELSE
         Q(JMUNU+IEKBIN)=BIG
      ENDIF
C
 9999 CONTINUE
      END
