*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GBREMM
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates discrete hard BREMSSTRAHLUNG by muons.              *
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GTMUON                                       *
C.    *       Author     L.Urban  *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gccuts.inc"
      DIMENSION RNDM(2)
C.
C.    ------------------------------------------------------------------
C.
C             Ensure cut-off avoids infra-red catastrophe.
C
      EEM1   = GETOT
      IF (GEKIN.LE.BCUTM) GO TO 900
      KCASE  = NAMEC(9)
C
      C1   = Q(JPROB+31)
      C2   = Q(JPROB+32)
      C3   = Q(JPROB+33)
      VC   = BCUTM/EEM1
      VM   = 1.-C3/EEM1
      IF(VM.LE.0.)GO TO 900
      CREJ = C2/EEM1
C
  50  CALL GRNDM(RNDM,2)
      R    = RNDM(1)
      V    = VC*(VM/VC)**R
C
C             Compute rejection function
C
      V1   = 1.-V
      F1   = C1-LOG(1.+CREJ*V/V1)
      SCREJ= (V1+0.75*V*V)*F1/C1
      IF(RNDM(2).GT.SCREJ)GO TO 50
C
C             Photon energy
C
      EGAMMA=EEM1*V
C
C        CUT ON ENERGY THRESHOLD ?
C
      IF((IBREM.NE.1).OR.(EGAMMA.LE.CUTGAM)) THEN
          NGKINE = 0
          DESTEP = DESTEP + EGAMMA
          GO TO 800
      ENDIF
C
C             Generate emitted photon angles with respect to a Z-axis
C             defined along parent track. PHI is generated isotropically
C             and THETA is assigned a universal angular distribution
C
      THETA  = GBTETH(EEM1, AMASS, V)*AMASS/EEM1
      SINTH  = SIN(THETA)
      COSTH  = COS(THETA)
      CALL GRNDM(RNDM,1)
      PHI    = TWOPI*RNDM(1)
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
C
C             Polar co-ordinates to momentum components.
C
      NGKINE = 1
      GKIN(1,1)=EGAMMA*SINTH*COSPHI
      GKIN(2,1)=EGAMMA*SINTH*SINPHI
      GKIN(3,1)=EGAMMA*COSTH
      GKIN(4,1)=EGAMMA
      GKIN(5,1)=1.
      TOFD(1)  =0.
      GPOS(1,1) = VECT(1)
      GPOS(2,1) = VECT(2)
      GPOS(3,1) = VECT(3)
C
C             Rotate photon into GEANT system
C
      CALL GVROT(VECT(4),GKIN)
C
C             Correct track for lost energy
C
  800 CONTINUE
      GEKIN = GEKIN - EGAMMA
      GETOT = GEKIN + AMASS
      VECT(7)=SQRT (ABS((GETOT+AMASS)*GEKIN))
      CALL GEKBIN
C
C             Update probabilities
C
  900 CALL GRNDM(RNDM,1)
      ZINTBR=-LOG(RNDM(1))
      SLBREM=SLENG
      STEPBR=BIG
C
      END
