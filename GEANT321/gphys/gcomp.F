*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GCOMP
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates photon-electron COMPTON scattering.                 *
C.    *                                                                *
C.    *  The scattered photon energy is sampled using the quantum-     *
C.    *  mechanical KLEIN-NISHINA formula. For this, the random       *
C.    *  number techniques of BUTCHER and MESSEL(NUC. PHYS.20(1960),   *
C.    *  15) are employed.                                             *
C.    *  NOTE :                                                        *
C.    *  (1) Effects due to binding of atomic electrons are            *
C.    *      ignored(recoil electron energy assumed large compared     *
C.    *      with binding energy).                                     *
C.    *                                                                *
C.    *    ==>Called by : GTGAMA                                       *
C.    *       Authors    G.Patrick, L.Urban  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"
      DIMENSION PGAM(3)
      REAL RNDM(5)
      LOGICAL ROTATE
C.
C.    ------------------------------------------------------------------
C.
      KCASE  = NAMEC(7)
      EGAM1  = VECT(7)
      EZERO  = EGAM1/EMASS
      EMINI  = 1.+2.*EZERO
      EMIN   = 1./EMINI
      DSIG1  = LOG(EMINI)
      DSIG2  = 0.5*(1.-EMIN*EMIN)
      DSIGT  = DSIG1+DSIG2
C
C             Decide which part of F(E)=1/E+E to sample from.
C
   10 CALL GRNDM(RNDM,3)
      IF (DSIG1.LT.DSIGT*RNDM(1))THEN
C
C             Sample from F2(E) distribution.
C
         BRD    = RNDM(2)
         CALL GRNDM(RNDM(4),1)
         IF (EZERO.GE.(EZERO+1.)*RNDM(4))THEN
            CALL GRNDM(RNDM(5),1)
            BRD    = MAX(BRD,RNDM(5))
         ENDIF
C
         BR     = EMIN+(1.-EMIN)*BRD
      ELSE
         BR     = EMIN*EXP(DSIG1*RNDM(2))
      ENDIF
C
C             Scattered photon energy.
C
      EGAM2  = BR*EGAM1
C
C             Calculate rejection function G(E).
C
      T      = EMASS*(1.-BR)/EGAM2
      SINTH  = MAX(0.,T*(2.-T))
      REJ    = 1.0-(BR*SINTH)/(1.+BR*BR)
      IF (RNDM(3).GT.REJ)                        GO TO 10
C
C             Successful sampling of scattered photon.
C
C             CUTS ON ENERGY THRESHOLDS ?
C
      TEL     = EGAM1-EGAM2
      IF((EGAM2.LE.CUTGAM).AND.(TEL.LE.CUTELE)) THEN
         ISTOP = 2
         VECT(7) = 0.
         GEKIN = 0.
         GETOT = 0.
         NGKINE = 0.
         DESTEP = DESTEP + EGAM1
         RETURN
      ENDIF
C
C             Generate photon angles with respect to a Z-axis
C             defined along the parent photon.
C             PHI is generated isotropically
C
      SINTH  = SQRT(SINTH)
      COSTH  = 1.-T
      CALL GRNDM(RNDM,1)
      PHI    = TWOPI*RNDM(1)
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
C
C             Polar co-ordinates to momentum components.
C
      PGAM(1) = EGAM2*SINTH*COSPHI
      PGAM(2) = EGAM2*SINTH*SINPHI
      PGAM(3) = EGAM2*COSTH
C
C             Momentum vector of recoil electron.
C
      NGKINE = 1
      EEL    = TEL + EMASS
      GKIN(1,1) = -PGAM(1)
      GKIN(2,1) = -PGAM(2)
      GKIN(3,1) = EGAM1-PGAM(3)
      GKIN(4,1)=EEL
      GKIN(5,1)=3
      TOFD(NGKINE)=0.
      GPOS(1,1) = VECT(1)
      GPOS(2,1) = VECT(2)
      GPOS(3,1) = VECT(3)
C
C             Rotate electron and scattered photon into GEANT system
C
      CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
      IF(ROTATE) THEN
         CALL GDROT(PGAM(1),COSTH,SINTH,COSPH,SINPH)
         CALL GDROT(GKIN,COSTH,SINTH,COSPH,SINPH)
      ENDIF
C
C             Correct photon for energy lost and scattered angle
C
      DO 60 I=1,3
   60 VECT(I+3) = PGAM(I)/EGAM2
      VECT(7) = EGAM2
      GETOT = EGAM2
      GEKIN = EGAM2
      CALL GEKBIN
C
C              Stop electron ?
C
      IF((ICOMP.NE.1).OR.(TEL.LE.CUTELE)) THEN
         NGKINE = 0
         DESTEP = DESTEP + TEL
      ENDIF
C
C             Update probabilities
C
      CALL GRNDM(RNDM,1)
      ZINTCO=-LOG(RNDM(1))
C
      END
