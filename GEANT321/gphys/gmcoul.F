*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:26  cernlib
* Geant
*
*
*    NEW
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMCOUL(OMEGA,DIN)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Generate NSCA single scatters in small angle approxima.  *
C.    *       for a particle with parameters VECT in common /GCTRAK/   *
C.    *                                                                *
C.    *       This subroutine must be called with the correct values   *
C.    *       of the constants OMC & CHC which depend of the medium    *
C.    *                                                                *
C.    *       OMC and CHC are computed at initialisation time (GMOLI)  *
C.    *       No lateral displacement of the particle with respect     *
C.    *       the incident direction is included.                      *
C.    *       No path length correction is included                    *
C.    *                                                                *
C.    *       Output angles overwrite VECT                             *
C.    *                                                                *
C.    *    ==>Called by : GMULTS                                       *
C.    *         Author G. Lynch (LBL)  (adapted by M. Maire  1.02.90)  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmulo.inc"
      PARAMETER (NSCMX = 50)
      DIMENSION DIN(3),RNDM(2*NSCMX)
      PARAMETER (OMCF=1.167)
*
*     ------------------------------------------------------------------
*
* *** Compute number of scatters (Poisson distr. with mean OMEGA0)
*
      OMEGA0 = OMCF*OMEGA
      CALL GPOISS (OMEGA0,NSCA,1)
      IF (NSCA.LE.0)THEN
         DIN(1) = 0.
         DIN(2) = 0.
         DIN(3) = 1.
         RETURN
      ENDIF
      NSCA = MIN(NSCA,NSCMX)
      CALL GRNDM (RNDM,2*NSCA)
*
* *** THMIN2 is the screening angle
      THMIN2 = (CHCMOL**2)/(OMCF*OMCMOL*(VECT(7)**2))
*
      SUMX = 0.
      SUMY = 0.
      DO 12 I=1,NSCA
         THET  = SQRT(THMIN2*((1./RNDM(I)) - 1.))
         PHI   = TWOPI*RNDM(NSCA+I)
         SUMX  = SUMX + THET*COS(PHI)
         SUMY  = SUMY + THET*SIN(PHI)
   12 CONTINUE
      THETA  = SQRT(SUMX**2 + SUMY**2)
      IF(THETA.NE.0.) THEN
      COSTH  = COS(THETA)
      SINTH  = SIN(THETA)
      CSPHI  = SUMX/THETA
      SNPHI  = SUMY/THETA
*
      DIN(1) = SINTH*CSPHI
      DIN(2) = SINTH*SNPHI
      DIN(3) = COSTH
      ELSE
         DIN(1) = 0.
         DIN(2) = 0.
         DIN(3) = 1.
      ENDIF
*
      END
