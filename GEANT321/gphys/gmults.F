*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMULTS
C
C     ******************************************************************
C     *                                                                *
C     *       Steering routine for the multiple scattering.            *
C     *       select Moliere theory , Gaussian approximation           *
C     *       or single Coulomb scattering depending of their range    *
C     *       of validity.                                             *
C     *                                                                *
C     *    ==>Called by : GTELEC , GTHADR , GTMUON                     *
C.    *       Author     M.Maire  *********                            *
C     *                                                                *
C     ******************************************************************
*
      PARAMETER (THRMOL=50.,THRGAU=0.01)
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
      DIMENSION DIN(3)
*
      BETA2 = (VECT(7)/GETOT)**2
      IF(BETA2.LE.0.) RETURN
      IF(IMULS.EQ.3) THEN
         CALL GMGAUS(BETA2,DIN)
*
      ELSE
         CHARG2=CHARGE**2
*
*     Here we decide whether we have to recalculate OMCMOL or not.
*     OMCMOL has been calculated in GMULOF according to a formula
*     which contains the term (1 + 3.34*(Alpha*z*Z/Beta)**2) where
*     z is the incident charge (CHARGE), setting Beta=1 and z=1.
*     We do not recalculate OMCMOL if:
*
*              3.34*(Z*Alpha)**2*((z/Beta)**2-1) << 1
*                   Z**2*(z**2-Beta**2)/Beta**2  << 3.34/Alpha**2 = 5500
*                   Z**2*(z**2-Beta**2)  <= 50.*Beta**2
*
*     We further multiply the first term for the number of elements
*     in the mixture because the approssimation is particularly bad
*     for mixtures so we make the condition more restrictive.
*     All this thanks to Gerry Lynch.
*
         IF(Z**2*(CHARG2-BETA2)*Q(JMA+11).GT.THRMOL*BETA2) THEN
            NLM=Q(JMA+11)
*
            IF(NLM.EQ.1)THEN
               CALL GMOLIO(A,Z,1.,1,DENS,BETA2,CHARG2,OMCMOL)
*
            ELSE
               CALL GMOLIO(Q(JMIXT+1),Q(JMIXT+NLM+1),Q(JMIXT+2*NLM+1),
     +         NLM,DENS,BETA2,CHARG2,OMCMOL)
*
            ENDIF
         ELSE
            JPROB = LQ(JMA-4)
            OMCMOL = Q(JPROB+21)*CHARG2
*
         ENDIF
*
         OMEGA = OMCMOL*STMULS/BETA2
*
         IF (OMEGA.LE.20.) THEN
            CALL GMCOUL(OMEGA,DIN)
*
         ELSE
            CALL GMOLIE(OMEGA,BETA2,DIN)
*
         ENDIF
      ENDIF
*
* *** Computes rotation matrix around particle direction
* *** Compute new direction cosines
*
      VMM = SQRT(VECT(4)*VECT(4)+VECT(5)*VECT(5))
      IF (VMM.NE.0.) THEN
         PD1=VECT(4)/VMM
         PD2=VECT(5)/VMM
         V4= PD1*VECT(6)*DIN(1) -PD2*DIN(2) +VECT(4)*DIN(3)
         V5= PD2*VECT(6)*DIN(1) +PD1*DIN(2) +VECT(5)*DIN(3)
         V6= -VMM*DIN(1) +VECT(6)*DIN(3)
      ELSE
         V4= DIN(1)
         V5= DIN(2)
         V6= DIN(3)*SIGN(1.,VECT(6))
      ENDIF
*
* *** Renormalize direction cosines
*
      VP = 1./SQRT(V4*V4+V5*V5+V6*V6)
      VECT(4) = V4*VP
      VECT(5) = V5*VP
      VECT(6) = V6*VP
*
      END
