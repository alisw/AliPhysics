*
* $Id$
*
* $Log$
* Revision 1.2  1996/02/27 10:08:11  ravndal
* Precision problem in cos(theta) solved
*
* Revision 1.1.1.1  1995/10/24 10:21:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GANNI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates positron annihilation                               *
C.    *                                                                *
C.    *    ==>Called by : GTELEC                                       *
C.    *       Author    L.Urban *********                              *
C.    *       10/06/93: modified by Georges Azuelos (Vancouver)        *
C     *                 to include 1-quantum annihilation              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmate.inc"
      DIMENSION PGAM(3),RNDM(2)
      LOGICAL ROTATE
      PARAMETER (ALFA=7.29735E-3)
C.
C.    ------------------------------------------------------------------
C.
      KCASE = NAMEC(11)
      IF((IANNI.NE.1).OR.((GETOT+EMASS).LE.CUTGAM)) THEN
         ISTOP = 2
         DESTEP = DESTEP + GETOT + EMASS
         GEKIN = 0.
         GETOT = EMASS
         VECT(7)= 0.
         GO TO 999
      ENDIF
C
      XE=GETOT
      GAM=XE/EMASS
      GAM2=GAM**2
      GAM1=MAX(GAM2-1.,0.)
      GAMP1=GAM+1.
      C=SQRT(GAM1)
*
      SIG=(GAM2+4.*GAM+1.)*LOG(GAM+C)/GAM1-(GAM+3.)/C
      SIG=0.5*Q(JPROB+17)*SIG/GAMP1
*
      BIND=0.5*(Z*ALFA)**2*EMASS
      E1Q=XE+EMASS-BIND
      IF(E1Q .GT. 0.)THEN
         GVE=VECT(7)/EMASS
         SIG1=GAM2+2.*(GAM+2.)/3.-(GAM+2.)/GVE*LOG(GAM+GVE)
         SIG1=2.*Q(JPROB+18)*SIG1/(GVE*GAMP1**2)
      ELSE
         SIG1=0.
      END IF
*
      SIG=SIG+SIG1
      CALL GRNDM(RNDM,1)
C
      IF(RNDM(1).GE.SIG1/SIG)THEN
         GAMP12=GAMP1**2
         P=VECT(7)
         E0=1./(GAMP1+C)
C
   10    CALL GRNDM(RNDM,2)
         E=E0*((1.-E0)/E0)**RNDM(1)
C
         SCREJ=(GAMP12+2.*GAMP1-2.-GAMP12*E-1./E)/(GAMP12-2.)
         IF(RNDM(2).GT.SCREJ) GOTO 10
C
         EPHOT1=(XE+EMASS)*E
C
         COSTH=(GEKIN+EMASS*(2.*E-1.)/E)/P
C
C restrict COSTH to [-1.,+1.]
C
         COSTH = MIN( 1. , MAX( -1. , COSTH ) )
         SINTH=SQRT((1.-COSTH)*(1.+COSTH))
         CALL GRNDM(RNDM,1)
         PHI = TWOPI * RNDM(1)
         COSPHI = COS(PHI)
         SINPHI = SIN(PHI)
C
         PGAM(1) = EPHOT1* SINTH * COSPHI
         PGAM(2) = EPHOT1* SINTH * SINPHI
         PGAM(3) = EPHOT1* COSTH
C
C             Rotate tracks into GEANT system and store.
C
         CALL GFANG(VECT(4),COSAL,SINAL,COSBT,SINBT,ROTATE)
C
C            Polar co-ordinates to momentum components.
C
         NGGAMM = 0
         IF(EPHOT1.GT.CUTGAM) THEN
            NGGAMM = NGGAMM + 1
            NGKINE = NGKINE + 1
            GKIN(1,NGKINE) = PGAM(1)
            GKIN(2,NGKINE) = PGAM(2)
            GKIN(3,NGKINE) = PGAM(3)
            GKIN(4,NGKINE) = EPHOT1
            GKIN(5,NGKINE) = 1
            TOFD(NGKINE)=0.
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
            IF(ROTATE)
     +      CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
         ELSE
            DESTEP = DESTEP + EPHOT1
         ENDIF
C
C             Momentum vector of second photon.
C
         EPHOT2 = GETOT + EMASS - EPHOT1
         IF(EPHOT2.GT.CUTGAM) THEN
            NGGAMM = NGGAMM + 1
            NGKINE = NGKINE + 1
            GKIN(1,NGKINE) = - PGAM(1)
            GKIN(2,NGKINE) = - PGAM(2)
            GKIN(3,NGKINE) = P - PGAM(3)
            GKIN(4,NGKINE) = EPHOT2
            GKIN(5,NGKINE) = 1
            TOFD(NGKINE)=0.
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
            IF(ROTATE)
     +      CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
         ELSE
            DESTEP = DESTEP + EPHOT2
         ENDIF
      ELSE
C 1-quantum annihilation
         P=VECT(7)
         EPHOT=E1Q
C Assume photon collinear with positron
         PGAM(1) = 0.
         PGAM(2) = 0.
         PGAM(3) = EPHOT
C
C             Rotate tracks into GEANT system and store.
C
         CALL GFANG(VECT(4),COSAL,SINAL,COSBT,SINBT,ROTATE)
C
C            Polar co-ordinates to momentum components.
C
         NGGAMM = 0
         IF(EPHOT.GT.CUTGAM) THEN
            NGGAMM = NGGAMM + 1
            NGKINE = NGKINE + 1
            GKIN(1,NGKINE) = PGAM(1)
            GKIN(2,NGKINE) = PGAM(2)
            GKIN(3,NGKINE) = PGAM(3)
            GKIN(4,NGKINE) = EPHOT
            GKIN(5,NGKINE) = 1
            TOFD(NGKINE)=0.
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
            IF(ROTATE)
     +      CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
         ELSE
            DESTEP = DESTEP + EPHOT
         ENDIF
      END IF
C
      IF(NGGAMM.GT.0) THEN
         ISTOP = 1
      ELSE
         ISTOP = 2
      ENDIF
C
 999  END
