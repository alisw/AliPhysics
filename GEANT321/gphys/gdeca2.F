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
      SUBROUTINE GDECA2(XM0,XM1,XM2,PCM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates two body decay process with isotropic angular       *
C.    *  distribution in CMS.                                          *
C.    *                                                                *
C.    *    ==>Called by : GDECAY                                       *
C.    *       Author    G.Patrick  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION PCM(4,2)
      DIMENSION RNDM(2)
C.
C.    ------------------------------------------------------------------
C.
C             Generate first decay product in CMS.
C
      E1=(XM0*XM0+XM1*XM1-XM2*XM2)/(2.*XM0)
      P1=SQRT(ABS((E1-XM1)*(E1+XM1)))
C
C             Isotropic decay angular distribution.
C
      CALL GRNDM(RNDM,2)
      COSTH=2.*RNDM(1)-1.
      IF(ABS(COSTH).GE.1.0) THEN
         COSTH=SIGN(1.,COSTH)
         SINTH=0.
      ELSE
         SINTH=SQRT((1.-COSTH)*(1.+COSTH))
      ENDIF
      PHI=TWOPI*RNDM(2)
C
C             Polar co-ordinates to momentum components.
C
      PCM(1,1)=P1*SINTH*COS(PHI)
      PCM(2,1)=P1*SINTH*SIN(PHI)
      PCM(3,1)=P1*COSTH
      PCM(4,1)=E1
C
C             Generate second decay product.
C
      PCM(1,2)=-PCM(1,1)
      PCM(2,2)=-PCM(2,1)
      PCM(3,2)=-PCM(3,1)
      PCM(4,2)=SQRT(PCM(1,2)**2+PCM(2,2)**2+PCM(3,2)**2+XM2*XM2)
C
      END
