*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.41  by  S.Giani
*-- Author :
*$ CREATE BETARN.FOR
*COPY BETARN
*
*=== betarn ===========================================================*
*
      FUNCTION BETARN(GAM,ETA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     New version:                                                     *
*     Created on 20 february 1991  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-feb-91     by    Alfredo Ferrari               *
*                                                                      *
*     Sampling from beta distribution in [0,1) :                       *
*                                                                      *
*     P(X) = X**(GAM-1.D0)*(1.D0-X)**(ETA-1)*GAMM(ETA+GAM)             *
*          / (GAMM(GAM*GAMM(ETA))                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      REAL RNDM(2)
*
      GAMI  = 1.D+00 / GAM
      ETAM1 = ETA - 1.D+00
      NTAM1 = NINT (ETA - 1.D+00)
*  +-------------------------------------------------------------------*
*  |
      IF ( ETAM1 - NTAM1 .NE. 0.D+00 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |  First sample from X**(gam-1) and then reject according to
*  |  |  (1-X)**(eta-1)
  100    CONTINUE
            CALL GRNDM(RNDM,2)
            BETARN = RNDM (1)**GAMI
            REJE   = ( 1.D+00 - BETARN )**ETAM1
         IF ( RNDM (2) .GE. REJE ) GO TO 100
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
*  |  +----------------------------------------------------------------*
*  |  |  First sample from X**(gam-1) and then reject according to
*  |  |  (1-X)**(eta-1)
  200    CONTINUE
            CALL GRNDM(RNDM,2)
            BETARN = RNDM (1)**GAMI
            REJE   = ( 1.D+00 - BETARN )**NTAM1
         IF ( RNDM (2) .GE. REJE ) GO TO 200
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*=== End of function betarn ===========================================*
      END
