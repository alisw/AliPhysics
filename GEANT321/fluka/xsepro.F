*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:05  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE XSEPRO.FOR
*COPY XSEPRO
*
*=== xsepro ===========================================================*
*
      FUNCTION XSEPRO ( E, ZTAR, ATAR )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 20 september 1991 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-sep-91     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/xsepar.inc"
       COMMON / FKNEGX / XSENEG
*
       IZ = NINT ( ZTAR )
       IF ( E .GT. 0.1D+00 ) THEN
          EMEV = E * 1.D+03
          XSEPRO = SITSAO ( EMEV, IZ, ATAR )
       ELSE
          X  = 1.D+01 * E
          X2 = X * X
          XSEPRO = AAPXSE (IZ) * ( X - ZZPXSE (IZ) ) / X * ( 1.D+00
     &           + BBPXSE (IZ) * X + CCPXSE (IZ) * X2 + DDPXSE (IZ)
     &           * X * X2 + EEPXSE (IZ) * X2 * X2 + FFPXSE (IZ) * X2
     &           * X * X2 )
          IF ( XSEPRO .LT. 0.D+00 ) THEN
             XSENEG = XSEPRO
             XSEPRO = 0.D+00
          END IF
       END IF
*=== End of function XSEPRO ===========================================*
       RETURN
       END
