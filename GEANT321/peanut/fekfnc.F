*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE FEKFNC.FOR
*COPY FEKFNC
*
*=== fekfnc ===========================================================*
*
      FUNCTION FEKFNC (PF,I)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
#include "geant321/nucdat.inc"
*
      IF ( PFRCEN (I) - PF .LE. ANGLGB ) THEN
         FEKFNC = EKFCEN (I)
      ELSE
         FEKFNC = SQRT ( AMNUSQ (I) + PF**2 ) - AMNUCL (I)
      END IF
      RETURN
      END
