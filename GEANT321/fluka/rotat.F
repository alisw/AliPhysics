*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE ROTAT.FOR
*COPY ROTAT
*
*=== rotat ============================================================*
*
      SUBROUTINE ROTAT(PX,PY,PZ,PXN,PYN,PZN,COTE,SITE,COPS,SIPS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      PXN=-PX*SIPS-PY*COTE*COPS+PZ*SITE*COPS
      PYN=PX*COPS-PY*COTE*SIPS+PZ*SITE*SIPS
      PZN=PY*SITE+PZ*COTE
      RETURN
      END
