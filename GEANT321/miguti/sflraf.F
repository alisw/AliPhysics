*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      FUNCTION SFLRAF(X)
C
#include "geant321/crandm.inc"
C
      SFLRAF = 2.0 * RANDC(ISEED)
      TEMP = 1.0 - SFLRAF
      IF(TEMP) 10 ,20 ,20
   10 SFLRAF = TEMP
   20 RETURN
      END
