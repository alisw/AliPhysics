*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
*
      SUBROUTINE GLIGHT(XXLL,YYLL,ZZLL,IILL)
*
#include "geant321/gcrayt.inc"
*
      IOMBRA=10
      XLPOS=XXLL
      YLPOS=YYLL
      ZLPOS=ZZLL
      INTEN=IILL
      IF(INTEN.EQ.0)IOMBRA=0
*
      END
