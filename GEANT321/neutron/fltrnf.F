*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
C*********************************************************************
      FUNCTION FLTRNF(X)
C*********************************************************************
#include "geant321/crandm.inc"
      FLTRNF = RANDC(ISEED)
      RETURN
      END
