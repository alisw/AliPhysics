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
      SUBROUTINE GTISO(U,V,W)
C
#include "geant321/crandm.inc"
C
   10 Z = RANDC(ISEED)
      X = 0.687368 * SFLRAF(Z)
      Y = 0.687368 * SFLRAF(Z)
      XSQ = X * X
      YSQ = Y * Y
      ZSQ = Z * Z
      D = XSQ + YSQ + ZSQ
      IF(D*D-Z) 20 ,20 ,10
   20 U = 2.0*X*Z/D
      V = 2.0*Y*Z/D
      W = (ZSQ-XSQ-YSQ)/D
      RETURN
      END
