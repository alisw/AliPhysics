*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE AZIRN(SIN,COS)
#include "geant321/crandm.inc"
C       THIS ROUTINE SELECTS THE AZIMUTHAL ANGLE UNIFORMLY IN THETA
   10 R1 = SFLRAF(DUM)
      R1SQ = R1 * R1
      R2 = RANDC(ISEED)
      R2SQ = R2 * R2
      RSQ = R1SQ + R2SQ
      IF(1.0-RSQ) 10 ,20 ,20
   20 SIN = 2.0 * R1 * R2 / RSQ
      COS = (R2SQ-R1SQ) / RSQ
      RETURN
      END
