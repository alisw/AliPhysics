*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE ALTRAF.FOR
*COPY ALTRAF
*
*=== altraf ===========================================================*
*
      SUBROUTINE ALTRAF(GA,BGA,CX,CY,CZ,COD,COF,SIF,PC,EC,P,PX,PY,PZ,E)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C
C--------------------------------------------------
C*** (S1 IN S2 / PARTICLE IN S1 / RESULT: PARTICLE IN S2)
C*** ARBITRARY LORENTZ TRANSFORM
C*** SI ARE THE DIFFERENT LORENTZ-SYSTEMS
C--------------------------------------------------
      BGX  = BGA * CX
      BGY  = BGA * CY
      BGZ  = BGA * CZ
      COD2 = COD * COD
      IF (COD2 .GT. 1.D0) COD2 = 1.D0
      SID  = SQRT(1.D0 - COD2) * PC
      PCX  = SID * COF
      PCY  = SID * SIF
      PCZ  = COD * PC
      EP = PCX * BGX + PCY * BGY + PCZ * BGZ
      E  = GA * EC + EP
      PE = EP / (GA + 1.D0) + EC
      PX = PCX + BGX * PE
      PY = PCY + BGY * PE
      PZ = PCZ + BGZ * PE
      P  = SQRT (PX * PX + PY * PY + PZ * PZ)
      PM = 1.D0 / P
      PX = PX * PM
      PY = PY * PM
      PZ = PZ * PM
      RETURN
      END
