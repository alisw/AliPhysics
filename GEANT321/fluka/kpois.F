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
*$ CREATE KPOIS.FOR
*COPY KPOIS
*
*=== kpois ============================================================*
*
      FUNCTION KPOIS(EXPLAM)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      REAL RNDM(1)
C*** EXPLAM IS EXP(-LAMBDA)
C*** SQMPLING FROM POISSON DISTRIBUTION, AVERAGE LAMBDA
C*** PK = EXP(-LAMBDA)*LAMBDA**K/(K-FACTORIAL)
      Y=1.D0
      KPOISS=0
5     CONTINUE
      CALL GRNDM(RNDM,1)
      Y=Y*RNDM(1)
      IF (Y.LT.EXPLAM) GO TO 10
      KPOISS=KPOISS+1
      GO TO 5
10    CONTINUE
      KPOIS=KPOISS
      RETURN
      END
