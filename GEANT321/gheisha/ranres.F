*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      FUNCTION RANRES(X)
C
C *** RESTRICTED RANDOM NUMBERS TO BE USED AS ARGUMENT IN LOG ETC... ***
C *** NVE 13-JUL-1988 CERN GENEVA ***
C
C NOTE : 0 < RANRES < 1
      DIMENSION RNDM(1)
C
 1    CONTINUE
      CALL GRNDM(RNDM,1)
      RANRES=RNDM(1)
      IF ((RANRES .LE. 0.) .OR. (RANRES .GE. 1.)) GO TO 1
      RETURN
      END
