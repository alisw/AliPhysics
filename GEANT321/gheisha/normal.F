*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE NORMAL(RAN)
C
C *** NVE 14-APR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (27-OCT-1983)
C
      DIMENSION RNDM(12)
      RAN=-6.
      CALL GRNDM(RNDM,12)
      DO 1 I=1,12
      RAN=RAN+RNDM(I)
 1    CONTINUE
      RETURN
      END
