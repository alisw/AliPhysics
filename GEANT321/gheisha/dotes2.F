*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.39  by  S.Giani
*-- Author :
      SUBROUTINE DOTES2(C,S,C2,S2,PR,I)
C
C *** DOUBLE PRECISION VERSION OF "ROTES2" ***
C *** NVE 29-MAR-1988 CERN GENEVA ***
C
C CALLED BY : PHASP
C ORIGIN : H.FESEFELDT (02-DEC-1986)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PR(50)
C
      K1 = 5*I - 4
      K2 = K1 + 1
      SA = PR(K1)
      SB = PR(K2)
      A      = SA*C - SB*S
      PR(K2) = SA*S + SB*C
      K2 = K2 + 1
      B = PR(K2)
      PR(K1) = A*C2 - B*S2
      PR(K2) = A*S2 + B*C2
      RETURN
      END
