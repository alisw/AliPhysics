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
      INTEGER FUNCTION NFAC(N)
C
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (27-OCT-1983)
C
      NFAC=1.
      M=N
      IF(M.LE.1) RETURN
      IF(M.GT.10) M=10
      DO 1 I=2,M
    1 NFAC=NFAC*I
      RETURN
      END
