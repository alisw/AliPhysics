*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 19/05/94  14.28.55  by  S.Ravndal
*-- Author :
      SUBROUTINE SVGAM7(A,X)
      DIMENSION U(2)
      REAL LA
      GA=A-1.
      LA=SQRT(2.*A-1.)
      EP=1.570796327+ATAN(GA/LA)
      RO=1.570796327-EP
    1 CALL GRNDM(U,2)
      X=GA+LA*TAN(EP*U(1)+RO)
      IF (X.EQ.0) GOTO 1
      Y=LOG(1.+((X-GA)/LA)**2) +GA*LOG(X/GA)-X+GA
      IF(LOG(U(2)).GT.Y) GOTO 1
      RETURN
      END
