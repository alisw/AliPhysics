*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE SVERL2(M,X)
      DIMENSION U(2)
      CALL GRNDM(U,2)
      X=0.62666*LOG((1.+U(1))/(1.-U(1)))
      IF(U(2).LT.0.5) X=-X
      X=M+X*SQRT(M*1.)
      RETURN
      END
