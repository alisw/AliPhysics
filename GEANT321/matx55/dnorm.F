*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:36  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.50  by  S.Giani
*-- Author :
      SUBROUTINE DNORM(X,AX)
*     ======================
*
*     Double Precision Vector Normalization
*
*     Author: V.Innocente (Naples)
*
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
#endif
      DIMENSION X(3)
*
      AX = SQRT(X(1)**2+X(2)**2+X(3)**2)
      AX1 = 1./AX
      X(1) = AX1*X(1)
      X(2) = AX1*X(2)
      X(3) = AX1*X(3)
*
      END
