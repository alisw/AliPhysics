*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:08  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GCUBS(X,Y,D1,D2,A)
C.
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Calculates a cubic through P1,(-X,Y1),(X,Y2),P2          *
C.    *       where Y2=-Y1                                             *
C.    *        Y=A(1)+A(2)*X+A(3)*X**2+A(4)*X**3                       *
C.    *        The coordinate system is assumed to be the cms system   *
C.    *        of P1,P2.                                               *
C.    *                                                                *
C.    *    ==>Called by : GIPLAN,GICYL                                 *
C.    *       Author    H.Boerner  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
      REAL X,Y,D1,D2,A(4)
C.
C.
C.    ------------------------------------------------------------------
C.
C.
      IF (X.EQ.0.)                               GO TO 10
C
C
      FACT   = (D1 - D2) * 0.25
      A(1)   = - 1. * FACT * X
      A(3)   = FACT / X
      A(2)   = (6. * Y - (D1 + D2) * X) / (4. * X)
      A(4)   = ((D1 + D2)*X - 2.*Y) / (4.*X**3)
      RETURN
C
  10  A(1)   = 0.
      A(2)   = 1.
      A(3)   = 0.
      A(4)   = 0.
      END
