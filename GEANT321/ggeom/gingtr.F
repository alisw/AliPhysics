*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINGTR(X,P,IYES)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *          SUBROUTINE GINGTR(X,P,IYES)                           *
C.    *      Routine performes GINME type function to determine        *
C.    *      whether the point X is in the GTRA shape with parameters  *
C.    *      P. IYES is returned 1 if it is 0 if not.                  *
C.    *          A.C.McPherson 22nd April 1985                         *
C.    *                                                                *
C.    *      Called by : GINME                                         *
C.    *                                                                *
C     ******************************************************************
C.
      DIMENSION X(3),P(30),XPT(4),YPT(4)
C.
C.            ----------------------------------------
C.
C             First the range in Z is checked.
C
      IYES=0
      IF(ABS(X(3)).GT.P(1)) GO TO 30
C
C             Next the four corners at the z of the input point (X(3))
C             are computed as:
C
      DO 10 I=1,4
      I0=I*4+11
      XPT(I) = P(I0)+P(I0+2)*X(3)
      YPT(I) = P(I0+1)+P(I0+3)*X(3)
   10 CONTINUE
C
C             Then the cross product between the vector from that
C             corner to the input point and the vector from that corner
C             to the next corner is computed and its sign is
C             interrogated. Since the quadralateral is convex, if any
C             one of these signs is wrong, the point is outside the
C             shape. This task is performed as follows:
C
      DO 20 I=1,4
      J = I+1
      IF(J.EQ.5) J=1
      CROSS = (X(1)-XPT(I))*(YPT(J)-YPT(I)) -
     +        (X(2)-YPT(I))*(XPT(J)-XPT(I))
      IF(CROSS.GT.0.0) GO TO 30
   20 CONTINUE
      IYES = 1
   30 CONTINUE
      RETURN
      END
