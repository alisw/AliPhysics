*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSAGTR(X,P,SAFETY,INSIDE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *              SUBROUTINE GSAGTR(X,P,SAFETY,INSIDE)              *
C.    *    Routine to cumpute the 'Safe' distance to nearest boundary  *
C.    *    of the general twisted trapezoid. Point is given in X(1-3), *
C.    *    parameters of trapezoid are in P, if INSIDE = 1 point is    *
C.    *    inside shape if 0 outside. 'SAFE' distance is returned in   *
C.    *    SAFETY.                                                     *
C.    *    I have not yet been able to come up with an exact           *
C.    *    computation of this. From the outside I use an exscribed    *
C.    *    cylinder with its axis as the line joining the centres of   *
C.    *    the trapezia at the ends. As far as I can see there is no   *
C.    *    straight forward way of finding a reliable conservative     *
C.    *    estimate from the inside; so I set SAFETY to 0.0 for all    *
C.    *    points inside the shape. The radius of the exscribed        *
C.    *    cylinder is given by the longest of the eight distances     *
C.    *    from the centre of the trapezium to each corner on each of  *
C.    *    the end faces.                                              *
C.    *        Called by : GSNGTR                                      *
C.    *            A.C.McPherson.      23rd April 1985.                *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),P(30)
C
      SAFETY=0.0
C
C                  Check if point is inside.
C
      IF(INSIDE.EQ.1) GO TO 999
C
C                  Compute radius of cylinder.
C
      RCYL=0.0
      DO 10 I=1,4
         I0=I*4+11
         RC2=(P(I0)+(P(I0+2)-P(13))*P(1))**2+
     +   (P(I0+1)+(P(I0+3)-P(14))*P(1))**2
         IF(RC2.GT.RCYL) RCYL=RC2
         RC2=(P(I0)-(P(I0+2)-P(13))*P(1))**2+
     +   (P(I0+1)-(P(I0+3)-P(14))*P(1))**2
         IF(RC2.GT.RCYL) RCYL=RC2
   10 CONTINUE
      IF(RC2.GT.0.0) RCYL=SQRT(RC2)
C
C                 The direction cosines of the axis of the cylinder
C                 are computed and the distance from the point to the
C                 axis is calculated from the cross product of these
C                 direction cosines with the vector from the origin to
C                 the point. The cross product of this cross product
C                 with the direction cosines gives the vector from the
C                 axis to the point. Subtracting this times the ratio
C                 (D-RCYL)/D where D is the distance of the point from
C                 the axis and RCYL is the radius of the cylinder gives
C                 the point on the cylinder nearest to the point. If
C                 this is outside the z range of the shape then the
C                 distance along the cylinder surface to the z limit is
C                 added in quadrature.
C
      IF(ABS(X(3)).GT.P(1)) SAFETY=ABS(X(3)-P(1))
      TTH2=P(13)**2+P(14)**2
      CTH2=1.0/(1.0+TTH2)
      DIR3=SQRT(CTH2)
      DIR1=P(13)*DIR3
      DIR2=P(14)*DIR3
      DX=DIR2*X(3)-DIR3*X(2)
      DY=DIR3*X(1)-DIR1*X(3)
      DZ=DIR1*X(2)-DIR2*X(1)
      D2=DX*DX+DY*DY+DZ*DZ
      IF(D2.LT.RCYL*RCYL) GO TO 999
C
C            Only Z component of vector is needed.
C
      DDZ=DX*DIR2-DY*DIR1
      D=SQRT(D2)
      Z=X(3)-DDZ*(D-RCYL)/D
      SAFETY=D-RCYL
      IF(ABS(Z).LT.P(1)) GO TO 999
      SAFETY=SQRT(SAFETY*SAFETY+(ABS(Z)-P(1))**2/CTH2)
  999 CONTINUE
      RETURN
      END
