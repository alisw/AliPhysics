*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNOTR1(P)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routine to extend  the  user supplied parameters P  to      *
C.    *    internal parameters   for the  trapezoidal shape TRAP       *
C.    *                                                                *
C.    *         Called by : GSPOSP, GSVOLU , GGDPAR                    *
C.    *         Author  R.Nierhaus  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
*     GNOTR1  computes the coordinates for the 8 vertices
*     X(8), Y(8) and Z(8).
*     Then GNOTR2 is called for each of the six boundary surfaces.
*
      DIMENSION P(35),X(8),Y(8),Z(8)
*
      HZ  =P(1)
      SX  =P(2)
      SY  =P(3)
      HY1 =P(4)
      HX11=P(5)
      HX12=P(6)
      SX1 =P(7)
      HY2 =P(8)
      HX21=P(9)
      HX22=P(10)
      SX2 =P(11)
      IF (HY1.LT.1.E-6) HY1=1.E-6
      IF (HX11.LT.1.E-6) HX11=1.E-6
      IF (HX12.LT.1.E-6) HX12=1.E-6
      IF (HY2.LT.1.E-6) HY2=1.E-6
      IF (HX21.LT.1.E-6) HX21=1.E-6
      IF (HX22.LT.1.E-6) HX22=1.E-6
      X(1)=-HZ*SX-HY1*SX1-HX11
      X(2)=-HZ*SX-HY1*SX1+HX11
      X(3)=-HZ*SX+HY1*SX1-HX12
      X(4)=-HZ*SX+HY1*SX1+HX12
      X(5)=+HZ*SX-HY2*SX2-HX21
      X(6)=+HZ*SX-HY2*SX2+HX21
      X(7)=+HZ*SX+HY2*SX2-HX22
      X(8)=+HZ*SX+HY2*SX2+HX22
      Y(1)=-HZ*SY-HY1
      Y(2)=-HZ*SY-HY1
      Y(3)=-HZ*SY+HY1
      Y(4)=-HZ*SY+HY1
      Y(5)=+HZ*SY-HY2
      Y(6)=+HZ*SY-HY2
      Y(7)=+HZ*SY+HY2
      Y(8)=+HZ*SY+HY2
      Z(1)=-HZ
      Z(2)=-HZ
      Z(3)=-HZ
      Z(4)=-HZ
      Z(5)=+HZ
      Z(6)=+HZ
      Z(7)=+HZ
      Z(8)=+HZ
      CALL GNOTR2(X,Y,Z,1,2,4,3,P(12))
      CALL GNOTR2(X,Y,Z,5,7,8,6,P(16))
      CALL GNOTR2(X,Y,Z,1,5,6,2,P(20))
      CALL GNOTR2(X,Y,Z,3,4,8,7,P(24))
      CALL GNOTR2(X,Y,Z,1,3,7,5,P(28))
      CALL GNOTR2(X,Y,Z,2,6,8,4,P(32))
*
      END
