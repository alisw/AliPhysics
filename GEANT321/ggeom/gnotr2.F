*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNOTR2(X,Y,Z,N1,N2,N3,N4,ABCD)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     GNOTR2 computes the coefficients of the implicit           *
C.    *     normalized plane equation.                                 *
C.    *     These are called ABCD(1),ABCD(2),ABCD(3) and ABCD(4) in    *
C.    *     GNOTR2. They are stored in the parameter array             *
C.    *     respectively in P(12),P(13),P(14),P(15) for the first      *
C.    *     surface, P(16),P(17),P(18),P(19) for the second surface    *
C.    *     and so on. The face is formed by the half-edges            *
C.    *     N1-N2,N2-N3,N3-N4,N4-N1, where N1,N2,N3,N4 are the         *
C.    *     indices of the corresponding vertices.                     *
C.    *         Called by : GNOTR1                                     *
C.    *         Author  R.Nierhaus  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(8),Y(8),Z(8),ABCD(4)
*
      CALL GNOTR9(X,Y,Z,N1,N2,N3,N4)
*
      A=((Y(N2)-Y(N1))*(Z(N1)+Z(N2))+(Y(N3)-Y(N2))*(Z(N2)+Z(N3))
     +  +(Y(N4)-Y(N3))*(Z(N3)+Z(N4))+(Y(N1)-Y(N4))*(Z(N4)+Z(N1)))
      B=((Z(N2)-Z(N1))*(X(N1)+X(N2))+(Z(N3)-Z(N2))*(X(N2)+X(N3))
     +  +(Z(N4)-Z(N3))*(X(N3)+X(N4))+(Z(N1)-Z(N4))*(X(N4)+X(N1)))
      C=((X(N2)-X(N1))*(Y(N1)+Y(N2))+(X(N3)-X(N2))*(Y(N2)+Y(N3))
     +  +(X(N4)-X(N3))*(Y(N3)+Y(N4))+(X(N1)-X(N4))*(Y(N4)+Y(N1)))
      S=SQRT(A**2+B**2+C**2)
      ABCD(1)=A/S
      ABCD(2)=B/S
      ABCD(3)=C/S
      ABCD(4)=-(ABCD(1)*X(N1)+ABCD(2)*Y(N1)+ABCD(3)*Z(N1))
*
      END
