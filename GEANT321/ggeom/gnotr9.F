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
      SUBROUTINE GNOTR9(X,Y,Z,N1,N2,N3,N4)
C.    ******************************************************************
C.    *                                                                *
C.    *     GNOTR9 checks coplanarity of the face described by the     *
C.    *     indices N1,N2,N3,N4 into the vertex array X,Y,Z.           *
C.    *     This check is perfomed since shape TRAP is described by    *
C.    *     11 values in the parameter array of which only 9 are       *
C.    *     independant.
C.    *         Called by : GNOTR2
C.    *         Author  R.Nierhaus  *********                          *
C.    *                                                                *
C.    ******************************************************************
*
#include "geant321/gcunit.inc"
*
      DIMENSION X(8),Y(8),Z(8)
*
      A=-((Y(N2)-Y(N1))*(Z(N3)-Z(N1))-(Y(N3)-Y(N1))*(Z(N2)-Z(N1)))
      B=-((Z(N2)-Z(N1))*(X(N3)-X(N1))-(Z(N3)-Z(N1))*(X(N2)-X(N1)))
      C=-((X(N2)-X(N1))*(Y(N3)-Y(N1))-(X(N3)-X(N1))*(Y(N2)-Y(N1)))
      V=(X(N4)-X(N1))*A+(Y(N4)-Y(N1))*B+(Z(N4)-Z(N1))*C
      IF (ABS(V)/SQRT(A**2+B**2+C**2).GT.1.E-3) THEN
         WRITE(CHMAIL,1001)
         CALL GMAIL(0,0)
 1001    FORMAT (' GNOTR9 ...TRAP shape :  surface not coplanar')
      END IF
*
      END
 
