*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DRKNYS(N,H,X,Y,YP,SUB,W)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='DRKNYS')
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RRKNYS(N,H,X,Y,YP,SUB,W)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RRKNYS')
#endif
      DIMENSION Y(*),YP(*),W(N,*)

      PARAMETER (Z1 = 1, R2 = Z1/2, R6 = Z1/6, R8 = Z1/8)

#if !defined(CERNLIB_DOUBLE)
      ENTRY RKNYS(N,H,X,Y,YP,SUB,W)
#endif

      IF(N .LE. 0 .OR. H .EQ. 0) RETURN
      H2=R2*H
      H6=R6*H
      HH2=H*H2
      HH6=H*H6
      HH8=R8*H*H
      XH=X+H
      XH2=X+H2
      CALL SUB(X,Y,YP,W(1,1))
      DO 1 J = 1,N
      W(J,4)=Y(J)+H2*YP(J)
      W(J,5)=W(J,4)+HH8*W(J,1)
    1 W(J,6)=YP(J)+H2*W(J,1)
      CALL SUB(XH2,W(1,5),W(1,6),W(1,2))
      DO 2 J = 1,N
      W(J,6)=YP(J)+H2*W(J,2)
      W(J,1)=W(J,1)+W(J,2)
    2 W(J,2)=W(J,1)+W(J,2)
      CALL SUB(XH2,W(1,5),W(1,6),W(1,3))
      DO 3 J = 1,N
      W(J,4)=W(J,4)+H2*YP(J)
      W(J,5)=W(J,4)+HH2*W(J,3)
      W(J,6)=YP(J)+H*W(J,3)
      W(J,1)=W(J,1)+W(J,3)
    3 W(J,2)=W(J,2)+2*W(J,3)
      CALL SUB(XH,W(1,5),W(1,6),W(1,3))
      DO 4 J = 1,N
      Y(J)=W(J,4)+HH6*W(J,1)
    4 YP(J)=YP(J)+H6*(W(J,2)+W(J,3))
      X=XH
      RETURN
      END
