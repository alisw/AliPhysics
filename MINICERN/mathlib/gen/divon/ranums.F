*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_VAX)
      OPTIONS/CHECK=NOOVERFLOW
#endif
      SUBROUTINE RANUMS (X,N)
      REAL X(N)
      INTEGER IA, IC, ITWO, IY, M2, M
      DOUBLE PRECISION HALFM
      DATA M2 / 0 /, ITWO / 2 /, IY /123456789/
      IF(M2.NE.0) GOTO 20
      M=1
 10   M2=M
      M=ITWO*M2
      IF(M.GT.M2) GOTO 10
      HALFM=M2
      IA=8*INT(HALFM*ATAN(1.0D0)/8.0D0)+5
      IC=2*INT(HALFM*(0.5D0-SQRT(3.0D0)/6.0D0))+1
      S=0.5/HALFM
 20   DO 30 I=1,N
      IY=IY*IA+IC
      IF(IY/2.GT.M2) IY=(IY-M2)-M2
      IF(IY.LT.0) IY=(IY+M2)+M2
 30   X(I)=IY*S
      RETURN
      END
