*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE NEWRO(M1,N,IBCX,DELTAX,POT,POTW,POTE,POTS,POTN)
      DIMENSION POT(M1,N)
      DIMENSION POTW(N),POTE(N),POTS(M1),POTN(M1)
      M=M1-1
      COEFX=2.*DELTAX
      GO TO (100,200,300,400),IBCX
 100  RETURN
 200  DO 1 J=1,N
      POT(1,J)=POT(1,J)+COEFX*POTW(J)
 1    POT(M1,J)=POT(M1,J)-COEFX*POTE(J)
      RETURN
 300  DO 2 J=1,N
      POT(1,J)=POTW(J)
      POT(M1,J)=POTE(J)
      POT(2,J)=POT(2,J)-POTW(J)
 2    POT(M,J)=POT(M,J)-POTE(J)
      RETURN
 400  DO 3 J=1,N
      POT(1,J)=POTW(J)
      POT(2,J)=POT(2,J)-POTW(J)
 3    POT(M1,J)=POT(M1,J)-COEFX*POTE(J)
      RETURN
      END
