*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ERPRIP(GA1,GA2,GA12,GB)
      DIMENSION GA1(9),GA2(9),GB(9),GA12(9)
      COMMON /SLATE/XQX(9), V(4), G1(9), G2(9), G(9)
      CALL MXMPY (G1, GA1,G,3,3,3)
      CALL MXMPY1 (G,G1,GB,3,3,3)
      CALL MXMPY (G2,GA2,G,3,3,3)
      CALL MXMAD1 (G,G2,GB,3,3,3)
      IF(GA12(1).EQ.0.) GOTO 1
      CALL MXMPY(G1,GA12,G,3,3,3)
      CALL MXMAD1(G,G2,GB,3,3,3)
      CALL MXMPY1(G2,GA12,G,3,3,3)
      CALL MXMAD1(G,G1,GB,3,3,3)
    1 G(3)=-1./V(4)**3
      G(1)=V(1)*G(3)
      G(2)=V(2)*G(3)
      G(3)=V(3)*G(3)
      G(8)=V(1)**2+V(2)**2
      G(6)=1./(V(4)**2*SQRT(G(8)))
      G(5)=-G(6)*V(3)
      G(4)=V(1)*G(5)
      G(5)=V(2)*G(5)
      G(6)=G(8)*G(6)
      G(8)=1./G(8)
      G(7)=-V(2)*G(8)
      G(8)=V(1)*G(8)
      G(9)=0.
      CALL MXMPY (G,GB,G1,3,3,3)
      CALL MXMPY1(G1,G,GB,3,3,3)
      RETURN
      END
