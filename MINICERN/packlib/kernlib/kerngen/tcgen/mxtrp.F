*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE MXTRP (A,B,I,J)
C
C CERN PROGLIB# F110    MXTRP           .VERSION KERNFOR  1.0   650809
C ORIG. 01/01/64 RKB
C
      DIMENSION A(*),B(*)
C
C
      IF((I.EQ.0).OR.(J.EQ.0))GO TO 100
      IB=1
      DO 1 K=1,J
      IA=K
      DO 1 L=1,I
      B(IB)=A(IA)
      IA=IA+J
    1 IB=IB+1
  100 RETURN
      END
