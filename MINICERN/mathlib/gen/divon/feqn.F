*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FEQN (NCUT,F,FGAM,FNLIN)
      INTEGER NCUT
      DOUBLE PRECISION FGAM
      DOUBLE PRECISION F(NCUT), FNLIN(NCUT)
      INTEGER I, NCUTM1
      IF(NCUT.EQ.1) GOTO 20
      NCUTM1=NCUT-1
      DO 10 I=1,NCUTM1
      FNLIN(I)=F(I)-F(I+1)
 10   CONTINUE
 20   FNLIN(NCUT)=F(1)-FGAM
      RETURN
      END
