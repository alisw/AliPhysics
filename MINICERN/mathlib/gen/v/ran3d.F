*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RAN3D (X,Y,Z,XLONG)
C         SUBROUTINE TO GENERATE RANDOM THREE-DIMENSIONAL VECTOR
C         OF LENGTH XLONG.       F. JAMES, MAY, 1976
      DIMENSION RBUF(30)
      SAVE IBUF,RBUF
      DATA IBUF/30/
    1 CONTINUE
      IF (IBUF .LT. 30)  GO TO 10
      CALL NRAN(RBUF,30)
      IBUF = 1
   10 CONTINUE
      A = RBUF(IBUF) - 0.5
      B = RBUF(IBUF+1) - 0.5
      C = RBUF(IBUF+2) - 0.5
      IBUF = IBUF + 3
      R2=A**2+B**2+C**2
      IF(R2.GT.0.25) GO TO 1
      RINV = XLONG/SQRT(R2)
      X=A*RINV
      Y=B*RINV
      Z=C*RINV
      RETURN
      END
