*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUHUNT (IT,IA,N,NW)
C
C CERN PROGLIB# V304    IUHUNT          .VERSION KERNFOR  2.09  751101
C ORIG. 01/09/72
C
      DIMENSION IA(9)
C
      IF (N.EQ.0)  GO TO 19
      IF (NW.EQ.0) GO TO 20
      DO 18 J=1,N,NW
      IF (IA(J).EQ.IT)  GO TO 21
   18 CONTINUE
   19 IUHUNT= 0
      RETURN
   20 IF (IT.NE.IA(1)) GO TO 19
      IUHUNT= 1
      RETURN
C
   21 IUHUNT= J
      RETURN
      END
