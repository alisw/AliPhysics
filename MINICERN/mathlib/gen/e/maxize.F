*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION MAXIZE(IA,N)
 
      DIMENSION IA(*),RA(*)
 
      J=1
      I=IA(1)
      DO 1 K = 2,N
      IF(IA(K) .LE. I) GO TO 1
      J=K
      I=IA(K)
    1 CONTINUE
      MAXIZE=J
      RETURN
 
      ENTRY MINIZE(IA,N)
 
      J=1
      I=IA(1)
      DO 2 K = 2,N
      IF(IA(K) .GE. I) GO TO 2
      J=K
      I=IA(K)
    2 CONTINUE
      MINIZE=J
      RETURN
 
      ENTRY MAXRZE(RA,N)
      ENTRY MAXFZE(RA,N)
 
      J=1
      F=RA(1)
      DO 3 K = 2, N
      IF(RA(K) .LE. F) GO TO 3
      J=K
      F=RA(K)
    3 CONTINUE
      MAXFZE=J
      RETURN
 
      ENTRY MINRZE(RA,N)
      ENTRY MINFZE(RA,N)
 
      J=1
      F=RA(1)
      DO 4 K = 2,N
      IF(RA(K) .GE. F) GO TO 4
      J=K
      F=RA(K)
    4 CONTINUE
      MINFZE=J
      RETURN
      END
