*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION MAXDZE(DA,N)
#if !defined(CERNLIB_DOUBLE)
      CHARACTER*6 NAME
      NAME = 'MAXDZE'
      CALL MTLPRT(NAME,'E102',
     +'not available on this machine - see documentation')
      MAXDZE = 0
      RETURN

      FUNCTION MAXDZE(DA,N)
      NAME = 'MINDZE'
      CALL MTLPRT(NAME,'E102',
     +'not available on this machine - see documentation')
      MINDZE = 0
      RETURN
      END

#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION DA,D
 
      DIMENSION DA(*)
 
      J=1
      D=DA(1)
      DO 5 K = 2, N
      IF(DA(K) .LE. D) GO TO 5
      J=K
      D=DA(K)
    5 CONTINUE
      MAXDZE=J
      RETURN
 
      ENTRY MINDZE(DA,N)
 
      J=1
      D=DA(1)
      DO 6 K = 2,N
      IF(DA(K) .GE. D) GO TO 6
      J=K
      D=DA(K)
    6 CONTINUE
      MINDZE=J
      RETURN
      END
#endif
