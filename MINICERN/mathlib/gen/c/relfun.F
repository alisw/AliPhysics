*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RELFUN(X,AK2,SN,CN,DN)
#include "gen/def64.inc"
     +   DSN,DCN,DDN,D

      ENTRY ELFUN(X,AK2,SN,CN,DN)
      SROUND(D)=D+(D-SNGL(D))
      CALL DELFUN(DBLE(X),DBLE(AK2),DSN,DCN,DDN)
      SN=SROUND(DSN)
      CN=SROUND(DCN)
      DN=SROUND(DDN)
      RETURN
      END
#endif
