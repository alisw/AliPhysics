*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DBZEJY(A,N,MODE,REL,X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DBZEJY')
      CALL MTLPRT(NAME,'C345.0','Not available on this machine')
      RETURN
      END
#endif
