*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DBESJ0(X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1, NAME2, NAME3, NAME4
      PARAMETER( NAME1='DBESJ0')
      PARAMETER( NAME2='DBESJ1')
      PARAMETER( NAME3='DBESY0')
      PARAMETER( NAME4='DBESY1')
      DBESJ0=0
      CALL MTLPRT(NAME1,'C312.0','Not available on this machine')
      RETURN

      ENTRY DBESJ1(X)
      DBESJ1=0
      CALL MTLPRT(NAME2,'C312.0','Not available on this machine')
      RETURN

      ENTRY DBESY0(X)
      DBESY0=0
      CALL MTLPRT(NAME3,'C312.0','Not available on this machine')
      RETURN

      ENTRY DBESY1(X)
      DBESY1=0
      CALL MTLPRT(NAME4,'C312.0','Not available on this machine')
      RETURN
      END
#endif
