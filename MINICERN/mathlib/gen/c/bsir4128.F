*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DBSIR4(X,NU)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1, NAME2, NAME3, NAME4
      PARAMETER( NAME1='DBSIR4')
      PARAMETER( NAME2='DEBIR4')
      PARAMETER( NAME3='DBSKR4')
      PARAMETER( NAME4='DEBKR4')
      DBSIR4=0
      CALL MTLPRT(NAME1,'C327.0','Not available on this machine')
      RETURN

      ENTRY DEBIR4(X,NU)
      DEBIR4=0
      CALL MTLPRT(NAME2,'C327.0','Not available on this machine')
      RETURN

      ENTRY DBSKR4(X,NU)
      DBSKR4=0
      CALL MTLPRT(NAME3,'C327.0','Not available on this machine')
      RETURN

      ENTRY DEBKR4(X,NU)
      DEBKR4=0
      CALL MTLPRT(NAME4,'C327.0','Not available on this machine')
      RETURN
      END
#endif
