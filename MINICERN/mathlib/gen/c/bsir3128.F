*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DBSIR3(X,NU)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1, NAME2, NAME3, NAME4
      PARAMETER( NAME1='DBSIR3')
      PARAMETER( NAME2='DEBIR3')
      PARAMETER( NAME3='DBSKR3')
      PARAMETER( NAME4='DEBKR3')
      DBSIR3=0
      CALL MTLPRT(NAME1,'C340.0','Not available on this machine')
      RETURN

      ENTRY DEBIR3(X,NU)
      DEBIR3=0
      CALL MTLPRT(NAME2,'C340.0','Not available on this machine')
      RETURN

      ENTRY DBSKR3(X,NU)
      DBSKR3=0
      CALL MTLPRT(NAME3,'C340.0','Not available on this machine')
      RETURN

      ENTRY DEBKR3(X,NU)
      DEBKR3=0
      CALL MTLPRT(NAME4,'C340.0','Not available on this machine')
      RETURN
      END
#endif
