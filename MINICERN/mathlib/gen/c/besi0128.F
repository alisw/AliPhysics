*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DBESI0(X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1, NAME2, NAME3, NAME4,
     +              NAME5, NAME6, NAME7, NAME8
      PARAMETER( NAME1='DBESI0')
      PARAMETER( NAME2='DEBSI0')
      PARAMETER( NAME3='DBESI1')
      PARAMETER( NAME4='DEBSI1')
      PARAMETER( NAME5='DBESK0')
      PARAMETER( NAME6='DEBSK0')
      PARAMETER( NAME7='DBESK1')
      PARAMETER( NAME8='DEBSK1')
      CALL MTLPRT(NAME1,'C313.0','Not available on this machine')
      RETURN

      ENTRY DEBSI0(X)
      DEBSI0=0
      CALL MTLPRT(NAME2,'C313.0','Not available on this machine')
      RETURN

      ENTRY DBESI1(X)
      DBESI1=0
      CALL MTLPRT(NAME3,'C313.0','Not available on this machine')
      RETURN

      ENTRY DEBSI1(X)
      DEBSI1=0
      CALL MTLPRT(NAME4,'C313.0','Not available on this machine')
      RETURN

      ENTRY DBESK0(X)
      DBESK0=0
      CALL MTLPRT(NAME5,'C313.0','Not available on this machine')
      RETURN

      ENTRY DEBSK0(X)
      DEBSK0=0
      CALL MTLPRT(NAME6,'C313.0','Not available on this machine')
      RETURN

      ENTRY DBESK1(X)
      DBESK1=0
      CALL MTLPRT(NAME7,'C313.0','Not available on this machine')
      RETURN

      ENTRY DEBSK1(X)
      DEBSK1=0
      CALL MTLPRT(NAME8,'C313.0','Not available on this machine')
      RETURN
      END
#endif
