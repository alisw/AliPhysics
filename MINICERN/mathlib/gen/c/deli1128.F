*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:10  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DELI1(X,AKP)
C
#include "gen/imp128.inc"
      DELI1=0
      CALL MTLPRT('DELI1','C346.0','Not available on this machine')
      RETURN

      ENTRY DELI2(X,AKP,A,B,MODE)
      DELI2=0
      CALL MTLPRT('DELI2','C346.0','Not available on this machine')
      RETURN

      ENTRY DELI3(X,AKP,P)
      DELI3=0
      CALL MTLPRT('DELI3','C346.0','Not available on this machine')
      RETURN
      END
#endif
