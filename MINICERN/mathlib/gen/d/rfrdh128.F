*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DFRDH1(F,G,M,TV,NGTV,WS,IDIM,N)
C
#include "gen/imp128.inc"
C
      CALL MTLPRT('DFRDH1','D601.0','Not available on this machine')
      RETURN
C
      ENTRY DFRDH2(F,G,X,WS,IDIM,N)
      CALL MTLPRT('DFRDH2','D601.0','Not available on this machine')
      RETURN
C
      ENTRY DFRDH3(H,WS,IDIM,N)
      CALL MTLPRT('DFRDH3','D601.0','Not available on this machine')
      RETURN
C
      END
#endif
