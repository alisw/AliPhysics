*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RGMLT1(FSUB1,A,B,NI,NG,X)

      ENTRY RGMLT2(FSUB2,A,B,NI,NG,X)
      ENTRY RGMLT3(FSUB3,A,B,NI,NG,X)
      ENTRY RGMLT4(FSUB4,A,B,NI,NG,X)
      ENTRY RGMLT5(FSUB5,A,B,NI,NG,X)
      ENTRY RGMLT6(FSUB6,A,B,NI,NG,X)

C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'RGMLT')
      CALL MTLPRT(NAME,'D110.0','Not available on this machine')
      RGMLT1=0
      RETURN
      END
#endif
