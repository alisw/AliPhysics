*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DGMLT1(FSUB1,A,B,NI,NG,X)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGMLT')
      DGMLT1=0
      ENTRY DGMLT2(FSUB2,A,B,NI,NG,X)
      DGMLT2=0
      ENTRY DGMLT3(FSUB3,A,B,NI,NG,X)
      DGMLT3=0
      ENTRY DGMLT4(FSUB4,A,B,NI,NG,X)
      DGMLT4=0
      ENTRY DGMLT5(FSUB5,A,B,NI,NG,X)
      DGMLT5=0
      ENTRY DGMLT6(FSUB6,A,B,NI,NG,X)
      DGMLT6=0
      CALL MTLPRT(NAME,'D110.0','Not available on this machine')
      RETURN
      END
#endif
