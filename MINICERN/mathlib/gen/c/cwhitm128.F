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
      FUNCTION WWHITM(Z,KA,MU)
C
#if (!defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX*16  WWHITM,Z,KA,MU
#endif
#if (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX(2)  WWHITM,Z,KA,MU
#endif
#if !defined(CERNLIB_DOUBLE)
C
      CHARACTER*(*) NAME
      PARAMETER( NAME='WWHITM')
      WWHITM=0
      CALL MTLPRT(NAME,'C328.0','Not available on this machine')
      RETURN
      END
#endif
#endif
