*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION WDIGAM(Z)
C
#if (!defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX*16 WDIGAM,Z
#endif
#if (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX(2) WDIGAM,Z
#endif
#if !defined(CERNLIB_DOUBLE)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'WDIGAM')
      WDIGAM=0
      CALL MTLPRT(NAME,'C307.0','Not available on this machine')
      RETURN
      END
#endif
#endif
