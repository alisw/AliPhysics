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
      FUNCTION WLGAMA(Z)
#include "gen/imp128.inc"
#if (!defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX*16 WLGAMA,WLOGAM
#endif
#if (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX(2) WLGAMA,WLOGAM
#endif
C
      CHARACTER*(*) NAME
      PARAMETER( NAME='WLGAMA')

      ENTRY WLOGAM(Z)

      WLGAMA=0
      CALL MTLPRT(NAME,'C306.0','Not available on this machine')
      RETURN
      END
#endif
