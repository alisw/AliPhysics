*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:04  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DFCONC(X,TAU,M)
#include "gen/imp128.inc"
#if (!defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX*16 WGAMMA,WLOGAM
#endif
#if (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX(2) WGAMMA,WLOGAM
#endif
#if !defined(CERNLIB_DOUBLE)
C
      CHARACTER*(*) NAME
      PARAMETER( NAME='DFCONC')
      DFCONC=0
      CALL MTLPRT(NAME,'C331.0','Not available on this machine')
      RETURN
      END
#endif
#endif
