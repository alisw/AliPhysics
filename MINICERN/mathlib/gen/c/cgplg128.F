*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION WGPLG(N,M,X)
C
#include "gen/imp128.inc"
#if (!defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX*16 WGPLG
#endif
#if (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      COMPLEX(2) WGPLG
#endif
#if !defined(CERNLIB_DOUBLE)
#include "gen/def128.inc"
     +     X
C
      CHARACTER*(*) NAME
      PARAMETER( NAME='WGPLG')
      WGPLG=0
      CALL MTLPRT(NAME,'C321.0','Not available on this machine')
      RETURN
      END
#endif
#endif
