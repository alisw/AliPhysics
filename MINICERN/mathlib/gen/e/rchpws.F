*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RCHPWS(N,C,A)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DCHPWS(N,C,A)
#endif
      CHARACTER NAMECP*(*),NAMEPC*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMECP = 'RCHPWS', NAMEPC = 'RPWCHS')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMECP = 'DCHPWS', NAMEPC = 'DPWCHS')
#endif

      CALL MTLPRT(NAMECP,'E408.0','NOT AVAILABLE IN THIS MACHINE')
      RETURN

#if defined(CERNLIB_DOUBLE)
      ENTRY RPWCHS(N,A,C)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY DPWCHS(N,A,C)
#endif

      CALL MTLPRT(NAMEPC,'E408.0','NOT AVAILABLE IN THIS MACHINE')
      RETURN
      END
