*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DBSKA(X,IA,JA,NL,B)
#include "gen/imp128.inc"
      CHARACTER*6 NAME
      NAME = 'DBSKA'
      CALL MTLPRT(NAME,'C341',
     +'not available on this machine - see documentation')
      RETURN

      ENTRY DEBSKA(X,IA,JA,NL,B)
      NAME = 'DEBSKA'
      CALL MTLPRT(NAME,'C341',
     +'not available on this machine - see documentation')
      RETURN
      END
#endif
