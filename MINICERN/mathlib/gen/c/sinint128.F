*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:06  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DSININ(X)
#include "gen/imp128.inc"
C
      CHARACTER*(*) NAME1, NAME2
      PARAMETER( NAME1='DSININ')
      PARAMETER( NAME2='DCOSIN')
      DSININ=0
      CALL MTLPRT(NAME1,'C336.0','Not available on this machine')
      RETURN

      ENTRY DCOSIN(X)
      DCOSIN=0
      CALL MTLPRT(NAME2,'C336.0','Not available on this machine')
      END
#endif
