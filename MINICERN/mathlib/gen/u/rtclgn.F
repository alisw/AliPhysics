*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RTCLGN(J1,J2,J3,M1,M2,M3,DNUM,DDEN,KPEX)
      CHARACTER*6 NAME
      NAME  ='RTCLGN'
      CALL MTLPRT(NAME,'U112',
     +'not available on this machine - see documentation')
      RETURN
      END
#endif
