*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RCHEBN(M,N,A,MDIM,B,TOL,RELERR,X,RESMAX,IRK,ITER,IOCD)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DCHEBN(M,N,A,MDIM,B,TOL,RELERR,X,RESMAX,IRK,ITER,IOCD)
#endif
      CHARACTER*6 NAME
#if defined(CERNLIB_DOUBLE)
      NAME = 'RCHEBN'
#endif
#if !defined(CERNLIB_DOUBLE)
      NAME = 'DCHEBN'
#endif
      CALL MTLPRT(NAME,'E222',
     +'not available on this machine - see documentation')
      RETURN
      END
