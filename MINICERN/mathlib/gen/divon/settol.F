*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SETTOL (N,FTOL,GTOL,DELTA,ETA,ALFMAX,MAXFUN,IPRINT)
      INTEGER N, MAXFUN, IPRINT
      DOUBLE PRECISION FTOL, GTOL, DELTA, ETA, ALFMAX
      FTOL=5.0D-2
      GTOL=1.0D-2
      DELTA=1.0D-5
      ETA=0.2D+0
      ALFMAX=1.0D+2
      MAXFUN=50*N
      IPRINT=-1
      RETURN
      END
