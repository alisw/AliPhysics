*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RLEN (LEN,VEC,VLEN)
      INTEGER   LEN
      DOUBLE PRECISION   VEC(LEN), VLEN
      INTEGER  I
      DOUBLE PRECISION   ABSV, RATIO, RMIN, TOL, TVMX, VMAX
      DATA TOL / 1.0D-20 /
#if defined(CERNLIB_IBM)||defined(CERNLIB_SINGLE)
      DATA RMIN / 1.0D-58 /
#endif
#if (!defined(CERNLIB_IBM))&&(defined(CERNLIB_DOUBLE))
      DATA RMIN / 1.0D-34/
#endif
      VMAX=0.0D+0
      DO 10 I=1,LEN
      ABSV=ABS(VEC(I))
      IF(ABSV.GT.VMAX) VMAX=ABSV
 10   CONTINUE
      VLEN=0.0D+0
      IF(VMAX.EQ.0.0D+0) RETURN
      TVMX=0.0D+0
      IF(VMAX.GT.RMIN) TVMX=TOL*VMAX
      DO 20 I=1,LEN
      IF(ABS(VEC(I)).LE.TVMX) GOTO 20
      RATIO=VEC(I)/VMAX
      VLEN=VLEN+RATIO*RATIO
 20   CONTINUE
      VLEN=VMAX*SQRT(VLEN)
      RETURN
      END
