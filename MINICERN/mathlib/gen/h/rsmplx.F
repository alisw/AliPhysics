*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RSMPLX(A,B,C,Z0,IDA,M,M1,N,N1,LW,IDW,W,X,Z,ITYPE)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DSMPLX(A,B,C,Z0,IDA,M,M1,N,N1,LW,IDW,W,X,Z,ITYPE)
#endif

      CHARACTER NAME*(*)
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RSMPLX')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DSMPLX')
#endif
      CALL MTLPRT(NAME,'H101.0','NOT AVAILABLE IN THIS MACHINE')
      RETURN
      END
