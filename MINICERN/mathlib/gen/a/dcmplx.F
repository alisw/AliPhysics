*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (defined(CERNLIB_QF2C))&&(!defined(CERNLIB_LINUX))
      SUBROUTINE DCMPLX(A,X,Y)
      DOUBLE PRECISION X,Y,A(2)
      A(1)=X
      A(2)=Y
      END
#endif
