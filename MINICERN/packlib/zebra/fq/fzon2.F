*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZFFNAT)
      SUBROUTINE FZON2 (V1,N1,V2,N2)

#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
#if !defined(CERNLIB_AUTODBL)
      INTEGER   V1(N1),V2(N2)
#endif
#if defined(CERNLIB_AUTODBL)
      REAL   V1(N1),V2(N2),RELIDX(2)
      EQUIVALENCE (RELIDX(1),IDX(1))
#endif

      IDX(1) = N1 + N2

#if !defined(CERNLIB_AUTODBL)
      WRITE (LUNX) IDX,V1,V2
#endif
#if defined(CERNLIB_AUTODBL)
      WRITE (LUNX) RELIDX,V1,V2
#endif
      IQ(KQSP+LQFX+20) = IQ(KQSP+LQFX+20) + IDX(1) + 3
      IQ(KQSP+LQFX+21) = IQ(KQSP+LQFX+21) + 1
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
