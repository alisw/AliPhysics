*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE CELFUN(Z,AK2,SN,CN,DN)
      COMPLEX Z,SN,CN,DN
#include "gen/defc64.inc"
     +    WSN,WCN,WDN
#include "gen/def64.inc"
     +       D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
#include "gen/gcmpfun.inc"
      SROUND(D)=D+(D-SNGL(D))
      CALL WELFUN
     1     (DCMPLX(DBLE(REAL(Z)),DBLE(AIMAG(Z))),DBLE(AK2),WSN,WCN,WDN)
      SN=CMPLX(SROUND(DREAL(WSN)),SROUND(DIMAG(WSN)))
      CN=CMPLX(SROUND(DREAL(WCN)),SROUND(DIMAG(WCN)))
      DN=CMPLX(SROUND(DREAL(WDN)),SROUND(DIMAG(WDN)))
      RETURN
      END
#endif
