*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MXSTEP (X,XLOW,XHI,P,STPPOS,STPNEG)
      DOUBLE PRECISION X, XLOW, XHI, P, STPPOS, STPNEG
      IF(P.EQ.0.0D+0) RETURN
      IF(P.GT.0.0D+0) GOTO 10
      STPPOS=-(X-XLOW)/P
      STPNEG=-(XHI-X)/P
      RETURN
 10   STPPOS=(XHI-X)/P
      STPNEG=(X-XLOW)/P
      RETURN
      END
