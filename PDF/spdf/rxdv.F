*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:15  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:29  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE RXDV (X,QQ,XDV)
C.......................................................................
C.
C.    DV
C.
C.......................................................................
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +       BBETA
      DATA XLAM /0.50D00/, Q0/2.1D0/, Q00/1.8D0/
      DATA ONED/1.D0/, THRED/3.D0/
      IF(QQ.LT.Q0) QQ=Q0
      S=LOG(QQ/XLAM**2)/LOG(Q00/XLAM**2)
      S=LOG(S)
      ETA3=0.85D00-0.24D00*S
      ETA4=3.35D00+0.816D00*S
      XDV=(X**ETA3)*((ONED-X)**ETA4)/BBETA(ETA3,ONED+ETA4)
      RETURN
      END
