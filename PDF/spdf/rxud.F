*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:18  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:29  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE RXUD (X,QQ,XUD)
C.......................................................................
C.
C.    F3   (GIVES VALENCE FOR U AND D)
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
      ETA1=0.7D00-0.176D00*S
      ETA2=2.6D00+0.8D00*S
      XUD=THRED*(X**ETA1)*((ONED-X)**ETA2)/BBETA(ETA1,ONED+ETA2)
      RETURN
      END
