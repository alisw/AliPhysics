*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:17  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:29  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE RXQBAR (X,QQ,XQBAR)
C.......................................................................
C.
C.    SEA
C.
C.......................................................................
#include "pdf/impdp.inc"
      DATA XLAM /0.50D00/, Q0/2.1D0/, Q00/1.8D0/
      DATA ONED/1.D0/, THRED/3.D0/
      IF(QQ.LT.Q0) QQ=Q0
      S=LOG(QQ/XLAM**2)/LOG(Q00/XLAM**2)
      S=LOG(S)
      AS=0.1467D00-0.1211D00*S+0.0274D00*S**2
      ETAS=7.0D00+0.0217D00*S+0.0037D00*S**2
      APS=0.1853D00*S-0.0608D00*S**2
      ETAPS=9.5041D00+1.0165D00*S-0.1049D00*S**2
      BS=0.1682D00*S+0.4473D00*S**2
      CS=25.8997D00+3.9572D00*S+1.6331D00*S**2
      XQBAR=AS*((ONED-X)**ETAS)+APS*((ONED-X)**ETAPS)+BS*EXP(-CS*X)
      RETURN
      END
