*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:41  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:05  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE BXQBAR (X,QQ,XQBAR)
C.......................................................................
C.
C.    SEA
C.
C.......................................................................
#include "pdf/impdp.inc"
      DATA XLAM /0.47D00/, Q02/5.D0/
      DATA ONED/1.D0/, FIVED2/5.2D0/
C
      S = LOG( LOG( MAX( FIVED2,QQ)/XLAM**2)/ LOG( Q02/XLAM**2))
      AA=0.99D00+0.72D00*S+0.96D00*S*S
      PP=8.10D00-1.49D00*S+5.10D00*S*S
      XQBAR=AA*(ONED-X)**PP
      RETURN
      END
