*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:16  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:29  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE RXGLUE (X,QQ,XGL)
C.......................................................................
C.
C.    GLUON
C.
C.......................................................................
#include "pdf/impdp.inc"
      DATA XLAM /0.50D00/, Q0/2.D0/, Q00/1.8D0/
      DATA ONED/1.D0/, THRED/3.D0/
      IF(QQ.LT.Q0) QQ=Q0
      S=LOG(QQ/XLAM**2)/LOG(Q00/XLAM**2)
      S=LOG(S)
      AG=2.412D00-1.9845D00*S+0.4443D00*S**2
      ETAG=5.0D00+1.5464D00*S-0.5287D00*S**2
      APG=3.6363D00*S-1.4146D00*S**2
      ETAPG=13.8237D00+0.7914D00*S-0.2873D00*S**2
      BG=7.6609D00*S-1.4595D00*S**2
      CG=36.7928D00+12.5884D00*S-1.1536D00*S**2
      XGL=AG*((ONED-X)**ETAG)+APG*((ONED-X)**ETAPG)+BG*EXP(-CG*X)
      RETURN
      END
