*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:40  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:05  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE BXGLUE (X,QQ,XGL)
C.......................................................................
C.
C.    GLUON
C.
C.......................................................................
#include "pdf/impdp.inc"
      DATA XLAM /0.47D00/, Q02/5.D0/
      DATA ONED/1.D0/, FIVED2/5.2D0/
C
      S = LOG( LOG( MAX( FIVED2,QQ)/XLAM**2)/ LOG( Q02/XLAM**2))
C     CGL=2.01D00-3.56D00*S+5.10D00*S**2
C ... modified according to Z. Phys. C2 (79) 265, HPB - 25.2.91
      CGL=2.01D00-3.56D00*S+1.98D00*S*S
      ETA1=-1.13D00*S+0.48D00*S*S
      ETA2=2.90D00+5.08D00*S*4.00D00/25.0D00
      XGL=CGL*(X**ETA1)*((ONED-X)**ETA2)
      RETURN
      END
