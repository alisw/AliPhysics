*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:39  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:04  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE BXF3 (X,QQ,XF3)
C.......................................................................
C.
C.    F3   (GIVES VALENCE FOR U AND D)
C.
C.......................................................................
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +       BBETA
      DATA XLAM /0.47D00/, Q02/5.D0/
      DATA ONED/1.D0/, FIVED2/5.2D0/
C
      S = LOG( LOG( MAX(FIVED2,QQ)/XLAM**2)/ LOG( Q02/XLAM**2))
      ETA1=0.56D00-0.92D00*S*4.00D00/25.00D00
      ETA2=2.71D00+5.08D00*S*4.00D00/25.00D00
      XF3=(X**ETA1)*((ONED-X)**ETA2)/BBETA(ETA1,ONED+ETA2)
      XF3=3.0D00*XF3
      RETURN
      END
