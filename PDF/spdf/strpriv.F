*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:49  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:05  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      SUBROUTINE STRPRIV(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C.......................................................................
C.
C :::::::::::: PRIVATE set of structure functions :::::::::::::::::
C.
C.......................................................................
#include "pdf/impdp.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50517.inc"
#include "pdf/w50516.inc"
      DATA ZEROD/0.D0/
      SAVE /W50516/
C
      IF (N6.LE.0) N6 = L6
      IF (FIRST) WRITE(N6,1000)
 1000 FORMAT(1H ,' This is a DUMMY private SF routine ')
      UPV  = ZEROD
      DNV  = ZEROD
      USEA = ZEROD
      DSEA = ZEROD
      STR  = ZEROD
      CHM  = ZEROD
      BOT  = ZEROD
      TOP  = ZEROD
      GL   = ZEROD
      RETURN
      END
