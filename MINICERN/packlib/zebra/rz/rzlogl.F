*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLOGL(LUN,LOGLEV)
*
************************************************************************
*
*        To set the LOG level
* Input:
*   LUN     Logical unit number for which the logging level has to be set
*   LOGLEV  Logging level
*           -3  Suppress all messages
*           -2  Error messages only
*           -1  Terse logging
*            0  Normal logging: RZFILE, RZEND
*            1  Log to watch rare events
*            2  Log to monitor calls
*            3  Short diagnostics
*               Short dumps to debug user-written output routines
*            4  Full diagnostics
*               Full dumps to debug user-written output routines
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 17.11.92 JDS - restrict loglevel to documented range
*                     so that only 3 bits are required to store it (15-17)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      LRZ=LQRS
  10  IF(LRZ.EQ.0)GO TO 99
      IF(IQ(KQSP+LRZ-5).NE.LUN)THEN
         LRZ=LQ(KQSP+LRZ)
         GO TO 10
      ENDIF
*
      LOGL = MIN(LOGLEV,4)
      LOGL = MAX(LOGL, -3) + 3
      CALL SBYT(LOGL,IQ(KQSP+LRZ),15,3)
      IF(LOGLEV.GE.0)WRITE(IQLOG,9022)LUN,LOGLEV
 9022 FORMAT(' RZLOGL.  File at LUN=',I6,', Diagnostic log level=',I2)
*
  99  RETURN
      END
