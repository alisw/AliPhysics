*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      LOGICAL FUNCTION RZSAME(IH1,IH2,N)
************************************************************************
*
*        To compare two hollerith strings
* Input:
*   IH1,IH2 two arrays of length N to be compared
*
* Called by RZDELT,RZEND,RZFDIR,RZMDIR
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 01.04.86
*
************************************************************************
      DIMENSION IH1(N),IH2(N)
*
*-----------------------------------------------------------------------
*
      IF(N.LE.0)GO TO 20
*
      DO 10 I=1,N
         IF(IH1(I).NE.IH2(I))GO TO 20
  10  CONTINUE
      RZSAME=.TRUE.
      GO TO 99
*
  20  RZSAME=.FALSE.
*
  99  RETURN
      END
