*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDFIX(N,NS)
C
C      SETS NS = 0 IF N=0
C                1 IF N IS FIXED POINT NUMBER
C                -1 IF N IS FLOATING POINT NUMBER
C      N.B. THIS VERSION ASSUMES FIXED POINT NUMBER IS KNOWN TO BE LESS
C      THAN 256
C
      IF(N)3,4,3
    4 NS=0
      RETURN
    3 IF(ABS(N)-255)5,5,6
    5 NS=1
      RETURN
    6 NS=-1
      RETURN
      END
