*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE EXITF (ICODE)
C
C CERN PROGLIB# Z035    EXITF           .VERSION KERNFOR  4.39  940228
C ORIG. 15/12/93  JS + JZ
C

      IC = IABS(ICODE)

#if defined(CERNLIB_QMVAX)

      CALL SYS$EXIT(%VAL(IC))

#endif
#if !defined(CERNLIB_QMVAX)

      IF (IC.EQ.0)   STOP

      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),IC
      STOP 255

    1 STOP 1
    2 STOP 2
    3 STOP 3
    4 STOP 4
    5 STOP 5
    6 STOP 6
    7 STOP 7
    8 STOP 8
    9 STOP 9
   10 STOP 10
   11 STOP 11
   12 STOP 12
   13 STOP 13
   14 STOP 14
   15 STOP 15
   16 STOP 16
   17 STOP 17
   18 STOP 18
   19 STOP 19
   20 STOP 20
#endif
      END
