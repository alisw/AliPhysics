*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TREVNI(FROM, LENGTH, TO, MAX, NEWLNG)
C
C   *INVERT* PUTS INTO 'TO' ALL NUMBERS 1 ... MAX WHICH ARE NOT
C   CONTAINED IN 'FROM(1)' ... 'FROM(LENGTH)'. 'NEWLNG' SPECIFIES THE
C   NUMBER OF ELEMENTS IN 'TO'.
C   (NEWLNG = MAX - LENGTH)
C
      INTEGER FROM(1), TO(1)
      NEW = 0
      DO 1 I = 1,MAX
      DO 2 J = 1,LENGTH
      IF(FROM(J) .EQ. I) GO TO 1
2     CONTINUE
      NEW = NEW + 1
      TO(NEW) = I
1     CONTINUE
      NEWLNG = NEW
      RETURN
      END
