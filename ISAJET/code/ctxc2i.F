#include "isajet/pilot.h"
      SUBROUTINE CTXC2I(CVAL,IVAL,NSIZE)
C-----------------------------------------------------------------------
C          Convert character variable CVAL to integer array IVAL
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      CHARACTER*(*) CVAL
      INTEGER I,NSIZE
      INTEGER IVAL(NSIZE)
C
      DO 100 I=1,NSIZE
100   IVAL(I)=ICHAR(CVAL(I:I))
C
      RETURN
      END
