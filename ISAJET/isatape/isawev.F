#include "isajet/pilot.h"
      SUBROUTINE ISAWEV
C
C          WRITE OUT MONTECARLO DATA IF EDIT IS TRUE
C
      LOGICAL EDIT
C
C        keep entry point WRTAPE for backward compatibility
      ENTRY WRTAPE
      IF(.NOT.EDIT(I)) RETURN
      CALL WGENS
      RETURN
      END
