*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:38  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CIOPEN (LUNDES, MODE, NAME, ISTAT)
C
C CERN PROGLIB# Z311    CIOPEN          .VERSION KERNFOR  4.31  911111
C ORIG. 12/10/91, JZ
C
      DIMENSION    LUNDES(9), ISTAT(9)
      CHARACTER    MODE*(*), NAME*(*)
      DIMENSION    NVMODE(4)
      CHARACTER    CHUSE*4


      LGN = LNBLNK (NAME)

      CHUSE = MODE
      CALL CUTOL (CHUSE)
      CALL VZERO (NVMODE,4)

      IF (INDEX(CHUSE,'a').NE.0)  NVMODE(1) = 2
      IF (INDEX(CHUSE,'w').NE.0)  THEN
          IF (NVMODE(1).NE.0)      GO TO 91
          NVMODE(1) = 1
        ENDIF
      IF (INDEX(CHUSE,'r').NE.0)  THEN
          IF (NVMODE(1).NE.0)      GO TO 91
        ENDIF

      IF (INDEX(CHUSE,'+').NE.0)  NVMODE(2) = 1

#if !defined(CERNLIB_QMVAX)
      CALL CIOPEI (LUNDES,NVMODE,NAME,ISTAT,LGN)
#endif
#if defined(CERNLIB_QMVAX)
      CALL CIOPEI (LUNDES,NVMODE,%REF(NAME),ISTAT,LGN)
#endif
      RETURN

   91 LUNDES(1) = 0
      ISTAT(1)  = -1
      RETURN
      END
