*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CFOPEN (LUNDES, MEDIUM, NWREC, MODE, NBUF, NAME, ISTAT)
C
C CERN PROGLIB# Z310    CFOPEN          .VERSION KERNFOR  4.27  910426
C ORIG. 12/01/91, JZ
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
      IF (INDEX(CHUSE,'l').NE.0)  NVMODE(3) = 1

#if !defined(CERNLIB_QMVAX)
      CALL CFOPEI (LUNDES,MEDIUM,NWREC,NVMODE,NBUF,NAME,ISTAT,LGN)
#endif
#if defined(CERNLIB_QMVAX)
      CALL CFOPEI (LUNDES,MEDIUM,NWREC,NVMODE,NBUF,%REF(NAME),ISTAT,LGN)
#endif
      RETURN

   91 LUNDES(1) = 0
      ISTAT(1)  = -1
      RETURN
      END
#ifdef CERNLIB_TCGEN_CFOPEN
#undef CERNLIB_TCGEN_CFOPEN
#endif
