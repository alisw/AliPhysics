*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      LOGICAL FUNCTION ACCESSF (FNAME, MODE)
C
C CERN PROGLIB# Z265    ACCESSF         .VERSION KERNFOR  4.37  930715
C ORIG. 06/10/92, RDM + JZ
C Fortran interface to access
C
      CHARACTER*(*)   FNAME
      INTEGER         ACCESI

      LGNAME = LNBLNK(FNAME)

#if defined(CERNLIB_QMVAX)
      I = ACCESI (%REF(FNAME), MODE, LGNAME)
#endif
#if !defined(CERNLIB_QMVAX)
      I = ACCESI (FNAME, MODE, LGNAME)
#endif

      IF (I .EQ. 0)  THEN
          ACCESSF = .TRUE.
          RETURN
        ENDIF

      ACCESSF = .FALSE.
      RETURN
      END
#ifdef CERNLIB_TCGEN_ACCESSF
#undef CERNLIB_TCGEN_ACCESSF
#endif
