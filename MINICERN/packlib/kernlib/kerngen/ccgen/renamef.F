*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/renamef.F"
#else
      INTEGER FUNCTION RENAMEF (FROM,TO)
C
C CERN PROGLIB# Z265    RENAMEF         .VERSION KERNFOR  4.26  910313
C ORIG. 22/02/91, JZ
C Fortran interface routine to rename
C
      CHARACTER*(*)  FROM, TO
      INTEGER        RENAMI

      LGFR = LNBLNK (FROM)
      LGTO = LNBLNK (TO)
      RENAMEF = RENAMI (FROM, TO, LGFR, LGTO)
      RETURN
      END
#endif
