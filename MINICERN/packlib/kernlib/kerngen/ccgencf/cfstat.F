*
* $Id$
*
* $Log$
* Revision 1.3  1996/10/17 16:16:30  cernlib
* Correct type of cfstati to be integer
*
* Revision 1.2  1996/10/17 09:23:15  cernlib
* add back in accidentally deleted declaration lines
*
* Revision 1.1  1996/10/16 12:57:34  cernlib
* Add cfstat. cfstati is used by cfstat. This uses rfio_stat if CERNLIB_SHIFT
* is set.
*
* Kernlib
*
*
#include "kerngen/pilot.h"
      INTEGER FUNCTION CFSTAT (NAME, INFO)
C
C CERN PROGLIB# Z310    CFSTAF
C ORIG. 11-Oct-96; stolen from ccgen/statf.F; GF.
C ORIG. 14/03/91, RDM
C Fortran interface routine to stat
C
C        NAME   the name of the file of which information is required
C     INFO(12)  buffer returning the information about "file"
C
C        A succesfull call returns 0. A failed call returns -1.
C
      CHARACTER*(*)  NAME
      INTEGER        INFO(12), CFSTATI

      LGN   = LNBLNK (NAME)
#if !defined(CERNLIB_QMVAX)
      CFSTAT = CFSTATI (NAME, INFO, LGN)
#else
      CFSTAT = CFSTATI (%REF(NAME), INFO, LGN)
#endif
      RETURN
      END
