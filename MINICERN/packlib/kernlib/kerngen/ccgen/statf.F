*
* $Id$
*
* $Log$
* Revision 1.3  1997/09/02 14:26:39  mclareni
* WINNT correction
*
* Revision 1.2  1997/02/04 17:34:44  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:29:44  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      INTEGER FUNCTION STATF (NAME, INFO)
C
C CERN PROGLIB# Z265    STATF           .VERSION KERNFOR  4.37  930715
C ORIG. 14/03/91, RDM
C Fortran interface routine to stat
C
C        NAME   the name of the file of which information is required
C     INFO(12)  buffer returning the information about "file"
C
C        A succesfull call returns 0. A failed call returns -1.
C
      CHARACTER*(*)  NAME
      INTEGER        INFO(12), STATI

      LGN   = LNBLNK (NAME)

#if defined(CERNLIB_QMVAX)
      STATF = STATI (%REF(NAME), INFO, LGN)
#endif
#if !defined(CERNLIB_QMVAX)
      STATF = STATI (NAME, INFO, LGN)
#endif

      RETURN
      END
