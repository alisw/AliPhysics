*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:24  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMOS9)
#include "os9gs/lstatf.F"
#else
      INTEGER FUNCTION LSTATF (NAME, INFO)
C
C CERN PROGLIB# Z265    LSTATF          .VERSION KERNFOR  4.36  930602
C ORIG. 24/03/93, RDM + JZ
C Fortran interface routine to lstat
C
C        NAME   the name of the file of which information is required
C     INFO(12)  buffer returning the information about "file"
C
C     A succesful call returns 0; a failed call returns -1.
C     also:  NDSLAT = 0  if NAME is a regular file
C            NESLAT = 0  if NAME is a symbolic link
C            NFSLAT = 0  if NAME is a directory
C
      CHARACTER*(*)  NAME
      INTEGER        INFO(12), LSTATI
      COMMON /SLATE/ISLATE(40)

      LGN    = LNBLNK (NAME)
      LSTATF = LSTATI (NAME, INFO, LGN, ISLATE)
      RETURN
      END
#endif
