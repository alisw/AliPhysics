*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:34:46  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:29:45  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/systemf.F"
#else
      INTEGER FUNCTION SYSTEMF (NAME)
C
C CERN PROGLIB# Z265    SYSTEMF         .VERSION KERNFOR  4.26  910313
C ORIG. 22/02/91, JZ
C Fortran interface routine to system
C
      CHARACTER    NAME*(*)
      INTEGER      SYSTEI

 
      LGNAME  = LNBLNK (NAME)
      SYSTEMF = SYSTEI (NAME,LGNAME)
      RETURN
      END
#endif
