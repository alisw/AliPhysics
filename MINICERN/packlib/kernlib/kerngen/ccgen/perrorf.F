*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE PERRORF (TEXT)
C
C CERN PROGLIB# Z265    PERRORF         .VERSION KERNFOR  4.26  910313
C ORIG. 22/02/91, JZ
C Fortran interface routine to perror
C
      CHARACTER    TEXT*(*)

#if !defined(CERNLIB_QMVAX)
      CALL PERROI (TEXT,LEN(TEXT))
#endif
#if defined(CERNLIB_QMVAX)
      CALL PERROI (%REF(TEXT),LEN(TEXT))
#endif
      RETURN
      END
