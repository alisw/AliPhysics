*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:22  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/getenvf.F"
#else
      SUBROUTINE GETENVF (NAME, TEXT)
C
C CERN PROGLIB# Z265    GETENVF         .VERSION KERNFOR  4.26  910313
C ORIG. 22/02/91, JZ
C Fortran interface routine to getenv
C
      COMMON /SLATE/ISL(40)
      CHARACTER*(*)  NAME, TEXT

      LGN = LNBLNK (NAME)
      NTX = LEN(TEXT)
      CALL GETENI (NAME,TEXT,NTX,LGN)
      ISL(1) = NTX
      RETURN
      END
#endif
