*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:46:57  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:10  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:48  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
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
