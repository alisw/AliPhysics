*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      INTEGER FUNCTION READLNF (PATH, TEXT)
C
C CERN PROGLIB# Z265    READLNF         .VERSION KERNFOR  4.36  930602
C ORIG. 24/03/93, JZ
C Fortran interface routine to readlink
C
      CHARACTER*(*)  PATH, TEXT
      INTEGER        READLNI

      LGN = LNBLNK (PATH)
      NTX = LEN(TEXT)
      READLNF = READLNI (PATH,TEXT,NTX,LGN)
      RETURN
      END
