*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QASCII)
      SUBROUTINE ITOCH (I,CH,*)
C
C CERN PROGLIB# M400    CHTOI           .VERSION KERNFOR  4.23  891215
C ORIG. 29/05/84 H.RENSHALL
C MOD.  06/11/89 JZ for UNIX & for speed
C
C  ITOCH - convert ASCII integer to character CH*1
C          error exit if 'unseen', setting '?'
C
      CHARACTER CH
C
      IF (I.GE.32)  THEN
          IF (I.LE.126)  THEN
              CH = CHAR (I)
              RETURN
            ENDIF
        ENDIF
C
      CH = '?'
      RETURN 1
      END
#endif
