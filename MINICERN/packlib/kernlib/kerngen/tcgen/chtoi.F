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
      SUBROUTINE CHTOI (CH,I,*)
C
C CERN PROGLIB# M400    CHTOI           .VERSION KERNFOR  4.23  891215
C ORIG. 29/05/84 H.RENSHALL
C MOD.  06/11/89 JZ for UNIX & for speed
C
C  CHTOI - convert CH*1 to its ASCII value
C          error exit if 'unseen' character
C
      CHARACTER CH*1
C
      I = ICHAR(CH)
      IF (I.GE.32)  THEN
          IF (I.LE.126)          RETURN
        ENDIF
C
      I = 0
      RETURN 1
      END
#endif
