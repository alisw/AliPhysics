*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      CHARACTER*(*) FUNCTION CFROMI(I)
C
C CERN PROGLIB# M431    CFROMI          .VERSION KERNFOR  4.18  880425
C ORIG. 10/02/86 M.METCALF
C
C   TO CONVERT 4 PACKED ASCII CODES TO CHARACTERS
C   ANSI FORTRAN
C
C   FOLLOWING STATEMENT IS ANSI EQUIVALENT OF M8 = ZFF
      PARAMETER (M8 = 255)
#if defined(CERNLIB_SGI)
       CHARACTER CHTMP
#endif
C
      CFROMI = ' '
      J = I
      LENCH = MIN(LEN(CFROMI), 4)
      JJ = LENCH
      DO 19  II = 4, 1, -1
         IF (II.LE.LENCH) THEN
#if defined(CERNLIB_SGI)
            CALL ITOCH(IAND(J, M8), CHTMP, *16)
   16       CFROMI(JJ:JJ)=CHTMP
            JJ = JJ-1
#else
            CALL ITOCH(IAND(J, M8), CFROMI(JJ:JJ), *16)
   16       JJ = JJ-1
#endif
         ENDIF
         J = ISHFT(J, -8)
   19 CONTINUE
      END
