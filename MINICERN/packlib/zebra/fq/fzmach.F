*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZMACH (JOPP)

C-    System operations  0 : REWIND   1 : UNLOAD,   2 : CLOSE

#include "zebra/zunit.inc"
#include "zebra/fzcf.inc"
C--------------    End CDE                             --------------


      IF (MEDIUF.GE.4)             GO TO 99
      IF (JOPP.GE.2)               GO TO 41

C----              Rewind

      IF (IFIFOF.EQ.2)             GO TO 99
#if defined(CERNLIB_FZLIBC)
      IF (IACMOF.EQ.2)  THEN
          IF (IADOPF.LT.0)         GO TO 99
          CALL CFREW (IADOPF, MEDIUF)
          GO TO 29
        ENDIF
#endif

C--       Non-Fortran REWIND

#include "fzrewex.inc"
* Ignoring t=pass

C--       Fortran REWIND

   28 REWIND LUNF
   29 CONTINUE
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.0)  WRITE (IQLOG,9029) LUNF
 9029 FORMAT (10X,'LUN=',I4,' Rewind')
#endif
      GO TO 99

C----              Close

   41 CONTINUE
#if defined(CERNLIB_FZLIBC)
      IF (IACMOF.EQ.2)  THEN
          IF (IADOPF.LT.0)         GO TO 99
          CALL CFCLOS (IADOPF, MEDIUF)
          IADOPF = -1
          GO TO 49
        ENDIF
#endif

C--       Non-Fortran CLOSE

#include "fzclosex.inc"
* Ignoring t=pass

C--       Fortran CLOSE

   48 CLOSE (LUNF)
   49 CONTINUE
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.0)  WRITE (IQLOG,9049) LUNF
 9049 FORMAT (10X,'LUN=',I4,' Close')
#endif
   99 RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
