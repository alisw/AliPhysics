*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIPRL (IHOW)

C-    Physical record length mis-match, reset if possible

C-    in :    IHOW = 0  reset if possible
C-                   1  reset only if start-of-run

C-    out :  IRETN = 0  reset done
C-                   1  reset fails, not start-of-run
C-                   2  reset fails, record larger than buffer

#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
#include "fzhci.inc"
C--------------    End CDE                             --------------


      IF (NWRI.GT.IQ(KQSP+LBPARI))    GO TO 42
      IRETN = 0
      IF (IHOW.NE.0)  THEN
          IF (JMODI.NE.1)          GO TO 41
          IF (NTLRI.NE.8)          GO TO 41
          L4CURI = L4STAI + 8
          IF (LQ(L4CURI+1).NE.1)   GO TO 41
          IF (LQ(L4CURI+2).LE.0)   GO TO 41
        ENDIF

#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.-1)  WRITE (IQLOG,9066) LUNI, NWRI
 9066 FORMAT (1X/' FZIPHR.  LUN=',I3,', Reset block size to',I6)
#endif

      MAXREI = NWRI
      NWMIN  = NWRI
#if defined(CERNLIB_FQNEEDPK)
      IF (IUPAKI.NE.0)             GO TO 38
#include "fznwmach.inc"
* Ignoring t=pass
#endif
   38 IQ(KQSP+LQFI+5)   = MAXREI
      IQ(KQSP+LBPARI+1) = NWMIN
      GO TO 48

   42 IRETN = 1
   41 IRETN = IRETN + 1
   48 IQUEST(1) = IRETN
      RETURN
      END
