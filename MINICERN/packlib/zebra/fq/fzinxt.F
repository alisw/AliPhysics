*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZDACC)
      SUBROUTINE FZINXT (LUNP,MDSA1,MDSA2)

C-    Set next d/s to be read in direct access mode, user called

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), MDSADR(2)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIN, 4HXT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZINXT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZINXT  ')
#endif

#include "zebra/qtrace.inc"

      LUNNI     = LUNP(1)
      MDSADR(1) = MDSA1
      MDSADR(2) = MDSA2
      IF (LUNNI.NE.LUNI)  CALL FZLOC (LUNNI,1)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.2)  WRITE (IQLOG,9111) LUNNI,MDSADR
 9111 FORMAT (' FZINXT-  called for LUN=',I4,' d/s adr=',2I7)
#endif
      IF (IFIFOI.NE.2)             GO TO 801

      IACTVI = 2
      IQ(KQSP+LQFI+2)  = IACTVI
      IQ(KQSP+LQFI+27) = 0
      IQ(KQSP+LQFI+31) = MDSADR(1)
      IQ(KQSP+LQFI+32) = MDSADR(2)
      IQ(KQSP+LQFI+33) = -7

      LBPARI = LQFI + INCBPI
      CALL VZERO (IQ(KQSP+LBPARI-9), 9)
#include "zebra/qtrace99.inc"
      RETURN

  801 CONTINUE
#if defined(CERNLIB_QDEBPRI)
      WRITE (IQLOG,9111) LUNNI,MDSADR
#endif
      CALL ZFATAM ('FZINXT - file not ready for direct access.')
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
