*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZODAT (LUNP,IXDIVP,LENTP)

C-    Write Direct Access Table

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), IXDIVP(9), LENTP(9), IUHEAD(8)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOD, 4HAT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZODAT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZODAT  ')
#endif
      DATA  IUHEAD / 8*0 /


#include "zebra/qtrace.inc"

      LUNNX   = LUNP(1)
      CALL FZLOC (LUNNX,2)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.2)  WRITE (IQLOG,9002) LUNX
 9002 FORMAT (' FZODAT-  called for LUN=',I4)
#endif

C----              Write the DAT forward reference record

      IF (LENTP(1).NE.0)           GO TO 24
      IPILX(3) = 2
      CALL FZOUT (LUNX, 0,0, 1, 'Z', 1,8,IUHEAD)
      IPILX(3) = 0
      IACTVX   = 14
      IQ(KQSP+LQFX+2) = IACTVX
      GO TO 999

C----              Write the d/a table record

   24 IPILX(3) = 1
      CALL FZOUT (LUNX, IXDIVP,LENTP, 1, 'S', 0,0,0)
      IPILX(3)  = 0
      IQ(KQSP+LQFX+34) = IQUEST(5)
      IQ(KQSP+LQFX+35) = IQUEST(6)

C--                Try to update the DAT forward reference record

#if defined(CERNLIB_FZDACC)
      IF (IFIFOX.NE.2)             GO TO 999
C--       flush the buffer
      CALL FZOUT (LUNX, -7,0, 13, 'FZEND', 0,0,0)
C--       update
      CALL FZUDAT (LUNP,0)
      LUNX = -1
#endif

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
