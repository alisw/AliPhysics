*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZINFO (LUNP)

C-    Locate FZ control-bank for unit LUN,
C-    copy control bank information to /FZSTAT/,
C-    signal failure by setting INFLUN = 0

#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzstat.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9)


      LUN  = LUNP(1)

      LQFF = LQFS
   14 IF (LQFF.NE.0)  THEN
          IF (IQ(KQSP+LQFF-5).EQ.LUN)  GO TO 21
          LQFF = LQ(KQSP+LQFF)
          GO TO 14
        ELSE
          INFLUN = 0
          RETURN
        ENDIF

C--                Bank for LUN found

   21 INFLUN = LUN
      INFSTA = IQ(KQSP+LQFF)
      CALL UCOPY (IQ(KQSP+LQFF+1), INFOFZ, 40)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
