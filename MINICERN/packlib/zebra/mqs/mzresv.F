*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:08  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZRESV

C-    Calculate free space for current division, system called
C-    setting JQMODE, JQKIND, JQSHAR, NQRESV

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


      JQDIVR = JQDIVI
      IF (JQDIVR.LT.3)             GO TO 41
      JQMODE = JBIT (IQMODE(KQT+JQDIVR), 1)
      JQKIND = JBYT (IQKIND(KQT+JQDIVR),21,4)
      IF (JQMODE.NE.0)             GO TO 31

C--                FORWARD DIVISION

      JQDIVN = JQDIVR + 1
      IF (JQDIVR.EQ.JQDVLL)        GO TO 24
      IF (JQDIVR.EQ.20)            GO TO 25

      IF (JBYT(IQKIND(KQT+JQDIVN),21,4).NE.JQKIND)  GO TO 25
      IF (JBIT(IQMODE(KQT+JQDIVN),1)   .EQ.JQMODE)  GO TO 25
      JQSHAR = JQDIVN
      JQSHR1 = JQDIVR
      JQSHR2 = JQDIVN
      NQRESV = LQSTA(KQT+JQDIVN) - LQEND(KQT+JQDIVR)
      RETURN

   24 JQDIVN = JQDVSY
   25 L      = MIN (LQSTA(KQT+JQDIVR)+NQDMAX(KQT+JQDIVR),
     +               LQSTA(KQT+JQDIVN) )
      NQRESV = L - LQEND(KQT+JQDIVR)
      JQSHAR = 0
      RETURN

C--                REVERSE DIVISION

   31 JQDIVN = JQDIVR - 1
      IF (JQDIVR.EQ.JQDVSY)        GO TO 34
      IF (JBYT(IQKIND(KQT+JQDIVN),21,4).NE.JQKIND)  GO TO 35
      IF (JBIT(IQMODE(KQT+JQDIVN),1)   .EQ.JQMODE)  GO TO 35
      JQSHAR = JQDIVN
      JQSHR1 = JQDIVN
      JQSHR2 = JQDIVR
      NQRESV = LQSTA(KQT+JQDIVR) - LQEND(KQT+JQDIVN)
      RETURN

   34 JQDIVN = JQDVLL
   35 L      = MAX (LQEND(KQT+JQDIVR)-NQDMAX(KQT+JQDIVR),
     +               LQEND(KQT+JQDIVN) )
      NQRESV = LQSTA(KQT+JQDIVR) - L
      JQSHAR = 0
      RETURN

C----              DIVISIONS 2 + 1

   41 JQKIND = 1
      JQSHR1 = 1
      JQSHR2 = 2
      NQRESV = LQSTA(KQT+2) - LQEND(KQT+1) - NQMINR
      IF (JQDIVR.EQ.1)             GO TO 44
      JQMODE = 1
      JQDIVN = 1
      JQSHAR = 1
      RETURN

   44 JQMODE = 0
      JQDIVN = 2
      JQSHAR = 2
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
