*
* $Id$
*
* $Log$
* Revision 1.3  1997/03/14 17:19:29  mclareni
* WNT mods
*
* Revision 1.2.2.1  1997/01/21 11:33:38  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.2  1996/04/24 17:26:20  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
      FUNCTION NZLEFT (IXDIV,IOPT)
#include "zebra/mqsys.inc"

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      CALL MZSDIV (IXDIV,1)


      JQDIVR = JQDIVI
      IF (JQDIVR.LT.3)                                     GO TO 41
      JQMODE = JBIT (IQMODE(KQT+JQDIVR), 1)
      JQKIND = JBYT (IQKIND(KQT+JQDIVR),21,4)
      JQSHAR = 0
      IF (JQMODE.NE.0)                                     GO TO 31


      JQDIVN = JQDIVR + 1
      IF (JQDIVR.EQ.JQDVLL)                                GO TO 24
      IF (JQDIVR.EQ.20)                                    GO TO 25

      IF (JBYT(IQKIND(KQT+JQDIVN),21,4).NE.JQKIND)         GO TO 25
      IF (JBYT(IQMODE(KQT+JQDIVN),21,4).EQ.JQMODE)         GO TO 25
      NZLEFT = LQSTA(KQT+JQDIVN) - LQEND(KQT+JQDIVR)
                                                           GO TO 100

   24 JQDIVN = JQDVSY
   25 L      = MIN(LQSTA(KQT+JQDIVR)+NQDMAX(KQT+JQDIVR),
     +               LQSTA(KQT+JQDIVN) )
      NZLEFT = L - LQEND(KQT+JQDIVR)
                                                           GO TO 100


   31 JQDIVN = JQDIVR - 1
      IF (JQDIVR.EQ.JQDVSY)                                GO TO 34
      IF (JBYT(IQKIND(KQT+JQDIVN),21,4).NE.JQKIND)         GO TO 35
      IF (JBYT(IQMODE(KQT+JQDIVN),21,4).EQ.JQMODE)         GO TO 35
      NZLEFT = LQSTA(KQT+JQDIVR) - LQEND(KQT+JQDIVN)
                                                           GO TO 100

   34 JQDIVN = JQDVLL
   35 NZLEFT = NQDMAX(KQT+JQDIVR)-LQEND(KQT+JQDIVR)+LQEND(KQT+JQDIVN)
                                                           GO TO 100


   41 NZLEFT = LQSTA(KQT+2) - LQEND(KQT+1) - NQMINR

  100 IQUEST(11) = LQEND(KQT+JQDIVI) - LQSTA(KQT+JQDIVI)

      END
