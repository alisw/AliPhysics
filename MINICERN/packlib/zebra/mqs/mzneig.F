*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION MZNEIG (JDIV)

C-    Find the neighbour to division JDIV :
C-    the division MZNEIG shares the reserve area with JDIV

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------


#include "zebra/q_jbit.inc"

      IF (JDIV.GE.3)               GO TO 21
      MZNEIG = 3 - JDIV
      RETURN

   21 MODE = JBIT (IQMODE(KQT+JDIV),1)
      NEIG = 1 - 2*MODE
      JN   = JDIV + NEIG
      IF (JN.GE.21)                GO TO 29
      IF (JN.EQ.JQDVSY-1)          GO TO 26
      IF (JN.NE.JQDVLL+1)          GO TO 27
      JN = JQDVSY
      GO TO 27

   26 JN = JQDVLL
   27 IF (JBIT(IQMODE(KQT+JN),1).EQ.MODE)  GO TO 29
      MZNEIG = JN
      RETURN

C--                No neighbour

   29 MZNEIG = JDIV
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
