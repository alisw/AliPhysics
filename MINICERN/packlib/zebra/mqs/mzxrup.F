*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:14  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZXRUP

C-    Derive IQRCU from IQRTO and IQRNO for current store

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_jbytet.inc"

C----              Loop all divisions

      JDVFR = 1
   21 IF (JDVFR.EQ.JQDVLL+1)  JDVFR = JQDVSY
      MASKTO = IQRTO(KQT+JDVFR)
      MASKNO = IQRNO(KQT+JDVFR)
      MASKCU = JBYT (MASKTO,1,20)

      IF (JBIT(MASKTO,24).EQ.1)  MASKCU = MSBIT1(MASKCU,JQDVSY)

      IF (JBYT(MASKTO,21,3).EQ.0)  GO TO 47

C--                Find divisions for generic

      JDVTO = 1
   31 IF (JDVTO.EQ.JQDVLL+1)  JDVTO = JQDVSY

      IF (JBIT(MASKNO,JDVTO).EQ.0)                   GO TO 37
      IF (JBYTET(IQKIND(KQT+JDVTO),MASKTO,1,23).EQ.0)   GO TO 37
      MASKCU = MSBIT1 (MASKCU,JDVTO)
   37 JDVTO = JDVTO + 1
      IF (JDVTO.LE.20)             GO TO 31

   47 IQRCU(KQT+JDVFR) = MASKCU
      JDVFR = JDVFR + 1
      IF (JDVFR.LE.20)             GO TO 21
C!    IX = MSBYT (0,JQSTOR,27,4)
C!    CALL MZXRPR (IX)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
