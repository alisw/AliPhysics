*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:31  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION LZBYT (IXSTOR,LGO,ITEXT,JTH,NBITS)

C-    FIND BANK CONTAINING ITEXT IN STATUS-BYTE  (JTH,NBITS)
C-    SEARCH LINEAR STRUCTURE STARTING WITH BANK AT LGO
C-    RETURN BANK-ADR AS VALUE

#include "zebra/q.inc"
C--------------    END CDE                             --------------

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      L = LGO
   21 IF (L.EQ.0)                              GO TO 29
      IF (JBYT(IQ(KQS+L),JTH,NBITS).EQ.ITEXT)  GO TO 29
      L = LQ(KQS+L)
      GO TO 21

   29 LZBYT = L
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
