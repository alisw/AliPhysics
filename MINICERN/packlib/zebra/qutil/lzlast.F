*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:37  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION LZLAST (IXSTOR,LGO)

C-    FIND LAST BANK IN THE LINEAR STRUCTURE STARTING WITH BANK AT LGO
C-    RETURN BANK-ADR AS VALUE

#include "zebra/q.inc"
C--------------    END CDE                             --------------

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      LF = LGO
      IF (LF.EQ.0)                 GO TO 29

   21 L = LQ(KQS+LF)
      IF (L.EQ.0)                  GO TO 29
      LF = L
      GO TO 21

   29 LZLAST = LF
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
