*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:57  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZPUDX (LP,NWP)

C-    Mark a dead area of NW words in the current store at adr L
C-    (slave to MZPUSH)

#include "zebra/zbcd.inc"
#include "zebra/mqsysh.inc"
C--------------    End CDE                             --------------
      DIMENSION    LP(9),NWP(9)

#include "zebra/q_sbit1.inc"
#include "zebra/q_sbyt.inc"

      L  = LP(1)
      NW = NWP(1)
      ND = NW - 10
      N  = MIN (10,NW)

      DO 12  J=0,N-1
   12 LQ(KQS+L+J) = 0

      IF (ND.GE.0)   THEN
          LQ(KQS+L)   = 12
          L = L + 9
          LQ(KQS+L-4) = IQLETT(4)
          LQ(KQS+L-1) = ND
        ELSE
          N = MSBYT (NW,N,17,6)
          LQ(KQS+L) = N
        ENDIF

      LQ(KQS+L) = MSBIT1 (LQ(KQS+L),IQDROP)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
