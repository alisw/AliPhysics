*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_QDEBPRI)
      SUBROUTINE FZIDUM (MREC,NW4USE)

C-    Dump start of or full physical record read

#include "zebra/zunit.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    MREC(3600)
#if defined(CERNLIB_HEX)
      PARAMETER    (NWPLIN=8)
#endif
#if !defined(CERNLIB_HEX)
      PARAMETER    (NWPLIN=10)
#endif

      N  = NW4USE
      IF (LOGLVI.LT.4)  N=16

      WRITE (IQLOG,9047) N
 9047 FORMAT (1X/' FZIPHx-  Dump first',I6,' words of block:')

      JE = 0
   48 JA = JE + 1
      JE = MIN (N,JA+NWPLIN-1)
      WRITE (IQLOG,9048) JA,(MREC(J),J=JA,JE)
      IF (JE.LT.N)                 GO TO 48
      RETURN
#if (defined(CERNLIB_HEX))&&(defined(CERNLIB_B32))
 9048 FORMAT (7X,I7,8(2X,Z8))
#endif
#if (defined(CERNLIB_HEX))&&(!defined(CERNLIB_B32))
 9048 FORMAT (7X,I7,8(1X,Z9))
#endif
#if !defined(CERNLIB_HEX)
 9048 FORMAT (7X,I7,5(2X,O12)/14X,5(2X,O12))
#endif
      END
#endif
