*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:42  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZRELX

C-    Execute normal relocation of links in banks and in link areas

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZRE, 4HLX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZRELX /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZRELX  ')
#endif


#include "zebra/qtrace.inc"

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.11)            GO TO 29
      WRITE (IQLOG,9821) LQMTA,LQMTE
 9821 FORMAT (1X/' DEVZE MZRELX,  Memory Occupation Table,'
     F,' LQMTA,LQMTE= ',2I8/16X,
     F' DIV ACT     NWS    LFBK  LLBK+1     JFE     JLE    FREE')

      WRITE (IQLOG,9822) (LQ(J),J=LQMTA,LQMTE-1)
 9822 FORMAT (16X,2I4,6I8)

      WRITE (IQLOG,9823) LQRTA,LQTA,LQTE,LQ(LQTA-1)
 9823 FORMAT (1H0,15X,'Link Relocation Table,'
     F,' LQRTA,LQTA,LQTE= ',3I8
     F/16X,'   LOC       L      LE    NREL  BG'
     F/30X,I8)

      IF (LQTE.LE.LQTA)            GO TO 25
      I  = LQRTA
      JA = LQTA - I
      JE = LQTE - I - 1
      WRITE (IQLOG,9824) (J,LQ(I+J),LQ(I+J+1),LQ(I+J+2),LQ(I+J+3),
     +                                           J=JA,JE,4)
 9824 FORMAT (16X,I6,3I8,I4)

   25 WRITE (IQLOG,9825) LQ(LQTE)
 9825 FORMAT (22X,I8)
   29 CONTINUE
#endif

      L = LQSYSS(KQT+1)
      IF (L.NE.0)  THEN
          IQ(KQS+L+3) = IQ(KQS+L+2) + NQLINK
          CALL MZRELL (IQ(KQS+L+1))
        ENDIF

      CALL MZRELB
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
