*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*------------------------------------------------------------
#include "zebra/pilot.h"
      SUBROUTINE ZPAGE(LUN,LINES)
#include "zebra/zheadp.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
      DIMENSION    LUNS(3)
      EQUIVALENCE (LUNS(1),IQPRNT)
      IF (LINES.LT.0)                                      GO TO 11
* For C printer
#if defined(CERNLIB_QMIBM)
    7 IF(NQUSED.LE.NQLMAX-2)                               GO TO 8
      NQUSED= NQUSED - NQLMAX + 2
#endif
#if !defined(CERNLIB_QMIBM)
    7 IF(NQUSED.LE.NQLMAX)                                 GO TO 8
      NQUSED= NQUSED - NQLMAX
#endif
                                                           GO TO 7
    8 IF (NQUSED+LINES.LE.NQLNOR-2)  RETURN
   11 JLUN = IUCOMP (LUN,LUNS,3) + 1
      IQPAGE = NQPAGE(JLUN) + 1
      NQPAGE(JLUN) = IQPAGE
      NQUSED = NQLPTH + 1
      WRITE (LUN,'(''1'',20A4,28X,I6,''/'',I4,''  PAGE'',I5)')
     W       IQHEAD,IQDATE,IQTIME,IQPAGE

  999 RETURN
      END
