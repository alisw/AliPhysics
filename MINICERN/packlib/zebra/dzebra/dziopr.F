*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZIOPR (IFIRST,ILAST)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzioc.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/questparq.inc"
#include "zebra/bkfoparq.inc"

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZIOPR')

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      IF (IFLOPT(MPOSZQ).NE.0) THEN
          MFO(1) = IFOBIQ
          MFO(2) = -1
      ELSE
          CALL MZIOCR(LQ(KQS+LN))

          IF (IQUEST(1).LT.0) THEN
              CALL DZBKDV(LS)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO,'(A,''/'',I8,'','',I8,''/'',Z8)')
     +         CQDIV,LS,LN,LQ(KQS+LN)
              CALL DZTEXT(MIOP1Q,CDUMMQ,0)
                                                           GO TO 999
          ELSE
              IQUEST(1) = 0
          ENDIF
      ENDIF

      IF (IFLOPT(MPOSDQ).NE.0) THEN
          CALL DZIOPD(IFIRST,ILAST)
      ELSE
          CALL DZIOPH(IFIRST,ILAST)
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
