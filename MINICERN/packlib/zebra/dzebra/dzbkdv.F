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
      SUBROUTINE DZBKDV(LBANK)
#include "zebra/questparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/dzc1.inc"

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZBKDV')

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      JQDIVI = 1
      IF (LBANK.LT.LQEND(KQT+JQDVLL))                      GO TO 100
      IF (LBANK.GE.LQEND(KQT+NDVMXQ)) THEN
          CALL DZTEXT(MBKD1Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF
      JQDIVI = JQDVSY

  100 IF (LBANK.LT.LQEND(KQT+JQDIVI))                      GO TO 200
      JQDIVI = JQDIVI + 1
                                                           GO TO 100

  200 IF (LBANK.LT.LQSTA(KQT+JQDIVI)) THEN
          CALL DZTEXT(MBKD2Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF
      CALL UHTOC(IQDN1(KQT+JQDIVI),4,CQDIV(1:),4)
      CALL UHTOC(IQDN2(KQT+JQDIVI),4,CQDIV(5:),4)

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
