*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE ZDVCOP (IXDVFR,LIN,IXDVTO,LOUT,*)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzct.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/storparq.inc"
      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'ZDVCOP')
#include "zebra/q_jbit.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CALL MZSDIV(IXDVFR,1)
      JSTOR1 = JQSTOR
      JDIV1  = JQDIVI
      CALL MZSDIV(IXDVTO,1)
      JSTOR2 = JQSTOR
      JDIV2  = JQDIVI
      IF (JSTOR1.NE.JSTOR2) THEN
          WRITE(IQPRNT,'('' ZDVCOP -- Stores different '',2I5)')
     W         JSTOR1,JSTOR2
                                                           GO TO 998
      ENDIF

      IF (JDIV1.EQ.JDIV2)                                  GO TO 999

      IF (JDIV1.GT.JQDVLL.AND.JDIV1.LT.JQDVSY) THEN
          WRITE(IQPRNT,'('' ZDVCOP -- Division 1 ID invalid '',I5)')
     W         JDIV1
                                                           GO TO 998
      ENDIF
      IF (JDIV2.GT.JQDVLL.AND.JDIV2.LT.JQDVSY) THEN
          WRITE(IQPRNT,'('' ZDVCOP -- Division 2 ID invalid '',I5)')
     W         JDIV2
                                                           GO TO 998

      ENDIF
      CALL MZGARB(IXDVFR,IXDVTO)
      JQDVM1 = 0
      JQDVM2 = 0
      NQDVMV = 0

      NWORDS = LQEND(KQT+JDIV1)-LQSTA(KQT+JDIV1)

      IF (JBIT(IQMODE(KQT+JDIV2),JDVBFQ).EQ.IDVFWQ)  THEN


          IF (JDIV2.EQ.JQDVLL) THEN
              JNEXT  = JQDVSY
          ELSE
              JNEXT  = JDIV2 + 1
          ENDIF
          IF (NWORDS.GT.LQSTA(KQT+JNEXT)-LQSTA(KQT+JDIV2))   THEN
              WRITE(IQPRNT,'('' ZDVCOP -- TARGET DIVISION TOO SMALL'')')
                                                           GO TO 998
          ENDIF
          CALL UCOPY
     U    (LQ(KQS+LQSTA(KQT+JDIV1)),LQ(KQS+LQSTA(KQT+JDIV2)),NWORDS)
          LQEND(KQT+JDIV2) = LQSTA(KQT+JDIV2)+NWORDS
          NMOVE  = LQSTA(KQT+JDIV2)-LQSTA(KQT+JDIV1)
      ELSE


          IF (NWORDS.GT.LQEND(KQT+JDIV2)-LQEND(KQT+JDIV2-1))   THEN
              WRITE(IQPRNT,'('' ZDVCOP -- TARGET DIVISION TOO SMALL'')')
                                                           GO TO 998
          ENDIF
          CALL UCOPY
     U    (LQ(KQS+LQSTA(KQT+JDIV1)),LQ(KQS+LQEND(KQT+JDIV2)-NWORDS)
     U     ,NWORDS                                                  )
          LQSTA(KQT+JDIV2) = LQEND(KQT+JDIV2)-NWORDS
          NMOVE  = LQEND(KQT+JDIV2)-LQSTA(KQT+JDIV1)-NWORDS
      ENDIF

      LQMTA  = NQOFFS(1) + LQEND(1)
      MQDVGA = 0
      MQDVWI = 0
      CALL SBIT1(MQDVGA,JDIV2)
      CALL MZTABM
      LQMTE = LQMTA + JDIV2*8
      DO 300 L=LQMTA,LQMTE-1,8
          IF(LQ(L).EQ.JDIV2)   THEN
              LQ(L+1) = 1
          ELSE
              LQ(L+1) = -1
          ENDIF
  300 CONTINUE
      LQRTA = LQMTE
      LQTA   = LQRTA + 2
      LQ(LQTA-1) = LQSTA(KQT+JDIV1)
      LQ(LQTA  ) = LQSTA(KQT+JDIV1)
      LQ(LQTA+1) = LQEND(KQT+JDIV1)
      LQ(LQTA+2) = NMOVE
      LQ(LQTA+3) = 0
      LQTE       = LQTA + 4
      LQ(LQTE  ) = LQEND(KQT+JDIV1)
      LSUP   = LQLORG(KQS+LIN)
      IF (LSUP.NE.0)    THEN
          LQLORG(KQS+LIN+NMOVE) = LSUP+LOCF(LOUT)-LOCF(LIN)
      ELSE
          LQLORG(KQS+LIN+NMOVE) = 0
      ENDIF
#if defined(CERNLIB_QDEVZE)
      IF(NQDEVZ.NE.0) CALL CQDTAB (0)
#endif
      CALL MZRELB
      LOUT = LIN + NMOVE

                                                           GO TO 999

  998 RETURN 1
  999 RETURN
      END
