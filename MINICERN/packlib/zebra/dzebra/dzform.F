*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:08  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZFORM(IXSTOR,LBANK)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zunit.inc"
#include "zebra/eqlqform.inc"
#include "zebra/dzc1.inc"
#include "zebra/questparq.inc"

      CHARACTER   CHROUT*(*)
      PARAMETER (CHROUT = 'DZFORM')

#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      IF (LBANK.EQ.0)    THEN


          IF (LQFORM.EQ.0)                                 GO TO 999

          CQMAP(1) = ' '
          CQMAP(2) =' DZFORM -- List of the IO characteristics '//
     X             'known to ZEBRA'
          CQMAP(3) = ' ========================================'//
     X             '=============='
          CALL DZTEXT(0,CDUMMQ,3)

          LIOD   = LQ(KQSP+LQFORM-2)
          LID    = LQFORM
          DO 200 I = 1,2
              LIX    = LQ(KQSP+LID-1)
              NWID   = IQ(KQSP+LID+1)
              IF (NWID.GT.0) THEN
                  CALL SORTZV
     X            (IQ(KQSP+LID+4),LQ(LQWKTB),MIN(NWID,NQWKTB),0,0,0)
              ENDIF
              DO 100 J = 0,NWID-1
                  CQLINE = ' '
                  INDXA  = LQ(LQWKTB+J)
                  IXIOD  = IQ(KQSP+LIX+INDXA)
                  NWIO   = JBYT(IQ(KQSP+LIOD+IXIOD+1),7,5) - 1
                  CALL MZIOCR(IQ(KQSP+LIOD+IXIOD+1))
                  IF(IQUEST(1).LT.0) THEN
                      WRITE(CQINFO,'(A4)') IQ(KQSP+LID+INDXA+3)
                      CALL DZTEXT(MFOR1Q,CDUMMQ,0)
                                                           GO TO 999
                  ELSE
                      IQUEST(1) = 0
                  ENDIF
                  WRITE(CQLINE(1:37),
     X            '('' Format for ID='',A4,'' is NWIO=1/2+'',I2,1X)')
     X             IQ(KQSP+LID+INDXA+3),NWIO
                  CALL DZFORP
  100         CONTINUE
              LID=LQ(KQSP+LID)
  200     CONTINUE

      ELSE


          IF (IXSTOR.NE.NCHEKQ) THEN
              CALL MZSDIV(IXSTOR,-1)
          ENDIF

          CALL MZCHLS(NCHEKQ,LBANK)

          IF(IQFOUL.NE.0) THEN
              CALL DZBKDV(LBANK)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO,'(A,''/'',I10,'','',I4)') CQDIV,LBANK,IQFOUL
              CALL DZTEXT(MFOR2Q,CDUMMQ,0)
                                                           GO TO 999
          ENDIF

          CALL MZIOCR(LQ(KQS+IQLN))

          IF(IQUEST(1).LT.0) THEN
              CALL DZBKDV(LBANK)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO,'(A,''/'',I10,'','',A4)') CQDIV,LBANK,IQID
              CALL DZTEXT(MFOR3Q,CDUMMQ,0)
                                                           GO TO 999
          ELSE
              IQUEST(1) = 0
          ENDIF

          CQLINE = ' '

          WRITE(CQLINE(1:37),
     X    '('' DZFORM for ID='',A4,'' is NWIO=1/2+'',I2,1X)') IQID,IQNIO

          CALL DZFORP
      ENDIF

  999 RETURN
      END
