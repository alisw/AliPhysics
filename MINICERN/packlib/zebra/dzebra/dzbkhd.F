*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:02  mclareni
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
      SUBROUTINE DZBKHD
#include "zebra/bankparq.inc"
#include "zebra/questparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"

      CHARACTER CHROUT*(*),CHSTAK*6,CLIVE*1,CLORIG*8
      PARAMETER (CHROUT = 'DZBKHD')

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IF (JBIT(IQ(KQS+LS),IQDROP).EQ.0)            THEN

          CALL DZBKDV(LS)
          IF (IQUEST(1).NE.0)                              GO TO 999


          CALL DZBKUP(LS)


          CALL DZBKXR(LS)

          CLIVE = ' '
      ELSE
          CLIVE = '('
      ENDIF


      LORIG  = LQLORG(KQS+LS)
      IF (LORIG.NE.LNULL) THEN
          IF (LORIG.LT.IQTABV(KQT+13).OR.LORIG.GT.IQTABV(KQT+14)) THEN
              CLORIG = '_'
          ELSE
              WRITE(CLORIG,'(I8)') LQ(KQS+LORIG)
          ENDIF
      ELSE
          CLORIG = '0'
      ENDIF


      WRITE(CQLINE,
     W '(1X,A1,A4,''.'',I6,I9,''('',A ,'') SY/US/IO '',Z4,''/'',Z5,
     W   ''/'',Z4,'' NL/NS/ND'',I5,''/'',I5,''/'',I8,'' N/U/O/@O'',
#if (!defined(CERNLIB_QMIBM))&&(!defined(CERNLIB_QMVAX))
     W   I8,''/'',I8,''/'',I8,''/'',A)',                     ERR=999)
#endif
#if defined(CERNLIB_QMVAX)
     W   I8,''/'',I8,''/'',I8,''/'',A)',                IOSTAT=ISTAT)
#endif
#if defined(CERNLIB_QMIBM)
     W   I8,''/'',I8,''/'',I8,''/'',A)')
#endif
     W   CLIVE, IQWIDH(KQS+LS),IQWIDN(KQS+LS),LS,CQDIV,
     W   JBYT(IQ(KQS+LS),JBKSYQ,NBKSYQ),JBYT(IQ(KQS+LS),JBKUSQ,NBKUSQ),
     W   JBYT(LQ(KQS+LN+MBKIOQ),JBKIOQ,NBKIOQ),NL,NS,ND,
     W   LQ(KQS+LS),LQLUP(KQS+LS),LORIG,CLORIG

      IF (IFLOPT(MPOSQQ).EQ.0) CALL DZTEXT(0,CDUMMQ,1)
      IQUEST(1) = 0

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
