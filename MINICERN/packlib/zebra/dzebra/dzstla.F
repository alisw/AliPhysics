*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:17  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZSTLA
      SAVE CLATYP
#include "zebra/divparq.inc"
#include "zebra/storparq.inc"
#include "zebra/questparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"

      CHARACTER CHROUT*(*),CHSTAK*6,            CLATYP(0:1)*9, CAKTIV*8
      PARAMETER (CHROUT = 'DZSTLA')
      DATA      CLATYP               /'PERMANENT','TEMPORARY'/

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      IF (LQSTA(KQT+JQDVSY).EQ.LQEND(KQT+JQDVSY))          GO TO 999

      LSYSB  = LQSYSS(KQT+MSYLAQ)
      CALL MZCHLS(NCHEKQ,LSYSB)
      IF (IQFOUL.NE.0)    THEN
          CALL DZTEXT(MARE1Q,CDUMMQ,0)
                                                           GO TO 999
      ENDIF


      NWTAB  = IQ(KQS+LSYSB+MLAUSQ)


      LENTRY = LSYSB + KQS + MLAUSQ

      CQMAP(1) = ' '
      CQMAP(2) = '  --- Link area parameters ---'
      CQMAP(3) = ' '
      CALL DZTEXT(0,CDUMMQ,3)

      DO 200 IENTRY = 1,(NWTAB-1)/NLAENQ

          LLAAR1 = IQ(LENTRY+MLAADQ)
          LLAARL = IQ(LENTRY+MLALTQ)
          NLANS  = IQ(LENTRY+MLANSQ)
          NAME1  = IQ(LENTRY+MLAN1Q)
          NAME2  = IQ(LENTRY+MLAN2Q)
          JTEMP  = JBIT(NLANS,JLATMQ)
          NTEMP  = NLATMQ*JTEMP
          NLANS  = JBYT(NLANS,JLANSQ,NLANSQ)-NTEMP
          IF                                                     (
     X      JTEMP.EQ.0
     X    .OR.
     X      (JTEMP.EQ.1.AND.LQ(KQS+LLAAR1+MLACTQ-1).NE.0)        )
     X                                                        THEN
              CAKTIV = '  ACTIVE'
          ELSE
              CAKTIV = 'INACTIVE'
          ENDIF

          LLAABS = LLAAR1 + LQSTOR
*                MAP addresses are in BYTES
#if !defined(CERNLIB_WORDMAP)
          LLAABS = LLAABS*4
#endif
          NDW    = LLAARL - LLAAR1  - NTEMP

          WRITE(CQLINE,'(1X,2A4,1X,A9,'' LIST AREA'',5X,
#if !defined(CERNLIB_OCTMAP)
     W    '' is at absolute '',Z8,'' NL/NS '',2I5,
#endif
#if defined(CERNLIB_OCTMAP)
     W    ''is at absolute '',O11,'' NL/NS '',2I5,
#endif
     W     5X,''status '',A8)'                    )
     W     NAME1,NAME2,CLATYP(JTEMP),LLAABS,NDW,NLANS,CAKTIV
          CALL DZTEXT(0,CDUMMQ,1)

          LENTRY = LENTRY + NLAENQ

  200 CONTINUE

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
