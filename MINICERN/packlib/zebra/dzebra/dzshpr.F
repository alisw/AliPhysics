*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:14  mclareni
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
      SUBROUTINE DZSHPR (LBANK,ILINK1,ILINK2,IDATA1,IDATA2)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZSHPR')

#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      CALL MZCHLS (NCHEKQ,LBANK)
      IF (IQFOUL.NE.0) THEN
          CALL DZBKDV(LBANK)
          IF (IQUEST(1).NE.0)                              GO TO 999
          WRITE(CQINFO,'(A,''/'',I8)') CQDIV,LBANK
          CALL DZTEXT(MSHP1Q,CDUMMQ,0)
                                                           GO TO 999
      ENDIF


      CALL DZBKDV(IQLS)
      IF (IQUEST(1).NE.0)                                  GO TO 999


      CALL DZBKXR(0)
      IF (IQUEST(1).NE.0)                                  GO TO 999

      LN = IQLN
      LS = IQLS
      NL = IQNL
      NS = IQNS
      ND = IQND
      IO = JBYT (LQ(KQS+IQLN+MBKIOQ),JBKIOQ,NBKIOQ)

      CALL DZBKHD
      IF (IQUEST(1).NE.0)                                  GO TO 999


      IF (ILINK1.GT.ILINK2.OR.NL.EQ.0)                     GO TO 100

      IL1 = ILINK1
      IL2 = ILINK2


      IF (IL1.LE.0.OR.IL1.GT.NL) IL1=1
      IF (IL2.LE.0.OR.IL2.GT.NL) IL2=NL

      IF (IFLOPT(MPOSDQ).NE.0)       THEN



          IBEGIN = MIN(NL,IL2)
          IF (IL1.LE.NS) THEN
              IBASE = NS
          ELSE
              IBASE = IL1-1
          ENDIF

          NDW    = IBEGIN-IBASE
          JDFD   = NDW + 1

          IF (NDW.GT.0)                                             THEN
              IF (NDW.GT.NQWKTB)         THEN
                  LBASE  = LS-IBEGIN-1
                  CALL DZDATA('REFERENCE links -- to be inverted --')
              ELSE
                  CALL UCOPIV (LQ(KQS+LS-IBEGIN),LQ(LQWKTB),NDW)
                  LBASE  = LQWKTB - 1 - KQS
                  CALL DZDATA('REFERENCE links')
              ENDIF
          ENDIF



          IBEGIN = MIN(NS,IL2)
          IBASE  = IL1-1

          NDW    = IBEGIN-IBASE
          JDFD   = NDW + 1

          IF (NDW.GT.0)                                             THEN
              IF (NDW.GE.NQWKTB)          THEN
                  LBASE  = LS-IBEGIN-1
                  CALL DZDATA('STRUCTURAL links -- to be inverted --')
              ELSE
                  CALL UCOPIV (LQ(KQS+LS-IBEGIN),LQ(LQWKTB),NDW)
                  LBASE  = LQWKTB - 1 - KQS
                  CALL DZDATA('STRUCTURAL links')
              ENDIF
          ENDIF

      ELSE


          CALL ZPAGE(IQPRNT,MIN((IL2-IL1-1)/10+2,5))
          CQLINE = ' --------  LINK part of bank  --------'
          CALL DZTEXT(0,CDUMMQ,1)
          DO 50 I=IL1,IL2,10
              WRITE (CQLINE,'(1X,I7,'' /'',10I12)')
     +        I,(LQ(KQS+LS-J),J=I,MIN(IL2,I+9))
              CALL DZTEXT(0,CDUMMQ,1)
   50     CONTINUE

      ENDIF


  100 IF (ND.EQ.0.OR.IDATA2.LT.IDATA1)                     GO TO 999

      ID1 = IDATA1
      ID2 = IDATA2


      IF (ID1.LE.0.OR.ID1.GT.ND) ID1=1
      IF (ID2.LE.0.OR.ID2.GT.ND) ID2=ND

      CALL DZIOPR(ID1,ID2)

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
