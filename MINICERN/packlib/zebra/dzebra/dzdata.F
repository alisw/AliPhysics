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
      SUBROUTINE DZDATA (CHTEXT)
      SAVE IHOLE
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"


      PARAMETER ( MLITXQ =  6 )
      PARAMETER ( NLINEQ = 24 )
      PARAMETER ( MLLKBQ = 17 , MLLKEQ = 24 )
      PARAMETER ( MLDRQ  = 10 )
      PARAMETER ( MLIDBQ = 11 , MLIDEQ = 14 )

      PARAMETER ( NENLNQ =  5 )
      PARAMETER ( NLNGRQ = 10 )

      CHARACTER*(*) CHTEXT

      CHARACTER CHROUT*(*),CHSTAK*6              ,CTEXT1*1,CSTART*3
      PARAMETER (CHROUT = 'DZDATA')
      DATA      IHOLE  /4H*HO*/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IF (NDW.LE.0)                                        GO TO 999

      CTEXT1 = CHTEXT(1:1)
      CSTART = ' '
      IF (LEN(CHTEXT).GT.1)   THEN
          IF (CTEXT1.EQ.'+') CSTART = CHTEXT(2:)
      ENDIF

    9 IF (CHTEXT.NE.CDUMMQ.AND.CTEXT1.NE.'+') THEN
          CQLINE      = ' '//CHTEXT
          CQLINE(60:) = '--------------------'
          CALL DZTEXT(0,CDUMMQ,1)
      ENDIF
      JDST = 1


   10 JDL  = JDST

   12 NSAME= IUSAME (LQ(LBASE+KQS+1),JDL,NDW,30,JSAME)

      JDL  = JSAME + NSAME - 1
      IF (JDL.NE.NDW)  JDL=10*(JDL/10)
      IF (JSAME.NE.JDST)                                   GO TO 20
   16 J = LBASE + JSAME
      N = JDL+1 - JSAME
      WRITE(CQLINE,'(T30,''====='',I9,'' WORDS from'',I9,
     X '' to'',I9,'' all contain'',Z18)') N,JSAME,JDL,LQ(J+KQS)
      CALL DZTEXT(0,CDUMMQ,1)
      IF (JDL.GE.NDW)                                      GO TO 999
      JDST= JDL + 1
                                                           GO TO 10

   20 JDE = JSAME - 1
      N   = JDE - JDST
      NPG    = N/NENLNQ + 1
      IF (NPG.GT.NLNGRQ) NPG = NLNGRQ
      NPG1= NPG + 1
      NGRV= N/(NPG*NENLNQ) + 1

      CALL ZPAGE (IQPRNT,NPG1)
      NGRP = (NQLNOR-NQUSED)/NPG1
      NGRV = MIN(NGRV,NGRP)

      NSTEP= NPG*NGRV
      JDE  = MIN(NDW,JDST-1+NENLNQ*NSTEP)
      IF (JDE.NE.NDW.AND.JDE.GE.JDL)                       GO TO 12
      JSAME= MAX(JDE+1,JSAME)
      NPG  = MIN(NPG,JDE+1-JDST)
      DO 300 JGROUP=1,NGRV
          IF (JGROUP.GT.1) CALL ZPAGE(IQPRNT,1)
          CQLINE = ' '//CSTART
          ILINE  = MLITXQ

          DO 200 JLINE=1,NPG
              JD = JDST

              DO 100 JWORD=1,NENLNQ
                  IF (JD.GT.JDE)                           GO TO 100
                  CALL DZTYP
                  IPBEG  = ILINE  +  1
                  IPEND  = ILINE  + NLINEQ
                  IF     (JTYP.LT.0)           THEN
                      WRITE(CQLINE(ILINE+1:ILINE+6),'(I6)') JD+IBASE
                      IF     (IQFOUL.EQ.0)              THEN
                          IF (JBIT(IQ(KQS+IQLS),IQDROP).EQ.1)
     I                    CQLINE(ILINE+MLDRQ:ILINE+MLDRQ) = '('
                          IF (IQND.LT.0) IQID = IHOLE
                          WRITE
     W                    (CQLINE(ILINE+MLIDBQ:ILINE+MLIDEQ),'(A4)')IQID
                      ELSEIF (IQFOUL.GT.0)              THEN
                          CQLINE(ILINE+MLDRQ:ILINE+MLDRQ)   = '?'
                          CQLINE(ILINE+MLIDBQ:ILINE+MLIDEQ) = '****'
                      ELSE
                          CQLINE(ILINE+MLDRQ:ILINE+MLDRQ)   = '?'
                          CQLINE(ILINE+MLIDBQ:ILINE+MLIDEQ) = '-'
                      ENDIF
                      WRITE
     W                (CQLINE(ILINE+MLLKBQ:ILINE+MLLKEQ),'(I8)')
     W                LQ(KQS+LBASE+JD)
                      ILINE  = ILINE + NLINEQ
                  ELSEIF (JTYP.EQ.1)           THEN
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMND3)||defined(CERNLIB_QMCRY)||defined(CERNLIB_QMDOS)
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,1X,''Z'',Z16)')
#endif
#if (!defined(CERNLIB_QMIBM))&&(!defined(CERNLIB_QMIBX))&&(!defined(CERNLIB_QMND3))&&(!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMDOS))
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,1X,''Z'',Z16.16)')
#endif
     W                      JD+IBASE,LQ(KQS+LBASE+JD)
                  ELSEIF (JTYP.EQ.2)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,I18)')
     W                      JD+IBASE,LQ(KQS+LBASE+JD)
                  ELSEIF (JTYP.EQ.3)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,8X,F10.1)')
     W                      JD+IBASE,Q(KQS+LBASE+JD-NOFLIQ)
                  ELSEIF (JTYP.EQ.4)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,6X,F12.3)')
     W                      JD+IBASE,Q(KQS+LBASE+JD-NOFLIQ)
                  ELSEIF (JTYP.EQ.5)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,2X,F16.7)')
     W                      JD+IBASE,Q(KQS+LBASE+JD-NOFLIQ)
                  ELSEIF (JTYP.EQ.6)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,2X,E16.7)')
     W                      JD+IBASE,Q(KQS+LBASE+JD-NOFLIQ)
                  ELSEIF (JTYP.EQ.7)           THEN
                      WRITE(CQLINE(IPBEG:IPEND),'(I6,13X,''"'',A4)')
     W                      JD+IBASE,LQ(KQS+LBASE+JD)
                  ENDIF
                  ILINE  = IPEND
                  JD = JD + NSTEP
  100         CONTINUE
              CALL DZTEXT(0,CDUMMQ,1)
              CQLINE = ' '//CSTART
              ILINE  = MLITXQ
  200     JDST= JDST + 1

  300 CONTINUE

      JDST= JDE + 1
      IF (JDST.LT.JSAME)                                   GO TO 20
      IF (JDST.LT.NDW)                                     GO TO 16

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
