*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:03  mclareni
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
      SUBROUTINE DZBKUP(LBANK)
      SAVE LASTAK,NLASTR
#include "zebra/bankparq.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/dzc1.inc"

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZBKUP')

      PARAMETER (NSTAKQ = 200)
      INTEGER LASTAK(NSTAKQ)

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      IF (LBANK.EQ.0) THEN


          LLASTR = 1
          NLASTR = 0

          IF (LQSTA(KQT+JQDVSY).EQ.LQEND(KQT+JQDVSY))      GO TO 999

          LSYSB  = LQSYSS(KQT+MSYLAQ)
          IF(LSYSB.GT.0)                   THEN
              CALL MZCHLS(NCHEKQ,LSYSB)
              IF (IQFOUL.NE.0)   THEN
                  WRITE(CQINFO,'(''Store = '',2A4)')
     X              NQSNAM(1),NQSNAM(2)
                  CALL DZTEXT(MARE1Q,CDUMMQ,0)
                  IQUEST(1) = 1
                                                           GO TO 999
              ENDIF


              NWTAB  = IQ(KQS+LSYSB+MLAUSQ)


              LENTRY = LSYSB + KQS + MLAUSQ

              DO 100 IENTRY = 1,(NWTAB-1)/NLAENQ

                  LLAAR1 = IQ(LENTRY+MLAADQ)
                  NLANS  = IQ(LENTRY+MLANSQ)
                  JTEMP  = JBIT(NLANS,JLATMQ)
                  NTEMP  = NLATMQ*JTEMP
                  NLANS  = JBYT(IQ(LENTRY+MLANSQ),JLANSQ,NLANSQ) - NTEMP
                  IF (NLANS.GT.0)   THEN
                      IF (NLASTR.GE.NSTAKQ/2) THEN
                         WRITE(CQINFO,'('' LLASTR NSTAKQ '',2I6)')
     X                     LLASTR,NSTAKQ
                         CALL DZTEXT(MBKU1Q,CDUMMQ,0)
                         IQUEST(1) = 1
                                                           GO TO 999
                      ENDIF
                      LASTAK(LLASTR  ) = LLAAR1+NTEMP
                      LASTAK(LLASTR+1) = LLAAR1+NTEMP+NLANS-1
                      NLASTR = NLASTR+1
                      LLASTR = LLASTR+2
                  ENDIF

                  LENTRY = LENTRY + NLAENQ
  100         CONTINUE
          ENDIF
                                                           GO TO 999
      ELSE


          LUP  = LQLUP(KQS+LBANK)
          IF (LUP.NE.0) THEN
              CALL MZCHLS(NCHEKQ,LUP)
              IF (IQFOUL.NE.0) THEN
                 CALL DZBKDV(LBANK)
                 IF (IQUEST(1).NE.0)                       GO TO 999
                 CQINFO = CQDIV//'/'
                 WRITE(CQINFO(10:),'(I8,'','',I8,'','',I4)')
     X           LBANK,LUP,IQFOUL
                 CALL DZTEXT(MBKU2Q,CDUMMQ,0)
                 IQUEST(1) = 1
                                                           GO TO 999
              ENDIF
          ENDIF


          LSUP = LQLORG(KQS+LBANK)
          IF (LSUP.EQ.0.AND.LUP.EQ.0)                      GO TO 999
          IF (LSUP.LT.LQSTA(KQT+1).OR.
     X        LSUP.GT.LQSTA(KQT+NDVMXQ+1)) THEN
              LOCSUP = LSUP
              DO 200 I = 1,NLASTR
                  IF (LOCSUP.GE.LASTAK(I*2-1)
     X            .AND.
     X                LOCSUP.LE.LASTAK(I*2))               GO TO 300
  200         CONTINUE
              CALL DZBKDV(LS)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO,'(A,''/'',I8,'','',I8,''('',Z8,'')'')')
     X        CQDIV,LS,LSUP,LSUP+LQSTOR
              CALL DZTEXT(MBKU3Q,CDUMMQ,0)
              IQUEST(1) = 1
                                                           GO TO 999
          ENDIF


  300     IF (LQ(LSUP+KQS).NE.LS) THEN
              WRITE(CQINFO,'(I8,''('',Z8,'')'',2I8)')
     X        LSUP,LSUP+LQSTOR,LQ(LSUP+KQS),LS
              CALL DZTEXT(MBKU4Q,CDUMMQ,0)
              IQUEST(1) = 1
                                                           GO TO 999
          ENDIF
                                                           GO TO 999
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
