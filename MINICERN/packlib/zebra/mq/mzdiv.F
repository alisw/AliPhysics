*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:22  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:18  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZDIV (IXSTOR,IXDIV,CHNAME,NWAP,NWMP,CHOPT)

C-    Create new division, user called

#include "zebra/zbcd.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIV(9),NWAP(9),NWMP(9)
      CHARACTER    *(*) CHNAME, CHOPT
      DIMENSION    NAME(2)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZDI, 4HV    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZDIV  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZDIV   ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_sbyt.inc"
#include "zebra/q_shiftl.inc"


#include "zebra/qtrace.inc"

      NWALLO = NWAP(1)
      NWMAX  = NWMP(1)
      CALL UOPTC (CHOPT,'RMLPC',IQUEST)
      MODE   = IQUEST(1) + 2*IQUEST(2)
      KIND   = MIN (2, IQUEST(3)+2*IQUEST(4)) + 1
      IOPTC  = IQUEST(5)

C--                Printing name of division

      NAME(1) = IQBLAN
      NAME(2) = IQBLAN
      N = MIN (8, LEN(CHNAME))
      IF (N.NE.0)  CALL UCTOH (CHNAME,NAME,4,N)

#include "zebra/qstore.inc"

C--       option M, match mode of new division to neighbour

      IF (MODE.LT.2)               GO TO 29
      MODE = 0
      IF (KIND.GE.2)               GO TO 24
      IF (JQDVLL.EQ.2)             GO TO 29
      IF (IQMODE(KQT+JQDVLL).NE.0)    GO TO 29
      GO TO 28

   24 IF (JQDVSY.EQ.20)            GO TO 28
      IF (IQMODE(KQT+JQDVSY+1).NE.0)  GO TO 29
   28 MODE = 1
   29 CONTINUE

#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.0)
     +WRITE (IQLOG,9029) NAME,JQSTOR,NWALLO,NWMAX,MODE,KIND
 9029 FORMAT (1X/' MZDIV.   Initialize Division  ',2A4,'  in Store',I3
     F/10X,'NW/NWMAX=',2I7,',  MODE/KIND=',2I3)
#endif
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

C--                Check parameters

      IF (NWALLO.LT.100)           GO TO 91
      IF (NWMAX .LT.NWALLO)        GO TO 91
      IF (MODE  .LT.0)             GO TO 91
      IF (MODE  .GE.2)             GO TO 91
      IF (KIND  .LT.0)             GO TO 91
      IF (KIND  .GE.4)             GO TO 91
      IF (JQDVLL+1.EQ.JQDVSY)      GO TO 92

C----              Push down low divisions, and maybe system division

      MQDVGA = 0
      MQDVWI = 0
      NQDVMV = -NWALLO
      IQTNMV = 0
      JQSTMV = JQSTOR
      JQDVM1 = 2
      IF (KIND.LT.2)  THEN
          KIND   = 1
          JQDVM2 = JQDVLL
        ELSE
          JQDVM2 = JQDVSY
          MQDVGA = MSBIT1 (MQDVGA,JQDVSY)
        ENDIF
      JQDIVI = JQDVM2
      NRESV1 = LQSTA(KQT+2) - LQEND(KQT+1) - NQMINR
      NRESV2 = LQEND(KQT+2) - LQ2END
      NQRESV = MIN (NRESV1,NRESV2) + NQDVMV
      IF (NQRESV.LT.0)             GO TO 93

      CALL MZTABM
      CALL MZTABS
      CALL MZTABR
      CALL MZTABX
      CALL MZTABF
      CALL MZRELX
      CALL MZMOVE

C--                Create high division, re-number system division

      IF (JQDIVI.NE.JQDVSY)        GO TO 61
      JQDVSY = JQDVSY - 1

       LQSTA(KQT+JQDVSY) =  LQSTA(KQT+JQDVSY+1)
       LQEND(KQT+JQDVSY) =  LQEND(KQT+JQDVSY+1)
      NQDMAX(KQT+JQDVSY) = NQDMAX(KQT+JQDVSY+1)
      IQMODE(KQT+JQDVSY) = IQMODE(KQT+JQDVSY+1)
      IQKIND(KQT+JQDVSY) = IQKIND(KQT+JQDVSY+1)
       IQRTO(KQT+JQDVSY) =  IQRTO(KQT+JQDVSY+1)
       IQRNO(KQT+JQDVSY) =  IQRNO(KQT+JQDVSY+1)
      NQDINI(KQT+JQDVSY) = NQDINI(KQT+JQDVSY+1)
      NQDWIP(KQT+JQDVSY) = NQDWIP(KQT+JQDVSY+1)
      NQDGAU(KQT+JQDVSY) = NQDGAU(KQT+JQDVSY+1)
      NQDGAF(KQT+JQDVSY) = NQDGAF(KQT+JQDVSY+1)
      NQDPSH(KQT+JQDVSY) = NQDPSH(KQT+JQDVSY+1)
      NQDRED(KQT+JQDVSY) = NQDRED(KQT+JQDVSY+1)
      NQDSIZ(KQT+JQDVSY) = NQDSIZ(KQT+JQDVSY+1)
       IQDN1(KQT+JQDVSY) =  IQDN1(KQT+JQDVSY+1)
       IQDN2(KQT+JQDVSY) =  IQDN2(KQT+JQDVSY+1)

      NQDWIP(KQT+JQDIVI) = 0
      NQDGAU(KQT+JQDIVI) = 0
      NQDGAF(KQT+JQDIVI) = 0
      NQDPSH(KQT+JQDIVI) = 0
      NQDRED(KQT+JQDIVI) = 0
      NQDSIZ(KQT+JQDIVI) = 0

      LSTA   = LQEND(KQT+JQDVSY)
      GO TO 64

C--                Create low division JQDVLL

   61 LSTA = LQEND(KQT+JQDVLL)
      MOLL = JBIT (IQMODE(KQT+JQDVLL),1)

      JQDVLL = JQDVLL + 1
      JQDIVI = JQDVLL
      IF (MOLL.EQ.0)  THEN
          LSTA = MAX (LSTA, LQSTA(KQT+JQDVLL-1)+IQTABV(KQT+15))
          LSTA = MIN (LSTA, LQSTA(KQT+JQDVSY)-NWALLO)
        ENDIF
      IQTABV(KQT+15) = NWALLO

C--                Update store tables

   64 IF (MODE.NE.0) LSTA=LSTA+NWALLO

      MKIND = MSBIT1 (0,JQDIVI)
      MKIND = MSBIT1 (MKIND,20+KIND)

      MREF  = 0
      IF (KIND.GE.3)               GO TO 65
      IF (IOPTC.NE.0)              GO TO 65
      MREF = ISHFTL (3, 20)

   65 LQSTA(KQT+JQDIVI)  = LSTA
      LQEND(KQT+JQDIVI)  = LSTA
      NQDMAX(KQT+JQDIVI) = NWMAX
      IQMODE(KQT+JQDIVI) = MODE
      IQKIND(KQT+JQDIVI) = MKIND
      IQRTO(KQT+JQDIVI)  = MREF
      IQRNO(KQT+JQDIVI)  = 9437183
      NQDINI(KQT+JQDIVI) = NWALLO
      IQDN1(KQT+JQDIVI)  = NAME(1)
      IQDN2(KQT+JQDIVI)  = NAME(2)

      IQTABV(KQT+8) = JQDVLL
      IQTABV(KQT+9) = JQDVSY
      CALL MZXRUP


C--                Return IXDIV

#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.0)
     +WRITE (IQLOG,9079) JQDIVI
 9079 FORMAT (10X,'Division',I3,' initialized.')
#endif

      IXDIV(1) = MSBYT (JQSTOR,JQDIVI,27,4)
#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

   93 NQCASE = 1
      NQFATA = 4
      IQUEST(18) = NQRESV
      IQUEST(19) = -NQDVMV
      IQUEST(20) = NRESV1
      IQUEST(21) = NRESV2

   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 7
      IQUEST(11) = NAME(1)
      IQUEST(12) = NAME(2)
      IQUEST(13) = NWALLO
      IQUEST(14) = NWMAX
      IQUEST(15) = MODE
      IQUEST(16) = KIND
      IQUEST(17) = IOPTC
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
