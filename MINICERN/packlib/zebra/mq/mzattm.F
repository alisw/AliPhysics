*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:16  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZATTM (IXSTOP,CHNAME
     +,                  MEMOR,LSTAP,NWMP,NWEXP,CHOPT,ITABLE)

C-    Attach flat memory as a Zebra store region, user called

#include "zebra/zbcd.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXSTOP(9),MEMOR(99),LSTAP(9),NWMP(9),NWEXP(9)
      DIMENSION    ITABLE(400)
      CHARACTER    *(*) CHNAME,CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZAT, 4HTM   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZATTM /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZATTM  ')
#endif

#include "zebra/q_sbit1.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_locf.inc"


#include "zebra/qtrace.inc"

      IXSTOR = IXSTOP(1)
      LSTA   = LSTAP(1)
      NWMEM  = NWMP(1)
      NWEX   = NWEXP(1)
      CALL UOPTC (CHOPT,'QIR',IQUEST)
      LOGQ = IQUEST(1)
      IREI = IQUEST(2)
      IONL = IQUEST(3)

      LASTOR = LOCF(MEMOR(1)) - 1
      LEND20 = 0

C------            Reset size of an existing read-only store

      IF (IXSTOR.EQ.0)             GO TO 21
      CALL MZSDIV (IXSTOR,0)
      IF (NQALLO(JQSTOR).GE.0)     GO TO  91
      IF (NQALLO(JQSTOR).LT.-3)    GO TO  92
      IF (IREI.NE.0)               GO TO 24
      IF (NQALLO(JQSTOR).EQ.-1)    GO TO  92
      GO TO 61

C------            Initialize an new read-only store

   21 JQSTOR = NQSTOR + 1
      IF   (NQSTOR)           93, 24, 22
   22 JQSTOR = IUFIND (-1,NQALLO(2),1,NQSTOR)
   24 CALL VZERO (KQT,25)

C--                Printing name of store

      NQSNAM(1) = IQBLAN
      NQSNAM(2) = IQBLAN
      N = MIN (8, LEN(CHNAME))
      IF (N.NE.0)  CALL UCTOH (CHNAME,NQSNAM,4,N)

C--                Set log level

      NQLOGL = NQLOGD
      IF (LOGQ.NE.0)  NQLOGL=-2

C--                Calculate table off-set

      LOCT = LOCF (ITABLE(1)) - 1
      KQT  = LOCT - LQATAB

#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.0)  THEN
          KQS = LASTOR - LQASTO
          WRITE (IQLOG,9021) JQSTOR,NQSNAM(1),NQSNAM(2)
     +,                  LASTOR,LOCT,LASTOR,LOCT,KQS,KQT,KQS,KQT
        ENDIF
 9021 FORMAT (1X/' MZATTM.  Attach Memory as Store',I3,'  in ',2A4,
     F/10X,'with Store/Table at absolute adrs',2I12
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))
     F/40X,'OCT',2(1X,O11)/40X,'OCT',2(1X,O11)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))
     F/40X,'HEX',2(1X,Z11)/40X,'HEX',2(1X,Z11)
#endif
#if defined(CERNLIB_QPRINT)
     F/30X,'relative adrs',2I12)
#endif

      IF (JQSTOR.GE.16)            GO TO  94

C----              Initialize divisions 1 + 2 + system

      NQOFFT(JQSTOR+1) = KQT
      CALL VZERO (IQTABV(KQT+1),NQTSYS)
      CALL VBLANK (IQDN1(KQT+1), 40)

      JQDVLL = 2
      JQDVSY = 20
      IQDN1(KQT+20) = IQLETT(19)
      IQDN1(KQT+2)  = IQLETT(4)

      IQKIND(KQT+1) = MSBIT1 (1,21)
      IQDN1(KQT+1)  = NQSNAM(1)
      IQDN2(KQT+1)  = NQSNAM(2)

      CALL UCOPY (IQCUR,IQTABV(KQT+1),16)

C--                Return IXSTOR
      IF (IXSTOR.EQ.0)  THEN
          IXSTOP(1) = ISHFTL (JQSTOR,26)
        ENDIF

      IF (JQSTOR.GT.NQSTOR)  NQSTOR = JQSTOR
      NQALLO(JQSTOR) = -3 + IONL

C------            Set size of the store

   61 LQSTOR = LASTOR
      KQS    = LQSTOR - LQASTO
      NQOFFS(JQSTOR+1) = KQS
      IQTABV(KQT+1)       = LQSTOR

      LEND1  = LSTA  + NWMEM
      LEND20 = LEND1 + NWEX

      NQDMAX(KQT+1) = NWMEM + NWEX
      LQSTA(KQT+1)  = LSTA
      LQEND(KQT+1)  = LEND1
      LQSTA(KQT+2)  = LEND20
      LQEND(KQT+2)  = LEND20
      LQSTA(KQT+20) = LEND20
      LQEND(KQT+20) = LEND20
      LQSTA(KQT+21) = LEND20
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.2)
     +  WRITE (IQLOG,9089) JQSTOR,NQSNAM(1),NQSNAM(2),LASTOR
     +,                    LSTA,NWMEM,NWEX
 9089 FORMAT (' MZATTM-  Store',I3,'  in ',2A4,' at adr',I12
     F/10X,'Memory starting at LSTA=',I8,' with',2I8,' words.')
#endif

#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

   94 NQCASE = NQCASE + 1
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = 4
      IQUEST(11) = NQSNAM(1)
      IQUEST(12) = NQSNAM(2)
      IQUEST(13) = LEND20
      IQUEST(14) = IXSTOR
      JQSTOR = -1
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
