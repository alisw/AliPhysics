*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:17  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:18  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZBOOK (IXP,LP,LSUPP,JBP, CHIDH,NL,NS,ND,NIOP,NZP)

C-    Book a bank, user called

#include "zebra/zbcd.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcl.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXP(9),LP(9),LSUPP(9),JBP(9),NIOP(9),NZP(9)
      CHARACTER    CHIDH*(*)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZBO, 4HOK   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZBOOK /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZBOOK  ')
#endif

#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      NQID = IQQUES
      NIO  = MIN (4, LEN(CHIDH))
      IF  (NIO.NE.0)  CALL UCTOH (CHIDH,NQID,4,NIO)
      NQNL  = NL
      NQNS  = NS
      NQND  = ND
      NQBIA = JBP(1)

      IODORG = NIOP(1)
      NIO = JBYT (IODORG,12,4)

      IF (NIO.EQ.0)  THEN
          NQIOCH(1) = IODORG
        ELSE
          CALL UCOPY (NIOP,NQIOCH,NIO+1)
          NQIOSV(1) = 0
        ENDIF

      CALL MZLIFT (IXP,LP,LSUPP,63, NQID, NZP)

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
