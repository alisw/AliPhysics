*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:06  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:16  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZINIT (IXPARA,CHPA1,CHPA2,IPA3,IPA4,IPA5,IPA6)

C-    Initialize JZ91 package


#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IPA3(7),IPA4(7),IPA5(7),IPA6(7)
      CHARACTER    CHPA1*(*), CHPA2*(*)
#if defined(CERNLIB_A4)
      CHARACTER    CHIAM*4
#endif
#if defined(CERNLIB_A8)
      CHARACTER    CHIAM*8
#endif
#if defined(CERNLIB_EQUHOLCH)
      EQUIVALENCE (CHIAM, NQME(1))
#endif

      DIMENSION    INIT(4),          MEXTR(4)
      EQUIVALENCE (INIT(1),JQTIME), (NACCE,NQME(6))
      DIMENSION    MMJZ91(5), MMCALL(5)
#if defined(CERNLIB_QDEBUG)
      DIMENSION    MMFL(5)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZIN, 4HIT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZINIT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZINIT  ')
#endif
#if defined(CERNLIB_QDEBUG)

      DATA  MEND   /4HEND /
      DATA  MMFL   /4HJZFL,  0,  0, -7, 2 /
#endif
      DATA  MMJZ91 /4HJZ91, -7, -7, -7, 2 /
      DATA  MMCALL /4HJZCA, -7,  0, -7, 0 /

      DATA    MJQAN/4HJZAN/, MJQFL/4HJZFL/

#include "zebra/q_shiftl.inc"

#include "zebra/qtraceq.inc"

      IF (JQLEV.NE.-1)       CALL ZFATAM ('re-init of JZ91.')
      CALL VZERO (JQTIME,16)

      CALL MZSDIV (IXPARA,-7)
      JQSTJZ = JQSTOR
      IXSTJZ = ISHFTL(JQSTJZ,26)

      JQMLEV = IPA3(1)
      JQCBNL = IPA4(1)
      JQCBND = IPA5(1)

      CALL UOPTC (CHPA2,'TQE',IQUEST)
      JQTIME = IQUEST(1)
      JQLLEV = 2
      IF (IQUEST(2).NE.0)  JQLLEV = 0
      IF (IQUEST(3).NE.0)  JQLLEV = 1

      CHIAM = CHPA1
#if !defined(CERNLIB_EQUHOLCH)
      CALL UCTOH (CHIAM, NQME(1),4,4)
#endif
      NQME(2) = 10
      NQME(3) = 10

      N = IPA6(1)
      IF (N.GE.0)  THEN
          N =  MIN (N,5)
          CALL UCOPY (IPA6(2),NQME(2),N)
        ENDIF

      IF (JQLLEV.GE.2)  WRITE (IQLOG,9001) JQSTJZ,INIT
 9001 FORMAT (/' JZINIT.  JZ91 in store',I2,
     F', IFTIME,MAXLEV,NL,ND =',7I6)

      MMCALL(2) = JQCBNL
      MMCALL(4) = JQCBND
      JQNACC   = 9 + NACCE
      MEXTR(1) = 3
      MEXTR(2) = NQME(3)
      MEXTR(3) = NQME(4)
      MEXTR(4) = NQME(5)

C--                Lift division JZ91

      NWINT = 1000
      NWMAX = LQSTA(KQT+21) / 2
      CALL MZDIV  (IXSTJZ,IXDVJZ,'JZ91',NWINT,NWMAX,'ML')
      CALL MZLINK (IXSTJZ,'JZ91',LQJZ,LQUP,LQAN)

C----              Lift main bank

      MMJZ91(2) = 2*JQMLEV + 6
      MMJZ91(3) =   JQMLEV + 5
      MMJZ91(4) = 2*JQMLEV + 2
      CALL MZLIFT (IXDVJZ,L,LQJZ,1,MMJZ91,0)
      IQ(KQS+LQJZ+1) = NACCE

C----              Lift the call banks

      DO 24 J=1,JQMLEV
      CALL MZLIFT (IXDVJZ,L,LQJZ,-J-5,MMCALL,0)
   24 CONTINUE

C----              Collect all JQAN titles

      CALL TZINQ (IXSTJZ,IXDVTT,LQAN,1)
      IF (LQAN.EQ.0)               GO TO 49
      CALL MZXREF (IXDVJZ,IXDVTT,'A')
      GO TO 33

   32 LQAN = LQ(KQS+LQAN)
   33 IF (LQAN.EQ.0)               GO TO 39
      IF (IQ(KQS+LQAN-4).NE.MJQAN)    GO TO 32

      L = LZFIND (IXSTJZ,LQ(KQS+LQJZ-1), IQ(KQS+LQAN+1),1)
      IF (L.NE.0)                  GO TO 34
      L    = LQAN
      LQAN = LQ(KQS+LQAN)
      CALL ZSHUNT (IXSTJZ,L,LQJZ,-1,0)
      GO TO 33

   34 CALL MZDROP (IXSTJZ,LQAN, '.')
      GO TO 32

   39 CALL ZTOPSY (IXSTJZ,LQ(KQS+LQJZ-1))

C----              Digest all JQFL titles

      LQAN = LQT(KQT+1)
      GO TO 43

   42 LQAN = LQ(KQS+LQAN)
   43 IF (LQAN.EQ.0)               GO TO 49
      IF (IQ(KQS+LQAN-4).NE.MJQFL)    GO TO 42

#if defined(CERNLIB_QDEBUG)
      JNEXT = 1
      JMAX  = IQ(KQS+LQAN-1)

   44 IF (JNEXT.GT.JMAX)           GO TO 47
      JID   = JNEXT
      JEND  = IUFIND (MEND,IQ(KQS+LQAN+1),JID+1,JMAX)
      JNEXT = JEND + 1
      ID = IQ(KQS+LQAN+JID)
      IF (LZFIND(IXSTJZ,LQ(KQS+LQJZ-2),ID,1).NE.0)   GO TO 44

      IF (JID .NE.1)               GO TO 45
      IF (JEND.LE.JMAX)            GO TO 45
      L    = LQAN
      LQAN = LQ(KQS+LQAN)
      CALL ZSHUNT (IXSTJZ,L, LQJZ,-2,0)
      GO TO 33

   45 MMFL(4) = MIN (JEND-JID,JQMFLW+1)
      CALL MZLIFT (IXDVTT,L,LQJZ,-2,MMFL,-1)
      CALL UCOPY (IQ(KQS+LQAN+JID),IQ(KQS+L+1),MMFL(4))
      GO TO 44

#endif
   47 CALL MZDROP (IXSTJZ,LQAN, '.')
      GO TO 42

   49 CONTINUE

C----              Lift the root level

#if defined(CERNLIB_JZTIME)
#include "zebra/jztimed.inc"
#endif

      JQREOD = 512
      JQEALL = -JQREOD

      CALL JZIN (CHIAM,1,NQME(2),MEXTR)

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
