      SUBROUTINE ZSORVH (IXSTOR,LGOP,JWORD,NWORDS)

C-    SORT BANKS AT LGO SUCH THAT THE 'NWORDS' LONG KEY STRINGS
C-    STARTING AT Q(L+JWORD) ARE IN INCREASING ORDER

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    JWORD(9), NWORDS(9), LGOP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZSOR, 4HVH   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZSORVH /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZSORVH  ')
#endif

#include "zebra/q_jbyt.inc"


      LGO = LGOP(1)
      IF (LGO.EQ.0)          RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstore.inc"

      JW  = JWORD(1)
      NCH = 4*NWORDS(1)
      KGO = LQSTA(KQT+2) - 1
      LQ(KQS+KGO) = LGO
      LL  = LGO
      IFL = 0

   11 LN = LQ(KQS+LL)
      IF (LN.EQ.0)                 GO TO 81
      IF (IUCOMH(IQ(KQS+LN+JW),IQ(KQS+LL+JW),NCH).LT.0)   GO TO 21
      LL = LN
      GO TO 11

C--                BANK  LN  OUT OF SEQUENCE

   21 LQ(KQS+LL)   = LQ(KQS+LN)
      IFL = 7
      K   = KGO

   24 L = LQ(KQS+K)
      IF (IUCOMH(IQ(KQS+LN+JW),IQ(KQS+L+JW),NCH).LT.0)   GO TO 29
      K = L
      GO TO 24

C--                PLACE FOR BANK  LN  FOUND

   29 LQ(KQS+LN)  = L
      LQ(KQS+K)   = LN
      GO TO 11

C----              FINISHED, CHAIN K-LINKS

   81 IF (IFL.EQ.0)                GO TO 999
      K = LQ(KQS+LGO+2)
      L = LQ(KQS+KGO)
      LQ(KQS+L+2) = K
      IF (K.NE.0)  LQ(KQS+K)=L
      LGOP(1) = L

   84 K = L
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 999
      LQ(KQS+L+2) = K
      K = L
      GO TO 84

#include "zebra/qtrace99.inc"
      RETURN
      END
