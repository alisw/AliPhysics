*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:47  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZFFNAT)
      SUBROUTINE FZOFFN (IUHEAD)

C-    Write operations for file format native,
C-    subsidiary to FZOUT

C-    Controlling parameter : IDX(2)
C-
C-    IDX(2)   = 1   write start/end-of-run
C-             > 1   write pilot for d/s
C-             = 0   write bank material for d/s
C-            (= -1  flush the buffer in exchange mode)
C-             = -2  End-of-File
C-             = -3  End-of-Data

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcwk.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcseg.inc"
C--------------    End CDE                             --------------
      DIMENSION    IUHEAD(99)
      DIMENSION    LV(6), NV(6)
      EQUIVALENCE (L1,LV(1)), (L2,LV(2)), (L3,LV(3))
     +,           (L4,LV(4)), (L5,LV(5)), (L6,LV(6))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOF, 4HFN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOFFN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOFFN  ')
#endif


#include "zebra/qtrace.inc"

      IF   (IDX(2))          401, 201, 101
  101 IF (IDX(2).EQ.1)             GO TO 301

C-----------------------------------------------------------
C------            WRITE PILOT INFORMATION, STARTING LR AND D/S
C-----------------------------------------------------------

      NWMAX  = MIN (MAXREX,1020)
      NWPILA = 10 + NWIOX
      NW     = NWPILA
      JSEND  = 0
      JDONE  = -1

C--                Do user header

      NW  = NW + NWUHX

C--                Do segment table

      IF (NWSEGX.NE.0)  THEN
          NW    = NW + NWSEGX
          JSEND = 2
          LV(1) = LOCF (IQSEGH(1,1)) - LQASTO
          LV(2) = LV(1) + 40
          NV(1) = 2*NQSEG
          NV(2) =   NQSEG
        ENDIF

C--                Do text vector

      IF (NWTXX.NE.0)  THEN
          LTX   = LTEXTX
          NW    = NW + NWTXX
          IF (NW.GT.NWMAX)         GO TO 127
          JSEND = JSEND + 1
          LV(JSEND) = KQSP+8 + LTX + 5
          NV(JSEND) = NWTXX
        ENDIF
      JDONE = 0

C--                Do early table

      IF (NWTABX.GE.41)            GO TO 127
      NW = NW + NWTABX
      IF (NW.GT.NWMAX)             GO TO 127

      IF (NWTABX.NE.0)  THEN
          JSEND = JSEND + 1
          LV(JSEND) = LQTA
          NV(JSEND) = NWTABX
        ENDIF
      JDONE = 1

C--                Transmit pilot record

  127 N  = NWPILA
      NU = NWUHX
      JSEND = JSEND + 1
      IF (NU.NE.0)  GO TO ( 140, 141, 142, 143, 144), JSEND
                    GO TO ( 130, 131, 132, 133, 134), JSEND

  130 CALL FZON1 (IPILX,N)
      GO TO 149

  131 CALL FZON2 (IPILX,N,LQ(L1),NV(1))
      GO TO 149

  132 CALL FZON3 (IPILX,N,LQ(L1),NV(1),LQ(L2),NV(2))
      GO TO 149

  133 CALL FZON4 (IPILX,N,LQ(L1),NV(1),LQ(L2),NV(2),LQ(L3),NV(3))
      GO TO 149

  134 CALL FZON5 (IPILX,N,LQ(L1),NV(1),LQ(L2),NV(2),LQ(L3),NV(3)
     +,                   LQ(L4),NV(4))
      GO TO 149

  140 CALL FZON2 (IPILX,N,IUHEAD,NU)
      GO TO 149

  141 CALL FZON3 (IPILX,N,IUHEAD,NU,LQ(L1),NV(1))
      GO TO 149

  142 CALL FZON4 (IPILX,N,IUHEAD,NU,LQ(L1),NV(1),LQ(L2),NV(2))
      GO TO 149

  143 CALL FZON5 (IPILX,N,IUHEAD,NU,LQ(L1),NV(1),LQ(L2),NV(2)
     +,                   LQ(L3),NV(3))
      GO TO 149

  144 CALL FZON6 (IPILX,N,IUHEAD,NU,LQ(L1),NV(1),LQ(L2),NV(2)
     +,                   LQ(L3),NV(3),LQ(L4),NV(4))

  149 IDX(2) = 4
      IQUEST(7) = JDONE
      IF   (JDONE)           171, 181, 999

C----              Pilot continuation: text/table

  171 IF (NWTABX.GE.41)            GO TO 174
      IF (NWTABX.EQ.0)             GO TO 174
      IF (NWTXX+NWTABX.GT.MAXREX)  GO TO 174
      CALL FZON2 (IQ(KQSP+LTX+5),NWTXX,LQ(LQTA),NWTABX)
      IQUEST(7) = 1
      GO TO 999

  174 CALL FZON1 (IQ(KQSP+LTX+5),NWTXX)

C--                Pilot continuation: table only  (only FZOUT)

  181 IF (NWTABX.EQ.0)             GO TO 999
      IF (ICOPYX.NE.0)             GO TO 999
      NT = NWTABX
      L  = LQTA
  182 N  = MIN (NT,MAXREX)
      CALL FZON1 (LQ(L),N)
      L  = L  + N
      NT = NT - N
      IF (NT.NE.0)                 GO TO 182

#include "zebra/qtrace99.inc"
      RETURN

C-----------------------------------------------------------
C--                WRITE BANK MATERIAL
C-----------------------------------------------------------

  201 MINREC = MAXREX/2
      LTEMPX = 0
      LTB    = LQTA
      JSEG   = 0
      NDOSG  = 0
      IF (NQSEG.EQ.0)  NDOSG=NWBKX
      IDX(2) = 7
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVX.GE.5)  WRITE (IQLOG,9801)
 9801 FORMAT (' FZOFFN-  Entered for bank material.')
#endif

  242 JSEND = 0
      NWS   = 0

C--                Load next sector

  243 L     = LQ(LTB)
      LE    = LQ(LTB+1)
      N     = LE - L
      JSEND = JSEND + 1
      LV(JSEND) = L
      NV(JSEND) = N

      NWS = NWS + N
      NOV = NWS - MAXREX
      LTB = LTB + 2

C--                Next segment ?

      IF (NDOSG.EQ.0)  THEN
          JSEG  = JSEG + 1
          NDOSG = IQSEGD(JSEG)
        ENDIF
      NDOSG = NDOSG - N

C--                Send ?

      IF (NOV.GT.0)                GO TO 261
      IF (LTB.EQ.LQTE)             GO TO 268
      IF (NDOSG.EQ.0)              GO TO 270
      IF (NOV.GE.-10)              GO TO 270
      IF (JSEND.LT.6)              GO TO 243
      IF (NWS.GE.MINREC)           GO TO 270

C--                6 sectors have less than MINREC words
C--                Compact to TEMP buffer

      IF (LTEMPX.EQ.0)  LTEMPX = LQWKFZ - KQS

      LOV = LTEMPX
      NOV = 0
      JGO = 1
      IF (LV(1).EQ.LTEMPX)  THEN
          NOV = NV(1)
          JGO = 2
        ENDIF

      DO 256  J=JGO,6
      L = LV(J)
      N = NV(J)
      CALL UCOPY (LQ(KQS+L),LQ(KQS+LOV+NOV),N)
  256 NOV = NOV + N
      LV(1) = LOV
      NV(1) = NOV

  257 JSEND = 1
      NWS   = NV(1)
      GO TO 243

C--                Last sector overflows MAXREX words

  261 N   = N - NOV
      LOV = L + N
      NV(JSEND) = N
      GO TO 270

C--                End of material, with overflow on last sector

  267 JSEND = 1

C--                End of material reached, send last record

  268 IDX(2) = 8

C------            Write 1 record

  270 GO TO ( 271, 272, 273, 274, 275, 276), JSEND

  271 CALL FZON1 (LQ(KQS+L1),NV(1))
      GO TO 279
  272 CALL FZON2 (LQ(KQS+L1),NV(1),LQ(KQS+L2),NV(2))
      GO TO 279
  273 CALL FZON3 (LQ(KQS+L1),NV(1),LQ(KQS+L2),NV(2),LQ(KQS+L3),NV(3))
      GO TO 279
  274 CALL FZON4 (LQ(KQS+L1),NV(1),LQ(KQS+L2),NV(2),LQ(KQS+L3),NV(3),
     +            LQ(KQS+L4),NV(4))
      GO TO 279
  275 CALL FZON5 (LQ(KQS+L1),NV(1),LQ(KQS+L2),NV(2),LQ(KQS+L3),NV(3),
     +            LQ(KQS+L4),NV(4),LQ(KQS+L5),NV(5))
      GO TO 279
  276 CALL FZON6 (LQ(KQS+L1),NV(1),LQ(KQS+L2),NV(2),LQ(KQS+L3),NV(3),
     +            LQ(KQS+L4),NV(4),LQ(KQS+L5),NV(5),LQ(KQS+L6),NV(6))

C----              Overflow material pending ?

  279 IF (IDX(2).EQ.8)             GO TO 999
      IF (NOV.LE.0)                GO TO 242
      LV(1) = LOV
      NV(1) = NOV
      NOV   = NOV - MAXREX
      IF (NOV.GT.0)                GO TO 282

C--                End of all material ?

      IF (LTB.EQ.LQTE)             GO TO 267
      IF (NOV.EQ.0)                GO TO 271

C--                End of segment ?

      IF (NDOSG.NE.0)              GO TO 257
      GO TO 271

C--                Overflow on overflow

  282 NV(1) = MAXREX
      LOV   = LOV + MAXREX
      GO TO 271

C-----------------------------------------------------------
C--                WRITE START/END-OF-RUN
C-----------------------------------------------------------

  301 JRUN   = IQUEST(11)
      NWUHU  = IDX(1) - 1
      IF (NWUHU.NE.0)  THEN
          CALL FZON2 (JRUN,1,IUHEAD,NWUHU)
        ELSE
          CALL FZON1 (JRUN,1)
        ENDIF
      GO TO 999

C-----------------------------------------------------------
C--                ENDFILE
C-----------------------------------------------------------

C-        NEOF = 1  EoF 1 only       IDX(1) = -2  EoF
C-               2  EOF 2 only                -3  EoD
C-               3  EOF 1 + 2

  401 NEOFU = IQUEST(11)
      NEOF  = IQUEST(12)
      IQUEST(11) = NEOFU
      NEOF = NEOFU
  412 ENDFILE LUNX
      NEOF = NEOF - 1
      IF (NEOF.GT.0)               GO TO 412
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)   WRITE (IQLOG,9414) LUNX,IQUEST(11)
 9414 FORMAT (' FZOFFN.  LUN=',I4,' Write',I2,' System EOF')
#endif
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
