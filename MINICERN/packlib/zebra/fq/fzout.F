*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:51  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOUT (LUNP,IXDIVP,LENTP,IEVP,CHOPT,NIOP,NUHP,IUHEAD)

C-    MAIN SEQUENTIAL OUTPUT ROUTINE, USER CALLED

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzct.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9),IXDIVP(9),LENTP(9),IEVP(9)
      DIMENSION    NIOP(9),NUHP(9),IUHEAD(99)
      CHARACTER    CHOPT*(*)
      DIMENSION    IUHDAT(3)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOU, 4HT    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOUT  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOUT   ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      LUNNX   = LUNP(1)
      IXDIVX  = IXDIVP(1)
      LENTRX  = LENTP(1)
      IEVFLX  = IEVP(1)
      IOCHX(1)= NIOP(1)
      NWUHOR  = MAX (NUHP(1),0)
      NWUHX   = MIN (NWUHOR,400)
      NWFILX  = 0
      NWMEMT  = 0
      ICOPYX  = 0

C--                Set current output unit

      IF (LUNNX.NE.LUNX)  CALL FZLOC (LUNNX,2)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.3)
     +   WRITE (IQLOG,9111) LUNNX,LENTRX,IEVFLX,IACTVX,CHOPT
 9111 FORMAT (1X/' FZOUT-   Enter for LUN=',I3,
     F' LENTRY,IEVFL,IACTV,OPT=',I8,2I3,1X,A)

#endif
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif
#if defined(CERNLIB_QPRINT)
      IF (NWUHOR.GT.NWUHX)  THEN
          IF (LOGLVX.GE.-2)  WRITE (IQLOG,9112) LUNX,NWUHOR
        ENDIF
 9112 FORMAT (1X/' FZOUT.   LUN=',I4,' Of ',I4,' user header words',
     F' only 400 are taken !!!')
#endif
#if defined(CERNLIB_FZCHANNEL)
      IF (IACMOX.EQ.3)  THEN
          IF (IADOPX.EQ.0)         GO TO 907
        ENDIF

#endif
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.EQ.3)  THEN
          IADOPX = IQ(KQSP+LQFX+8)
          IF (IADOPX.EQ.0)         GO TO 907
          IQ(KQSP+LQFX+1) = IADOPX
        ENDIF

#endif
C-      from FZRUN to 401
      IF (IEVFLX.LT.0)             GO TO 401
C-      from FZENDO to 411
      IF (IEVFLX.GE.2)             GO TO 411
      IF (IACTVX.GE.16)            GO TO 901

   21 NWTXX   = 0
      NWSEGX  = 0
      NWTABX  = 0
      NWBKX   = 0
      NWUHCX  = 0
      NWIOX   = 0

C--                Ready I/O characteristic

      IF (NWUHX.EQ.0)              GO TO 39

      IF   (IOCHX(1))        34, 32, 33
   32 IOCHX(1) = 3
   33 NWIOX = 1
      IF (IOCHX(1).LT.8)           GO TO 38
   34 NWIOX = JBYT (IOCHX(1), 7,5)
      J     = JBYT (IOCHX(1),12,5)
      IF     (JBYT (IOCHX(1), 1,6).NE.1)  GO TO 903
      IF (NWIOX.GT.16)             GO TO 903
      IF (NWIOX.NE.J+1)            GO TO 903

      IF (NWIOX.GT.1)  CALL UCOPY (NIOP,IOCHX,NWIOX)
   38 NWUHCX = NWUHX + NWIOX
   39 CONTINUE

C----              Ready text vector

      LTEXTX = LQ(KQSP+LQFX-2)
      IF (LTEXTX.NE.0)  CALL FZOTXT

C----              Construct table of material to be written

      CALL UOPTC (CHOPT,'DISZMLNP',IOPTXD)
      IOPTXN = IOPTXN + IOPTXS + IOPTXZ

      IF (IOPTXZ.NE.0)             GO TO 121

      CALL MZSDIV (IXDIVX,0)
      JQSTMV = -1
      MODTBX = 0
      JFLGAX = 0
      CALL FZOTAB
      IF (IQUEST(1).NE.0)          GO TO 999
      IF (IOPTXN.NE.0)  NWTABX=0

C------            Output of pilot records

  121 IDX(2) = 3
      IF (IEVFLX.NE.0)  IDX(2)=2
      NWMEMT = 20 + NWUHCX + NWSEGX + NWTXX + NWTABX + NWBKX

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.EQ.2)
     +WRITE (IQLOG,9121) LUNX,JQSTOR,LENTRX,NWTABX,NWBKX,IEVFLX,CHOPT
 9121 FORMAT (' FZOUT-   LUN=',I4,' Store/LENTRY=',I3,I9,
     F' NWTAB,NWBANK=',I5,I7,' Evfl/Opt=',I3,1X,A)

#endif
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (IUHEAD)
          GO TO 124
        ENDIF
#endif
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.EQ.3)  THEN
          IF (NWMEMT.GT.IQ(KQSP+LQFX+9))      GO TO 909
        ENDIF
#endif
      CALL FZOFFX (IUHEAD)
  124 IQ(KQSP+LQFX+15) = IQ(KQSP+LQFX+15) + 1

      IACTVX  = 12
      IF (NWBKX.EQ.0)              GO TO 190

C------            Write out material according to table

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9140) NWTABX,NWBKX
 9140 FORMAT (10X,'NWTAB/NWBK =',2I9)
#endif
      IQ(KQSP+LQFX+16) = IQ(KQSP+LQFX+16) + 1

      IDX(2) = 0
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (0)
          GO TO 191
        ENDIF
#endif
      CALL FZOFFX (0)
      GO TO 191

C----              Test for pseudo end-of-tape

  190 IQ(KQSP+LQFX+17) = IQ(KQSP+LQFX+17) + 1
  191 NUM1 = IQ(KQSP+LQFX+19)
      NUM2 = IQ(KQSP+LQFX+20)
  192 IF (NUM2.GE.1000000)  THEN
          NUM1 = NUM1 + 1
          NUM2 = NUM2 - 1000000
          IQ(KQSP+LQFX+19) = NUM1
          IQ(KQSP+LQFX+20) = NUM2
          GO TO 192
        ENDIF

      LIM1 = IQ(KQSP+LQFX+37)
      LIM2 = IQ(KQSP+LQFX+38)

      IF (LIM1+LIM2.EQ.0)          GO TO 991
      IF (NUM1-LIM1)         991, 196, 197
  196 IF (NUM2.LT.LIM2)            GO TO 991

  197 IQUEST(1) = 1
      GO TO 992

C-------------------------------------------------
C-                 Write start-of-run
C-                 Write end-of-run, end-of-file
C-------------------------------------------------

C-      from FZRUN to 401
  401 JRUN   = LENTRX
      JRUNCR = IQ(KQSP+LQFX+29)
      IF (JRUN.LT.0)               GO TO 410
      IF (IACTVX.LT.11)            GO TO 406
      IF (IACTVX.GE.14)            GO TO 406
      JTRA  = 1
      NWUHU = 0
C-      precede start-of-run by end-of-run
      GO TO 422

C--                Write start-of-run

  406 IACTVX = 11
      JRUNX  = 1
      JRUN   = LENTRX
      IF (JRUN.EQ.0)  JRUN=JRUNCR+1
      IQ(KQSP+LQFX+29) = JRUN
      IQ(KQSP+LQFX+14) = IQ(KQSP+LQFX+14) + 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)  WRITE (IQLOG,9406) LUNX,JRUN
 9406 FORMAT (1X/' FZOUT.   LUN=',I4,' Write Start-of-Run ',I6)
#endif
      JTRA  = 2
      NWUHU = NWUHX
      GO TO 424

C-  IACTV : last action            IEVFL : action requested
C-
C-     11   StoR written
C-     12   d/s written
C-     13   buffer flushed            13   flush buffer
C-     14   EOR written               14   write EOR
C-     15   first  EOF written        15   write EoF
C-     16   second EOF written        16   write EoD
C-     17   attempted write after eod

  410 IEVFLX = 14
C-      from FZENDO to 411
  411 IF (IEVFLX.LT.13)            GO TO 902
      IF (IEVFLX.GE.17)            GO TO 902
      IF (IACTVX.LT.11)            GO TO 991
      IF (IEVFLX.LE.IACTVX)        GO TO 991

      IF (IEVFLX.GE.14)            GO TO 421
      IACTVX = 13
      JTRA   = 3
      GO TO 441

C----              Write end-of-run

  421 IF (IACTVX.GE.14)            GO TO 461
      JRUNCR = IQ(KQSP+LQFX+29)
      JTRA   = 4
      NWUHU  = NWUHX

  422 IQ(KQSP+LQFX+13) = IQ(KQSP+LQFX+13) + 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)  WRITE (IQLOG,9421) LUNX,JRUNCR
 9421 FORMAT (' FZOUT.   LUN=',I4,' Write End-of-Run ',I8)
#endif
C--       flush the buffer before EoR
      IF   (IACTVX.LT.13)  THEN
        IF (IFIFOX.NE.0)   THEN
            IDX(2) = -1
            CALL FZOFFX (IUHEAD)
          ENDIF
        ENDIF
      IACTVX = 14
      JRUNX  = 0
      JRUN   = 0

C--                Write StoR/EoR record

  424 IDX(1) = 1 + NWUHU
      IDX(2) = 1
      IQUEST(11) = JRUN
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (IUHEAD)
          GO TO 426
        ENDIF
#endif
      CALL FZOFFX (IUHEAD)
  426 CONTINUE

      IF (JBIT(MSTATX,15).NE.0)    GO TO 441
      IF (JTRA.EQ.2)               GO TO 991

C--                Flush the buffer after operation

  441 IF (IFIFOX.NE.0)  THEN
          IDX(2) = -1
          CALL FZOFFX (0)
        ENDIF

      GO TO ( 406, 991, 991, 461, 471), JTRA

C----              Write end-of-file / end-of-data

  461 IF (IEVFLX.LT.15)            GO TO 991
      IF (IACTVX.EQ.15)            GO TO 471
      NWUHU = 0
      JTRA  = 5
      IQ(KQSP+LQFX+12) = IQ(KQSP+LQFX+12) + 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)  WRITE (IQLOG,9464) LUNX
 9464 FORMAT (' FZOUT.   LUN=',I4,' Write Zebra EoF')
#endif

      JRUNX  = 0
      IDX(1) = 4
      IDX(2) = 1
      IUHDAT(1)  = IQ(KQSP+LQFX+34)
      IUHDAT(2)  = IQ(KQSP+LQFX+35)
      IUHDAT(3)  = IQ(KQSP+LQFX+33) + 1
      IQUEST(11) = -1
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (IUHDAT)
          GO TO 426
        ENDIF
#endif
      CALL FZOFFX (IUHDAT)
      GO TO 426

C--                Write true EOF

  471 NEOF  = JBYT (MSTATX,13,2)
      IF (IEVFLX.EQ.15)  THEN
          NEOFU = MOD (NEOF,2)
        ELSE
          NEOFU = MAX (0,NEOF-1)
        ENDIF
      IF (NEOFU.EQ.0)              GO TO 474

      IQUEST(11) = NEOFU
      IQUEST(12) = NEOF
      IDX(2)     = 13 - IEVFLX
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (0)
          GO TO 473
        ENDIF
#endif
      CALL FZOFFX (0)
  473 IQ(KQSP+LQFX+11) = IQ(KQSP+LQFX+11) + IQUEST(11)

  474 IACTVX = IEVFLX
#if defined(CERNLIB_QPRINT)
      IF   (LOGLVX.GE.0)  THEN
         IF (IEVFLX.NE.16)  THEN
            WRITE (IQLOG,9474) LUNX
          ELSE
            WRITE (IQLOG,9475) LUNX
          ENDIF
        ENDIF
 9474 FORMAT (' FZOUT.   LUN=',I4,' End-of-File')
 9475 FORMAT (' FZOUT.   LUN=',I4,' End-of-Data')
#endif
      GO TO 991

C-------------------------------------------------
C-                 ERROR HANDLING
C-------------------------------------------------

  901 IF (IACTVX.EQ.18)            GO TO 21
      IF (IACTVX.EQ.17)   CALL ZFATAM ('FZOUT - Going beyond EoD.')
      IACTVX = 17
      IQUEST(1) = -1
      GO TO 992

  902 IQUEST(11) = IEVFLX
      CALL ZFATAM ('FZOUT - Faulty parameter IEVENT.')

  903 IQUEST(11) = IOCHX(1)
      CALL ZFATAM ('FZOUT - IOCH invalid.')

#if defined(CERNLIB_FZCHANNEL)||defined(CERNLIB_FZMEMORY)
  907 IQUEST(1)  = -2
      IQUEST(11) = LUNX
      IF (IEVFLX.GE.2)             GO TO 992
      CALL ZFATAM ('FZOUT - User routine / buffer not connected.')
#endif
#if defined(CERNLIB_FZMEMORY)
  909 IQUEST(1)  = -2
      IQUEST(2)  = 14
      IQUEST(8)  = IQ(KQSP+LQFX+9)
      IQUEST(9)  = NWMEMT
      IF (IOPTXP.EQ.0)  CALL ZTELL (14,1)
      GO TO 999
#endif
C-------------------------------------------------
C-                 COMMON EXIT
C-------------------------------------------------

  991 IQUEST(1)  = 0
  992 IQUEST(2)  = 0
      IQUEST(5)  = IQ(KQSP+LQFX+31)
      IQUEST(6)  = IQ(KQSP+LQFX+32)
      IQ(KQSP+LQFX+2)  = IACTVX
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.EQ.3)
     +IQUEST(9)  = IQ(KQSP+LQFX+1) - IQ(KQSP+LQFX+8)
#endif
      IQUEST(10) = NWMEMT
      IQUEST(11) = NWBKX
      IQUEST(12) = NWTABX
      IQUEST(13) = IQ(KQSP+LQFX+15)
      IQUEST(14) = IQ(KQSP+LQFX+19)
      IQUEST(15) = IQ(KQSP+LQFX+20)
      IQUEST(16) = IQ(KQSP+LQFX+21)
      IQUEST(17) = IQ(KQSP+LQFX+22)
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
