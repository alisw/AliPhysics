*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:45  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZOASC (LUNP,LBUFAP,LBUFEP,JFAST,LWORKA,LWORKE)

C-    Dump one physical record to cards in ALFA mode,
C-    service routine to FZOUT.

C-    The record is stored in words LQ(LBUFA) to LQ(LBUFE-1)
C-    It is a 'fast' record if JFAST .NE. 0
C-    The working space available to FZOASC
C-       is  LQ(LWORKA) to LQ(LWORKE)  if LWORKE .NE. 0
C-       or  LQ(LWORKA) to LQ(LBUFC)   if LWORKE .EQ. 0

#include "zebra/zbcd.inc"
#include "zebra/zkrakc.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), LBUFAP(9), LBUFEP(9)

      DIMENSION MTYPTR(7), ICHSUM(2)
      EQUIVALENCE  (LUN ,IQUEST(90))
      EQUIVALENCE  (JPUT,IQUEST(91)), (IFLEND,IQUEST(92))
      EQUIVALENCE (LUPKA,IQUEST(93)), (LUPKE, IQUEST(94))

#if defined(CERNLIB_QMVDS)
      SAVE         MTYPTR
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOA, 4HSC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOASC  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOASC  ')
#endif
      DATA  MTYPTR / 0, 0, 5, 10, 15, 20, 36 /

#include "zebra/q_jbyt.inc"

#include "zebra/qtraceq.inc"

      LUN   = LUNP(1)
      LBUFA = LBUFAP(1)
      LBUFE = LBUFEP(1)
      LWKA  = LWORKA
      LWKE  = LWORKE

      LUPKA = LWKA
      LUPKE = LWKE
      NWWKU = LWKE - LWKA
      IF (LUPKE.EQ.0)  THEN
          LUPKE = LBUFA
          NWWKU = 0
        ENDIF

      IQCETK(1) = IQGREA
      IQCETK(2) = 58
      IF (JFAST.EQ.0)  THEN
          IQCETK(3) = 27
          LBUFC = LBUFA + 4
        ELSE
          IQCETK(3) = 26
          LBUFC = LBUFA
        ENDIF
      JPUT   = 3
      IFLEND = 0
      ITYPE  = 0
      ITYPC  = -1

      ICHSUM(1) = 0
      ICHSUM(2) = 0
      GO TO 40

C--------          Do next word

   31 IWORD = LQ(LBUFC)
   32 JTYPE = LQ(LUPK)
      IF (JTYPE.LT.0)              GO TO 61
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.9)  THEN
          WRITE (IQLOG,9832) LBUFC+1-LBUFA,IWORD,(LQ(LUPK+J),J=0,9)
          DO  33  J=0,7
          IF (IABS(LQ(LUPK+J)).GT.99)  WRITE (IQLOG,9833) LQ(LUPK+J)
   33     CONTINUE
        ENDIF

#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_B32))
 9832 FORMAT (' DEVZE FZOASC, do word',I4,Z9,' Upk=',8I3,' Next',2I5)
 9833 FORMAT (22X,Z9)
#endif
#if (defined(CERNLIB_QDEVZE))&&(!defined(CERNLIB_B32))
 9832 FORMAT (' DEVZE FZOASC, do word',I4,Z17,' Upk=',8I3,' Next',2I5)
 9833 FORMAT (22X,Z17)

#endif
      JTYPS = LQ(LUPK+1)
      JTYPC = MTYPTR(JTYPE+1) + JTYPS

C--       same type running
      IF (JTYPC.EQ.ITYPC)          GO TO 36

C--                M-type 0 : integers 0 -> 9

      IF (JTYPE.NE.0)              GO TO 41
      IF (ITYPC.EQ.24)             GO TO 36
      JBA = 0

C------            Store next word

C--       terminate running type by ']', if any
   34 IF (ITYPE.NE.0) THEN
          JPUT = JPUT + 1
          IQCETK(JPUT) = 53
          ITYPE = 0
          ITYPC = -1
        ENDIF

C--       place type-code
   35 JPUT = JPUT + 1
      IQCETK(JPUT) = JTYPC
      IF (JBA.EQ.0)                GO TO 38

C--       place significant bytes
   36 DO 37  J=JBA,JBE
      JPUT = JPUT + 1
   37 IQCETK(JPUT) = LQ(LUPK+J)
   38 LUPK = LUPK + 8

C--       put blank separator for tests
#if defined(CERNLIB_FQABLANK)
      JPUT = JPUT + 1
      IQCETK(JPUT) = 44
#endif

C--                Cumulate check-sum

      ICHSUM(1) = ICHSUM(1) + JBYT (IWORD,17,16)
      ICHSUM(2) = ICHSUM(2) + JBYT (IWORD, 1,16)

C--                Write line if full

      IF (JPUT.GE.80)   CALL FZOALN

C--                End-of-record ?

      LBUFC = LBUFC + 1
      IF (LBUFC.EQ.LBUFE)          GO TO 81

C--                Unpack the next lot, if necc.

      IF (LUPK.LT.LUPKE)           GO TO 31

      IF (NWWKU.NE.0)  THEN
          IF (NWWKU.GE.LBUFC-LBUFA) THEN
              LUPKE = LWKE
            ELSE
              LUPKA = LBUFA
              LUPKE = LBUFC
            ENDIF
        ELSE
          LUPKE = LBUFC
        ENDIF

   40 LUPKEL = LUPKE
      CALL FZOAPK (LBUFC,LBUFE)
      LUPK = LUPKA
      GO TO 31

C------            Type analysis
C-                 come to here if the type of the current word
C-                 is not the one of the running set, and if
C-                 the current M-type is not zero


   41 NTYPE = LQ(LUPK+8)
      NTYPS = LQ(LUPK+9)
      NTYPC = -2
      IF (NTYPE.GE.0)  NTYPC = MTYPTR(NTYPE+1) + NTYPS

C--                Check next word after current has again ITYPC
C-                 and ITYPC covers JTYPC with 1 or 2 zero bytes

      IF (JTYPE.NE.ITYPE)          GO TO 51
      IF (NTYPC.NE.ITYPC)          GO TO 51
      N = JTYPC - ITYPC
      IF (N.LT.0)                  GO TO 51
      IF (N.LT.3)                  GO TO 36

C----              The running set is definitely not continued

   51 IF (JTYPE.LT.5)  THEN
          JBA = 2
          JBE = 7 - JTYPS
        ELSE
          JBA = 3 + JTYPS
          JBE = 7
        ENDIF

C----              Start of a new set if this and the next
C-                 two numbers are of the same type

      IF (NTYPC.NE.JTYPC)          GO TO 34
      IF (LUPK+16.GE.LUPKE)        GO TO 34
      IF (LQ(LUPK+16).NE.JTYPE)    GO TO 34
      IF (LQ(LUPK+17).NE.JTYPS)    GO TO 34

      ITYPE = JTYPE
      ITYPC = JTYPC
      JPUT  = JPUT + 1
      IQCETK(JPUT) = 57
      GO TO 35

C------            Repeat last number N+1 times

   61 CONTINUE
#if defined(CERNLIB_QDEVZE)
      IF (JTYPE.NE.-43)  CALL ZFATAM ('FZOASC - trouble.')
#endif
      JPUT = JPUT + 1
      IQCETK(JPUT) = 43

      LUPK  = LUPK + 1
      IWORD = LQ(LUPK)
      LBUFC = LBUFC + IWORD
#if defined(CERNLIB_QDEVZE)
      LUPKEL = LUPKEL - 4
      LUPK   = LUPKEL
      LQ(LUPK) = IWORD
#endif

      LUPKE = LUPKEL
      CALL FZOAPK (LUPK,LUPK+1)
   68 LUPK  = LUPKA
      ITYPE = 0
      ITYPC = -1
      GO TO 32

C--------          End of record

   81 IF (IFLEND.NE.0)             GO TO 84
      JPUT = JPUT + 1
      IQCETK(JPUT) = 59
      LST   = LUPKA
      LUPKA = LUPKA + 2
      LUPKE = LUPKEL

      LQ(LST)   = ICHSUM(1)
      LQ(LST+1) = ICHSUM(2)
      CALL FZOAPK (LST,LST+2)

      IFLEND = -1
      LBUFC  = LBUFE - 1
      IWORD  = ICHSUM(1)
      GO TO 68

   84 IF (IFLEND.NE.-1)            GO TO 86
      IFLEND = -2
      LBUFC  = LBUFE - 1
      IWORD  = LQ(LST+1)
      GO TO 32

   86 IFLEND = 7
      CALL FZOALN

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
