*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:34  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIMTB

C-    ready memory occupation table for input
C-    ready space for the relocation table
C-    called from FZIN

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcseg.inc"
#include "zebra/fzcocc.inc"
C--------------    End CDE                             --------------

      DIMENSION    ITOSOR(20), ISORDV(20), ISORSP(20)
      DIMENSION    LSTAV(20),             LENDV(20)
      EQUIVALENCE (LSTAV(1),IQUEST(60)), (LENDV(1),IQUEST(80))

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIM, 4HTB   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIMTB /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIMTB  ')
#endif

#include "zebra/q_shiftl.inc"
#include "zebra/q_jbyt.inc"

C----              Input information

C-    if NQSEG = 0 :   read non-segmented

C-    if NQSEG > 0 :   read segmented according to :
C-       IQSEGD        segment selection table set by the user
C-       IQ(LFISEG+1)  copy of the segmentation table from the pilot

C----              Output information

C-    The memory occupation table at LQ(LQMTA) ready for the material
C-    to be accepted :
C-        skip :   LQ(LMT)   = 0
C-                 LQ(LMT+1) = 0
C-                 LQ(LMT+3) = -NW words to be skipped

C-        read :   LQ(LMT)   = JDIV
C-                 LQ(LMT+1) = 1
C-                 LQ(LMT+3) = LSTA
C-                 LQ(LMT+4) = LEND

C-    The space at LQ(LQTA) ready (and big enough) to receive
C-    the relocation table.

C-    The common /FZOCC/ indicates the amount of space IQOCSP(J)
C-    reserved at the end/start of division IQOCDV(J) for J=1,NQOCC
C-    and blocked by a dummy bank.


#include "zebra/qtrace.inc"

      IFLGAR = 0
      IF (NQSEG.LE.0)  THEN

C--                Single segment

          NQSEG = 1
          NSOR  = 1
          NOCC  = 1
          ITOSOR(1) = 1
          IQSEGD(1) = JQDIVI
          ISORDV(1) = JQDIVI
          IQOCDV(1) = JQDIVI
          ISORSP(1) = NWBKI
          IQOCSP(1) = NWBKI
          GO TO 41
        ENDIF

C--                Multiple segments, tidy division numbers

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)   WRITE (IQLOG,9016)  (J,
     +                   IQSEGH(1,J),IQSEGH(2,J),IQSEGD(J),J=1,NQSEG)
 9016 FORMAT (1X/' FZIMTB-  Segment Selection Table as set by the user'
#endif
#if (defined(CERNLIB_QDEBPRI))&&(!defined(CERNLIB_HEX))
     F/(10X,I2,1X,2A4,O23))
#endif
#if (defined(CERNLIB_QDEBPRI))&&(defined(CERNLIB_HEX))
     F/(10X,I2,1X,2A4,Z17))
#endif

      LFISEG = LQFI + JAUSEG
      IF (3*NQSEG.NE.IQ(KQSP+LFISEG))  GO TO 715
      LSPACE = KQSP + LFISEG + 2*NQSEG

      DO 27  JS=1,NQSEG
      IXDIV = IQSEGD(JS)
      IF   (IXDIV)           22, 23, 24
   22 IF (IXDIV.LT.-7)             GO TO 714
      ITOSOR(JS) = -IQ(LSPACE+JS)
      GO TO 27

   23 JDIV = JQDIVI
      GO TO 25

   24 JDIV = JBYT (IXDIV,1,26)
      IF (JDIV.GT.20)              GO TO 714
      JSTO = JBYT (IXDIV,27,4)
      IF   (JSTO.NE.0)  THEN
          IF (JSTO.NE.JQSTOR)      GO TO 714
        ENDIF

      IF (JDIV.EQ.0)               GO TO 23
      IF   (JDIV.GT.JQDVLL)  THEN
          IF (JDIV.LT.JQDVSY)      GO TO 714
        ENDIF

   25 IQSEGD(JS) = JDIV
      ITOSOR(JS) = 0
   27 CONTINUE

C----              Construct tables sorted by division number

C-             ISOR  table :  one entry per segment
C-             IQOC  table :  one entry per target division

      NSOR  = 0
      NOCC  = 0
      JANX  = 1
      JENX  = NQSEG

C--                Find largest division number not yet done

   31 JDVBIG = 0
      JA   = JANX
      JANX = 0
      JE   = JENX
      JENX = 0
      DO 35  JS=JA,JE
      IF (ITOSOR(JS).NE.0)         GO TO 35
      JENX = JS
      IF (JANX.EQ.0)  JANX=JS
      JDIV = IQSEGD(JS)
      IF (JDIV.LE.JDVBIG)          GO TO 35
      JDVBIG = JDIV
      JSBIG  = JS
   35 CONTINUE
      IF (JDVBIG.EQ.0)             GO TO 41

      NSOR = NSOR + 1
      ITOSOR(JSBIG) = NSOR
      ISORDV(NSOR)  = JDVBIG
      ISORSP(NSOR)  = IQ(LSPACE+JSBIG)

      NOCC = NOCC + 1
      IQOCDV(NOCC) = JDVBIG
      IQOCSP(NOCC) = IQ(LSPACE+JSBIG)

C--                Find further segments for the same division

      IF (JSBIG.EQ.JENX)           GO TO 31
      DO 37  JS=JSBIG+1,JENX
      IF (ITOSOR(JS).NE.0)         GO TO 37
      IF (IQSEGD(JS).NE.JDVBIG)    GO TO 37
      NSOR = NSOR + 1
      ITOSOR(JS)   = NSOR
      ISORDV(NSOR) = JDVBIG
      ISORSP(NSOR) = IQ(LSPACE+JS)
      IQOCSP(NOCC) = IQOCSP(NOCC) + IQ(LSPACE+JS)
   37 CONTINUE
      GO TO 31

C----              Reserve space in divisions

   41 IF (NOCC.EQ.0)               GO TO 81
      JOCC   = 0
   42 JOCC   = JOCC + 1
      JQDIVI = IQOCDV(JOCC)
      NW     = IQOCSP(JOCC)
      CALL MZRESV
      NQRESV = NQRESV - NW
      IF (NQRESV.LT.0)       CALL MZGAR1

C--                Forward division

      IF (JQMODE.EQ.0)  THEN
          IQLN = LQEND(KQT+JQDIVI)
          IQNX = IQLN + NW
          LQEND(KQT+JQDIVI) = IQNX

C--                Reverse division

        ELSE
          IQNX = LQSTA(KQT+JQDIVI)
          IQLN = IQNX - NW
          LQSTA(KQT+JQDIVI) = IQLN
        ENDIF
      NQOCC = JOCC

C--                Fill space reserved with dummy bank

      LQ(KQS+IQLN)   = 12
      LQ(KQS+IQLN+1) = 0
      LQ(KQS+IQLN+2) = 0
      LQ(KQS+IQLN+3) = 0
      LQ(KQS+IQLN+5) = IQLETT(1)
      LQ(KQS+IQLN+6) = 0
      LQ(KQS+IQLN+7) = 0
      LQ(KQS+IQLN+8) = NW - 10
      LQ(KQS+IQLN+9) = 0
      IF (JOCC.NE.NOCC)            GO TO 42

C----              Allocate space for the relocation tables

   46 NWTR = 2*NWTABI + 2
      NWTM = 8*NQSEG
      IF (NWTR+NWTM.LT.NQWKTB)  THEN
C--                Both tables in the normal work area
          LQMTA = LQWKTB
          LQRTA = LQMTA + NWTM
        ELSE
          JQSTMV = -1
          CALL MZFGAP
          IF (NQGAPN.EQ.0)         GO TO 61
          IF (IQGAP(1,1).LT.NWTR)  THEN
C--                MO table in the (smaller) gap, LR table in area
              IF (NQWKTB.LT.NWTR)  GO TO 61
              LQMTA = IQGAP(2,1)
              LQRTA = LQWKTB
            ELSE
C--                MO table in area, LR table in (bigger) gap
              LQMTA = LQWKTB
              LQRTA = IQGAP(2,1)
            ENDIF
        ENDIF

      LQMTE = LQMTA + NWTM
      LQTA  = LQRTA + 1
      LQTE  = LQTA  + 2*NWTABI
      LQRTE = LQTE
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.4)  WRITE (IQLOG,9141) LQMTA,LQMTE,LQTA,LQTE
 9141 FORMAT (' FZIMTB-  Tables : LQMTA,LQMTE,LQTA,LQTE=',4I9)
#endif

C----              Calculate start adrs on the sorted table

      JSOR   = 1
      JOCC   = 1
   52 JQDIVI = ISORDV(JSOR)

      IF (IQMODE(KQT+JQDIVI).EQ.0)  THEN
          LSTA = LQEND(KQT+JQDIVI) - IQOCSP(JOCC)
        ELSE
          LSTA = LQSTA(KQT+JQDIVI)
        ENDIF

      LEND = LSTA + ISORSP(JSOR)

      LENDV(JSOR) = LEND
      LSTAV(JSOR) = LSTA
      JOCC = JOCC + 1

   54 IF (JSOR.EQ.NSOR)            GO TO 55
      JSOR = JSOR + 1
      IF (ISORDV(JSOR).NE.JQDIVI)  GO TO 52
      LSTA = LEND
      LEND = LSTA + ISORSP(JSOR)
      LENDV(JSOR) = LEND
      LSTAV(JSOR) = LSTA
      GO TO 54

C----              Transfer the information from the sorted table
C-                 to the memory occupation table

   55 LMT = LQMTA

      DO  59  JS=1,NQSEG
      JSOR = ITOSOR(JS)
      IF (JSOR.GE.0)               GO TO 57
      LQ(LMT)   = 0
      LQ(LMT+1) = 0
      LQ(LMT+2) = 0
      LQ(LMT+3) = JSOR
      LQ(LMT+4) = JSOR
      LQ(LMT+5) = 0
      LQ(LMT+6) = 0
      LQ(LMT+7) = 0
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)   WRITE (IQLOG,9055) JS, -JSOR
 9055 FORMAT (' FZIMTB-  skip segment',I3,I9,' WORDS')
#endif
      GO TO 59

   57 LQ(LMT)   = ISORDV(JSOR)
      LQ(LMT+1) = 1
      LQ(LMT+2) = 0
      LQ(LMT+3) = LSTAV(JSOR)
      LQ(LMT+4) = LENDV(JSOR)
      LQ(LMT+5) = 0
      LQ(LMT+6) = 0
      LQ(LMT+7) = 0
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  THEN
          WRITE (IQLOG,9058)  JS,LQ(LMT),LQ(LMT+3),LQ(LMT+4)
        ENDIF
 9058 FORMAT (' FZIMTB-  read segment',I3,' into division/from/to'
     F,I3,2I9)
#endif
   59 LMT  = LMT + 8
#include "zebra/qtrace99.inc"
      RETURN

C------            Garbage collection to make room for the table

   61 IF (IFLGAR.GE.2)             GO TO 721
      IXSTOR = ISHFTL (JQSTOR,26)
      IF (IFLGAR.NE.0)             GO TO 63

      IXSTOR = MZIXCO (IXSTOR+21,22,23,24)
      CALL MZGARB (IXSTOR, 0)
      IFLGAR = 1
      IF (JQSTOR.NE.0)             GO TO 46
      IFLGAR = 2
      GO TO 46

C--                Collect also the primary store

   63 IFLGAR = 2
      J = MZIXCO (21,22,23,24)
      CALL MZGARB (J, 0)
      CALL MZSDIV (IXSTOR,-7)
      GO TO 46

C----              All segments to be skipped

   81 NWBKI  = 0
      JRETCD = -4
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)   WRITE (IQLOG,9081)
 9081 FORMAT (' FZIMTB-  skip all segments')
#endif
      GO TO 999

C-------------------------------------------------
C-                 ERROR CONDITIONS
C-------------------------------------------------

C----              User error

C-    JERROR = 15  NQSEG has been changed by the user
  715 JERROR = 15
      IQUEST(14)= NQSEG
      IQUEST(15)= IQ(KQSP+LFISEG)
      NWERR = 2
      GO TO 719

C-    JERROR = 14  invalid division number in segment table
  714 JERROR = 14
      IQUEST(14)= JS
      IQUEST(15)= 0
      IQUEST(16)= IXDIV
      NWERR = 3

  719 JRETCD = 4
      GO TO 999

C----              Not enough space

C-    JERROR = 21  not enough space
  721 JERROR = 21

      JRETCD = 3
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
