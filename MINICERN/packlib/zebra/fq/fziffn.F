*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:31  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZFFNAT)
      SUBROUTINE FZIFFN (JSTAGE)

C-    Read operations for file format native
C-    subsidiary to FZIN

C-    Controlling parameter : JSTAGE
C-
C-    JSTAGE =  1 :  read the first pilot record
C-              2 :  read the tables into memory
C-              3 :  read the bank material

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcwk.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcseg.inc"
#include "zebra/fzcocc.inc"
C--------------    End CDE                             --------------
      DIMENSION    MPILOT(10)
      REAL         CHDATA
#if defined(CERNLIB_QMVDS)
      SAVE         CHDATA
#endif
      EQUIVALENCE (MPILOT(1),IPILI(1))
      EQUIVALENCE (LRTYP,IDI(2)), (ICHDAT,CHDATA)
#if defined(CERNLIB_QREADFULL)
      COMMON /SLATE/ NRSLAT, DUMMY(39)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIF, 4HFN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIFFN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIFFN  ')
#endif
      DATA  CHDATA / 12345.0 /

#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      GO TO (101,201,301), JSTAGE

C-----------------------------------------------------
C-                 obtain and digest next pilot record
C-----------------------------------------------------

  101 LIN  = LQWKFZ

  102 NWR = NQWKTT
      CALL XINBS (LUNI,LRTYP,1,LQ(LIN),NWR)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9102) NRECAI+1,NWR,LRTYP
 9102 FORMAT (' FZIFFN-  hunt for pilot, seen LR #',I7,
     F' with NWRL/LRTYP=',I7,I3)
#endif
      IF (NWR.EQ.0)                GO TO 412
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + NWR + 3
      IF (LRTYP.LE.0)              GO TO 741
      IF (LRTYP.GT.9)              GO TO 741
      IF (LRTYP.GE.4)              GO TO 102
      IF (LRTYP.EQ.1)              GO TO 427

      IACTVI = 2
      IEVFLI = 3 - LRTYP
      IQ(KQSP+LQFI+15) = IQ(KQSP+LQFI+15) + 1
      IF (IOPTIR.NE.0)             GO TO 102
      IF (IEVFLI.LT.IOPTIE)        GO TO 102

      CALL UCOPY (LQ(LIN),IPILI,26)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9107,ERR=107) MPILOT,MPILOT
  107 CONTINUE
 9107 FORMAT (10X,'the 10 pilot control words : ',
#endif
#if (defined(CERNLIB_QDEBPRI))&&(!defined(CERNLIB_HEX))
     F/1X,4O23/1X,4O23/1X,2O23
     F/1H0,9X,F13.1,3(9X,I14)/1X,4(9X,I14)/1X,2(9X,I14))
#endif
#if (defined(CERNLIB_QDEBPRI))&&(defined(CERNLIB_HEX))
     F/10X,4Z17/10X,4Z17/10X,2Z17
     F/1H0,9X,F17.1,3I17/10X,4I17/10X,2I17)
#endif
      IF (IPILI(1).NE.ICHDAT)      GO TO 742

C--                Copy user header

      NWDONE = NWUHCI + 10
      NWIOI  = 0
      NWUHI  = 0
      LFIIOC = LQFI + JAUIOC
      IQ(KQSP+LFIIOC) = 0
      IF (NWUMXI.LE.0)             GO TO 121

      IF (NWUHCI.EQ.0)   THEN
          NWUHI = 0
          GO TO 121
        ENDIF

      NWIOI = JBYT (IOCHI(1),7,5)
      IF (NWIOI.EQ.0)   THEN
          NWIOI = 1
          IF (IOCHI(1).GE.8)       GO TO 743
          IF (IOCHI(1).LT.0)       GO TO 743
          IF (IOCHI(1).EQ.0)  IOCHI(1)=1
        ELSE
          J = JBYT (IOCHI(1),12,5)
          IF (J+1.NE.NWIOI)        GO TO 743
        ENDIF
      IQ(KQSP+LFIIOC) = NWIOI
      CALL UCOPY (IOCHI,IQ(KQSP+LFIIOC+1),NWIOI)

      NWUHI  = NWUHCI - NWIOI
      LUHEAI = LIN + 10 + NWIOI
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  THEN
          N = MIN (8,NWUHI)
          IF (LOGLVI.GE.4)  THEN
              WRITE (IQLOG,9113) NWIOI
              WRITE (IQLOG,9115) (IOCHI(J),J=1,NWIOI)
              N = NWUHI
            ENDIF
          WRITE (IQLOG,9114) NWUHI,N
          WRITE (IQLOG,9115) (LQ(LUHEAI+J),J=0,N-1)
        ENDIF
 9113 FORMAT (10X,I2,' words I/O characteristic for UHV =')
 9114 FORMAT (10X,I4,' words of User Header accepted, dump the first'
     F,I5)
#endif
#if (defined(CERNLIB_QDEBPRI))&&(!defined(CERNLIB_HEX))
 9115 FORMAT (1X,4O23)
#endif
#if (defined(CERNLIB_QDEBPRI))&&(defined(CERNLIB_HEX))
 9115 FORMAT (10X,4Z17)
#endif

C----              Save the segment table

  121 LFISEG = LQFI + JAUSEG
      IF (NWSEGI.EQ.0)             GO TO 124
      IF (NWSEGI.GE.61)            GO TO 744
      IF (NWSEGI.LT.0)             GO TO 744
      IF (NWTABI.EQ.0)  THEN
          NWDONE = NWDONE + NWSEGI
          NWSEGI = 0
        ELSE
          CALL UCOPY (LQ(LIN+NWDONE),IQ(KQSP+LFISEG+1),NWSEGI)
          NWDONE = NWDONE + NWSEGI
        ENDIF
  124 IQ(KQSP+LFISEG) = NWSEGI

C----              Save the text vector

      LTEXT = LQ(KQSP+LQFI-2)
      IF (LTEXT.NE.0)  IQ(KQSP+LTEXT+1)=0
      IF (NWTXI.EQ.0)              GO TO 141

C--                Increase the size of the text buffer if necessary

      IF (LTEXT.EQ.0)              GO TO 131
      NINC = NWTXI + 4 - IQ(KQSP+LTEXT-1)
      IF (NINC.LE.0)               GO TO 131

      NQWKTB = NQWKTT
      CALL MZPUSH (JQPDVS,LTEXT,0,NINC,'.')

C--                Read the text vector

  131 IF (NWDONE.LT.NWR)           GO TO 136
  132 LIN = LQWKTB
      NWR = NQWKTT
      CALL XINBS (LUNI,MRTYP,1,LQ(LIN),NWR)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9132) NRECAI+1,NWR,MRTYP
 9132 FORMAT (10X,'expect pilot continuation, seen LR #',I7,
     F' with NWRL/LRTYP=',I7,I3)
#endif
      IF (NWR.EQ.0)                GO TO 411
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + NWR + 3

      IF (MRTYP-4)           134, 135, 133
  133 IF (MRTYP.LT.7)              GO TO 132
  134 LRTYP = MRTYP
      IF (LRTYP.EQ.1)              GO TO 424
      GO TO 731

  135 NWDONE = 0

C--                Store into the text buffer

  136 IF (NWR-NWDONE.LT.NWTXI)     GO TO 745
      IF (LTEXT.NE.0)  THEN
          IQ(KQSP+LTEXT+1) = NWTXI
          CALL UCOPY (LQ(LIN+NWDONE),IQ(KQSP+LTEXT+5),NWTXI)
        ENDIF
      NWDONE = NWDONE + NWTXI

C----              Save early table words

  141 NTBE   = NWR - NWDONE
      LFIEAR = LQFI + JAUEAR
      IQ(KQSP+LFIEAR) = NTBE

      IF (NTBE.NE.0)  THEN
          IF (NTBE.NE.NWTABI)      GO TO 746
          IF (NTBE.GE.41)          GO TO 746
          CALL UCOPY (LQ(LIN+NWDONE),IQ(KQSP+LFIEAR+1),NTBE)
        ENDIF
      GO TO 999

C-------------------------------------------------
C-                 read table
C-------------------------------------------------

  201 LIN = LQTA + NWTABI
  204 NWR = LQTE - LIN
      CALL XINBS (LUNI,LRTYP,1,LQ(LIN),NWR)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9204) NRECAI+1,NWR,LRTYP
 9204 FORMAT (' FZIFFN-  expect table, seen LR #',I7,
     F' with NWRL/LRTYP=',I7,I3)
#endif
      IF (NWR.EQ.0)                GO TO 411
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + NWR + 3
      IF (LRTYP.EQ.1)              GO TO 424
      IF (LRTYP.LT.4)              GO TO 732
      IF (LRTYP.GE.7)              GO TO 732
      IF (LRTYP.GE.5)              GO TO 204
      LIN = LIN + NWR
      IF   (LIN-LQTE)        204, 999, 747

C-------------------------------------------------
C-                 read the data
C-------------------------------------------------

  301 LMT  = LQMTA
  302 IF (LQ(LMT+1).NE.0)          GO TO 311

C--                Skip segment to be ignored

      NWSK = LQ(LMT+3)
      LIN  = LQWKFZ
#if defined(CERNLIB_QREADFULL)
  304 NWR  = NQWKTT
      CALL XINBS (LUNI,LRTYP,1,LQ(LIN+2),NWR)
      IDI(1) = NWR
      NWRU   = NWR
#endif
#if !defined(CERNLIB_QREADFULL)
  304 NWR  = 3
      CALL XINBF (LUNI,LQ(LIN),NWR)
      IDI(1) = LQ(LIN)
      IDI(2) = LQ(LIN+1)
      NWRU   = 1
#endif
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9314) NRECAI+1,IDI
#endif
      IF (NWR.EQ.0)                GO TO 411
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + IDI(1) + 3
      IF (LRTYP.EQ.1)              GO TO 421
      IF (LRTYP.LT.5)              GO TO 733
      IF (LRTYP.LT.7)              GO TO 304
      IF (LRTYP.GE.9)              GO TO 733
      NWSK = NWSK + IDI(1)
      IF (NWSK.GE.0)               GO TO 307
      IF (LRTYP.EQ.7)              GO TO 304
      GO TO 734

  307 IF (NWSK.EQ.0)               GO TO 318
      GO TO 734

C--                Read segment to accept

  311 LSTA = LQ(LMT+3)
      LEND = LQ(LMT+4)

      LIN = LSTA
  314 NWR = LEND - LIN
      CALL XINBS (LUNI,LRTYP,1,LQ(KQS+LIN),NWR)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9314) NRECAI+1,NWR,LRTYP
 9314 FORMAT (' FZIFFN-  expect bank material, seen LR #',I7,
     F' with NWRL/LRTYP=',I7,I3)
#endif
      IF (NWR.EQ.0)                GO TO 411
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + NWR + 3
      IF (LRTYP.EQ.1)              GO TO 424
      IF (LRTYP.LT.5)              GO TO 733
      IF (LRTYP.LT.7)              GO TO 314
      IF (LRTYP.GE.9)              GO TO 733
      LIN = LIN + NWR
      IF (LIN.GE.LEND)             GO TO 317
      IF (LRTYP.EQ.7)              GO TO 314
      GO TO 735

  317 IF (LIN.GT.LEND)             GO TO 735
  318 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 302
      IF (LRTYP.NE.8)              GO TO 736
      GO TO 999

C-------------------------------------------------
C-                 end-of-file / end-of-run
C-------------------------------------------------

C--                Unexpected end-of-file

  411 JRETCD = 1
      GO TO 999

C--                Normal EOF

  412 JRETCD = -1
      GO TO 999

C------            Unexpected start/end of run

  421 LIN  = LIN  + 2
      NWR  = NWRU - 2

  424 JRETCD = 2
      NWRU   = MIN (78, NWR)
      LFIIOC = LQFI + JAUIOC
      IQ(KQSP+LFIIOC) = NWRU
      CALL UCOPY (LQ(LIN),IQ(KQSP+LFIIOC+1),NWRU)
      GO TO 999

C----              Normal S/E-OF-RUN

  427 JRETCD = -2
      NWUHI  = NWR
      LUHEAI = LIN
      GO TO 999

C-------------------------------------------------
C-                 ERROR CONDITIONS
C-------------------------------------------------

C--                BAD DATA

C-    JERROR = 137  emergency stop record seen
  737 JERROR = 137
      JRETCD = 8
      GO TO 780

C-    JERROR = 136  last bank material record needed is not type 8
  736 JERROR = 136
      GO TO 739

C-    JERROR = 135  end of segm read does not coincide with LR
  735 JERROR = 135
      IQUEST(14)= (LMT-LQMTA)/8 + 1
      IQUEST(15)= LEND - LSTA
      IQUEST(16)= LIN - LEND
      NWERR  = 3
      GO TO 739

C-    JERROR = 134  end of segm skipped does not coincide with LR
  734 JERROR = 134
      IQUEST(14)= (LMT-LQMTA)/8 + 1
      IQUEST(15)= -LQ(LMT+3)
      IQUEST(16)= NWSK
      NWERR  = 3
      GO TO 739

C-                  record of unexpected record type read
C-    JERROR = 133  expect bank material
  733 JERROR = 1

C-    JERROR = 132  expect pilot continuation for table
  732 JERROR = JERROR + 1

C-    JERROR = 131  expect pilot continuation for text vector
  731 JERROR = 131 + JERROR
      IF (LRTYP.EQ.9)              GO TO 737

  739 JRETCD = 5
      GO TO 780

C--                BAD CONSTRUCTION

C-    JERROR = 147  table end does not coincide with LR
  747 JERROR = 147
      IQUEST(14)= NWTABI
      IQUEST(15)= LIN - LQTE
      NWERR  = 2
      GO TO 749

C-    JERROR = 146  number of early table words wrong
  746 JERROR = 146
      IQUEST(14)= NWTABI
      IQUEST(15)= NTBE
      NWERR  = 2
      GO TO 749

C-    JERROR = 145  text vector NWTXI words longer than its record
  745 JERROR = 145
      IQUEST(14)= NWTXI
      IQUEST(15)= NWR - NWDONE
      IQUEST(16)= NWDONE
      NWERR  = 3
      GO TO 749

C-    JERROR = 144  NWSEGI wrong
  744 JERROR = 144
      IQUEST(14)= NWSEGI
      NWERR  = 1
      GO TO 749

C-    JERROR = 143  control wd of I/O char. for user header invalid
  743 JERROR = 143
      IQUEST(14)= NWUHCI
      IQUEST(15)= 0
      IQUEST(16)= IOCHI(1)
      NWERR  = 3
      GO TO 749

C-    JERROR = 142  check-word which should be 12345.0 is wrong
  742 JERROR = 142
      IQUEST(14)= 0
      IQUEST(15)= 0
      IQUEST(16)= IPILI(1)
      IQUEST(17)= ICHDAT
      NWERR  = 4
      GO TO 749

C-    JERROR = 141  LRTYP invalid when hunting for pilot
  741 JERROR = 141

  749 JRETCD = 6
      GO TO 780

C--                READ ERROR

  751 JERROR = 151
      IQUEST(14) = -NWR
      NWERR  = 1
      JRETCD = 7

  780 CONTINUE
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
