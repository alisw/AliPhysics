*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:33  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIFFX (JSTAGE)

C-    Service routine to FZIN for exchange mode

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/mzioc.inc"
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
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIF, 4HFX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIFFX /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIFFX  ')
#endif
      DATA  CHDATA / 12345.0 /

#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      GO TO (101,201,301), JSTAGE

C-----------------------------------------------------
C-                 OBTAIN AND DIGEST NEXT PILOT RECORD
C-----------------------------------------------------

  100 CONTINUE
  101 ISTENI = 1
      LBPARI = LQFI + INCBPI
#if defined(CERNLIB_FZCHANNEL)
      IF (IACMOI.EQ.3)  THEN
          IF (IADOPI.EQ.0)         GO TO 740
        ENDIF
#endif
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOI.EQ.3)  THEN
          IADOPI = IQ(KQSP+LQFI+8)
          IQ(KQSP+LQFI+1) = IADOPI
          IF (IADOPI.EQ.0)         GO TO 740
          CALL VZERO (IQ(KQSP+LBPARI-9),9)
#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_FZMEMORY))
          IF (LOGLVI.GE.3)  WRITE (IQLOG,9142) IADOPI
 9142         FORMAT (' DEVZE FZIN, relative buffer adr =',I9)
#endif
#if defined(CERNLIB_FZMEMORY)
        ENDIF
#endif
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9101) IOPTIR,IOPTIE
 9101 FORMAT (' FZIFFX-  Hunt for pilot, options RUN/EVENT =',2I2)
#endif

      IFLAGI = 1
      CALL FZIREC
      IF (IFLAGI.NE.0)             GO TO 999
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9102) NRECAI,IDI
 9102 FORMAT (' FZIFFX-  Seen LR #',I7,' with NWRL/LRTYP=',I7,I3)
#endif
      ISTENI = 0
      IF (LRTYP.GE.4)              GO TO 100
      IF (LRTYP.EQ.1)              GO TO 422

#if defined(CERNLIB_FQXISN)
      CALL FZITRN (IPILI,10)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (IPILI,10)
        ELSE
          MFO(1) = 3
          MFO(2) = 1
          MFO(3) = 1
          MFO(4) = -1
          JFOEND = 4
          CALL FZITRX (IPILI,10)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999

      IACTVI = 2
      IEVFLI = 3 - LRTYP
      IQ(KQSP+LQFI+15) = IQ(KQSP+LQFI+15) + 1

C--       check values in pilot head

      IF (IPILI(1).NE.ICHDAT)      GO TO 741
      IF (NWTXI .LT.0)             GO TO 743
      IF (NWSEGI.GE.61)            GO TO 743
      IF (NWSEGI.LT.0)             GO TO 743
      IF (NWTABI.LT.0)             GO TO 744
      IF (NWBKI .LT.0)             GO TO 744
      IF (NWUHCI.LT.0)             GO TO 745

C--       skip to start/end of run
      IF (IOPTIR.NE.0)  THEN
          IQ(KQSP+LQFI+26) = 0
          GO TO 100
        ENDIF

      IF (IEVFLI.LT.IOPTIE)        GO TO 100

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9107,ERR=107) MPILOT,MPILOT
  107 CONTINUE
 9107 FORMAT (10X,'The 10 pilot control words : ',
#endif
#if (defined(CERNLIB_QDEBPRI))&&(!defined(CERNLIB_HEX))
     F/1X,4O23/1X,4O23/1X,2O23
     F/1H0,9X,F13.1,3(9X,I14)/1X,4(9X,I14)/1X,2(9X,I14))
#endif
#if (defined(CERNLIB_QDEBPRI))&&(defined(CERNLIB_HEX))
     F/10X,4Z17/10X,4Z17/10X,2Z17
     F/1H0,9X,F17.1,3I17/10X,4I17/10X,2I17)
#endif

C--                Read I/O char. for User Header Vector

      NWIOI  = 0
      NWUHI  = 0
      LFIIOC = LQFI + JAUIOC
      IQ(KQSP+LFIIOC) = 0
      IF (NWUHCI.EQ.0)             GO TO 121

#if defined(CERNLIB_FQXISN)
      CALL FZITRN (IOCHI,1)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (IOCHI,1)
        ELSE
          MFO(1) = 1
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (IOCHI,1)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999

      NWIOI = JBYT (IOCHI(1),7,5)
      IF (NWIOI.EQ.0)   THEN
          NWIOI = 1
          IF (IOCHI(1).GE.8)       GO TO 742
          IF (IOCHI(1).LT.0)       GO TO 742
          IF (IOCHI(1).EQ.0)  IOCHI(1)=1
        ELSE
          J = JBYT (IOCHI(1),12,5)
          IF (J+1.NE.NWIOI)        GO TO 742
        ENDIF
      IF (NWIOI.LT.2)              GO TO 116

#if defined(CERNLIB_FQXISN)
      CALL FZITRN (IQ(KQSP+LFIIOC+2),NWIOI-1)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (IQ(KQSP+LFIIOC+2),NWIOI-1)
        ELSE
          MFO(1) = 1
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (IOCHI(2),NWIOI-1)
          CALL UCOPY (IOCHI(2),IQ(KQSP+LFIIOC+2),NWIOI-1)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999

  116 IQ(KQSP+LFIIOC)   = NWIOI
      IQ(KQSP+LFIIOC+1) = IOCHI(1)

C--                Read User Header Vector

      LUHEAI = LQWKFZ
      NWUHA  = NWUHCI - NWIOI
      NWUHI  = MIN (NWUHA,NWUMXI)
#if defined(CERNLIB_FQXISN)
      CALL FZITRN (LQ(LUHEAI),NWUHI)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (LQ(LUHEAI),NWUHI)
        ELSE
          CALL MZIOCR (IOCHI)
          CALL FZITRX (LQ(LUHEAI),NWUHI)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999
      N4SKII = NWUHA - NWUHI

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


C----              Read the Segment Table

  121 LFISEG = LQFI + JAUSEG
      IF (NWSEGI.EQ.0)             GO TO 124

#if defined(CERNLIB_FQXISN)
      CALL FZITRN (IQ(KQSP+LFISEG+1),NWSEGI)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (IQ(KQSP+LFISEG+1),NWSEGI)
        ELSE
          NSEG = NWSEGI / 3
          MFO(1) = 5
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (IQ(KQSP+LFISEG+1),2*NSEG)
          IF (IFLAGI.NE.0)         GO TO 999
          MFO(1) = 1
          MFO(2) = -1
          CALL FZITRX (IQ(KQSP+LFISEG+1+2*NSEG),NSEG)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999

  124 IQ(KQSP+LFISEG) = NWSEGI

C----              Read the Text Vector

      LTEXT = LQ(KQSP+LQFI-2)
      IF (LTEXT.NE.0)  IQ(KQSP+LTEXT+1)=0
      IF (NWTXI.EQ.0)              GO TO 141

      IF (LTEXT.EQ.0)  THEN
          N4SKII = N4SKII + NWTXI
          GO TO 141
        ENDIF

C--       increase the size of the text buffer if necessary

      NINC = NWTXI + 4 - IQ(KQSP+LTEXT-1)
      IF (NINC.GT.0)  THEN
          CALL FZIREC
          CALL MZPUSH (JQPDVS,LTEXT,0,NINC,'.')
          L4STAI = 0
          CALL FZIREC
        ENDIF

C--       transmit to text buffer

#if defined(CERNLIB_FQXISN)
      CALL FZITRN (IQ(KQSP+LTEXT+5),NWTXI)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (IQ(KQSP+LTEXT+5),NWTXI)
        ELSE
          MFO(1) = 5
          MFO(2) = 0
          JFOEND = 2
          JFOREP = 0
          CALL FZITRX (IQ(KQSP+LTEXT+5),NWTXI)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999
      IQ(KQSP+LTEXT+1) = NWTXI

C----              No early table words

  141 LFIEAR = LQFI + JAUEAR
      IQ(KQSP+LFIEAR) = 0

      CALL FZIREC
      GO TO 999

C-------------------------------------------------
C-                 READ TABLE
C-------------------------------------------------

  201 L4STAI = 0
      IFLAGI = 0
      CALL FZIREC
      LIN = LQTA + NWTABI
#if defined(CERNLIB_FQXISN)
      CALL FZITRN (LQ(LIN),NWTABI)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (LQ(LIN),NWTABI)
        ELSE
          MFO(1) = 1
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (LQ(LIN),NWTABI)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999
      CALL FZIREC
      GO TO 999

C-------------------------------------------------
C-                 READ THE DATA
C-------------------------------------------------

  301 L4STAI = 0
      IFLAGI = 0
      CALL FZIREC
      LMT  = LQMTA
  302 IF (LQ(LMT+1).NE.0)          GO TO 311

C----              Skip segment to be ignored

      N4SKII = N4SKII - LQ(LMT+3)
      GO TO 348

C--------          Read segment to accept

  311 LSTA = LQ(LMT+3)
      LEND = LQ(LMT+4)
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)             GO TO 341

C------            Read segment in exchange mode, bank-by-bank

      LIN = LSTA

C-       Next bank : first word
  322 MFO(1) = 1
      MFO(2) = -1
      JFOEND = 2
      CALL FZITRX (LQ(KQS+LIN),1)
      IF (IFLAGI.NE.0)             GO TO 999
      IWD  = LQ(KQS+LIN)
      NST  = JBYT (IWD,1,16) - 12
#endif
#if (!defined(CERNLIB_FQXISN))&&(defined(CERNLIB_QDEBPRI))
      IF (LOGLVI.GE.4)  WRITE (IQLOG,9322) LIN,L4CURI-L4STAI-1,NST
 9322 FORMAT (' FZIFFX-  Next bank : Lst, NWdone, NIO+NL =',I9,2I6)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (NST.LT.0)                GO TO 326

C--                True bank

      IQLN  = LIN
      IQLS  = LIN + NST + 1
      IF (IQLS+8.GE.LEND)          GO TO 752

C-       I/O words, links, next-link, up-link
      MFO(1) = 1
      MFO(2) = NST + 2

C-       origin, numeric ID
      MFO(3) = 2
      MFO(4) = 2

C-       Hollerith ID
      MFO(5) = 5
      MFO(6) = 1

C-       NL, NS, ND, status
      MFO(7) = 1
      MFO(8) = -1
      JFOEND = 8
      CALL FZITRX (LQ(KQS+LIN+1), NST+9)
      IF (IFLAGI.NE.0)             GO TO 999
      IQNIO = JBYT (IQ(KQS+IQLS),19,4)
      IQNL  = IQ(KQS+IQLS-3)
      IQND  = IQ(KQS+IQLS-1)
#endif
#if (!defined(CERNLIB_FQXISN))&&(defined(CERNLIB_QDEBPRI))
      IF (LOGLVI.GE.4)  THEN
          IQID  = IQ(KQS+IQLS-4)
          WRITE (IQLOG,9323) IQID,IQNL,IQND,IQNIO
        ENDIF
 9323 FORMAT (' FZIFFX-       ID,NL,ND,NIO = ',A4,2I8,I4)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IQNIO+IQNL.NE.NST)       GO TO 751
      IF (IQND.LT.0)               GO TO 751

C-       data words
      LIN   = IQLS + IQND + 9
      IF (IQND.EQ.0)               GO TO 324
      IF (LIN.GT.LEND)             GO TO 753

      CALL MZIOCR (LQ(KQS+IQLN))
      CALL FZITRX (IQ(KQS+IQLS+1), IQND)
      IF (IFLAGI.NE.0)             GO TO 999
  324 IF (LIN.LT.LEND)             GO TO 322
      GO TO 348

C--                Short dead region

  326 NWD = JBYT (IWD,17,IQDROP-17)
      IF (NWD.EQ.0)                GO TO 751
      IF (NWD.NE.NST+12)           GO TO 751
      IF (JBYT(IWD,IQDROP,IQBITW-IQDROP).NE.1)  GO TO 751
      IF (NWD.GT.1)  THEN
          IF (LIN+NWD.GT.LEND)     GO TO 754
          MFO(1) = 0
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (LQ(KQS+LIN+1), NWD-1)
          IF (IFLAGI.NE.0)         GO TO 999
        ENDIF
      LIN = LIN + NWD
      IF (LIN.LT.LEND)             GO TO 322
      GO TO 348
#endif

C------            Read segment in native mode

  341 CALL FZITRN (LQ(KQS+LSTA),LEND-LSTA)
      IF (IFLAGI.NE.0)             GO TO 999

C--                Verify bank chaining

      LIN  = LSTA
  344 IWD  = LQ(KQS+LIN)
      NST  = JBYT (IWD,1,16) - 12
      IF (NST.LT.0)                GO TO 346

C--                True bank

      IQLN  = LIN
      IQLS  = LIN + NST + 1
      IF (IQLS+8.GE.LEND)          GO TO 752

      IQNIO = JBYT (IQ(KQS+IQLS),19,4)
      IQNL  = IQ(KQS+IQLS-3)
      IQND  = IQ(KQS+IQLS-1)
      IF (IQNIO+IQNL.NE.NST)       GO TO 751
      IF (IQND.LT.0)               GO TO 751

      LIN   = IQLS + IQND + 9
      IF   (LEND-LIN)        753, 348, 344

C--                Short dead region

  346 NWD = JBYT (IWD,17,IQDROP-17)
      IF (NWD.EQ.0)                GO TO 751
      IF (NWD.NE.NST+12)           GO TO 751
      IF (JBYT(IWD,IQDROP,IQBITW-IQDROP).NE.1)  GO TO 751
      LIN = LIN + NWD
      IF   (LEND-LIN)        754, 348, 344

C----              End of segment

  348 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 302

      IF (N4SKII.NE.0)  THEN
#if defined(CERNLIB_FQXISN)
          CALL FZITRN (IPILI,0)
#endif
#if !defined(CERNLIB_FQXISN)
          IF (IDAFOI.EQ.0)  THEN
              CALL FZITRN (IPILI,0)
            ELSE
              CALL FZITRX (IPILI,0)
            ENDIF
#endif
          IF (IFLAGI.NE.0)         GO TO 999
        ENDIF
      CALL FZIREC
      IF (N4RESI.NE.0)             GO TO 755
      GO TO 999

C-------------------------------------------------
C-                 START / END-OF-RUN
C-------------------------------------------------

  422 LUHEAI = LQWKFZ
      NWUHI  = MIN (IDI(1), NQWKTT)
#if defined(CERNLIB_FQXISN)
      CALL FZITRN (LQ(LUHEAI),NWUHI)
#endif
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOI.EQ.0)  THEN
          CALL FZITRN (LQ(LUHEAI),NWUHI)
        ELSE
          MFO(1) = 2
          MFO(2) = -1
          JFOEND = 2
          CALL FZITRX (LQ(LUHEAI),NWUHI)
        ENDIF
#endif
      IF (IFLAGI.NE.0)             GO TO 999

      CALL FZIREC
      JRETCD = -2
      IQ(KQSP+LBPARI-7) = 0
      GO TO 999

C-------------------------------------------------
C-                 ERROR CONDITIONS
C-------------------------------------------------

#if defined(CERNLIB_FZMEMORY)||defined(CERNLIB_FZCHANNEL)
C-    JERROR = 240  user routine or buffer not connected for C/M mode
  740 JERROR = 240
      JRETCD = 4
      GO TO 999
#endif

C--                Bad construction

C-    JERROR = 241  check-word which should be 12345.0 is wrong
  741 JERROR = 241
      IQUEST(14)= 0
      IQUEST(15)= 0
      IQUEST(16)= IPILI(1)
      IQUEST(17)= ICHDAT
      NWERR  = 4
      GO TO 749

C-    JERROR = 242  control wd of I/O char. for user header invalid
  742 JERROR = 242
      IQUEST(14)= NWUHCI
      IQUEST(15)= 0
      IQUEST(16)= IOCHI(1)
      NWERR  = 3
      GO TO 749

C-    JERROR = 243  NWSEGI or NWTXI wrong
  743 JERROR = 243
      IQUEST(14)= NWSEGI
      IQUEST(15)= NWTXI
      NWERR  = 2
      GO TO 749

C-    JERROR = 244  NWTABI or NWBKI wrong
  744 JERROR = 244
      IQUEST(14)= NWTABI
      IQUEST(15)= NWBKI
      NWERR  = 2
      GO TO 749

C-    JERROR = 245  NWUHCI wrong
  745 JERROR = 245
      IQUEST(14)= NWUHCI
      NWERR  = 2

  749 JRETCD = 6
      GO TO 780

C--                Bad data

C-    JERROR = 251  inconsistent bank parameters
  751 JERROR = 251
      GO TO 759

C-    JERROR = 252  link part of bank overshoots segment end
  752 JERROR = 252
      GO TO 759

C-    JERROR = 253  data part of bank overshoots segment end
  753 JERROR = 253
      GO TO 759

C-    JERROR = 254  short dead region overshoots segment end
  754 JERROR = 254
      GO TO 759

C-    JERROR = 255  bank material does not end exactly with LR
  755 JERROR = 255

  759 JRETCD = 5

  780 IQ(KQSP+LBPARI-9) = -3
      IQ(KQSP+LBPARI-1) = 0
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
