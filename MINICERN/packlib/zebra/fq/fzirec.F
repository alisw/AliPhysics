*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIREC

C-    Logical record controls, exchange file format
C-    service routine to FZIFFX of FZIN

C-    This routine is called with
C-       IFLAGI > 0  from FZIFFX to initiate the next d/s making sure
C-                   that the block containing its beginning is
C-                   in the buffer, skipping if nec. trailing records
C-                   of the previous d/s
C-              = 0  from FZIFFX to copy the buffer control params.
C-                   between the control-bank and /FZCI/ (for speed)
C-              < 0  from FZITR to read continuation blocks (phR)
C-                   into the buffer

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
#include "fzntolds.inc"
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIR, 4HEC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIREC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIREC  ')
#endif


#include "zebra/qtrace.inc"


C-    buffer parameters :
C-
C-          -9         non-zero : last LR abended
C-          -7         expected next Ph record number
C-          -6         =0 if current block is steering
C-                     =1 if current block is last in burst
C-                     =2 if cur. block is last-but-one in burst, etc
C-          -5         # of fast records still pending
C-
C-          -4 N4SKII, # of words to be skipped
C-                                      before next transmission
C-          -3 N4RESI, # of words in LR still to be done
C-          -2 N4DONI, # of words already out of buffer
C-          -1 N4ENDI, # of words in buffer before start of next LR
C-                     if =0 :      buffer empty
C-                     if =MAXREI : LR continues in next block
C-   LBPARI +0         maximum size of buffer, words
C-          +1         expected size of PhR, local machine words
C-          +2 INCBUF  step to buffer
C-          +3         off-set from start-of-buffer for reading
C-                                  the packed record
C-          +4         (off-set for output)
C-          -1         space for left double-precision word saved
C-   L4STAI +0         first word of normal buffer

      IFLIN  = IFLAGI
      IFLAGI = 0
      LBPARI = LQFI + INCBPI
      IF   (IFLIN)                 61, 71, 21

C-----------------------------------------------------------
C------            IFLAGI > 0 :  start new d/s
C-----------------------------------------------------------

   21 ICARRL = 0
      IFLRST = IQ(KQSP+LBPARI-9)
      INCBUF = IQ(KQSP+LBPARI+2)
      L4STAI = KQSP+8 + LBPARI + INCBUF

      NLRPAD = 0
      JFAST  = IQ(KQSP+LBPARI-6)
      N4SKII = 0
      N4ENDI = IQ(KQSP+LBPARI-1)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +    WRITE (IQLOG,9022) IFLIN,IFLRST,JFAST,N4ENDI,MAXREI
 9022 FORMAT (1X/' FZIREC-  Going for next LR, Buffer Status :'
     F/10X,'IFLAGI, Restart?, Fast?, NWtoLR, NWbuf =',6I6)
#endif

C--                Re-start LR after error

      IF (IFLRST.NE.0)  THEN
          IQ(KQSP+LBPARI-9) = 0
          IF (N4ENDI.GT.0)  THEN
              IF (N4ENDI.LT.MAXREI-1)  GO TO 51
            ENDIF
          IFLAGI = -1
          GO TO 27
        ENDIF

C--                Start afresh

      IF (N4ENDI.NE.0)             GO TO 31
      IFLAGI = -2
   27 CALL VZERO (IQ(KQSP+LBPARI-8),8)
      N4RESI = 0
      N4DONI = 0
      GO TO 46

C----              Normal continuation on current buffer

   31 N4RESI = IQ(KQSP+LBPARI-3)
      N4DONI = IQ(KQSP+LBPARI-2)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9031) N4RESI,N4DONI,N4ENDI
 9031 FORMAT (10X,'Status of last LR : NWrest, NWdone, NWend ='
     F,3I5)
#endif

C--                Skip trailing words of last LR

      IF (N4RESI.EQ.0)             GO TO 41
   34 NWNEW  = MIN (N4ENDI-N4DONI, N4RESI)
      N4DONI = N4DONI + NWNEW
      N4RESI = N4RESI - NWNEW
      IF (N4RESI.EQ.0)             GO TO 41
      IF (N4ENDI.NE.MAXREI)        GO TO 821
      N4SKII = N4RESI
      NWFAST = IQ(KQSP+LBPARI-5) * MAXREI
      IF (NWFAST.LE.N4SKII)        GO TO 46
      N4SKII = NWFAST
      N4RESI = NWFAST
      IFLAGI = -2
      GO TO 46

C--                New Ph record if buffer exhausted

   41 IF (N4DONI.NE.N4ENDI)        GO TO 822
      ICARRL = 0
      IF (N4ENDI.LT.MAXREI-1)      GO TO 51
      IF (N4ENDI.GT.MAXREI-1)      GO TO 46
      MCARRL = LQ(L4STAI+MAXREI-1)
      IF (MCARRL.EQ.0)             GO TO 46
      ICARRL = 1
      JRECNO = IQ(KQSP+LQFI+33)
      JDSADR = N4ENDI
      JFAST  = IQ(KQSP+LBPARI-6)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9043) MCARRL
 9043 FORMAT (' FZIREC-  Saved LR length',I6,' across to next PhR')
#endif

   46 IF (IFIFOI.LE.1)  THEN
          CALL FZIPHR
#if defined(CERNLIB_FZDACC)
        ELSEIF (IFIFOI.EQ.2)  THEN
          CALL FZIPHD
#endif
#if defined(CERNLIB_FZMEMORY)
        ELSEIF (IFIFOI.EQ.3) THEN
          CALL FZIPHM
#endif
#if defined(CERNLIB_FZALFA)
        ELSEIF (IFIFOI.EQ.4)  THEN
          CALL FZIPHA
#endif
        ENDIF

C--       read failure ?
      IF (IFLAGI.NE.0)  THEN
          IF (IFLIN.GT.0)  THEN
              IF (JRETCD.EQ.1)  JRETCD=-1
            ENDIF
          GO TO 999
        ENDIF

C--       LR length last word of previous block ?
      IF (ICARRL.NE.0)  THEN
          N4ENDI = N4DONI - 1
          LQ(L4STAI+N4ENDI) = MCARRL
        ENDIF

C--       still discarding previous d/s ?
      IF (IFLIN.GT.0)  THEN
          IF (N4RESI.GT.0)         GO TO 34
        ENDIF

C----              Start new logical record

   51 N4DONI = N4ENDI
      N4SKII = 0
      IF (N4ENDI.EQ.MAXREI-1)      GO TO 41
      IF (ICARRL.EQ.0)  THEN
          JRECNO = IQ(KQSP+LQFI+33)
          JDSADR = N4DONI
          JFAST  = IQ(KQSP+LBPARI-6)
        ENDIF
      ICARRL = 0
      L4CURI = L4STAI + N4DONI + 2
      IDI(1) = LQ(L4CURI-2)
      IDI(2) = LQ(L4CURI-1)
      IF (IDI(1).EQ.0)  IDI(2)=5
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9053) JFAST,N4DONI,IDI
 9053 FORMAT (' FZIREC-  Start LR : Fast?, NWdone, Length, Type='
     F,I5,I6,I8,I3)
#endif
      N4RESI = IDI(1)
      IF (N4RESI.LT.0)             GO TO 826
      IF (N4RESI.GT.NTOLDS)        GO TO 826
      IF (IDI(2).GE.5)             GO TO 57
      IF (IDI(2).EQ.4)             GO TO 53

C--                LR type 1, 2, 3

      IF (IDI(2).LE.0)             GO TO 825
      IF (IFLIN.LT.0)              GO TO 824
      IF (JFAST.NE.0)              GO TO 827
      IQ(KQSP+LQFI+31) = JRECNO
      IQ(KQSP+LQFI+32) = JDSADR

   53 N4DONI = N4DONI + 2
      N4ENDI = MIN (MAXREI, N4DONI+N4RESI)
      L4ENDI = L4STAI + N4ENDI
      NRECAI = NRECAI + 1
      IF (IFLIN.GE.0)              GO TO 77
      N4SKII = NSKISV
      GO TO 77

C--                Skip padding records

   57 IF (IDI(2).GE.7)             GO TO 825
      IF (N4RESI.GT.MAXREI)        GO TO 828
      NLRPAD = NLRPAD + 1
      IF (NLRPAD.GE.5)             GO TO 829
      N4DONI = N4DONI + 1
      N4ENDI = MIN (MAXREI, N4DONI+N4RESI)

      NWNEW  = MIN (N4ENDI-N4DONI, N4RESI)
      N4DONI = N4DONI + NWNEW
      N4RESI = N4RESI - NWNEW
      GO TO 41

C-----------------------------------------------------------
C------            IFLAGI = -1 : ready continuation Ph / L record
C-----------------------------------------------------------

   61 ICARRY = 0
      N4DONI = L4CURI - L4STAI
      NWNEW  = N4DONI - IQ(KQSP+LBPARI-2)
      N4RESI = N4RESI - NWNEW

      IF (N4DONI.EQ.MAXREI+1)  THEN
          ICARRY = 1
          MCARRY = LQ(L4CURI-2)
          N4RESI = N4RESI + 1
        ENDIF

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9061) N4RESI,N4DONI,MAXREI,ICARRY
 9061 FORMAT (1X/' FZIREC-  Need continuation PhR/LR, Buffer Status :'
     F/10X,'Rest LR, NWdone, NWbuf, NWcarry =',3I7,I4)
#endif

C--                Start continuation LR

      IF (N4RESI.EQ.0)  THEN
          NLRPAD = 0
          NSKISV = N4SKII
          N4SKII = 0
          GO TO 41
        ENDIF

C--                Get continuation PhR

      IF (IFIFOI.LE.1)  THEN
          CALL FZIPHR
#if defined(CERNLIB_FZDACC)
        ELSEIF (IFIFOI.EQ.2)  THEN
          CALL FZIPHD
#endif
#if defined(CERNLIB_FZMEMORY)
        ELSEIF (IFIFOI.EQ.3) THEN
          CALL FZIPHM
#endif
#if defined(CERNLIB_FZALFA)
        ELSEIF (IFIFOI.EQ.4)  THEN
          CALL FZIPHA
#endif
        ENDIF
      IF (IFLAGI.NE.0)             GO TO 999

      L4CURI = L4STAI + N4DONI
      L4ENDI = L4STAI + N4ENDI

      IF (ICARRY.NE.0)  THEN
          N4RESI = N4RESI + 1
          N4DONI = N4DONI - 1
          L4CURI = L4CURI - 1
          LQ(L4CURI) = MCARRY
        ENDIF
      GO TO 77

C-----------------------------------------------------------
C------            IFLAGI = 0 :  set buffer parameters
C-----------------------------------------------------------

   71 IF (L4STAI.NE.0)             GO TO 74

C--                Copy parameters from control bank to COMMON

      N4SKII = IQ(KQSP+LBPARI-4)
      N4RESI = IQ(KQSP+LBPARI-3)
      N4DONI = IQ(KQSP+LBPARI-2)
      N4ENDI = IQ(KQSP+LBPARI-1)
      INCBUF = IQ(KQSP+LBPARI+2)
      L4STAI = KQSP+8 + LBPARI + INCBUF
      L4CURI = L4STAI + N4DONI
      L4ENDI = L4STAI + N4ENDI
      GO TO 999

C--                Up-date parameters in control bank

   74 N4DONI = L4CURI - L4STAI
      NWNEW  = N4DONI - IQ(KQSP+LBPARI-2)
      N4RESI = N4RESI - NWNEW

   77 IQ(KQSP+LBPARI-4) = N4SKII
      IQ(KQSP+LBPARI-3) = N4RESI
      IQ(KQSP+LBPARI-2) = N4DONI
      IQ(KQSP+LBPARI-1) = N4ENDI

#include "zebra/qtrace99.inc"
      RETURN

C-----------------------------------------------------------
C-                 ERROR CONDITIONS
C-----------------------------------------------------------

C-    JERROR = 221  LR overshoots physical record control
  821 JERROR = 221
      IQUEST(14) = N4DONI
      IQUEST(15) = N4RESI
      NWERR = 2
      GO TO 871

C-    JERROR = 222  LR undershoots physical record control
  822 JERROR = 222
      IQUEST(14) = N4DONI
      IQUEST(15) = N4ENDI
      NWERR = 2
      GO TO 871

C-    JERROR = 223  LR expected to start at start of new PhR
C!823 JERROR = 223
C!    IQUEST(14) = N4ENDI
C!    NWERR = 1
C!    GO TO 871

C-    JERROR = 224  LR type 1,2,3 seen when 4 expected
  824 JERROR = 224
      IQUEST(14) = N4ENDI
      NWERR = 1
      GO TO 871


C-    JERROR = 225  Faulty LR type
  825 JERROR = -1

C-    JERROR = 226  Faulty LR length
  826 JERROR = JERROR - 1

C-    JERROR = 227  LR of type 1,2,3 must start on steering block
  827 JERROR = JERROR - 1

C-    JERROR = 228  Padding record longer than one physical record
  828 JERROR = JERROR - 1

C-    JERROR = 229  More than 4 padding records in a row
  829 JERROR = 229 + JERROR

      IQUEST(14) = N4DONI
      IQUEST(15) = N4RESI
      IQUEST(16) = IDI(2)
      NWERR  = 3
      N4ENDI = 0

  871 IQ(KQSP+LBPARI-9) = -2
      IFLAGI = -2
      JRETCD = 6
      GO TO 77
      END
*      ==================================================
#include "zebra/qcardl.inc"
