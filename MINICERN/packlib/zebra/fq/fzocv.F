*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:46  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if !defined(CERNLIB_FQXISN)
      SUBROUTINE FZOCV (MS,MT)

C-    Convert for output with copy
C-    from source in native to target in exchange data format

#include "zebra/quest.inc"
#include "zebra/mzioc.inc"
C--------------    End CDE                             --------------
      DIMENSION    MS(99), MT(99)

      DOUBLE PRECISION   THDB
      DIMENSION    THIS(2)
      EQUIVALENCE  (THDB,THIS)
      EQUIVALENCE  (ITHA,THA,THIS(1)), (ITHB,THB,THIS(2))

C----            Conversion Control in /MZIOC/ :

C-    for a given call translation source MS -> target MT is done
C-        either  for a complete batch of NWFOTT words
C-                      if NWFODN.EQ.0  and  NWFOAV.GE.NWFOTT

C-            or  for the first instalment of a batch of NWFOTT words
C-                      if NWFODN.EQ.0  and  NWFOAV.LT.NWFOTT

C-            or  for a new instalment of N=MIN(NWFOAV,NWFOTT-NWFODN)
C-                      words, if  NWFODN.NE.0

C-    *  marks words to be initialized by the caller
C-       only for the call at the beginning of a new batch
C-
C-         NWFOAV  number of words available in the buffer to receive
C-                                                         the result
C-                    set by the caller whenever a new lot of data
C-                    becomes available for conversion (new buffer)
C-                    counted down by FZOCV
C-       * NWFOTT  total number of words in the batch to be done,
C-                    maybe in several instalments
C-       * NWFODN  number of words in the batch already done
C-                    set to zero by the caller at start of batch
C-                    (in fact MZIOCR sets it to zero)
C-         NWFORE  n.w. remaining to be done for the pending batch
C-                    set by FZOCV, zero if end of batch
C-
C-         IFOCON  remembers the last conversion problem
C-            (1)  error code if -ve, expected type if +ve
C-            (2)  location of the word
C-            (3)  content  of the word

C-         MFOSAV
C-          (1+2)  saves type and word-count for re-entry

C-       * JFOEND  position of last sector plus 1
C-       * JFOREP  position of repeat sector descr.
C-         JFOCUR  position of current sector description

C-   * MFO(JFO+1)  t= sector type as in format
C-   *        +2)  c= word count  as in format
C-                    c > 0 :  no. of words
C-                    c = 0 :  dynamic sector
C-                    c < 0 :  indefinite sector, rest of the bank

C-            JMT  # of words done so far for the current call
C-            JMS  # of words done so far for the current batch

#include "fzocvd1.inc"
#include "fzocvd2.inc"


      JMT = 0
      IF (NWFODN.NE.0)             GO TO 30

      NWFORE = NWFOTT
      JMTEX  = MIN (NWFORE,NWFOAV)

      JMS    = 0
      JFOCUR = 0
      IFOCON(1) = 0

C------            Start next sector

   21 ITYPE  = MFO(JFOCUR+1)
      IF (ITYPE.EQ.7)              GO TO 24
      NWSEC = MFO(JFOCUR+2)
      IF   (NWSEC)           22, 23, 31

C--                Rest of the bank

   22 NWSEC = NWFORE
      GO TO 31

C--                Dynamic sector

   23 IWORD = MS(JMS+1)
      NWSEC = IWORD
      GO TO 25

C--                Self-describing sector

   24 IWORD = MS(JMS+1)
      ITYPE = MOD (IWORD,16)
      NWSEC = IWORD/16

   25 MT(JMT+1) = IWORD
      JMT = JMT + 1
      JMS = JMS + 1
      NWFORE = NWFORE - 1

      IF (ITYPE.GE.8)              GO TO 27
      IF (NWSEC.EQ.0)              GO TO 29
      IF (NWSEC.GT.0)              GO TO 31

C--                Faulty sector control word

   27 IFOCON(1) = -1
      IFOCON(2) = JMS
      IFOCON(3) = IWORD

C--                Rest of the bank is unused

   29 ITYPE = 0
      NWSEC = NWFORE
      GO TO 31

C--                RE-ENTRY TO CONTINUE

   30 JMTEX  = MIN (NWFORE,NWFOAV)
      JMS    = NWFODN
      ITYPE  = MFOSAV(1)
      NWSEC  = MFOSAV(2)

C------            CONVERSION LOOPS

   31 NWDO   = MIN (NWSEC,JMTEX-JMT)
      IF (NWDO.EQ.0)               GO TO 801
      IF (ITYPE.LE.0)              GO TO 91
      GO TO (101,201,301,401,501,101,101), ITYPE

C--                Rest of the bank unused

   91 CALL VZERO (MT(JMT+1),NWDO)
      JMS = JMS + NWDO
      JMT = JMT + NWDO
      GO TO 801

C--                B - bit strings

C--                I - integers

#include "fzocvfi.inc"
* Ignoring t=pass

C--                F - floating

#include "fzocvff.inc"
* Ignoring t=pass

C--                D - double precision

  401 NDPN   = (NWDO+1) / 2
      NWDODB =  NDPN * 2
#include "fzocvfd.inc"
      IF (NWDODB.EQ.NWDO)          GO TO 801
      IF (NWDODB.GT.NWSEC)         GO TO 471
      IF (NWDODB.GT.NWFORE)        GO TO 471
      NWDO = NWDODB
      GO TO 801

C--                Error : odd number of double-precision words

  471 JMS = JMS - 1
      JMT = JMT - 1
      IFOCON(1) = -2
      IFOCON(2) = JMS
      IFOCON(3) = NWDO
      GO TO 801

C--                H - hollerith

  501 CONTINUE
#include "fzocvfh.inc"
* Ignoring t=pass

C----              COPY AS IS

#include "fzocvjf.inc"
  201 CONTINUE
  101 CONTINUE
#include "fzocvfai.inc"

C------            END OF SECTOR

  801 NWFORE = NWFOTT - JMS
      IF (JMT.GE.JMTEX)            GO TO 804
      JFOCUR = JFOCUR + 2
      IF (JFOCUR.LT.JFOEND)        GO TO 21
      JFOCUR = JFOREP
      GO TO 21

C--                Data or buffer exhausted

  804 IQUEST(1) = JMT
      NWFOAV = NWFOAV - JMT
      IF (NWFORE.EQ.0)             RETURN

C--                Ready for re-entry

      NWFODN    = JMS
      MFOSAV(1) = ITYPE
      MFOSAV(2) = NWSEC - NWDO
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
