*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:50  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOTAB

C-    Construct table of material to be output,
C-    called from FZOUT, RZOUT, MZCOPY

#include "zebra/zlimit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcseg.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LMT,LQMTB), (LS,IQLS), (LNX,IQNX)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOT, 4HAB   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOTAB /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOTAB  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit0.inc"


#include "zebra/qtrace.inc"

  101 MQDVGA = 0
      MQDVWI = 0
      NWBKX  = 0
      NQSEG  = 0

      IF (MODTBX.NE.0)  THEN
          NSTEP = 4
        ELSE
          NSTEP = 2
        ENDIF

      CALL MZCHLS (-7,LENTRX)
      IF (IQFOUL.NE.0)                    GO TO 911
      IF (JBIT(IQ(KQS+LENTRX),IQDROP).NE.0)  GO TO 911

      IF (IOPTXD.NE.0)             GO TO 301
      IF (IOPTXS.EQ.0)             GO TO 201

C----              Single bank output

      CALL MZTABM
      LQTE = LQTA + NSTEP

      LQ(LQTA)   = IQLN
      LQ(LQTA+1) = IQNX
      NWBKX      = IQNX - IQLN

      IOPTXL = 0
      IOPTXD = 0
      GO TO 361

C-------------------------------------------------
C-                 d/s defined by LENTRX
C-------------------------------------------------

  201 IF (IOPTXM.EQ.0)  THEN
          IF (IOPTXL.EQ.0)  THEN
              CALL MZMARK (IXDIVX,LENTRX,'.',0,0)
            ELSE
              CALL MZMARK (IXDIVX,LENTRX,'L',0,0)
            ENDIF
        ENDIF

      LIMLOW = LQLIML - IQ(KQS+LQLIML-3) - JBYT(IQ(KQS+LQLIML),19,4) - 1
      LIMHI  = LQLIMH + IQ(KQS+LQLIMH-1) + 9

C----              Memory occupation table

      CALL MZTABM

C--                First division concerned

  226 IF (LIMLOW.LT.LQ(LQMTA+4))   GO TO 227
      LQMTA = LQMTA + 8
      GO TO 226

C--                Last division concerned

  227 LMT = LQMTA
      LQ(LMT+3) = LIMLOW
  228 IF (LIMHI.LE.LQ(LMT+4))      GO TO 229
      LMT = LMT + 8
      GO TO 228

  229 LQ(LMT+4) = LIMHI
      LQMTE = LMT + 8

C----              Relocation table

      LMT   = LQMTA

C--                Next division

  232 IF (LQ(LMT+1).LT.0)          GO TO 249
      NWDIV  = 0
      JDIV   = LQ(LMT)
      LNX    = LQ(LMT+3)
      LQMTC2 = LQ(LMT+4)
      NEW = 0

C--                Next bank

  234 MODE = NEW
  235 LN   = LNX
      IF (LN.GE.LQMTC2)            GO TO 247
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLN (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 912
#endif
#if !defined(CERNLIB_QDEBUG)
      NST = JBYT(LQ(KQS+LN),1,16) - 11
      IF (NST.GT.0)  THEN
          LS  = LN + NST
          LNX = LS + IQ(KQS+LS-1) + 9
        ELSE
          LS  = LN - 8
          LNX = LN + NST + 11
        ENDIF
#endif
      NEW = JBIT(IQ(KQS+LS),IQMARK)
      IF (NEW.EQ.1)  IQ(KQS+LS)= MSBIT0 (IQ(KQS+LS),IQMARK)
      IF (NEW.EQ.MODE)             GO TO 235
      IF (NEW.EQ.0)                GO TO 244

C--                Start of a live group

      LQ(LQTE) = LN
      IF (LQTE.LT.LQRTE)           GO TO 234
      CALL MZTABH
      IF (IQPART.LE.0)             GO TO 234
      GO TO 261

C--                Start of a dead group

  244 NWDIV = NWDIV + (LN - LQ(LQTE))
      LQ(LQTE+1) = LN
      LQTE = LQTE + NSTEP
      GO TO 234

C--                End of division

  247 IF (MODE.EQ.0)               GO TO 248
      NWDIV = NWDIV + (LN - LQ(LQTE))
      LQ(LQTE+1) = LN
      LQTE = LQTE + NSTEP

C--                Add to segment table

  248 IF (NWDIV.EQ.0)              GO TO 249
      NWBKX = NWBKX + NWDIV
      NQSEG = NQSEG + 1
      IQSEGH(1,NQSEG) = IQDN1(KQT+JDIV)
      IQSEGH(2,NQSEG) = IQDN2(KQT+JDIV)
      IQSEGD(NQSEG)   = NWDIV

C--                Step to next division

  249 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 232
      GO TO 361

C------            Not enough space, garbage collection

  261 IF (JFLGAX.NE.0)             GO TO 264
      JFLGAX = 1
      LQSYSR(KQT+1) = LENTRX
      IXGARB = MZIXCO (IXDIVX,21,22,23)
      IXGARB = MZIXCO (IXGARB,24,0,0)
      CALL MZGARB (IXGARB,0)
      LENTRX = LQSYSR(KQT+1)
      GO TO 269

  264 IF (JFLGAX.GE.2)             GO TO 913
      JFLGAX = 2
      IF (JQSTOR.EQ.0)             GO TO 913
      IXGARB = MZIXCO (21,22,23,24)
      CALL MZGARB (IXGARB,0)
      CALL MZSDIV (IXDIVX,0)
  269 IF (MODTBX.EQ.0)             GO TO 101
      IQUEST(1) = 1
      GO TO 101

C-------------------------------------------------
C-                 complete divisions
C-------------------------------------------------

  301 MQDVGA = MZDVAC (IXDIVX)
      CALL MZTABM
      IOPTXL = 1

      IF (IOPTXI.EQ.0)             GO TO 321

C----              Immediate dump of divisions

  305 IF (LQ(LMT+1).NE.3)          GO TO 307
      NWDIV = LQ(LMT+4) - LQ(LMT+3)
      NWBKX = NWBKX + NWDIV
      LQ(LQTE)   = LQ(LMT+3)
      LQ(LQTE+1) = LQ(LMT+4)
      LQTE = LQTE + NSTEP

C--                Add to segment table

      JDIV  = LQ(LMT)
      NQSEG = NQSEG + 1
      IQSEGH(1,NQSEG) = IQDN1(KQT+JDIV)
      IQSEGH(2,NQSEG) = IQDN2(KQT+JDIV)
      IQSEGD(NQSEG)   = NWDIV
  307 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 305
      GO TO 361

C----              Filtered dump of divisions

  321 CONTINUE

C--                Next division

  332 IF (LQ(LMT+1).NE.3)          GO TO 349
      NWDIV  = 0
      JDIV   = LQ(LMT)
      LNX    = LQ(LMT+3)
      LQMTC2 = LQ(LMT+4)
      NEW = 1

C--                Next bank

  334 MODE = NEW
  335 LN   = LNX
      IF (LN.GE.LQMTC2)            GO TO 347
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLN (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 912
#endif
#if !defined(CERNLIB_QDEBUG)
      NST = JBYT(LQ(KQS+LN),1,16) - 11
      IF (NST.GT.0)  THEN
          LS  = LN + NST
          LNX = LS + IQ(KQS+LS-1) + 9
        ELSE
          LS  = LN - 8
          LNX = LN + NST + 11
        ENDIF
#endif
      NEW = JBIT(IQ(KQS+LS),IQDROP)
      IF (NEW.EQ.0)                GO TO 336
      IF (MODE.NE.0)               GO TO 335
      GO TO 344

C--                Inspect structural links of live bank

  336 K = LS - IQ(KQS+LS-2) - 1
  337 K = K + 1
      IF (K.GT.LS)                 GO TO 340
      L = LQ(KQS+K)
      IF (L.EQ.0)                    GO TO 337
      IF (JBIT(IQ(KQS+L),IQDROP).EQ.0)  GO TO 337
      KD = LQ(KQS+L+2)
  338 L  = LQ(KQS+L)
      IF (L.EQ.0)                    GO TO 339
      IF (JBIT(IQ(KQS+L),IQDROP).NE.0)  GO TO 338
      LQ(KQS+K) = L
      IF (KD.NE.K)                 GO TO 337
      LQ(KQS+L+2) = K
      GO TO 337

  339 LQ(KQS+K) = 0
      GO TO 337

  340 IF (MODE.EQ.0)               GO TO 335

C--                Start of a live group

      LQ(LQTE) = LN
      IF (LQTE.LT.LQRTE)           GO TO 334
      CALL MZTABH
      IF (IQPART.LE.0)             GO TO 334
      GO TO 261

C--                Start of a dead group

  344 NWDIV = NWDIV + (LN - LQ(LQTE))
      LQ(LQTE+1) = LN
      LQTE = LQTE + NSTEP
      GO TO 334

C--                End of division

  347 IF (MODE.NE.0)               GO TO 348
      NWDIV = NWDIV + (LN - LQ(LQTE))
      LQ(LQTE+1) = LN
      LQTE = LQTE + NSTEP

C--                Add to segment table

  348 IF (NWDIV.EQ.0)              GO TO 349
      NWBKX = NWBKX + NWDIV
      NQSEG = NQSEG + 1
      IQSEGH(1,NQSEG) = IQDN1(KQT+JDIV)
      IQSEGH(2,NQSEG) = IQDN2(KQT+JDIV)
      IQSEGD(NQSEG)   = NWDIV

C--                Step to next division

  349 LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 332

  361 IQUEST(1) = 0
      NWTABX    = LQTE - LQTA
      IF (NQSEG.EQ.1)  NQSEG=0
      NWSEGX = 3*NQSEG
      GO TO 999

C-------------------------------------------------
C-                 ERROR HANDLING
C-------------------------------------------------

  911 IQUEST(2)  = 11
      IQUEST(11) = LENTRX
      GO TO 971

  912 IQUEST(2)  = 12
      IQUEST(11) = JQSTOR
      IQUEST(12) = JDIV
      GO TO 971

  913 IQUEST(2)  = 13

  971 IF (MODTBX.EQ.0)  THEN
          IF (IOPTXP.EQ.0)       CALL ZTELL (IQUEST(2),1)
        ENDIF
      IQUEST(1) = -2

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
