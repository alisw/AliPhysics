*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:19  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZCOPY (IXDVFR,LENTP,IXDVTO,LSUPP,JBIASP,CHOPT)

C-    Copy a data-structure, User called

#include "zebra/zbcd.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcn.inc"
#include "zebra/fzcx.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXDVFR(9),LENTP(9),IXDVTO(9),LSUPP(9),JBIASP(9)
      CHARACTER    CHOPT*(*)
      DIMENSION    LADESV(6)
#if defined(CERNLIB_QMVDS)
      SAVE         LADESV
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZCO, 4HPY   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZCOPY /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZCOPY  ')
#endif
      DATA  LADESV / 6, 5*0 /

#include "zebra/q_locf.inc"

#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      LENTRX = LENTP(1)
      LSUP   = LSUPP(1)
      NWBKFR = IQUEST(1)
      NOFFFR = IQUEST(2)
      NWBKMX = IQUEST(3)
      NOFFTO = IQUEST(4)
      CALL UOPTC (CHOPT,'DISZMLNPFT',IOPTXD)
      JFLGAX = 0
      IHANDL = 0
      NWBKU  = -7
      LOGLEV = NQLOGD

C----              TO option  yes / no

      IF (IOPTXT.NE.0)  THEN
          LOCTO  = LOCF(IXDVTO(1)) - 1 - NOFFTO
          KSQTO  = LOCTO - LQASTO
          JSTOTO = -2
        ELSE
          IXDIVI = IXDVTO(1)
          CALL MZSDIV (IXDIVI,7)
          IF (JQDIVI.EQ.0)         GO TO 97
          KSQTO  = KQS
          JSTOTO = JQSTOR
          JDIVTO = JQDIVI
          JDMODE = IQMODE(KQT+JQDIVI)
          NWBKMX = NQDMAX(KQT+JQDIVI)
          LOGLEV = MAX (LOGLEV,NQLOGL)
          LQSYSR(KQT+2) = LSUP
        ENDIF

C----              FROM option  yes / no

      IF (IOPTXF.NE.0)  THEN
          LOCFR  = LOCF(IXDVFR(1)) - 1 - NOFFFR
          KSQFR  = LOCFR - LQASTO
          JSTOFR = -1
          NWBKX  = NWBKFR
          IOPTXD = 1
          IOPTXL = 1
        ELSE
          IXDIVX = IXDVFR(1)
          CALL MZSDIV (IXDIVX,0)
          KSQFR  = KQS
          JSTOFR = JQSTOR
          LOGLEV = MAX (LOGLEV,NQLOGL)
          LQSYSR(KQT+1) = LENTRX
          IF ((JSTOFR.EQ.JSTOTO).AND.(IOPTXZ.EQ.0))  IHANDL=-1
        ENDIF
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLEV.GE.2)
     + WRITE (IQLOG,9009) JSTOFR,JSTOTO,CHOPT
 9009 FORMAT (' MZCOPY-  Store From/To =',2I3,' Options = ',A)
#endif

C----------        Table building         ----------------------

   11 MODTBX = 1
      IQPART = 0

C----              TO option  yes / no

      JQSTMV = -1
      IF (IOPTXT.NE.0)             GO TO 21

C--                TO option  no

      IF (NWBKU.GE.0)       GO TO 21
      JQSTMV = JSTOTO
      JQDVM1 = JDIVTO
      JQDVM2 = JDIVTO
      NQDVMV = 1 - 2*JDMODE
      IQTNMV = -7

C----              FROM option  yes / no

   21 IF (IOPTXF.EQ.0)             GO TO 24

      JQGAPM = 0
      JQGAPR = 0
      LQRTA  = LQWKTB
      LQTA   = LQRTA + 2
      LQTE   = LQTA  + 4
      LQRTE  = LQTE  + 1
      LQMTA  = LQRTE + 2

      LQ(LQTA)   = NOFFFR + 1
      LQ(LQTA+1) = LQ(LQTA) + NWBKX
      GO TO 33

C--                FROM option  no

   24 IF (JQSTOR.NE.JSTOFR)  CALL MZSDIV (IXDIVX,0)

      CALL FZOTAB
      IF (IQUEST(1).EQ.0)          GO TO 31
      IF (IQUEST(1).EQ.1)          GO TO 11
      IF (IQUEST(2).NE.13)         GO TO  91
      IF (JFLGAX.GE.3)             GO TO  93

C--                Collect garbage in the TO space, if no room

      JFLGAX = 3
      IF (JSTOTO.EQ.JQSTOR)        GO TO  93
      IF (JSTOTO.LE.0)             GO TO  93

      IXGARB = MZIXCO (IXDIVI,21,22,23)
      IXGARB = MZIXCO (IXGARB,24,0,0)
      CALL MZGARB (IXGARB,0)
      GO TO 11

C----------        Reserve target space     --------------------

   31 CONTINUE
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.7)             GO TO 33
      WRITE (IQLOG,9831) LQMTA,LQMTE
 9831 FORMAT (1X/' DEVZE MZCOPY,  Memory Occupation Table,'
     F,' LQMTA,LQMTE= ',2I8/16X,
     F' DIV ACT     NWS    LFBK  LLBK+1     JFE     JLE    FREE')

      WRITE (IQLOG,9832) (LQ(J),J=LQMTA,LQMTE-1)
 9832 FORMAT (16X,2I4,6I8)
#endif
   33 IF (NWBKX.GT.NWBKMX)         GO TO  94
      IF (NWBKX.EQ.0)              GO TO  95
      IF (IOPTXT.NE.0)             GO TO 41
      IF (JQSTOR.NE.JSTOTO)  CALL MZSDIV (IXDIVI,0)
      IF (NWBKU.GE.0)       GO TO 37

      JQDIVI = JDIVTO
      CALL MZRESV
      NQRESV = NQRESV - NWBKX
      IF (NQRESV.LT.0)  THEN
          CALL MZGAR1
          IQPART = -7
          IF (JSTOTO.EQ.JSTOFR)  LENTRX=LQSYSR(KQT+1)
        ENDIF

      IF (JDMODE.EQ.0)  THEN
          LTO  = LQEND(KQT+JDIVTO)
          LTOE = LTO + NWBKX
          LQEND(KQT+JDIVTO) = LTOE
        ELSE
          LTOE = LQSTA(KQT+JDIVTO)
          LTO  = LTOE - NWBKX
          LQSTA(KQT+JDIVTO) = LTO
        ENDIF
      IF (IQPART.EQ.0)             GO TO 42

      LQ(KQS+LTO)   = 12
      LQ(KQS+LTO+1) = 0
      LQ(KQS+LTO+2) = 0
      LQ(KQS+LTO+3) = 0
      LQ(KQS+LTO+4) = 0
      LQ(KQS+LTO+5) = IQLETT(1)
      LQ(KQS+LTO+6) = 0
      LQ(KQS+LTO+7) = 0
      LQ(KQS+LTO+8) = NWBKX - 10
      LQ(KQS+LTO+9) = 0
      NWBKU = NWBKX
      GO TO 11

C----              Target space already reserved, get it

   37 NDIF  = NWBKX - NWBKU
      IF (JDMODE.EQ.0)  THEN
          LTOE = LQEND(KQT+JDIVTO)
          LTO  = LTOE - NWBKU
          IF (NDIF.EQ.0)           GO TO 42
          LTOE = LTOE + NDIF
          LQEND(KQT+JDIVTO) = LTOE
        ELSE
          LTO  = LQSTA(KQT+JDIVTO)
          LTOE = LTO + NWBKU
          IF (NDIF.EQ.0)           GO TO 42
          LTO = LTO - NDIF
          LQSTA(KQT+JDIVTO) = LTO
        ENDIF
      GO TO 42

C----------        Move material to destination    -------------

C--                TO option yes

   41 LTO    = NOFFTO + 1
      LTOE   = LTO + NWBKX
      JQSTOR = -1
      LQSTOR = LOCTO
      KQS    = KSQTO
      KQT    = KQFT

      LQFSTA(1)  = LTO
      LQFSTA(21) = LTOE
      GO TO 44

C--                TO option no

   42 NWBKU = NWBKX
      LSUP  = LQSYSR(KQT+2)

   44 LTOA = LTO
      NREL = LTO
      LFRE = 0
      LTB  = LQTA

   46 LFR  = LQ(LTB)
      NREL = NREL - (LFR-LFRE)
      LFRE = LQ(LTB+1)

      LQ(LTB+2) = NREL
      LQ(LTB+3) = IHANDL

      NW  = LFRE - LFR
      CALL UCOPY (LQ(KSQFR+LFR), LQ(KSQTO+LTO), NW)
      LTO = LTO + NW
      LTB = LTB + 4
      IF (LTB.LT.LQTE)             GO TO 46

C----------        Relocation        ---------------------------

      LQ(LQMTA+1) = 1
      LQ(LQMTA+3) = LTOA
      LQ(LQMTA+4) = LTOE
      LQMTE = LQMTA + 8
      IF (IOPTXN.NE.0)             GO TO 61

      LQ(LQTE)   = LQ(LQTE-3)
      LQ(LQTA-1) = LQ(LQTA)

      IF (IHANDL.LT.0)  THEN
C--          retain links pointing outside the d/s
          IQFLIO = 0
        ELSE
C--          zero links pointing outside
          IQFLIO = 7
        ENDIF
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.7)             GO TO 56
      LQ(LQMTA)   = 0
      LQ(LQMTA+2) = 0
      LQ(LQMTA+5) = 0
      LQ(LQMTA+6) = 0
      LQ(LQMTA+7) = 0

      WRITE (IQLOG,9853) LQRTA,LQTA,LQTE,LQ(LQTA-1)
 9853 FORMAT (1X/' DEVZE MZCOPY,  Link Relocation Table,'
     F,' LQRTA,LQTA,LQTE= ',3I8
     F/16X,'   LOC       L      LE    NREL  BG'
     F/28X,I10)

      IF (LQTE.LE.LQTA)            GO TO 55
      I  = LQRTA
      JA = LQTA - I
      JE = LQTE - I - 1
      WRITE (IQLOG,9854) (J,LQ(I+J),LQ(I+J+1),LQ(I+J+2),LQ(I+J+3),
     +                                           J=JA,JE,4)
 9854 FORMAT (16X,I6,3I8,I4)

   55 WRITE (IQLOG,9855) LQ(LQTE)
 9855 FORMAT (20X,I10)
   56 CONTINUE
#endif
      CALL MZRELB
      IF (IQFLIO.LT.0)             GO TO  96

C--                Relocate the entry link

      LADESV(2) = LOCF(LENTRX) - LQSTOR
      LADESV(3) = LADESV(2) + 1
      LADESV(5) = IQLETT(9)
      LADESV(6) = IQLETT(15)
      CALL MZRELL (LADESV)

      IF (IOPTXL.EQ.0)  LQ(KQS+LENTRX)=0
      LQ(KQS+LENTRX+1) = 0
      LQ(KQS+LENTRX+2) = 0
      GO TO 64

C--                No-link option

   61 CALL FZILIN
      IF (IQFOUL.NE.0)             GO TO  96
      LENTRX = IQUEST(1)

C----              Connect d/s

   64 IF (IOPTXT.NE.0)  THEN
          LSUPP(1)  = LENTRX
          JBIASP(1) = NWBKX
        ELSE
          JB = JBIASP(1)
          LSUPP(1) = LQSYSR(KQT+2)
          CALL ZSHUNT (IXDIVI,LENTRX,LSUPP,JB,1)
        ENDIF

      IQUEST(1) = 0
      IQUEST(2) = NWBKX
#include "zebra/qtrace99.inc"
      RETURN

C----------        Error conditions      -----------------------

C--                LENTRX invalid / Bank chaning clobbered in FZOTAB

   91 IQUEST(2) = IQUEST(2) - 10
      GO TO 98

C--                Not enough table space

   93 IQUEST(2) = 3
      GO TO 98

C--                D/s larger than the target space

   94 IQUEST(2) = 4
      IQUEST(11) = NWBKX
      IQUEST(12)= NWBKMX
      GO TO 98

C--                D/s empty

   95 IQUEST(2) = 5
      GO TO 98

C--                MZRELB or FZILIN find bank chaining clobbered

   96 IQUEST(2) = 6
      GO TO 98

C--                Target division not specified

   97 IQUEST(2) = 7

   98 CONTINUE
      IF (NWBKU.GE.0)  THEN
          IF (JQSTOR.NE.JSTOTO)  CALL MZSDIV (IXDIVI,0)
          IF (JDMODE.EQ.0)  THEN
              LQEND(KQT+JDIVTO) = LQEND(KQT+JDIVTO) - NWBKU
            ELSE
              LQSTA(KQT+JDIVTO) = LQSTA(KQT+JDIVTO) + NWBKU
            ENDIF
        ENDIF

      IF (IOPTXP.EQ.0)  CALL ZTELL (15,1)
      IQUEST(1) = IQUEST(2)
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
