*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:36  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIN (LUNP,IXDIVP,LSUPP,JBIASP,CHOPT,NUHP,IUHEAD)

C-    MAIN SEQUENTIAL INPUT ROUTINE, USER CALLED

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcseg.inc"
#include "zebra/fzcocc.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9),IXDIVP(9),LSUPP(9),JBIASP(9)
      DIMENSION    NUHP(9),IUHEAD(99)
      DIMENSION    MOPTV(6)
      EQUIVALENCE (MOPTV(1),IOPTIE)
      EQUIVALENCE (LRTYP,IDI(2))
      CHARACTER    CHOPT*(*)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIN, 4H     /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIN   /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIN    ')
#endif


#include "zebra/qtrace.inc"

      LUNNI  = LUNP(1)
      IXDIVI = IXDIVP(1)
      CALL UOPTC (CHOPT,'ERSATDFGH234',IOPTIE)
      IOPTIF = IOPTIF + IOPTIG + IOPTIH
      IF (IOPTI2(1).NE.0)  IOPTIR = 2
      IF (IOPTI2(2).NE.0)  IOPTIR = 3
      IF (IOPTI2(3).NE.0)  IOPTIR = 4

C--                Set current input unit

      IF (LUNNI.NE.LUNI)  CALL FZLOC (LUNNI,1)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.2)   THEN
          IF (LOGLVI.GE.3)  WRITE (IQLOG,9110)
          WRITE (IQLOG,9111) LUNNI,NUHP(1),MOPTV
        ENDIF
 9110 FORMAT (1X)
 9111 FORMAT (' FZIN-    Entered for LUN=',I4,' NUH=',I5,
     F' E/R/S/A/T/D= ',6I1)
#endif

      NWRDAI = IQ(KQSP+LQFI+20)
      NRECAI = IQ(KQSP+LQFI+21)
      JRETCD = 0
      JERROR = 0
      NWERR  = 0
      NQOCC  = 0
      JOPT   = IOPTIT + IOPTIA + IOPTID
      IF (JOPT.EQ.0)               GO TO 141

C--                Re-entry for pending d/s

      IF (IOPTIR.NE.0)             GO TO 714
      LRTYP  = IQ(KQSP+LQFI+27)
      NWTABI = IQ(KQSP+LQFI+41)
      NWBKI  = IQ(KQSP+LQFI+42)
      LENTRI = IQ(KQSP+LQFI+43)
      IF (JOPT.GE.2)               GO TO 711
      IF (LRTYP.GE.5)              GO TO 712
      IF (LRTYP.LE.1)              GO TO 712
      IF (IOPTIT.NE.0)             GO TO 121
      IEVFLI = 3 - LRTYP
      NWUMXI = 0
      IF (IOPTIA.NE.0)             GO TO 154
      GO TO 151

C----              Ready segment table for the user

  121 LFISEG = LQFI + JAUSEG
      NQSEG  = IQ(KQSP+LFISEG) / 3
      IF (NQSEG.NE.0)  THEN
          N = 2*NQSEG
          CALL UCOPY (IQ(KQSP+LFISEG+1),IQSEGH,N)
          CALL VZERO (IQSEGD,NQSEG)
        ENDIF
      IQSGLU = LUNI
      IQUEST(1) = 0
      GO TO 999

C-----------------------------------------------------
C-                 obtain and digest next pilot record
C-----------------------------------------------------

  141 IF (IACTVI.GE.6)             GO TO 701
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)  THEN
          ISV = IQVID(2)
          CALL ZVAUTX
          IF (IQVID(2).NE.ISV)  THEN
              IQVREM(1,3) = IQVID(1)
              IQVREM(2,3) = IQVID(2)
            ENDIF
        ENDIF
#endif
  142 NWUMXI = NUHP(1)
      NWUHI  = 0

      JPENDG = IQ(KQSP+LQFI+30)
      IQ(KQSP+LQFI+30) = 0
      IF (JPENDG.EQ.1)             GO TO 301
      IF (JPENDG.EQ.2)             GO TO 311

#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOI.EQ.0)  THEN
          CALL FZIFFN (1)
          GO TO 144
        ENDIF
#endif
      CALL FZIFFX (1)
  144 IF (JRETCD.NE.0)             GO TO 391

      IF (NWBKI.EQ.0)  IQ(KQSP+LQFI+17) = IQ(KQSP+LQFI+17) + 1

C--                check DAT record wanted / needed

      IF (IOPTIF+IPILI(3).EQ.0)    GO TO 145

      IF (IOPTIH.NE.0)  THEN
          IF (IPILI(3).NE.1)       GO TO 142
          GO TO 145
        ENDIF

      IF (IPILI(3).NE.0)  THEN
          IF (IOPTIF.EQ.0)         GO TO 142
        ENDIF

C--                Copy header vector to user

  145 IF (NWUMXI.LE.0)             GO TO 146
      NWUHI   = MIN (NWUHI, NWUMXI)
      NUHP(1) = NWUHI
      IF (NWUHI.GT.0)  CALL UCOPY (LQ(LUHEAI),IUHEAD,NWUHI)

C--                Return if selective call

  146 IF (IOPTIS.EQ.0)             GO TO 154
      IQ(KQSP+LQFI+41) = NWTABI
      IQ(KQSP+LQFI+42) = NWBKI
      IQ(KQSP+LQFI+43) = LENTRI
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9146)
 9146 FORMAT (' FZIN-    Exit for S option.')
#endif
      GO TO 991

C-------------------------------------------------
C-                 read table and d/s, relocate
C-------------------------------------------------

C--                Re-entry with D-option

  151 IF (IQSGLU.NE.LUNI)          GO TO 713
      GO TO 155

C--                Re-entry with A-option

  154 NQSEG = 0
  155 IF (NWBKI.NE.0)              GO TO 157

C--                Empty d/s

      LRTYP  = 0
      LENTRI = 0
      GO TO 189

C--                Get store / division

  157 CALL MZSDIV (IXDIVI,7)
      IF (JQDIVI.EQ.0)  JQDIVI=2
      LQSYSR(KQT+1) = LSUPP(1)

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +  WRITE (IQLOG,9157) JQSTOR,JQDIVI,NWTABI,NWBKI
 9157 FORMAT (' FZIN-    accept d/s into STORE/def.DIV.=',2I4,
     F' need memory for NWTB/NWBK=',I5,I7)
#endif

C--                Ready memory occupation table, reserve space

      CALL FZIMTB
      IF (JRETCD.NE.0)             GO TO 390

C--                Early table

      LFIEAR = LQFI + JAUEAR
      NTBE   = IQ(KQSP+LFIEAR)
      IF (NTBE.NE.0)  THEN
          LTBR  = LQTA + NWTABI
          CALL UCOPY (IQ(KQSP+LFIEAR+1),LQ(LTBR),NWTABI)
          GO TO 159
        ENDIF

C--                Read the long table

      IF (NWTABI.EQ.0)             GO TO 159
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOI.EQ.0)  THEN
          CALL FZIFFN (2)
          GO TO 158
        ENDIF
#endif
      CALL FZIFFX (2)
  158 IF (JRETCD.NE.0)             GO TO 391

C--                Read the data

  159 CONTINUE
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOI.EQ.0)  THEN
          CALL FZIFFN (3)
          GO TO 160
        ENDIF
#endif
      CALL FZIFFX (3)
  160 IF (JRETCD.NE.0)             GO TO 391

C--                Relocation

      CALL FZIREL
      IF (JRETCD.NE.0)             GO TO 391

C-------------------------------------------------------------
C-              done : link and return
C-------------------------------------------------------------

  181 JB = JBIASP(1)
      IF (JB.GE.2)  THEN
          LSUPP(1) = LENTRI
        ELSE
          LSUPP(1) = LQSYSR(KQT+1)
          CALL ZSHUNT (IXDIVI,LENTRI,LSUPP,JB,1)
        ENDIF

      IQ(KQSP+LQFI+16) = IQ(KQSP+LQFI+16) + 1

  189 LRTYP  = 8
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9189) NWTABI,NWBKI,IEVFLI
 9189 FORMAT (' FZIN-    exit for d/s with NWTB,NWBK,IEVENT=',
     F 2I6,I9,I2)
#endif
      GO TO 991

C-------------------------------------------------------------
C-              Exit to exceptions
C-------------------------------------------------------------

C--                Pending EoF

  301 JRETCD = -1
      GO TO 781

C--                Pending start/end of run

  311 LFIIOC = LQFI + JAUIOC
      NWUHI  = IQ(KQSP+LFIIOC)
      LUHEAI = LQWKFZ
      CALL UCOPY (IQ(KQSP+LFIIOC+1),LQ(LUHEAI),NWUHI)
      JRETCD = -2
      GO TO 421

C--                Side exceptions

  390 IF (JRETCD.LT.0)             GO TO 155
  391 IF (JRETCD.NE.-2)            GO TO 781

C-------------------------------------------------
C-                 end-of-file / end-of-run
C-------------------------------------------------

C----              Normal S/E-OF-RUN

  421 CALL FZIDIA
      IF (LQ(LUHEAI).GT.0)         GO TO 424

C--                End of run

      IF (IOPTIR.EQ.1)             GO TO 142
      IF (IOPTIR.GE.IACTVI)        GO TO 142
      GO TO 427

C--                Start of run

  424 IF (IOPTIR.GE.2)             GO TO 142
  427 IF (NWUMXI.LE.0)             GO TO 997
      NWUHI = MIN (NWUHI-1, NWUMXI)
      NUHP(1) = NWUHI
      IF (NWUHI.GT.0)  CALL UCOPY (LQ(LUHEAI+1),IUHEAD(1),NWUHI)
      GO TO 997

C-------------------------------------------------
C-                 ERROR CONDITIONS
C-------------------------------------------------

C--                Reading beyond end-of-data

  701 IF (IACTVI.EQ.8)             GO TO 142
      IF (IACTVI.EQ.18)            GO TO 142
      JRETCD = -3
      GO TO 781

C----              User error

C-    JERROR = 14  options (R,2,3,4) not allowed with (T,A,D)
  714 JERROR = 1

C-    JERROR = 13  no segment table for entry with D option
  713 JERROR = JERROR + 1

C-    JERROR = 12  no pending d/s for entry with T A D options
  712 JERROR = JERROR + 1

C-    JERROR = 11  multiple options T A D not allowed
  711 JERROR = 11 + JERROR
      JRETCD = 4

C------            Print error message

  781 CALL FZIDIA

C--                Abandon reserved space, if any

      IF (NQOCC.EQ.0)              GO TO 997
      DO 784  J=1,NQOCC
      JDIV = IQOCDV(J)
      NW   = IQOCSP(J)
      IF (IQMODE(KQT+JDIV).EQ.0)  THEN
          LQEND(KQT+JDIV) = LQEND(KQT+JDIV) - NW
        ELSE
          LQSTA(KQT+JDIV) = LQSTA(KQT+JDIV) + NW
        ENDIF
  784 CONTINUE
      GO TO 997

C--                Common exit

  991 IQUEST(1)  = 0
      IQUEST(5)  = IQ(KQSP+LQFI+31)
      IQUEST(6)  = IQ(KQSP+LQFI+32)
      IQUEST(11) = IEVFLI
      IQUEST(12) = IPILI(3)
      IQUEST(13) = LENTRI
      IQUEST(14) = NWBKI
      LFIIOC = LQFI + JAUIOC
      NWIOI  = IQ(KQSP+LFIIOC)
      CALL UCOPY (IQ(KQSP+LFIIOC),IQUEST(20),NWIOI+1)

      IF (NWRDAI.GE.1000000)  THEN
          IQ(KQSP+LQFI+19) = IQ(KQSP+LQFI+19) + 1
          NWRDAI = NWRDAI - 1000000
        ENDIF

  997 IQ(KQSP+LQFI+26) = IQUEST(1)
      IQ(KQSP+LQFI+27) = LRTYP
      IQ(KQSP+LQFI+2)  = IACTVI
      IQ(KQSP+LQFI+20) = NWRDAI
      IQ(KQSP+LQFI+21) = NRECAI
      IQUEST(2) = NRECAI
      IQUEST(3) = IQ(KQSP+LQFI+22)
      IQSGLU = 0

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
