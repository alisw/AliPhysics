*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:18  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZCOPY (LUNIN,LUNOUT,IEVP,CHOPT,NIOP,NUHP,IUHEAD)

C-    Control routine to copy d/s without expansion, user called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcwk.inc"
#include "zebra/mzct.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcseg.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNIN(9),LUNOUT(9),IEVP(9)
      DIMENSION    NIOP(9),NUHP(9),IUHEAD(99)
      DIMENSION    NEWOPT(3)
      CHARACTER    CHOPT*(*)
      EQUIVALENCE (LRTYP,IDI(2))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZCO, 4HPY   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZCOPY  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZCOPY   ')
#endif

#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      LUNNI   = LUNIN(1)
      LUNNX   = LUNOUT(1)
      IEVFLX  = IEVP(1)
      IOCHX(1)= NIOP(1)
      NWUHOR  = MAX (NUHP(1),0)
      NWUHX   = MIN (NWUHOR,400)
      NWFILX  = 0
      NWMEMT  = 0
      IFSENT  = 0
      ICOPYX  = 7
      CALL UOPTC (CHOPT,'ITP',NEWOPT)

C----              Set current output unit

      IF (LUNNX.NE.LUNX)  CALL FZLOC (LUNNX,2)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.2)  THEN
          IF (LOGLVX.GE.3)   WRITE (IQLOG,9110)
          WRITE (IQLOG,9111) LUNNI,LUNNX,IEVFLX,IACTVX,CHOPT
        ENDIF
 9110 FORMAT (1X)
 9111 FORMAT (' FZCOPY-  LUNin/out=',2I3,' IEVFL,IACTV,OPT=',2I3,1X,A)

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
      IF (IACTVX.GE.16)            GO TO 901

C----              Set current input unit

      IF (LUNNI.NE.LUNI)  CALL FZLOC (LUNNI,1)

C--                Check compatibility

      IF (IDAFOX.NE.IDAFOI)        GO TO 941
      IF (IFIFOX.EQ.4)             GO TO 944
      IF (IFIFOI.EQ.4)             GO TO 944
      IF (IFIFOI.EQ.0)  THEN
          IF (MAXREI.GE.NQWKTB)    GO TO 942
        ENDIF

C--                Get the parameters of the pending d/s

      LFIIOC = LQFI + JAUIOC
      LFISEG = LQFI + JAUSEG

      NWRDAI = IQ(KQSP+LQFI+20)
      NRECAI = IQ(KQSP+LQFI+21)
      LRTYP  = IQ(KQSP+LQFI+27)
      NWTABI = IQ(KQSP+LQFI+41)
      NWBKI  = IQ(KQSP+LQFI+42)
      NWSEGI = IQ(KQSP+LFISEG)
      LENTRI = IQ(KQSP+LQFI+43)
      IF (LRTYP.GE.5)              GO TO 949
      IF (LRTYP.LE.1)              GO TO 949

      JRETCD = 0
      JERROR = 0
      NWERR  = 0
#if defined(CERNLIB_FZCHANNEL)
      IF (IACMOI.EQ.3)  THEN
          IF (IACMOX.EQ.3)         GO TO 943
        ENDIF
#endif

C------            Transmit the pilot

      NWSEGX = NWSEGI
      NWTABX = NWTABI
      NWBKX  = NWBKI
      LENTRX = LENTRI

      NWTXX   = 0
      NWUHCX  = 0
      NWIOX   = 0

C--                Ready I/O characteristic

      IF (NWUHX.EQ.0)              GO TO 39
      IF (NEWOPT(1).EQ.0)  THEN
          NWIOX  = IQ(KQSP+LFIIOC)
          CALL UCOPY (IQ(KQSP+LFIIOC+1),IOCHX,NWIOX)
          GO TO 38
        ENDIF

      IF   (IOCHX(1))              34, 32, 33
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

C--                Ready text vector

      IF (NEWOPT(2).EQ.0)  THEN
          LTEXTX = LQ(KQSP+LQFI-2)
        ELSE
          LTEXTX = LQ(KQSP+LQFX-2)
        ENDIF
      IF (LTEXTX.NE.0)  CALL FZOTXT

C--                Ready segment table

      NQSEG  = NWSEGX / 3
      IF (NQSEG.NE.0)  THEN
          N = 2*NQSEG
          CALL UCOPY (IQ(KQSP+LFISEG+1),   IQSEGH,N)
          CALL UCOPY (IQ(KQSP+LFISEG+1+N), IQSEGD,NQSEG)
        ENDIF
      IQSGLU = LUNI

C--                Ready early table

C-    ISTTAB : transmission status of table words
C-      -1 not yet read,  0 read into LQ(LQTA),  +1 written

      LQTA = LQWKTB
      LQTE = LQTA + 2*NWTABI

      ISTTAB = 1
      IF (NWTABX.EQ.0)             GO TO 121
      ISTTAB = -1
#if defined(CERNLIB_FZFFNAT)
      IF (NWTABX.GE.41)            GO TO 121

      IF (IFIFOI.EQ.0)  THEN
          LFIEAR = LQFI + JAUEAR
          NTBE   = IQ(KQSP+LFIEAR)
          IF (NTBE.NE.0)  THEN
              CALL UCOPY (IQ(KQSP+LFIEAR+1),LQ(LQTA),NWTABI)
              ISTTAB = 0
              GO TO 121
            ENDIF
        ENDIF

      IF (IFIFOX.NE.0)             GO TO 121
      IF (IFIFOI.EQ.0)  THEN
          CALL FZIFFN (2)
        ELSE
          CALL FZIFFX (2)
        ENDIF
      IF (JRETCD.NE.0)             GO TO 971
      LQTA   = LQTA + NWTABI
      ISTTAB = 0
#endif

C------            Output the pilot record

  121 IDX(2) = 3
      IF (IEVFLX.NE.0)  IDX(2)=2
      NWMEMT = 20 + NWUHCX + NWSEGX + NWTXX + NWTABX + NWBKX

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.2)  WRITE (IQLOG,9121) LENTRX,NWTABX,NWBKX
 9121 FORMAT (10X,' LENTRY=',I9,' NWTAB,NWBANK=',I5,I7)

#endif
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZOFFN (IUHEAD)
          IF (IQUEST(7).EQ.1)  ISTTAB = 1
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
      IFSENT  = 7

      IACTVX  = 12
      IF (NWBKX.EQ.0)              GO TO 190

C------            Copy table and bank material

      IQ(KQSP+LQFI+16) = IQ(KQSP+LQFI+16) + 1
      IQ(KQSP+LQFX+16) = IQ(KQSP+LQFX+16) + 1

      IQUEST(7) = ISTTAB
      IDX(2) = 0
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOI.EQ.0)  THEN
          CALL FZCFFN
          GO TO 189
        ENDIF
#endif
      CALL FZCFFX
  189 IF (JRETCD.NE.0)             GO TO 971
      GO TO 191

C----              Test for pseudo end-of-tape

  190 IQ(KQSP+LQFX+17) = IQ(KQSP+LQFX+17) + 1
      LRTYP = 0
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
      IF    (NUM1-LIM1)            991, 196, 197
  196 IF (NUM2.LT.LIM2)            GO TO 991

  197 IQUEST(1) = 1
      GO TO 992

C-------------------------------------------------
C-                 ERROR HANDLING
C-------------------------------------------------

  901 IF (IACTVX.EQ.17)   CALL ZFATAM ('FZCOPY - Going beyond EoD.')
      IACTVX = 17
      GO TO 929

  903 IQUEST(11) = IOCHX(1)
      CALL ZFATAM ('FZCOPY - IOCH invalid.')

#if defined(CERNLIB_FZCHANNEL)||defined(CERNLIB_FZMEMORY)
  907 IQUEST(11) = LUNX
      CALL ZFATAM ('FZCOPY - User routine / buffer not connected.')

#endif
#if defined(CERNLIB_FZMEMORY)
  909 IQUEST(2)  = 14
      IQUEST(8)  = IQ(KQSP+LQFX+9)
      IQUEST(9)  = NWMEMT
      IQUEST(11) = LUNX
      IF (NEWOPT(3).EQ.0)  CALL ZTELL (14,1)
#endif
  929 IQUEST(1) = 2
      GO TO 998

C-    JERROR = 409  input not positioned
  949 JERROR = 4

C-    JERROR = 404  Alfa mode not allowed
  944 JERROR = JERROR + 1

C-    JERROR = 403  input/output both channel mode
  943 JERROR = JERROR + 1

C-    JERROR = 402  native input record length too long
  942 JERROR = JERROR + 1

C-    JERROR = 401  input/output different data format
  941 JERROR = 401 + JERROR
      IF (NEWOPT(3).EQ.0)  NEWOPT(3)=-1
      JRETCD = 4

C-------           Input errors

C--       write emergency stop signal

  971 IF (IFSENT.EQ.0)             GO TO 974
      IDX(2) = 9
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)  THEN
          CALL FZON1 (9,1)
          GO TO 972
        ENDIF
#endif
      CALL FZOREC
  972 IACTVX = 12

  974 CALL FZIDIA
      IQ(KQSP+LQFI+26) = IQUEST(1)
      IF (NEWOPT(3).LT.0)  CALL ZFATAL
      GO TO 997

C-------------------------------------------------
C-                 COMMON EXIT
C-------------------------------------------------

  991 IQUEST(1)  = 0
  992 IQUEST(2)  = 0
      IQUEST(5)  = IQ(KQSP+LQFX+31)
      IQUEST(6)  = IQ(KQSP+LQFX+32)
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

C--       update the input control bank
      IQ(KQSP+LQFI+26) = 0
  997 IQ(KQSP+LQFI+2)  = IACTVI
      IQ(KQSP+LQFI+20) = NWRDAI
      IQ(KQSP+LQFI+21) = NRECAI
      IQ(KQSP+LQFI+27) = LRTYP
  998 IQ(KQSP+LQFX+2)  = IACTVX
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
