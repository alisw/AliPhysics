*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:14:52  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZFREE

C--   Read 1 title bank in free field format for TZINIT

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
#include "zebra/tzuc.inc"
#include "zebra/tzc1.inc"
C--------------    END CDE                             -----------------
      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT,NUM(2),DUMMY(34)
      CHARACTER    COL(LGL)*1
      EQUIVALENCE (COL(1), LINE(1:1))

      CHARACTER    CHTYP*(*), FAULT*20
      PARAMETER   (CHTYP = '#''":')


      IF (JCOLA.EQ.0)  JCOLA = 1
      IF (JCOLE.EQ.0)  JCOLE = LGL
      IF (NCHPW.EQ.0)  NCHPW = 4
      IPRHEA = IFLLOG
      NREADY = 7
      MULT   = 0
      LPUTX  = LPUTA

C--       Print delimitation, unless full line

      IF (IFLPRI.EQ.0)             GO TO 21
      IF (JCOLE.EQ.LGL)  THEN
          IF (JCOLA.EQ.1)          GO TO 21
          WRITE (IQLOG,9016)  BLANK(1:JCOLA)
        ELSE
          NB = JCOLE-1 - JCOLA
          IF (NB.GE.3)  THEN
              WRITE (IQLOG,9016)  BLANK(1:JCOLA), BLANK(1:NB-2)
            ELSE
              WRITE (IQLOG,9017)  BLANK(1:JCOLA), BLANK(1:NB)
            ENDIF
        ENDIF

 9016 FORMAT (5X,A,'<-',A,'->')
 9017 FORMAT (5X,A,'<',A,'>')

C----              Read line by line

   21 CONTINUE
#include "zebra/tzread1.inc"

      IF (COL(1).EQ.'*')           RETURN

      IF (IFLPRI.NE.0)  WRITE (IQLOG,9001) LINE(1:NCHORG)
 9001 FORMAT (6X,A)


C------            Crack field by field

      JSTATE = 0
      JTKEND = MIN (NCHORG,JCOLE)
      JTKE   = JCOLA - 1
   24 J  = ICNEXT (LINE,JTKE+1,JTKEND)
      IF (J.GT.JTKEND)  THEN
          IF (MULT.NE.0)           GO TO 82
          GO TO 21
        ENDIF

      JTKA = J
      JTKE = NESLAT - 1
      NTK  = NDSLAT

      JTYP = INDEX (CHTYP,COL(JTKA)) + 1
      GO TO (61, 31, 41, 41, 43), JTYP
C-            1   2   3   4   5
C-      numeric   #   '   "   :

C----              # item

   31 IF (NTK.LT.2)                GO TO 81
      JIT = INDEX ('.ADNadnOo0Xx', COL(JTKA+1))
C-                  123456789012
      IF (JIT.EQ.0)                GO TO 81
      IF (JIT.GE.8)                GO TO 61
      IF (JIT.GE.5)  JIT = JIT - 3
      IF (MULT.NE.0)               GO TO 82

C--       #. comment #
      IF (JIT.EQ.1)  THEN
          JTKE = ICFIND ('#', LINE,JTKA+2,JTKEND)
          GO TO 24
        ENDIF

C--       #AnCW  control item
      IF (JIT.EQ.2)  THEN
          CALL TZACW (LINE(JTKA+1:JTKE))
          IF (NCHPW.GE.0)          GO TO 24
          NCHPW = 4
          GO TO 81
        ENDIF

C--       #Double control item
      IF (JIT.EQ.3)  THEN
          IFLDBL = 1
          GO TO 24
        ENDIF

C--       #Normal control item
      IFLDBL = 0
      GO TO 24

C----              Handling Hollerith

   41 JTKE = ICFIND (COL(JTKA), LINE,JTKA+1,JTKEND)
      IF (JTKE.GT.JTKEND)          GO TO 83
      NTK  = JTKE - JTKA

   43 JTKA = JTKA + 1
      NTK  = NTK  - 1

      NWDH = (NTK-1) / NCHPW + 1
      NWDS = NWDH + IFLHC + IFLHW
      NWDT = NWDS
      IF (MULT.NE.0)  NWDT = NWDT * MULT
      IF (LPUTX+NWDT.GT.LPUTE)     GO TO 84

      JST = LPUTX
      IF (IFLHC.NE.0)  THEN
          LQ(JST) = NTK
          JST = JST + 1
        ENDIF

      IF (IFLHW.NE.0)  THEN
          LQ(JST) = NWDH
          JST = JST + 1
        ENDIF

      CALL UCTOH (LINE(JTKA:JTKE), LQ(JST),NCHPW, NTK)

      JST   = LPUTX
      LPUTX = LPUTX + NWDT
      IF (MULT.EQ.0)               GO TO 24

      MULT = MULT - 1
      CALL UCOCOP (LQ(JST),LQ(JST+NWDS),MULT,NWDS,0,NWDS)
      MULT = 0
      GO TO 24

C----              Handle numeric

   61 CALL CKRACK (LINE,JTKA,JTKE,IFLDBL)
      IF (NFSLAT.LE.0)             GO TO 86
      IF (NGSLAT.NE.0)             GO TO 67
      NWDS = MAX (1, NFSLAT-2)
      NWDT = NWDS
      IF (MULT.NE.0)  NWDT = NWDT * MULT
      IF (LPUTX+NWDT.GT.LPUTE)     GO TO 84

      JST   = LPUTX
      LPUTX = LPUTX + NWDT
      IF (MULT.NE.0)               GO TO 64

      LQ(JST) = NUM(1)
      IF (NWDS.EQ.1)               GO TO 24
      LQ(JST+1) = NUM(2)
      GO TO 24

   64 IF (NWDS.EQ.1)  THEN
          CALL VFILL (LQ(JST),MULT,NUM(1))
        ELSE
          CALL UCOCOP (NUM,LQ(JST),MULT,NWDS,0,NWDS)
        ENDIF
      MULT = 0
      GO TO 24

C--       handle repeat*

   67 JE = NESLAT
      IF (COL(JE).NE.'*')          GO TO 86
      IF (NFSLAT.GE.3)             GO TO 81
      IF (MULT.NE.0)               GO TO 82
      IF (NUM(1).LE.1)             GO TO 85
      MULT = NUM(1)
      JTKE = JE
      GO TO 24

C----              EoF seen

   78 NREADY = -7
      RETURN

C----              Error handling

C-        81 :  nothing special
C-        82 :  pending repeat not allowed
C-        83 :  missing terminator " or '
C-        84 :  more data than expected
C-        85 :  invalid repeat count
C-        86 :  invalid numeric

   86 JSTATE = 1
      JTKA = NESLAT
   85 JSTATE = JSTATE + 1
   84 JSTATE = JSTATE + 1
   83 JSTATE = JSTATE + 1
   82 JSTATE = JSTATE + 1
   81 JSTATE = JSTATE + 1

      IF (IFLPRI.EQ.0)  THEN
          IF (IPRHEA.EQ.0)  THEN
              WRITE (IQLOG,9012) LHEAD(1:NHEAD)
              IPRHEA = 7
            ENDIF
          WRITE (IQLOG,9001) LINE(1:JTKEND)
        ENDIF

      IF (JSTATE.EQ.1)  THEN
           FAULT = 'fault'
        ELSEIF (JSTATE.EQ.2)  THEN
           FAULT = 'pending repeat'
        ELSEIF (JSTATE.EQ.3)  THEN
           FAULT = 'missing terminator'
        ELSEIF (JSTATE.EQ.4)  THEN
           FAULT = 'too much data'
        ELSEIF (JSTATE.EQ.5)  THEN
           FAULT = 'invalid repeat count'
        ELSEIF (JSTATE.EQ.6)  THEN
           FAULT = 'invalid'
        ENDIF

      WRITE (IQLOG,9091) BLANK(1:JTKA),FAULT
      NFAULT = NFAULT + 1
      MULT   = 0
      IF (JSTATE.NE.4)             GO TO 21
      NREADY = 0
      RETURN

 9012 FORMAT (3X,' > ',A)
 9091 FORMAT ('  !!f',A,'^-> !!! ',A)
      END
*      ==================================================
#include "zebra/qcardl.inc"
