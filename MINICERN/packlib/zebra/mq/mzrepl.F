*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:04  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZREPL (IXDIV,LIXP,CHOPT)

C-    Link replacement banks and relocate, user called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqmst.inc"
#include "zebra/mzcl.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIV(9),LIXP(9)
      CHARACTER    CHOPT*(*)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZRE, 4HPL   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZREPL /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZREPL  ')
#endif

#include "zebra/q_sbit1.inc"
#include "zebra/q_locf.inc"

#include "zebra/qtrace.inc"

      CALL MZSDIV (IXDIV,7)

      CALL UOPTC (CHOPT,'KI',IQUEST)
      IFLAG = IQUEST(1)
      ISOLA = IQUEST(2)

C-         IFLAG = 0  drop old and index
C-                 1  keep old and index

      IXGARB = 0
      LIXO  = LIXP(1)
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.LT.2)             GO TO 12
      WRITE (IQLOG,9001) JQSTOR,JQDIVI,LIXO,CHOPT
 9001 FORMAT (1X/' MZREPL-  Store/Div/Lix/CHOPT :',2I3,I8,1X,A)

#endif
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.3)             GO TO 12
      WRITE (IQLOG,9811)
 9811 FORMAT (1X/20X,'Lix   link 1        link 2')
 9824 FORMAT (5X,3I9,1X,A4,I9,1X,A4)
#endif

   12 IF (LIXO.EQ.0)               GO TO 999

C--       Find division

      IF (JQDIVI.EQ.0)  THEN
          JQDIVI = MZFDIV (-7, LQ(KQS+LIXO-1))
          IF (JQDIVI.EQ.0)         GO TO 91
        ENDIF

      LIMDLO = LQSTA(KQT+JQDIVI)
      LIMDUP = LQEND(KQT+JQDIVI)

C--       Get the number of index banks

      LIX  = LIXO
      NIX  = 0
   24 NIX  = NIX + 1
#if defined(CERNLIB_QDEBUG)

C--       check index bank valid
      CALL MZCHLS (-7,LIX)
      IF (IQFOUL.NE.0)             GO TO 93
      LOLD = LQ(KQS+LIX-1)
      LNEW = LQ(KQS+LIX-2)

      IF (IQ(KQS+LIX-3).LT.2)         GO TO 94
      IF (IFLAG.NE.0)  THEN
          IF (IQ(KQS+LIX-1).EQ.0)     GO TO 94
        ENDIF

C--       check old/new all in same division
      IF (LOLD.LT.LIMDLO)          GO TO 95
      IF (LOLD.GE.LIMDUP)          GO TO 95
      IF (LNEW.LT.LIMDLO)          GO TO 95
      IF (LNEW.GE.LIMDUP)          GO TO 95

C--       check old/new valid banks
      CALL MZCHLS (-7,LOLD)
      IF (IQFOUL.NE.0)             GO TO 96

      CALL MZCHLS (-7,LNEW)
      IF (IQFOUL.NE.0)             GO TO 97

C--       check origin-link consistent
      K = LQ(KQS+LOLD+2)
      IF (K.NE.0)  THEN
          IF (LQ(KQS+K).NE.LOLD)      GO TO 98
        ENDIF

C--       check next-link consistent
      L = LQ(KQS+LOLD)
      IF (L.NE.0)  THEN
          IF (LQ(KQS+L+2).NE.LOLD)    GO TO 99
        ENDIF

#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.3)  WRITE (IQLOG,9824) NIX,LIX,
     +                     LOLD,IQ(KQS+LOLD-4),LNEW,IQ(KQS+LNEW-4)
#endif
      LIX  = LQ(KQS+LIX)
      IF (LIX.NE.0)                GO TO 24
      IF (ISOLA.NE.0)              GO TO 61

C--       Re-link the index structure to MZ working link

      CALL MZCHNB (LIXP)
      LQSYSR(KQT+1)  = LQ(KQS+LIXO+1)
      LQSYSR(KQT+2)  = LQ(KQS+LIXO+2)
      LQMST(KQT+1)   = LIXO
      LQ(KQS+LIXO+2) = LOCF (LQMST(KQT+1)) - LQSTOR
      IF (NIX.EQ.1)                GO TO 31

C--       Sort the index banks for increasing LOLD

      LIXN = LQ(KQS+LIXO)
      IF (LQ(KQS+LIXN-1).LT.LQ(KQS+LIXO-1))  THEN
          CALL ZTOPSY (IXDIV,LIXO)
          LIXO  = LQMST(KQT+1)
        ENDIF

      CALL ZSORTI (IXDIV,LIXO,-9)
      LIXO  = LQMST(KQT+1)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.3)             GO TO 31
      WRITE (IQLOG,9811)
      LIX  = LIXO
      JIX  = 1
   27 LOLD = LQ(KQS+LIX-1)
      LNEW = LQ(KQS+LIX-2)
      WRITE (IQLOG,9824) JIX,LIX,
     +                   LOLD,IQ(KQS+LOLD-4),LNEW,IQ(KQS+LNEW-4)
      LIX  = LQ(KQS+LIX)
      JIX  = JIX + 1
      IF (LIX.NE.0)                GO TO 27
#endif

C------            Relocation

   31 MQDVGA = 0
      MQDVWI = 0
      JQSTMV = -1

      CALL MZTABM

      LMT  = LQMTA - 8
   34 LMT  = LMT + 8
      IF (LQ(LMT).NE.JQDIVI)       GO TO 34
      LQ(LMT+1) = 2

      CALL MZTABX
      LQMTE = LQMTLU

C--       construct the link relocation table

      IFIRST = 7
      LIX  = LIXO
   37 LOLD = LQ(KQS+LIX-1)
      LNEW = LQ(KQS+LIX-2)

#if !defined(CERNLIB_QDEBUG)
      IF (LOLD.LT.LIMDLO)          GO TO 95
      IF (LOLD.GE.LIMDUP)          GO TO 95
      IF (LNEW.LT.LIMDLO)          GO TO 95
      IF (LNEW.GE.LIMDUP)          GO TO 95
#endif
      NL  = IQ(KQS+LOLD-3)
      ND  = IQ(KQS+LOLD-1)
      NLC = MIN (IQ(KQS+LNEW-3), NL)
      NDC = MIN (IQ(KQS+LNEW-1), ND)

      IF (IFIRST.NE.0)  THEN
          LQ(LQTA-1) = LOLD - NLC
          LQTE  = LQTA
          IFIRST = 0
        ELSE
          LQ(LQTE)   = LQ(LQTE-3)
          LQ(LQTE+1) = LOLD - NLC
          LQ(LQTE+2) = 0
          LQ(LQTE+3) = 0
          LQTE = LQTE + 4
        ENDIF

      LQ(LQTE)   = LOLD - NLC
      LQ(LQTE+1) = LOLD + NDC + 9
      LQ(LQTE+2) = LNEW - LOLD
      LQ(LQTE+3) = 0

      LQTE = LQTE + 4
      IF (LQTE.GE.LQRTE)  THEN
          CALL MZTABH
          IF (IQPART.NE.0)         GO TO 51
        ENDIF
      LIX  = LQ(KQS+LIX)
      IF (LIX.NE.0)                GO TO 37

      LQ(LQTE) = LQ(LQTE-3)

C--       Structural replacement of old by new

      LIX  = LIXO
   42 LOLD = LQ(KQS+LIX-1)
      LNEW = LQ(KQS+LIX-2)
      LQ(KQS+LNEW)   = LQ(KQS+LOLD)
      LQ(KQS+LNEW+1) = LQ(KQS+LOLD+1)
      LQ(KQS+LNEW+2) = LQ(KQS+LOLD+2)
      LQ(KQS+LOLD)   = 0
      LQ(KQS+LOLD+1) = LIX
      LQ(KQS+LOLD+2) = LIX - 1
      IQ(KQS+LOLD-2) = 0
      IQ(KQS+LIX-2)  = 1
      IF (IFLAG.NE.0)  THEN
C--       the old banks to remain accessible
          IQ(KQS+LIX+1) = LQ(KQS+LIX-1)
        ELSE
          IQ(KQS+LOLD) = MSBIT1(IQ(KQS+LOLD),IQDROP)
        ENDIF
      LQ(KQS+LIX-1)  = 0
      LIX  = LQ(KQS+LIX)
      IF (LIX.NE.0)                GO TO 42

C--                Relocate

      CALL MZRELX

C------            Finished, reset LIX

      LIXP(1)     = LIXO
      LQ(KQS+LIXO+1) = LQSYSR(KQT+1)
      LQ(KQS+LIXO+2) = LQSYSR(KQT+2)
      LQMST(KQT+1) = 0

      IF (IFLAG.EQ.0)  THEN
          CALL MZDROP (IXDIV, LIXO, 'L')
          GO TO 999
        ENDIF

C--       restore LOLD's in the index structure

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.3)  THEN
          WRITE (IQLOG,9811)
          JIX  = 1
        ENDIF
#endif
      LIX  = LIXO
   47 LOLD = IQ(KQS+LIX+1)
      LQ(KQS+LIX-1) = LOLD
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.3)  THEN
          LNEW = LQ(KQS+LIX-2)
          WRITE (IQLOG,9824) JIX,LIX,
     +                   LOLD,IQ(KQS+LOLD-4),LNEW,IQ(KQS+LNEW-4)
          JIX  = JIX + 1
        ENDIF
#endif
      LIX  = LQ(KQS+LIX)
      IF (LIX.NE.0)                GO TO 47
#include "zebra/qtrace99.inc"
      RETURN

C----     Not enough table space, collect garbage

   51 IF (IXGARB.NE.0)             GO TO 92

      JDIVSV = JQDIVI

      IXGARB = MZIXCO (IXDIV, 21, 22, 24)
      CALL MZGARB (IXGARB, 0)

      LIXO   = LQMST(KQT+1)
      JQDIVI = JDIVSV
#if !defined(CERNLIB_QDEBUG)
      LIMDLO = LQSTA(KQT+JQDIVI)
      LIMDUP = LQEND(KQT+JQDIVI)
#endif
      GO TO 31

C--------          Isolated case, update structural links  ---------

   61 LIX  = LIXO
   62 LOLD = LQ(KQS+LIX-1)
      LNEW = LQ(KQS+LIX-2)
      LQ(KQS+LNEW)   = LQ(KQS+LOLD)
      LQ(KQS+LNEW+1) = LQ(KQS+LOLD+1)
      LQ(KQS+LNEW+2) = LQ(KQS+LOLD+2)

C----              Update according to origin-link

      K = LQ(KQS+LNEW+2)
      IF (K.NE.0)  LQ(KQS+K) = LNEW

C----              Update according to next-link

      L = LQ(KQS+LNEW)
      IF (L.NE.0)  LQ(KQS+L+2) = LNEW

C----              Update k- and up-link in vertically dependent banks

      JBIAS = IQ(KQS+LNEW-2) + 1
   64 JBIAS = JBIAS - 1
      IF (JBIAS.LE.0)              GO TO 68
      KO = LOLD - JBIAS
      KN = LNEW - JBIAS
      L  = LQ(KQS+KN)
      IF (L.EQ.0)                  GO TO 64
      IF (LQ(KQS+L+2).NE.KO)          GO TO 64
      LQ(KQS+L+2) = KN

C--                          and its linear structure

      LF = L
   67 LQ(KQS+L+1) = LNEW
      L = LQ(KQS+L)
      IF (L.EQ.LF)                 GO TO 64
      IF (L.NE.0)                  GO TO 67
      GO TO 64

   68 LQ(KQS+LOLD)   = 0
      IQ(KQS+LOLD-2) = 0
      IQ(KQS+LIX-2)  = 1
      IF (IFLAG.NE.0)  THEN
          LQ(KQS+LOLD+1) = LIX
          LQ(KQS+LOLD+2) = LIX - 1
        ELSE
          LQ(KQS+LOLD+1) = 0
          LQ(KQS+LOLD+2) = 0
          IQ(KQS+LOLD)   = MSBIT1(IQ(KQS+LOLD),IQDROP)
          LQ(KQS+LIX-1)  = 0
        ENDIF
      LIX  = LQ(KQS+LIX)
      IF (LIX.NE.0)                GO TO 62

      IF (IFLAG.EQ.0)  CALL MZDROP (IXDIV, LIXO, 'L')
      GO TO 999

C-----------       Error conditions

   99 NQCASE = 1
   98 NQCASE = NQCASE + 5
      NQFATA = 4
      IQUEST(15) = LNEW
      IQUEST(16) = LOLD
      IQUEST(17) = LQ(KQS+LOLD+2)
      IQUEST(18) = LQ(KQS+LOLD)
      GO TO 93

   97 NQCASE = 1
   96 NQCASE = NQCASE + 1
   95 NQCASE = NQCASE + 1
   94 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 5
      IQUEST(15) = IQ(KQS+LIX-3)
      IQUEST(16) = IQ(KQS+LIX-2)
      IQUEST(17) = IQ(KQS+LIX-1)
      IQUEST(18) = LOLD
      IQUEST(19) = LNEW

   93 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 3
      IQUEST(12) = LIMDLO
      IQUEST(13) = LIMDUP
      IQUEST(14) = LIX

   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(11) = JQDIVI
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
