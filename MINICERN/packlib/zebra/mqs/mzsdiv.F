*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:09  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZSDIV (IXDIVP,IFLAGP)

C-    Set current store JQSTOR and maybe division JQDIVI from IXDIVP
C-    IFLAG = -ve :  store only
C-           else :  set also JQDIVI:
C-                   to the division if a specific single division,
C-                   or JQDIVI = 0 if no such given, 
C-                with these restrictions:
C-            > 0 :  IXDIVP may not be a compound or generic index
C-            = 4 :  specific single division required, =0 not allowed
C-    system called, could be user called

#include "zebra/zstate.inc"
#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXDIVP(9), IFLAGP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZSD, 4HIV   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZSDIV /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZSDIV  ')
#endif

#include "zebra/q_jbyt.inc"


      IXIN  = IXDIVP(1)
      IFLAG = IFLAGP(1)
      JSTO  = JBYT (IXIN,27,4)
      IF (JSTO.NE.JQSTOR)          GO TO 41
      IF (IFLAG.LT.0)              GO TO 48

   21 JDIV = JBYT (IXIN,1,26)
#if defined(CERNLIB_B32)
      JCOM = JBYT (IXIN,31,2)
#endif
#if defined(CERNLIB_B36)
      JCOM = JBYT (IXIN,31,6)
#endif
#if defined(CERNLIB_B48)
      JCOM = JBYT (IXIN,31,18)
#endif
#if defined(CERNLIB_B60)
      JCOM = JBYT (IXIN,31,30)
#endif
#if defined(CERNLIB_B64)
      JCOM = JBYT (IXIN,31,34)
#endif
      IF   (JCOM-1)          22, 31, 91

C--                SINGLE DIVISION INDEX

   22 IF (JDIV.GE.25)              GO TO 92
      IF (JDIV.GE.21)              GO TO 24
      IF (JDIV.GT.JQDVLL)  THEN
          IF (JDIV.LT.JQDVSY)      GO TO 92
        ENDIF
      IF (JDIV.EQ.0)  THEN
          IF (IFLAG.EQ.4)          GO TO 94
        ENDIF
      JQDIVI = JDIV
      RETURN

   24 IF (JDIV.EQ.24)              GO TO 26
      IF (IFLAG.GT.0)              GO TO 93
      JQDIVI = 0
      RETURN

   26 JQDIVI = JQDVSY
      RETURN

C--                COMPOUND INDEX

   31 IF (IFLAG.GT.0)              GO TO 93
      IF (JDIV.GE.16777216)        GO TO 92
      JQDIVI = 0
      RETURN

C----              SWITCH ZEBRA STORE

   41 IF (JSTO.GT.NQSTOR)          GO TO 91
      JQSTOR = JSTO
      JQDIVR = 0
      KQT = NQOFFT(JQSTOR+1)
      KQS = NQOFFS(JQSTOR+1)

      DO 44  J=1,12
   44 IQCUR(J) = IQTABV(KQT+J)
      NQLOGM = NQLOGL
      IF (IFLAG.GE.0)              GO TO 21
   48 JQDIVI = 0
      RETURN

C------            ERROR CONDITIONS

   94 NQCASE = 1
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
      NQFATA = 1
      IQUEST(14) = JDIV
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 3
      IQUEST(11) = IXIN
      IQUEST(12) = IFLAG
      IQUEST(13) = JSTO
#include "zebra/qtrace.inc"
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
