*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:07  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION MZIXCO (IXAA,IXBB,IXCC,IXDD)

C-    join IXAA, IXBB, ... into composite division index
C-    ignore zero,  user called

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXAA(9), IXBB(9), IXCC(9), IXDD(9), IXV(4)
      EQUIVALENCE (IXV(1),IQUEST(11))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZIX, 4HCO   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZIXCO /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZIXCO  ')
#endif

#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_sbyt.inc"
#include "zebra/q_mbytor.inc"


      IXV(1) = IXAA(1)
      IXV(2) = IXBB(1)
      IXV(3) = IXCC(1)
      IXV(4) = IXDD(1)
      IXCOMP = 0

      DO  49  JL=1,4
      IXIN = IXV(JL)
      IF (IXIN.EQ.0)               GO TO 49
      JDV = JBYT (IXIN,1,26)
      JST = JBYT (IXIN,27,6)
      IF (JST.LT.16)               GO TO 31

C--                COMPOSITE INDEX

      JST = JST - 16
      IF (JST.GT.NQSTOR)           GO TO 91
      IF (JDV.GE.16777216)         GO TO 92
      IF (JL.NE.1)                 GO TO 24
      IXCOMP = IXIN
      JSTORU = JST
      GO TO 49

   24 IF (JST.NE.JSTORU)           GO TO 93
      IXCOMP = MBYTOR (JDV,IXCOMP,1,26)
      GO TO 49

C--                SINGLE DIVISION INDEX

   31 IF (JST.GT.NQSTOR)           GO TO 91
      IF (JDV.GE.25)               GO TO 92
      IF (JDV.EQ.0)                GO TO 92
      IF (JL.NE.1)                 GO TO 34
      IXCOMP = MSBYT (JST+16,IXCOMP,27,5)
      JSTORU = JST
      GO TO 47

   34 IF (JST.EQ.JSTORU)           GO TO 47
      IF (JST.NE.0)                GO TO 93
      IF (JDV.LT.3)                GO TO 47
      IF (JDV.LT.21)               GO TO 93

   47 IXCOMP = MSBIT1 (IXCOMP,JDV)
   49 CONTINUE

   59 MZIXCO = IXCOMP
      RETURN

C------            ERROR CONDITIONS

   93 NQCASE = 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = 7
      IQUEST(15) = JL
      IQUEST(16) = JST
      IQUEST(17) = JDV
#include "zebra/qtrace.inc"
#include "zebra/qtofatal.inc"
      MZIXCO = 0
      END
*      ==================================================
#include "zebra/qcardl.inc"
