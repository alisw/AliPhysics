*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:07  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZLOG (CHOPT)

C-    Set log options

C-    Structure of JQLLEV :
C-        bits  1-4  =  0      Q, suppress all messages
C-                      1      E, error messages only
C-                      2      normal logging

C-        bits  5-8  =  1      T, JZTELL printing

C-        bits  9-12 :  bit 1  A, JZIN logging
C                           2  B, and dump call bank
C-                          3  C, and dump AN parameters

C-        bits 13-16 :  bit 1  X, JZOUT logging
C-                          2  Y, and dump call bank

#include "zebra/zunit.inc"
#include "zebra/quest.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      CHARACTER    CHOPT*(*)

#include "zebra/q_or.inc"
#include "zebra/q_shiftl.inc"

#if defined(CERNLIB_QPRINT)
      IF (JQLLEV.GE.2)  WRITE (IQLOG,9007) CHOPT
 9007 FORMAT (' JZLOG.   option = ',A)
#endif

      CALL UOPTC (CHOPT,'QETABCXY',IQUEST)
C-                       _:.=+=.:
      JQLLEV = 2
      IF (IQUEST(1).NE.0)  JQLLEV = 0
      IF (IQUEST(2).NE.0)  JQLLEV = 1
      IF (JQLLEV.LT.2)             RETURN

      IF (IQUEST(3).NE.0)  JQLLEV = 18

      J = IQUEST(5) + 2*IQUEST(6)
      IF (J.EQ.0)  THEN
          J = IQUEST(4)
        ELSE
          J = 2*J + 1
        ENDIF
      IF (J.NE.0)  JQLLEV = IOR (JQLLEV, ISHFTL(J,8))

      J = IQUEST(8)
      IF (J.EQ.0)  THEN
          J = IQUEST(7)
        ELSE
          J = 2*J + 1
        ENDIF
      IF (J.NE.0)  JQLLEV = IOR (JQLLEV, ISHFTL(J,12))
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
