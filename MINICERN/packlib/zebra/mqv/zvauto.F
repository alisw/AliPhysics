*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:23  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVAUTO (MSTORE,IDVTH1,IDVTH2,LOGLEV)

C-    Initiate automatic verification

#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/zvfaut.inc"
#include "zebra/zvfwkc.inc"
C--------------    END CDE                             --------------

#include "zebra/q_jbyt.inc"

      IF (IQVID(2).NE.0)           RETURN

      IQVSTA = -1
      MASKST = MSTORE
      IQVLOG = LOGLEV
      IQVTHR(1) = IDVTH1
      IQVTHR(2) = IDVTH2

      IF (IQVTHR(1).LT.0)  IQVSTA = 1
      IF (IQVTHR(2).LT.0)  IQVSTA = 1

      IF (MSTORE.LE.0)  THEN
          IQUEST(1) = -1
          MASKST = JBYT (IQUEST(1),2,16)
        ENDIF

      WRITE (IQLOG,9001) MASKST,IQVTHR,IQVLOG
 9001 FORMAT (1X/' ZVAUTO.  Automatic verification initiated',
     F' with stores mask= ',Z4,','
     F/10X,'Threshold VID=',I5,I10,', Log level=',I3)

      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
