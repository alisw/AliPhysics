*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:25  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZEBRA (LIST)

C-    Initialise memory manager,  user called

#include "zebra/zbcd.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zceta.inc"
#include "zebra/zheadp.inc"
#include "zebra/zmach.inc"
#include "zebra/znatur.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/mzcl.inc"
#include "zebra/jzuc.inc"
#include "zebra/rzcount.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzcstr.inc"
C--------------    End CDE                             --------------
      DIMENSION    LIST(9), INKEYS(3)
      CHARACTER    VIDQQ*(*)

#if defined(CERNLIB_QMVDS)
      SAVE         INIT, INKEYS
#endif
#if defined(CERNLIB_PY_VS5)
*      PARAMETER (VIDQQ= '@(#)? >')
#endif
#if !defined(CERNLIB_PY_VS5)
      PARAMETER   (VIDQQ = '@(#)' //
     +  'ZEMQ     3.77 '
     + // '>')
#endif
#include "zebra/zfatalch.inc"
#if defined(CERNLIB_QHOLL)
      DATA  INKEYS / 4HEBRA, 4HINIT, 4HDONE /
#endif
#if !defined(CERNLIB_QHOLL)
      DATA  INKEYS / 2345432, 3456543, 4567654 /
#endif


C--                Check re-entry to MZEBRA
C
C     DO 11  J=1,3
C     IF (MQKEYS(J).NE.INKEYS(J))  GO TO 12
C  11 CONTINUE
C     RETURN
      CALL CLEFT (VIDQQ,1,0)

   12 NQSTOR = -1
      JQSTOR = -99
      CALL VZERO (NQOFFT,66)

C--                Initialize all passive COMMONs

      CALL MZINCO (LIST)
      NQDCUT = 201
      NQWCUT = 500

      CALL UCOPY (INKEYS,MQKEYS,3)

C----              Clear /MZCL/, /JZUC/, /RZCOUNT/, /RZCSTR/

      CALL VZERO (NQLN, 28)
      CALL VZERO (LQJZ, 16)
      CALL VZERO (NSTRIP, MAXFILES)
      JQLEV = -1

      RZXIO(1) = 0.
      RZXIO(2) = 0.

      IMODEH   = 0

      CALL VFILL (IQFENC,4,IQNIL)


C----              Test validity of ZFATAL

      NQINIT = INIT
      CALL ZFATAL
      CALL ZFATAM ('.')
      NQINIT = -1
C-        force loading for "what"
      IF (NQSTOR.EQ.99)  CALL MZVERS
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
