*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:47  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZEND(CHDIR)
*
************************************************************************
*
*        To close all transactions with file CHDIR.
*          Corresponding directories are dropped
* Input:
*   CHDIR   Character variable specifying  the name of the  top directory
*           to be associated with unit LUN.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 07.12.92 Loglevel manipulation
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/zstate.inc"
#include "zebra/rzcl.inc"
      CHARACTER  CHDIR*(*)
      DIMENSION IHDIR(4)
      LOGICAL RZSAME
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
*           Save current directory
*
      CALL RZSAVE
*
      NCHD=LEN(CHDIR)
      IF(NCHD.GT.16)NCHD=16
      CALL VBLANK(IHDIR,4)
      CALL UCTOH(CHDIR,IHDIR,4,NCHD)
      CALL ZHTOI(IHDIR,IHDIR,4)
*
      LRZ=LQRS
  10  IF(LRZ.NE.0)THEN
         IF(.NOT.RZSAME(IHDIR,IQ(KQSP+LRZ+1),4))THEN
            LRZ=LQ(KQSP+LRZ)
            GO TO 10
         ENDIF
*
      LTOP=LRZ
      LOGLV = JBYT(IQ(KQSP+LTOP),15,3)-3
      IF(LOGLV.GE.0) WRITE(IQLOG,9019) CHDIR
 9019 FORMAT(' RZEND. called for ',A)
*
*           Release lock if update mode
*
         IF(JBIT(IQ(KQSP+LTOP),3).NE.0)THEN
            LCDIR=LTOP
            CALL RZFREE('RZFILE')
         ENDIF
*
*           Drop control bank
*
         CALL MZDROP(JQPDVS,LTOP,' ')
         LTOP = 0
         LCDIR= 0
*
*           Error
*
      ELSEIF(NQLOGD.GE.-2)THEN
         WRITE(IQLOG,1000) CHDIR
 1000    FORMAT(' RZEND. Unknown directory ',A)
      ENDIF
*
  99  RETURN
      END
