*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:03  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZMODS(CHROUT,IFLAG)
*
************************************************************************
*
*           Check if directory can be modified
*           Enter directory in the list of modified directories
*  Input:
*    CHROUT  Name of the calling routine
* Output:
*   IFLAG   Status flag =0 if OK  1 otherwise
*
* Called by RZCOPY,RZDELK,RZDELT,RZMDIR,RZOUT,RZPASS,RZPURG,RZQUOT,RZVOUT
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.04.86
*  Last mod: 04.10.90
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
      CHARACTER*255 CHPATH
      CHARACTER*6 CHROUT
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

*
*            Check permission
*
      IF(JBIT(IQ(KQSP+LCDIR),1).NE.0)GO TO 90
      IF(IFLAG.EQ.0)THEN
         IF(LFREE.EQ.0)GO TO 90
         IF(IQ(KQSP+LFREE+1).EQ.0)GO TO 90
      ELSE
         IFLAG=0
      ENDIF
*
      CALL SBIT1(IQ(KQSP+LCDIR),2)
      CALL SBIT1(IQ(KQSP+LTOP),2)
      GO TO 99
*
  90  CONTINUE
*
*     Get current directory name
*
      CALL RZPAFF(CHPAT,NLPAT,CHPATH)
      LPATH = LENOCC(CHPATH)
      IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) THEN
         WRITE(IQLOG,1000)CHROUT
         WRITE(IQLOG,*) '(',CHPATH(1:LPATH),')'
      ENDIF
      IFLAG=1
      IQUEST(1)=3
 1000 FORMAT(1X,A,'. No authorisation to write in that directory')
*
  99  RETURN
      END
