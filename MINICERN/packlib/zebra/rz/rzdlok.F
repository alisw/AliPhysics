*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZDLOK
*
************************************************************************
*
*      To reset the LOCKing word in first record.
*       This is usefull, in case a previous job could not finish
*       in a clean way (time limit, system crash,CRTL/Y,etc
*
* Called by <USER>,RZFILE
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.10.86
*  Last mod: 02.04.90
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      IF(LTOP.EQ.0)GO TO 99
*
      IF(LRIN.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LRIN,LTOP,-7,'RZIN',0,0,LREC+1,2,-1)
         IQ(KQSP+LRIN-5)=IQ(KQSP+LTOP-5)
      ENDIF
      NWL =50
  10  CALL RZIODO(LUN,NWL,1,IQ(KQSP+LRIN+1),1)
      IF(IQ(KQSP+LRIN+2).GT.NWL.AND.IQUEST(1).EQ.0)THEN
         NWL=IQ(KQSP+LRIN+2)
         GO TO 10
      ENDIF
      NWL=IQ(KQSP+LRIN+2)
      IF(LCDIR.NE.0)IQ(KQSP+LTOP+KIRIN)=0
      IQ(KQSP+LRIN+3)=0
      CALL RZIODO(LUN,MAX(NWL,50),1,IQ(KQSP+LRIN+1),2)
*
  99  RETURN
      END
