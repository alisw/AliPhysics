*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZHOOK(LUNIN,CHDIR,TARGET,LRECL,CHOPT)
*
*     Connect user routine for RZFILE on unit LUNIN
*
      CHARACTER*(*) CHOPT,CHDIR
      CHARACTER*20  CHOP
 
      EXTERNAL   TARGET
#include "zebra/quest.inc"
 
      IQUEST(10) = LRECL
      IQUEST(11) = JUMPAD(TARGET)
 
      LCH        = LENOCC(CHOPT)
      CHOP       = CHOPT(1:LCH)
      IF(INDEX(CHOP(1:LCH),'H').EQ.0) THEN
         LCH = LCH + 1
         CHOP(LCH:) = 'H'
      ENDIF
 
      CALL RZFILE(LUNIN,CHDIR,CHOP(1:LCH))
*
*     RZIODO calls user routine via
*        call jumpst(iquest(11))
*        call jumpx2(ibuf,ioway)
*
*     where IOWAY = 0: read buffer
*           IOWAY = 1: write buffer
*
*     and   IQUEST(1) = LUN; on exit: status flag (0 = success)
*           IQUEST(2) = number of machine words for transmission/transmitted
*           IQUEST(4) = record number
*
 
      END
