*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVERS(CHTITL,CHCOMP,CHOPT)
*
*     Print or return RZ version number
*
      CHARACTER*(*) CHTITL,CHCOMP,CHOPT
      CHARACTER*80  CTITLE,CCOMPL
 
#include "zebra/zunit.inc"
 
      IDATQQ = 0
      ITIMQQ = 0
*
*     Clear title and compilation date/time
*
      CTITLE = ' '
      CCOMPL = ' '
 
      CTITLE =
     +  'RZEBRA   3.01/08'
      WRITE(CCOMPL,9001) IDATQQ,ITIMQQ
9001  FORMAT(' This version created on ',I6,' at ',I4)
 
      LTITLE = LENOCC(CTITLE)
      LCOMPL = LENOCC(CCOMPL)
 
      IF(INDEX(CHOPT,'P').NE.0) THEN
 
         WRITE(IQPRNT,*) CTITLE(1:LTITLE)
         WRITE(IQPRNT,*) CCOMPL(1:LCOMPL)
 
      ELSE
 
         CHTITL = CTITLE(1:LTITLE)
         CHCOMP = CCOMPL(1:LCOMPL)
 
      ENDIF
 
      END
