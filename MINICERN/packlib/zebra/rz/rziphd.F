*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZIPHD(LUN,JREC,IREC,IBUF,IRC)
      DIMENSION IBUF(JREC)
      IRC = 0
      READ(LUN,REC=IREC,ERR=99,IOSTAT=ISTAT) IBUF
      RETURN
 99   IRC = ISTAT
      END
