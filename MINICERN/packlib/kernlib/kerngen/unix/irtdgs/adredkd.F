*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:06  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
@PROCESS AUTODBL(NONE)
      SUBROUTINE ADREDKD(LUN,IA,ND,IREC,ISW)
      INTEGER IA(2*ND)
      INTEGER IREC,ISW
*
      READ(LUN,REC=IREC,IOSTAT=ISW) (IA(J),J=1,2*ND)
      END
