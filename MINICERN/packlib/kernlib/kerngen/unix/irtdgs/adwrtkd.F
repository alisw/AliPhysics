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
      SUBROUTINE ADWRTKD(LUN,IA,ND,IREC,ISW)
      INTEGER IA(2*ND)
      INTEGER ISW
*
      WRITE(LUN,REC=IREC,IOSTAT=ISW) (IA(J),J=1,2*ND)
      END
