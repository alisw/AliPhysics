*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION IGET(MEM, IPOS)
      DIMENSION MEM(1)
      COMMON /BITSXB/ NBITPW,NBYTPW
C
C   *IGET* RETRIEVES THE VALUES STORED BY *TUP*.
C
      LW=(IPOS+NBYTPW-1)/NBYTPW
      LA=8*(IPOS-NBYTPW*(LW-1))-7
      IGET=JBYT(MEM(LW),LA,8)
      RETURN
      END
