*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:46:57  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:10  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:48  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/abend.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/abend.F"
#else
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END
#endif
