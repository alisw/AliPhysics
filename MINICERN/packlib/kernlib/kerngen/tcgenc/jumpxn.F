*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/jumpx2.F"
#elif defined(CERNLIB_QMMPW)
#include "mpwgs/jumpxn.F"
#else
      SUBROUTINE JUMPXN
C CERN PROGLIB# Z042    JUMPXN  dummy   .VERSION KERNFOR  4.23  891215

      ENTRY JUMPAD (TARGET)
      ENTRY JUMPST (TARGET)
      ENTRY JUMPX0
      ENTRY JUMPX1 (PAR1)
      ENTRY JUMPX2 (PAR1,PAR2)
      ENTRY JUMPX3 (PAR1,PAR2,PAR3)
      ENTRY JUMPX4 (PAR1,PAR2,PAR3,PAR4)

      PRINT 9001
      CALL ABEND

 9001 FORMAT (' **** Dummy JUMPAD / JUMPXN reached ****')
      END
#endif
