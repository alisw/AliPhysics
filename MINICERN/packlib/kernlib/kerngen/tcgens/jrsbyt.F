*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:12  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/jrsbyt.F"
#else
      FUNCTION JRSBYT (VAL,X,JTH,NBITS)
C
C CERN PROGLIB# M421    JRSBYT          .VERSION KERNFOR  4.07  830624
C ORIG. 03/05/72 JZ
C
      INTEGER   VAL, X
C
#include "kerngen/q_jbit.inc"
C
      JRSBYT= JBYT (X,JTH,NBITS)
      CALL SBYT (VAL,X,JTH,NBITS)
      RETURN
      END
#endif
