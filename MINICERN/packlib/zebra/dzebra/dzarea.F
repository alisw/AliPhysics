*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:05  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZAREA  (CHTEXT,IXSTOR,CLA,LLA,CHOPT)
#include "zebra/quest.inc"
#include "zebra/dzc1.inc"
      CHARACTER *(*) CLA,CHOPT,CHTEXT

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZAREA')

      ENTRY      DZLINK  (CHTEXT,IXSTOR,CLA,LLA,CHOPT)
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL DZOPT(CHOPT)

      CALL MZSDIV(IXSTOR,-1)

      CALL DZARE1(CHTEXT,CLA,LLA,CHOPT)

  999 RETURN
      END
