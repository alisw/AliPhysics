*
* $Id$
*
* $Log$
* Revision 1.3  1997/11/05 16:00:44  mclareni
* Correct typo
*
* Revision 1.2  1997/11/05 10:13:26  mclareni
* Include in WNT
*
* Revision 1.1.1.1  1996/02/15 17:50:12  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if (defined(CERNLIB_QMSUN))&&(!defined(CERNLIB_SOLARIS))
#include "sungs/jobnam.F"
#else
      SUBROUTINE JOBNAM (ID)
C
C CERN PROGLIB# Z100    JOBNAM          .VERSION KERNFOR  4.23  891215
C ORIG. 06/12/89, JZ+FCA
C
      DIMENSION ID(2)
      CHARACTER    JOBNN*8
#if defined(CERNLIB_QMHPX)
      PARAMETER (JOBNN='HPUnix')
#elif defined(CERNLIB_QMSGI)
      PARAMETER (JOBNN='SILGRAPH')
#elif defined(CERNLIB_QMOS9)
      PARAMETER (JOBNN='MWos9')
#elif defined(CERNLIB_QMMPW)
      PARAMETER (JOBNN='MAC II  ')
#elif defined(CERNLIB_QMMPW)
      PARAMETER (JOBNN='MAC II  ')
#elif defined(CERNLIB_WINNT)
      PARAMETER (JOBNN='Windows ')
#else
      PARAMETER    (JOBNN = 'NotBatch')
#endif

      CALL UCTOH (JOBNN,ID,8,8)
      RETURN
      END
#endif
