*
* $Id$
*
* $Log$
* Revision 1.2  1997/01/07 10:25:42  cernlib
* Remove #ifdef CERNLIB_MAIN; this shall be done via Imakefile.
*
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      PROGRAM GXINT
*
*     GEANT main program. To link with the MOTIF user interface
*     the routine GPAWPP(NWGEAN,NWPAW) should be called, whereas
*     the routine GPAW(NWGEAN,NWPAW) gives access to the basic
*     graphics version.
*
#if !defined(CERNLIB_IBM)
      PARAMETER (NWGEAN=3000000,NWPAW=1000000)
#endif
#if defined(CERNLIB_IBM)
      PARAMETER (NWGEAN=1000000,NWPAW=500000)
#endif
      COMMON/GCBANK/GEANT(NWGEAN)
      COMMON/PAWC/PAW(NWPAW)
#if defined(CERNLIB_HPUX)
      ON REAL UNDERFLOW IGNORE
#endif
*
#if defined(CERNLIB_IBM)
      CALL INITC
      CALL ERRSET(151,999,-1)
#endif
#if (defined(CERNLIB_MOTIF))&&(!defined(CERNLIB_IBM))
      CALL GPAWPP(NWGEAN,NWPAW)
#endif
#if !defined(CERNLIB_MOTIF)||defined(CERNLIB_IBM)
      CALL GPAW(NWGEAN,NWPAW)
#endif
*
      END
      SUBROUTINE QNEXT
      END
#if !defined(CERNLIB_CZ)
      SUBROUTINE CZOPEN
      END
      SUBROUTINE CZTCP
      END
      SUBROUTINE CZCLOS
      END
      SUBROUTINE CZPUTA
      END
#endif
#if defined(CERNLIB_IBM)
      FUNCTION IOSCLR()
      IOSCLR=0
      END
#endif
