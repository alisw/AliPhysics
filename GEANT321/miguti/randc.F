*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION FUNCTION RANDC(IDUMM)
#endif
#if defined(CERNLIB_SINGLE)
      FUNCTION RANDC(IDUMM)
#endif
C
      DIMENSION RND1(1)
C
      CALL GRNDM(RND1,1)
      RANDC = RND1(1)
      RETURN
      END
