*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE GUOUT
*
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called at the end of each event.            *
C.    *                                                                *
C.    ******************************************************************
C.
C.
#include "geant321/gcomis.inc"
C.    ------------------------------------------------------------------
C.
#if defined(CERNLIB_COMIS)
      IF(JUOUT.NE.0)THEN
         CALL CSJCAL(JUOUT,0,X,X,X,X,X,X,X,X,X,X)
         GO TO 99
      ENDIF
#endif
*
  99  END
C
