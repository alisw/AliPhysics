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
      SUBROUTINE UGLAST
*
************************************************************************
*                                                                      *
*            Termination routine to print histograms and statistics    *
*                                                                      *
************************************************************************
#include "geant321/gcomis.inc"
*
*     -----------------------------------------------------------------
*
#if defined(CERNLIB_COMIS)
      IF(JULAST.NE.0)THEN
         CALL CSJCAL(JULAST,0,X,X,X,X,X,X,X,X,X,X)
         GO TO 99
      ENDIF
#endif
      CALL GLAST
*
*             Close HIGZ
*
      CALL IGEND
*
  99  END
 
