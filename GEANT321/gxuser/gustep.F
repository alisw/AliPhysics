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
      SUBROUTINE GUSTEP
*
************************************************************************
*                                                                      *
*             User routine called at the end of each tracking step     *
*             MEC   is the mechanism origin of the step                *
*             INWVOL is different from 0 when the track has reached    *
*                    a volume boundary                                 *
*             ISTOP is different from 0 if the track has stopped       *
*                                                                      *
************************************************************************
*
#include "geant321/gcking.inc"
#include "geant321/gcomis.inc"
*
*     -----------------------------------------------------------------
*
#if defined(CERNLIB_COMIS)
      IF(JUSTEP.NE.0)THEN
         CALL CSJCAL(JUSTEP,0,X,X,X,X,X,X,X,X,X,X)
         GO TO 99
      ENDIF
#endif
      IF(NGKINE.GT.0) THEN
        DO 5 I=1,NGKINE
          ITYPA  = GKIN(5,I)
          IF(ITYPA.NE.4)  CALL GSKING(I)
   5    CONTINUE
      ENDIF
      CALL GDEBUG
  99  END
