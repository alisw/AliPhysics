*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUSTEP
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called at the end of each tracking step     *
C.    *       INWVOL is different from 0 when the track has reached    *
C.    *              a volume boundary                                 *
C.    *       ISTOP is different from 0 if the track has stopped       *
C.    *                                                                *
C.    *    ==>Called by : GTRACK                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
C.
C.    ------------------------------------------------------------------
C.
      END
