*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUPARA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called every time a particle falls below    *
C.    *       parametrization threshold. This routine should create    *
C.    *       the parametrization stack, and, when this is full,       *
C.    *       parametrize the shower and track the geantinos.          *
C.    *                                                                *
C.    *    ==>Called by : GTRACK                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcparm.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
C.
C.    ------------------------------------------------------------------
C.
*      CALL GSPARA
*
      END
