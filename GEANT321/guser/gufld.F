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
      SUBROUTINE GUFLD(X,F)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *                                                                *
C.    *    ==>Called by : GHELIX,GRKUTA                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),F(3)
C.
C.    ------------------------------------------------------------------
C.
      F(1)=0.
      F(2)=0.
      F(3)=0.
C
      END
